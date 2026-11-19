! Numerical core of the perfectly-plastic reconstruction: the Godunov update of the
! static eikonal |grad u| = F, fast sweeping, the right-hand side of the w = H^2
! equation in :flat and :full mode, the Godunov update of the surface equation
! (mode "surface", constant tau_b or the Coulomb yield law), and the Godunov residuals.
!
! All routines are pure array kernels on explicit-shape (nx,ny) arrays, i = x, j = y.
! They carry no state, so several solver instances can run concurrently.
module plastic_eikonal

    use plastic_kinds, only: dp, BIG

    implicit none

    private
    public :: plastic_par_t, plastic_coulomb_law_t
    public :: plastic_flotation_thickness
    public :: godunov_eikonal, godunov_surface, godunov_coulomb, coulomb_neff
    public :: sweep_once, sweep_eikonal, update_at
    public :: KERNEL_EIKONAL, KERNEL_SURFACE, KERNEL_COULOMB, MARGIN_DMIN
    public :: grad_central
    public :: flotation_bc, rhs_flat, rhs_full
    public :: godunov_residual, godunov_residual_surface

    integer,  parameter :: KERNEL_EIKONAL = 1, KERNEL_SURFACE = 2, KERNEL_COULOMB = 3

    ! Floor of the spacing factor to an ice-free neighbour, (1 - margin_frac) >= MARGIN_DMIN:
    ! keeps the update well-conditioned as the margin approaches the ice-cell centre.
    real(dp), parameter :: MARGIN_DMIN = 0.01_dp

    ! Physical constants and solver controls.
    type plastic_par_t
        real(dp) :: rho_i  = 917.0_dp       ! ice density [kg m-3]
        real(dp) :: rho_w  = 1025.0_dp      ! sea-water density [kg m-3]
        real(dp) :: g      = 9.80665_dp     ! gravity [m s-2]
        real(dp) :: H_min  = 1.0_dp         ! thickness floor in the :full RHS [m]
        character(len=8) :: mode = "surface"   ! "surface" | "flat" | "full"
        integer  :: max_sweeps = 200        ! inner fast-sweeping cap
        integer  :: n_outer    = 100        ! :full outer fixed-point cap
        real(dp) :: tol        = 1.0e-6_dp  ! inner early exit, on max decrease of w [m2] (z_s [m] in "surface")
        real(dp) :: outer_tol  = 1.0e-3_dp  ! :full outer exit, on max change of H [m]
        real(dp) :: relax      = 0.5_dp     ! :full under-relaxation in (0,1]
        ! Margin position between an ice cell and an ice-free neighbour, without margin_frac:
        ! "face" (half way, margin_frac = 0.5) or "centre" (at the ice-free cell centre,
        ! margin_frac = 0, as PlasticIceSheet.jl before v0.2).
        character(len=8) :: margin = "face"
        ! Sweep ordering: "serial" (row by row) or "wavefront" (the cells of each
        ! anti-diagonal in parallel with OpenMP; the same Gauss-Seidel dependencies, hence
        ! bit-identical results for any thread count).
        character(len=10) :: sweep = "serial"
        ! Warm start: fall back to a cold start when more than this fraction of the
        ! ice cells is new to the mask (see plastic_solve).
        real(dp) :: warm_frac_max = 0.1_dp
    end type plastic_par_t

    ! The Coulomb yield law of plastic_solve_coulomb, tau_b = mu N_eff (coulomb_neff):
    ! N = rho_i g H max((1 - H_f/H)^p, delta) (Leguy et al. 2014 with a floor, as Bueler and
    ! van Pelt 2015), 1/N_eff = 1/N + 1/(rho_i g H_N).
    type plastic_coulomb_law_t
        real(dp) :: p     = 1.0_dp          ! exponent of the ocean connection, in [0, 1]
        real(dp) :: H_N   = 0.0_dp          ! [m] cap of N_eff/(rho_i g) (required, > 0)
        real(dp) :: delta = 0.05_dp         ! floor of N as a fraction of the overburden, in (0, 1]
    end type plastic_coulomb_law_t

contains

    ! Flotation thickness at a marine margin, (rho_w/rho_i) max(z_sl - z_bed, 0).
    elemental function plastic_flotation_thickness(z_bed, z_sl, par) result(H_flt)
        real(dp),            intent(in) :: z_bed, z_sl
        type(plastic_par_t), intent(in) :: par
        real(dp) :: H_flt
        H_flt = (par%rho_w / par%rho_i) * max(z_sl - z_bed, 0.0_dp)
    end function plastic_flotation_thickness

    ! Smallest u satisfying ((u-a)+/da)^2 + ((u-b)+/db)^2 = F^2, where a, b are the
    ! upwind (smaller) neighbour values along x and y. A neighbour >= BIG carries no
    ! information, and the update falls back to the one-sided form.
    pure function godunov_eikonal(a, da, b, db, F) result(u)
        real(dp), intent(in) :: a, da, b, db, F
        real(dp) :: u
        real(dp) :: ida2, idb2, qa, qb, qc, disc

        if (a >= BIG .and. b >= BIG) then
            u = BIG
            return
        else if (a >= BIG) then
            u = b + F * db
            return
        else if (b >= BIG) then
            u = a + F * da
            return
        end if

        ! Two-sided update: qa u^2 + qb u + qc = 0, accepted if both axes are upwind.
        ida2 = 1.0_dp / (da * da)
        idb2 = 1.0_dp / (db * db)
        qa   = ida2 + idb2
        qb   = -2.0_dp * (a * ida2 + b * idb2)
        qc   = a * a * ida2 + b * b * idb2 - F * F
        disc = qb * qb - 4.0_dp * qa * qc
        if (disc >= 0.0_dp) then
            u = (-qb + sqrt(disc)) / (2.0_dp * qa)
            if (u >= a .and. u >= b) return
        end if

        ! Otherwise the smaller single-axis update governs.
        u = min(a + F * da, b + F * db)
    end function godunov_eikonal

    ! Godunov update of the plastic surface equation |grad z_s| = k / H, H = z_s - zb,
    ! k = tau_b / (rho_i g), for the surface u = z_s of a cell with bed zb, from the
    ! upwind (smaller) neighbour surfaces a (x, spacing da) and b (y, spacing db).
    !
    ! Each axis contributes the difference (u - a)+ times the edge-mean thickness
    ! (H + (a - zb)+)/2, with the neighbour thickness measured from this cell's bed:
    !
    !   P_a(H) = (H - alpha)+ (H + alpha+) / 2,  alpha = a - zb,
    !   (P_a/da)^2 + (P_b/db)^2 = k^2.
    !
    ! P_a is increasing in H (so the Hamiltonian is proper) and nonincreasing in a
    ! (monotone scheme): fast sweeping applies, with no outer iteration. The edge mean
    ! integrates the margin singularity exactly on a flat bed (it reproduces the Nye
    ! profile, like the flat w = H^2 update); a cell-centred k/H would put the first
    ! cell at sqrt(k dx) instead of sqrt(2 k dx). Where both neighbours stand above zb
    ! this is the w-eikonal with w_a = alpha^2 and F = 2k (closed form, for the increment
    ! over the smaller w_a to avoid cancellation); otherwise the
    ! two-sided root is found by Newton from the smaller one-sided root (P^2 is convex
    ! and increasing there, so Newton decreases monotonically to the root).
    pure function godunov_surface(a, da, b, db, zb, k) result(u)
        real(dp), intent(in) :: a, da, b, db, zb, k
        real(dp) :: u
        real(dp) :: al, be, wa, wb, wm, dm, dother, c, dw, qa
        real(dp) :: hA, hB, h0, other, h, pa, pb, dpa, dpb, phi, dphi, dh
        integer  :: it

        if (a >= BIG .and. b >= BIG) then
            u = BIG
            return
        end if
        if (.not. (k > 0.0_dp)) then
            u = max(zb, min(a, b))
            return
        end if

        al = 0.0_dp
        be = 0.0_dp
        if (a < BIG) al = a - zb
        if (b < BIG) be = b - zb

        ! Both neighbour surfaces at or above this bed: the w-eikonal (w = H^2, F = 2k)
        ! in closed form, solved for the increment over the smaller neighbour w
        ! (well-conditioned: no cancellation against w itself).
        if (al >= 0.0_dp .and. be >= 0.0_dp) then
            if (a >= BIG) then
                u = zb + sqrt(be * be + 2.0_dp * k * db)
            else if (b >= BIG) then
                u = zb + sqrt(al * al + 2.0_dp * k * da)
            else
                wa = al * al
                wb = be * be
                if (wa <= wb) then
                    wm = wa; dm = da; c = wb - wa; dother = db
                else
                    wm = wb; dm = db; c = wa - wb; dother = da
                end if
                dw = 2.0_dp * k * dm                      ! one-sided in the smaller neighbour
                if (dw > c) then                          ! both axes upwind: two-sided root in (c, dw)
                    qa = 1.0_dp / (dm * dm) + 1.0_dp / (dother * dother)
                    dw = (c / (dother * dother) + sqrt(max(qa * (2.0_dp * k)**2 &
                          - (c / (dm * dother))**2, 0.0_dp))) / qa
                end if
                u = zb + sqrt(wm + dw)
            end if
            return
        end if

        ! One-sided roots P(H) = k d; the smaller one bounds the two-sided root.
        hA = BIG
        hB = BIG
        if (a < BIG) hA = one_sided(al, k * da)
        if (b < BIG) hB = one_sided(be, k * db)
        if (hA <= hB) then
            h0 = hA
            other = b
        else
            h0 = hB
            other = a
        end if
        if (other >= BIG .or. zb + h0 <= other) then
            u = zb + h0
            return
        end if

        ! Two-sided: Newton on phi(H) = (P_a/da)^2 + (P_b/db)^2 - k^2 from H = h0.
        h = h0
        do it = 1, 60
            pa  = (h - al) * (h + max(al, 0.0_dp)) / 2.0_dp
            pb  = (h - be) * (h + max(be, 0.0_dp)) / 2.0_dp
            dpa = (2.0_dp * h + max(al, 0.0_dp) - al) / 2.0_dp
            dpb = (2.0_dp * h + max(be, 0.0_dp) - be) / 2.0_dp
            phi  = (pa / da)**2 + (pb / db)**2 - k * k
            dphi = 2.0_dp * (pa * dpa / (da * da) + pb * dpb / (db * db))
            if (.not. (phi > 0.0_dp .and. dphi > 0.0_dp)) exit
            dh = phi / dphi
            h  = h - dh
            if (dh <= 1.0e-14_dp * h) exit
        end do
        u = zb + h

    contains

        ! Root H > 0 of (H - al)(H + al+)/2 = kd.
        pure function one_sided(al, kd) result(hh)
            real(dp), intent(in) :: al, kd
            real(dp) :: hh
            if (al >= 0.0_dp) then
                hh = sqrt(al * al + 2.0_dp * kd)
            else
                hh = 4.0_dp * kd / (sqrt(al * al + 8.0_dp * kd) - al)
            end if
        end function one_sided

    end function godunov_surface

    ! Effective-pressure head of the Coulomb yield law, N_eff/(rho_i g) [m], at thickness H
    ! with flotation thickness hf: N = rho_i g H max((1 - hf/H)^p, delta) (Leguy et al. 2014,
    ! (1 - hf/H)^p = 0 at or below flotation, with the floor delta of the overburden of
    ! Bueler and van Pelt 2015), capped smoothly by N_max = rho_i g H_N: 1/N_eff = 1/N + 1/N_max.
    ! On land (hf = 0) N = rho_i g H. 0 for H <= 0.
    elemental function coulomb_neff(H, hf, law) result(ne)
        real(dp),                    intent(in) :: H, hf
        type(plastic_coulomb_law_t), intent(in) :: law
        real(dp) :: ne, n
        if (.not. H > 0.0_dp) then
            ne = 0.0_dp
            return
        end if
        if (hf > 0.0_dp) then
            n = 0.0_dp
            if (H > hf) then
                if (law%p == 1.0_dp) then
                    n = H - hf
                else
                    n = H * (1.0_dp - hf / H)**law%p
                end if
            end if
            n = max(n, law%delta * H)
        else
            n = H
        end if
        ne = n * law%H_N / (n + law%H_N)
    end function coulomb_neff

    ! Godunov update of the plastic surface equation with the Coulomb yield law,
    ! tau_b = mu N_eff(H), i.e. the slope |grad z_s| = s(H) = mu N_eff(H) / (rho_i g H)
    ! (coulomb_neff; mu = c_f lambda of the cell), for the surface u = z_s of a cell with bed
    ! zb and flotation thickness hf, from the upwind neighbour surfaces a (x, spacing da)
    ! and b (y, spacing db); law the exponent p, cap H_N and floor delta.
    !
    ! The upwind structure of godunov_surface: each axis carries the slope at its edge-mean
    ! thickness h_a = (H + (a - zb)+)/2 (the neighbour thickness measured from this bed),
    !
    !   T_a(u) = (u - a)+ / (da s(h_a)),   T_a^2 + T_b^2 = 1.
    !
    ! For s = k/h (constant tau_b = rho_i g k) this is godunov_surface's equation; for a
    ! constant slope (H << H_N on land) the eikonal with F = mu. On a marine cell (hf > 0)
    ! a neighbour below this cell's flotation surface is raised to it, a <- max(a, zb + hf)
    ! (grounded ice is at least at flotation: the step down to a lower neighbour is a
    ! grounding-line step without basal stress), so u >= zb + hf.
    !
    ! Uniqueness and monotonicity of the update: with 0 <= p <= 1, 0 < delta <= 1 and (after
    ! that clamp) every upwind a >= zb + hf, each T_a is strictly increasing in u > a and
    ! nonincreasing in a. Off the floor, d ln s/dh <= p hf/(h (h - hf)) and h - (a - zb) <= h - hf
    ! give (u - a) s'/(2 s) <= p hf/h < 1; on the floor (N = delta rho_i g h) s = mu delta H_N/
    ! (delta h + H_N) is nonincreasing; -d ln s/dh < 1/h in both gives the bound in a (a < zb on
    ! land: h does not depend on a). So sum T^2 - 1 has one root above the upwind values, and
    ! the update is nondecreasing in a and b (a monotone scheme; fast sweeping applies). On
    ! land (hf = 0) s = mu H_N/(h + H_N) is nonincreasing, which needs no clamp.
    !
    ! Uniqueness of the solution (causality): the floor bounds the slope from below on every
    ! grounded cell, s(h) >= s_min(h) = mu delta H_N/(delta h + H_N) > 0 (mu > 0), so the update
    ! lies above the smaller upwind value by at least d s_min/sqrt(2) (one axis has T >= 1/sqrt(2)).
    ! Every fixed point is then ordered from the margin inward (the smallest ice-cell value takes
    ! only margin values), i.e. unique, and continuous in the bed and sea level. Without the
    ! floor (delta = 0, N = 0 at flotation) the slope vanishes at flotation: cells at their
    ! flotation surfaces can support each other with zero increments (u = a after the clamp),
    ! the discrete problem has several fixed points (a marine trough at flotation or grounded),
    ! and the sweep's choice flips on decimetre changes of the bed.
    !
    ! The root is bracketed by lo = the smaller (clamped) neighbour (or zb + hf) and hi = min
    ! over the axes of a + mu d (s <= mu, so T >= 1 there), and found by Newton safeguarded by
    ! bisection (rtsafe).
    pure function godunov_coulomb(a, da, b, db, zb, mu, hf, law) result(u)
        real(dp),                    intent(in) :: a, da, b, db, zb, mu, hf
        type(plastic_coulomb_law_t), intent(in) :: law
        real(dp) :: u
        real(dp) :: ac, bc, lo, hi, phi, dphi, un, tol, dxo, dxs
        logical  :: ha, hb
        integer  :: it

        ha = a < BIG
        hb = b < BIG
        if (.not. (ha .or. hb)) then
            u = BIG
            return
        end if
        ac = a
        bc = b
        if (hf > 0.0_dp) then
            if (ha) ac = max(a, zb + hf)
            if (hb) bc = max(b, zb + hf)
        end if
        lo = max(min(ac, bc), zb + hf)
        if (.not. (mu > 0.0_dp)) then
            u = lo
            return
        end if
        hi = BIG
        if (ha) hi = min(hi, ac + mu * da)
        if (hb) hi = min(hi, bc + mu * db)
        call eval(lo, phi, dphi)
        if (hi <= lo .or. phi >= 0.0_dp) then
            u = lo
            return
        end if

        ! Newton from the upper bracket, bisection when a step leaves (lo, hi) or does not halve
        ! the step before last (rtsafe; the kink of s at the floor can make plain Newton cycle).
        u = hi
        call eval(u, phi, dphi)
        dxo = hi - lo
        dxs = dxo
        do it = 1, 100
            if (phi >= 0.0_dp) then
                hi = u
            else
                lo = u
            end if
            un = 0.5_dp * (lo + hi)
            if (dphi > 0.0_dp) then
                if (u - phi / dphi > lo .and. u - phi / dphi < hi .and. abs(2.0_dp * phi) <= abs(dxo * dphi)) &
                    un = u - phi / dphi
            end if
            dxo = dxs
            dxs = abs(un - u)
            tol = 1.0e-12_dp * max(1.0_dp, u - zb)
            if (abs(un - u) <= tol .or. hi - lo <= tol) then
                u = un
                exit
            end if
            u = un
            call eval(u, phi, dphi)
        end do

    contains

        ! phi(uu) = T_a^2 + T_b^2 - 1 and its derivative.
        pure subroutine eval(uu, f, df)
            real(dp), intent(in)  :: uu
            real(dp), intent(out) :: f, df
            f  = -1.0_dp
            df = 0.0_dp
            if (ha) call axis(uu, ac, da, f, df)
            if (hb) call axis(uu, bc, db, f, df)
        end subroutine eval

        pure subroutine axis(uu, aa, dd, f, df)
            real(dp), intent(in)    :: uu, aa, dd
            real(dp), intent(inout) :: f, df
            real(dp) :: h, q, dlq, t
            if (.not. uu > aa) return
            h = 0.5_dp * (uu - zb + max(aa - zb, 0.0_dp))
            call ratio(h, q, dlq)
            if (.not. q > 0.0_dp) then
                f = BIG
                return
            end if
            t  = (uu - aa) / (dd * mu * q)
            f  = f + t * t
            df = df + 2.0_dp * t * (1.0_dp - 0.5_dp * (uu - aa) * dlq) / (dd * mu * q)
        end subroutine axis

        ! q = s/mu = N_eff/(rho_i g h) at h, and d ln q/dh (on the floor r = delta is constant:
        ! d ln q/dh = (c - 1)/h).
        pure subroutine ratio(h, q, dlq)
            real(dp), intent(in)  :: h
            real(dp), intent(out) :: q, dlq
            real(dp) :: r, n, c
            if (hf > 0.0_dp) then
                r = 0.0_dp
                if (h > hf) then
                    if (law%p == 1.0_dp) then
                        r = 1.0_dp - hf / h
                    else
                        r = (1.0_dp - hf / h)**law%p
                    end if
                end if
                if (r < law%delta) then
                    n   = h * law%delta
                    c   = law%H_N / (n + law%H_N)
                    q   = law%delta * c
                    dlq = (c - 1.0_dp) / h
                else
                    n   = h * r
                    c   = law%H_N / (n + law%H_N)
                    q   = r * c
                    dlq = ((1.0_dp - law%p) / h + law%p / (h - hf)) * c - 1.0_dp / h
                end if
            else
                q   = law%H_N / (h + law%H_N)
                dlq = -1.0_dp / (h + law%H_N)
            end if
        end subroutine ratio

    end function godunov_coulomb

    ! The update of kernel `kind`: KERNEL_EIKONAL (F the slowness), KERNEL_SURFACE (F =
    ! tau_b/(rho_i g), zb the cell's bed) or KERNEL_COULOMB (F = mu, hf the cell's flotation
    ! thickness, law the Coulomb law).
    pure function kernel_update(kind, a, da, b, db, zb, F, hf, law) result(u)
        integer,                     intent(in) :: kind
        real(dp),                    intent(in) :: a, da, b, db, zb, F, hf
        type(plastic_coulomb_law_t), intent(in) :: law
        real(dp) :: u
        select case (kind)
        case (KERNEL_EIKONAL)
            u = godunov_eikonal(a, da, b, db, F)
        case (KERNEL_SURFACE)
            u = godunov_surface(a, da, b, db, zb, F)
        case default
            u = godunov_coulomb(a, da, b, db, zb, F, hf, law)
        end select
    end function kernel_update

    ! Godunov update at ice cell (i,j) from its four neighbours (clamp-to-edge), with the
    ! eikonal (kind KERNEL_EIKONAL, F = slowness), the surface kernel (KERNEL_SURFACE,
    ! F = tau_b/(rho_i g), zbij the cell's bed) or the Coulomb surface kernel
    ! (KERNEL_COULOMB, F = mu, hfij the cell's flotation thickness, law the Coulomb law).
    ! The spacing to a neighbour is d dfac(nb):
    ! dfac = 1 on ice cells, < 1 on ice-free cells whose margin lies between the two
    ! centres (sub-cell margin position). With equal spacings on both sides of an axis
    ! the upwind neighbour is the smaller one; otherwise the axis term is the Godunov
    ! max_s (u - a_s)+ / d_s over the two sides, i.e. the update is the minimum over the
    ! side combinations. Each candidate is monotone, so the minimum is.
    pure function update_at(kind, nx, ny, i, j, u, F, zbij, dfac, dx, dy, hfij, law) result(unew)
        integer,                     intent(in) :: kind, nx, ny, i, j
        real(dp),                    intent(in) :: u(nx, ny), F, zbij, dfac(nx, ny), dx, dy, hfij
        type(plastic_coulomb_law_t), intent(in) :: law
        real(dp) :: unew, ax(2), dxs(2), ay(2), dys(2), c
        integer  :: im, ip, jm, jp, nxc, nyc, m, q

        im = max(i - 1, 1)
        ip = min(i + 1, nx)
        jm = max(j - 1, 1)
        jp = min(j + 1, ny)
        if (dfac(im, j) == dfac(ip, j)) then
            nxc = 1
            ax(1)  = min(u(im, j), u(ip, j))
            dxs(1) = dx * dfac(im, j)
        else
            nxc = 2
            ax(1) = u(im, j); dxs(1) = dx * dfac(im, j)
            ax(2) = u(ip, j); dxs(2) = dx * dfac(ip, j)
        end if
        if (dfac(i, jm) == dfac(i, jp)) then
            nyc = 1
            ay(1)  = min(u(i, jm), u(i, jp))
            dys(1) = dy * dfac(i, jm)
        else
            nyc = 2
            ay(1) = u(i, jm); dys(1) = dy * dfac(i, jm)
            ay(2) = u(i, jp); dys(2) = dy * dfac(i, jp)
        end if

        unew = BIG
        do q = 1, nyc
            do m = 1, nxc
                c = kernel_update(kind, ax(m), dxs(m), ay(q), dys(q), zbij, F, hfij, law)
                unew = min(unew, c)
            end do
        end do
    end function update_at

    ! One Gauss-Seidel sweep in the ordering (xdir, ydir) = (+-1, +-1). Updates ice
    ! cells only, keeps the monotone min, and returns the largest decrease of u. Cells
    ! with equal spacing on both sides of each axis (all but those at a sub-cell margin)
    ! take the single-candidate update inline; the others go through update_at.
    function sweep_once(kind, nx, ny, u, F, zb, mask, dfac, dx, dy, xdir, ydir, hf, law) result(change)
        integer,  intent(in)    :: kind, nx, ny
        real(dp), intent(inout) :: u(nx, ny)
        real(dp), intent(in)    :: F(nx, ny), zb(nx, ny), dfac(nx, ny), hf(nx, ny)
        type(plastic_coulomb_law_t), intent(in) :: law
        logical,  intent(in)    :: mask(nx, ny)
        real(dp), intent(in)    :: dx, dy
        integer,  intent(in)    :: xdir, ydir
        real(dp) :: change
        integer  :: i, j, i0, i1, j0, j1, im, ip, jm, jp
        real(dp) :: unew, a, b

        if (xdir > 0) then
            i0 = 1;  i1 = nx
        else
            i0 = nx; i1 = 1
        end if
        if (ydir > 0) then
            j0 = 1;  j1 = ny
        else
            j0 = ny; j1 = 1
        end if

        change = 0.0_dp
        do j = j0, j1, ydir
            jm = max(j - 1, 1)
            jp = min(j + 1, ny)
            do i = i0, i1, xdir
                if (.not. mask(i, j)) cycle
                im = max(i - 1, 1)
                ip = min(i + 1, nx)
                if (dfac(im, j) == dfac(ip, j) .and. dfac(i, jm) == dfac(i, jp)) then
                    a = min(u(im, j), u(ip, j))
                    b = min(u(i, jm), u(i, jp))
                    unew = kernel_update(kind, a, dx * dfac(im, j), b, dy * dfac(i, jm), zb(i, j), F(i, j), &
                                         hf(i, j), law)
                else
                    unew = update_at(kind, nx, ny, i, j, u, F(i, j), zb(i, j), dfac, dx, dy, hf(i, j), law)
                end if
                if (unew < u(i, j)) then
                    change  = max(change, u(i, j) - unew)
                    u(i, j) = unew
                end if
            end do
        end do
    end function sweep_once

    ! sweep_once in wavefront order: the anti-diagonals p + q = const of the sweep
    ! ordering, one after the other, the cells of a diagonal in parallel. A cell reads its
    ! four neighbours: (m-1,q), (m,q-1) on the previous diagonal (updated) and (m+1,q),
    ! (m,q+1) on the next (not yet), exactly as in the row-by-row sweep.
    function sweep_once_wavefront(kind, nx, ny, u, F, zb, mask, dfac, dx, dy, xdir, ydir, hf, law) result(change)
        integer,  intent(in)    :: kind, nx, ny
        real(dp), intent(inout) :: u(nx, ny)
        real(dp), intent(in)    :: F(nx, ny), zb(nx, ny), dfac(nx, ny), hf(nx, ny)
        type(plastic_coulomb_law_t), intent(in) :: law
        logical,  intent(in)    :: mask(nx, ny)
        real(dp), intent(in)    :: dx, dy
        integer,  intent(in)    :: xdir, ydir
        real(dp) :: change
        integer  :: d, m, q, i, j, i0, j0, im, ip, jm, jp
        real(dp) :: unew, a, b

        i0 = 1; if (xdir < 0) i0 = nx
        j0 = 1; if (ydir < 0) j0 = ny
        change = 0.0_dp
        !$omp parallel private(d, m, q, i, j, im, ip, jm, jp, unew, a, b) reduction(max:change)
        do d = 2, nx + ny
            !$omp do schedule(static)
            do q = max(1, d - nx), min(ny, d - 1)
                m = d - q
                i = i0 + (m - 1) * xdir
                j = j0 + (q - 1) * ydir
                if (.not. mask(i, j)) cycle
                im = max(i - 1, 1); ip = min(i + 1, nx)
                jm = max(j - 1, 1); jp = min(j + 1, ny)
                if (dfac(im, j) == dfac(ip, j) .and. dfac(i, jm) == dfac(i, jp)) then
                    a = min(u(im, j), u(ip, j))
                    b = min(u(i, jm), u(i, jp))
                    unew = kernel_update(kind, a, dx * dfac(im, j), b, dy * dfac(i, jm), zb(i, j), F(i, j), &
                                         hf(i, j), law)
                else
                    unew = update_at(kind, nx, ny, i, j, u, F(i, j), zb(i, j), dfac, dx, dy, hf(i, j), law)
                end if
                if (unew < u(i, j)) then
                    change  = max(change, u(i, j) - unew)
                    u(i, j) = unew
                end if
            end do
            !$omp end do
        end do
        !$omp end parallel
    end function sweep_once_wavefront

    ! Solve |grad u| = F on the ice cells (mask), with Dirichlet values bc on the
    ! ice-free cells, placed at d dfac from the neighbouring ice-cell centres. u starts
    ! from BIG on ice cells (the monotone scheme only decreases u). The four diagonal
    ! orderings alternate; the early exit needs a full cycle of four. Returns the number
    ! of sweeps and whether tol was reached. With zb present, solves the surface
    ! equation instead (u = z_s, F = tau_b/(rho_i g), see godunov_surface); with zb, hf
    ! and law present, the Coulomb surface equation (F = mu, hf the flotation thickness,
    ! see godunov_coulomb).
    subroutine sweep_eikonal(nx, ny, u, F, mask, bc, dfac, dx, dy, max_sweeps, tol, n_sweeps, ok, zb, wavefront, &
                             hf, law)
        integer,  intent(in)  :: nx, ny
        real(dp), intent(out) :: u(nx, ny)
        real(dp), intent(in)  :: F(nx, ny), bc(nx, ny), dfac(nx, ny)
        logical,  intent(in)  :: mask(nx, ny)
        real(dp), intent(in)  :: dx, dy, tol
        integer,  intent(in)  :: max_sweeps
        integer,  intent(out) :: n_sweeps
        logical,  intent(out) :: ok
        real(dp), optional, intent(in) :: zb(nx, ny), hf(nx, ny)
        logical,  optional, intent(in) :: wavefront
        type(plastic_coulomb_law_t), optional, intent(in) :: law
        integer, parameter :: xd(4) = [1, -1, -1, 1]
        integer, parameter :: yd(4) = [1, 1, -1, -1]
        integer  :: i, j, s, k, kind
        real(dp) :: change
        type(plastic_coulomb_law_t) :: lw
        logical  :: wave

        do j = 1, ny
            do i = 1, nx
                if (mask(i, j)) then
                    u(i, j) = BIG
                else
                    u(i, j) = bc(i, j)
                end if
            end do
        end do

        kind = KERNEL_EIKONAL
        if (present(zb)) kind = KERNEL_SURFACE
        if (present(zb) .and. present(hf) .and. present(law)) then
            kind = KERNEL_COULOMB
            lw   = law
        end if
        wave = .false.
        if (present(wavefront)) wave = wavefront
        ok = .false.
        n_sweeps = 0
        do s = 1, max_sweeps
            k = mod(s - 1, 4) + 1
            select case (kind)
            case (KERNEL_COULOMB)
                if (wave) then
                    change = sweep_once_wavefront(kind, nx, ny, u, F, zb, mask, dfac, dx, dy, xd(k), yd(k), hf, lw)
                else
                    change = sweep_once(kind, nx, ny, u, F, zb, mask, dfac, dx, dy, xd(k), yd(k), hf, lw)
                end if
            case (KERNEL_SURFACE)
                if (wave) then
                    change = sweep_once_wavefront(kind, nx, ny, u, F, zb, mask, dfac, dx, dy, xd(k), yd(k), F, lw)
                else
                    change = sweep_once(kind, nx, ny, u, F, zb, mask, dfac, dx, dy, xd(k), yd(k), F, lw)
                end if
            case default
                if (wave) then
                    change = sweep_once_wavefront(kind, nx, ny, u, F, F, mask, dfac, dx, dy, xd(k), yd(k), F, lw)
                else
                    change = sweep_once(kind, nx, ny, u, F, F, mask, dfac, dx, dy, xd(k), yd(k), F, lw)
                end if
            end select
            n_sweeps = s
            if (s >= 4 .and. change < tol) then
                ok = .true.
                exit
            end if
        end do
    end subroutine sweep_eikonal

    ! Central-difference gradient with clamp-to-edge (one-sided at the array edges).
    subroutine grad_central(nx, ny, f, dx, dy, fx, fy)
        integer,  intent(in)  :: nx, ny
        real(dp), intent(in)  :: f(nx, ny), dx, dy
        real(dp), intent(out) :: fx(nx, ny), fy(nx, ny)
        integer :: i, j, im, ip, jm, jp

        do j = 1, ny
            jm = max(j - 1, 1)
            jp = min(j + 1, ny)
            do i = 1, nx
                im = max(i - 1, 1)
                ip = min(i + 1, nx)
                if (ip > im) then
                    fx(i, j) = (f(ip, j) - f(im, j)) / (real(ip - im, dp) * dx)
                else
                    fx(i, j) = 0.0_dp
                end if
                if (jp > jm) then
                    fy(i, j) = (f(i, jp) - f(i, jm)) / (real(jp - jm, dp) * dy)
                else
                    fy(i, j) = 0.0_dp
                end if
            end do
        end do
    end subroutine grad_central

    ! Margin Dirichlet value of w = H^2: the squared flotation thickness (0 on land).
    subroutine flotation_bc(nx, ny, z_bed, z_sl, par, bc)
        integer,             intent(in)  :: nx, ny
        real(dp),            intent(in)  :: z_bed(nx, ny), z_sl(nx, ny)
        type(plastic_par_t), intent(in)  :: par
        real(dp),            intent(out) :: bc(nx, ny)
        integer  :: i, j
        real(dp) :: h

        do j = 1, ny
            do i = 1, nx
                h = plastic_flotation_thickness(z_bed(i, j), z_sl(i, j), par)
                bc(i, j) = h * h
            end do
        end do
    end subroutine flotation_bc

    ! :flat right-hand side, G = 2 tau_b / (rho_i g) (independent of w).
    subroutine rhs_flat(nx, ny, tau_b, par, G)
        integer,             intent(in)  :: nx, ny
        real(dp),            intent(in)  :: tau_b(nx, ny)
        type(plastic_par_t), intent(in)  :: par
        real(dp),            intent(out) :: G(nx, ny)
        integer  :: i, j
        real(dp) :: rg

        rg = par%rho_i * par%g
        do j = 1, ny
            do i = 1, nx
                G(i, j) = 2.0_dp * tau_b(i, j) / rg
            end do
        end do
    end subroutine rhs_flat

    ! :full right-hand side, F = sqrt(G^2 - 4 sqrt(w) (grad z_b . grad w) - 4 w |grad z_b|^2),
    ! with sqrt(w) floored at H_min and F^2 floored at (1e-3 G)^2. zbx, zby is the
    ! (precomputed) bed gradient; grad w is taken here, central with clamp-to-edge.
    subroutine rhs_full(nx, ny, w, G, zbx, zby, dx, dy, H_min, F)
        integer,  intent(in)  :: nx, ny
        real(dp), intent(in)  :: w(nx, ny), G(nx, ny), zbx(nx, ny), zby(nx, ny)
        real(dp), intent(in)  :: dx, dy, H_min
        real(dp), intent(out) :: F(nx, ny)
        integer  :: i, j, im, ip, jm, jp
        real(dp) :: wx, wy, sqrtw, grad_zb2, dot_zbw, Gij, hmin2

        hmin2 = H_min**2
        do j = 1, ny
            jm = max(j - 1, 1)
            jp = min(j + 1, ny)
            do i = 1, nx
                im = max(i - 1, 1)
                ip = min(i + 1, nx)
                if (ip > im) then
                    wx = (w(ip, j) - w(im, j)) / (real(ip - im, dp) * dx)
                else
                    wx = 0.0_dp
                end if
                if (jp > jm) then
                    wy = (w(i, jp) - w(i, jm)) / (real(jp - jm, dp) * dy)
                else
                    wy = 0.0_dp
                end if
                sqrtw    = sqrt(max(w(i, j), hmin2))
                grad_zb2 = zbx(i, j) * zbx(i, j) + zby(i, j) * zby(i, j)
                dot_zbw  = zbx(i, j) * wx + zby(i, j) * wy
                Gij      = G(i, j)
                F(i, j)  = sqrt(max(Gij * Gij - 4.0_dp * sqrtw * dot_zbw &
                                    - 4.0_dp * w(i, j) * grad_zb2, (1.0e-3_dp * Gij)**2))
            end do
        end do
    end subroutine rhs_full

    ! Godunov residual of a field w for the right-hand side F: w - (Godunov update) on
    ! ice cells, w - bc on ice-free cells. Zero at a converged solution.
    subroutine godunov_residual(nx, ny, w, F, mask, bc, dfac, dx, dy, r)
        integer,  intent(in)  :: nx, ny
        real(dp), intent(in)  :: w(nx, ny), F(nx, ny), bc(nx, ny), dfac(nx, ny), dx, dy
        logical,  intent(in)  :: mask(nx, ny)
        real(dp), intent(out) :: r(nx, ny)
        integer :: i, j

        do j = 1, ny
            do i = 1, nx
                if (mask(i, j)) then
                    r(i, j) = w(i, j) - update_at(KERNEL_EIKONAL, nx, ny, i, j, w, F(i, j), 0.0_dp, &
                                                  dfac, dx, dy, 0.0_dp, plastic_coulomb_law_t())
                else
                    r(i, j) = w(i, j) - bc(i, j)
                end if
            end do
        end do
    end subroutine godunov_residual

    ! Godunov residual of the surface equation: u - (surface update) on ice cells,
    ! u - bc on ice-free cells (u = z_s, K = tau_b/(rho_i g), bc = z_b + H_flt).
    subroutine godunov_residual_surface(nx, ny, u, K, zb, mask, bc, dfac, dx, dy, r)
        integer,  intent(in)  :: nx, ny
        real(dp), intent(in)  :: u(nx, ny), K(nx, ny), zb(nx, ny), bc(nx, ny), dfac(nx, ny), dx, dy
        logical,  intent(in)  :: mask(nx, ny)
        real(dp), intent(out) :: r(nx, ny)
        integer :: i, j

        do j = 1, ny
            do i = 1, nx
                if (mask(i, j)) then
                    r(i, j) = u(i, j) - update_at(KERNEL_SURFACE, nx, ny, i, j, u, K(i, j), zb(i, j), &
                                                  dfac, dx, dy, 0.0_dp, plastic_coulomb_law_t())
                else
                    r(i, j) = u(i, j) - bc(i, j)
                end if
            end do
        end do
    end subroutine godunov_residual_surface

end module plastic_eikonal
