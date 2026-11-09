! Numerical core of the perfectly-plastic reconstruction: the Godunov update of the
! static eikonal |grad u| = F, fast sweeping, the right-hand side of the w = H^2
! equation in :flat and :full mode, the Godunov update of the surface equation
! (mode "surface"), and the Godunov residuals.
!
! All routines are pure array kernels on explicit-shape (nx,ny) arrays, i = x, j = y.
! They carry no state, so several solver instances can run concurrently.
module plastic_eikonal

    use plastic_kinds, only: dp, BIG

    implicit none

    private
    public :: plastic_par_t
    public :: plastic_flotation_thickness
    public :: godunov_eikonal, godunov_surface
    public :: sweep_once, sweep_eikonal, update_at
    public :: KERNEL_EIKONAL, KERNEL_SURFACE, MARGIN_DMIN
    public :: grad_central
    public :: flotation_bc, rhs_flat, rhs_full
    public :: godunov_residual, godunov_residual_surface

    integer,  parameter :: KERNEL_EIKONAL = 1, KERNEL_SURFACE = 2

    ! Floor of the spacing factor to an ice-free neighbour, (1 - margin_frac) >= MARGIN_DMIN:
    ! keeps the update well-conditioned as the margin approaches the ice-cell centre.
    real(dp), parameter :: MARGIN_DMIN = 0.01_dp

    ! Physical constants and solver controls.
    type plastic_par_t
        real(dp) :: rho_i  = 917.0_dp       ! ice density [kg m-3]
        real(dp) :: rho_w  = 1025.0_dp      ! sea-water density [kg m-3]
        real(dp) :: g      = 9.80665_dp     ! gravity [m s-2]
        real(dp) :: H_min  = 1.0_dp         ! thickness floor in the :full RHS [m]
        character(len=8) :: mode = "full"   ! "full" | "flat" | "surface"
        integer  :: max_sweeps = 200        ! inner fast-sweeping cap
        integer  :: n_outer    = 100        ! :full outer fixed-point cap
        real(dp) :: tol        = 1.0e-6_dp  ! inner early exit, on max decrease of w [m2] (z_s [m] in "surface")
        real(dp) :: outer_tol  = 1.0e-3_dp  ! :full outer exit, on max change of H [m]
        real(dp) :: relax      = 0.5_dp     ! :full under-relaxation in (0,1]
        ! Margin position between an ice cell and an ice-free neighbour, without margin_frac:
        ! "face" (half way, margin_frac = 0.5) or "centre" (at the ice-free cell centre,
        ! margin_frac = 0, as PlasticIceSheet.jl).
        character(len=8) :: margin = "face"
        ! Warm start: fall back to a cold start when more than this fraction of the
        ! ice cells is new to the mask (see plastic_solve).
        real(dp) :: warm_frac_max = 0.1_dp
    end type plastic_par_t

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

    ! Godunov update at ice cell (i,j) from its four neighbours (clamp-to-edge), with the
    ! eikonal (kind KERNEL_EIKONAL, F = slowness) or the surface kernel (KERNEL_SURFACE,
    ! F = tau_b/(rho_i g), zbij the cell's bed). The spacing to a neighbour is d dfac(nb):
    ! dfac = 1 on ice cells, < 1 on ice-free cells whose margin lies between the two
    ! centres (sub-cell margin position). With equal spacings on both sides of an axis
    ! the upwind neighbour is the smaller one; otherwise the axis term is the Godunov
    ! max_s (u - a_s)+ / d_s over the two sides, i.e. the update is the minimum over the
    ! side combinations. Each candidate is monotone, so the minimum is.
    pure function update_at(kind, nx, ny, i, j, u, F, zbij, dfac, dx, dy) result(unew)
        integer,  intent(in) :: kind, nx, ny, i, j
        real(dp), intent(in) :: u(nx, ny), F, zbij, dfac(nx, ny), dx, dy
        real(dp) :: unew, ax(2), dxs(2), ay(2), dys(2), c
        integer  :: im, ip, jm, jp, nxc, nyc, p, q

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
            do p = 1, nxc
                if (kind == KERNEL_EIKONAL) then
                    c = godunov_eikonal(ax(p), dxs(p), ay(q), dys(q), F)
                else
                    c = godunov_surface(ax(p), dxs(p), ay(q), dys(q), zbij, F)
                end if
                unew = min(unew, c)
            end do
        end do
    end function update_at

    ! One Gauss-Seidel sweep in the ordering (xdir, ydir) = (+-1, +-1). Updates ice
    ! cells only, keeps the monotone min, and returns the largest decrease of u. Cells
    ! with equal spacing on both sides of each axis (all but those at a sub-cell margin)
    ! take the single-candidate update inline; the others go through update_at.
    function sweep_once(kind, nx, ny, u, F, zb, mask, dfac, dx, dy, xdir, ydir) result(change)
        integer,  intent(in)    :: kind, nx, ny
        real(dp), intent(inout) :: u(nx, ny)
        real(dp), intent(in)    :: F(nx, ny), zb(nx, ny), dfac(nx, ny)
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
                    if (kind == KERNEL_EIKONAL) then
                        unew = godunov_eikonal(a, dx * dfac(im, j), b, dy * dfac(i, jm), F(i, j))
                    else
                        unew = godunov_surface(a, dx * dfac(im, j), b, dy * dfac(i, jm), zb(i, j), F(i, j))
                    end if
                else
                    unew = update_at(kind, nx, ny, i, j, u, F(i, j), zb(i, j), dfac, dx, dy)
                end if
                if (unew < u(i, j)) then
                    change  = max(change, u(i, j) - unew)
                    u(i, j) = unew
                end if
            end do
        end do
    end function sweep_once

    ! Solve |grad u| = F on the ice cells (mask), with Dirichlet values bc on the
    ! ice-free cells, placed at d dfac from the neighbouring ice-cell centres. u starts
    ! from BIG on ice cells (the monotone scheme only decreases u). The four diagonal
    ! orderings alternate; the early exit needs a full cycle of four. Returns the number
    ! of sweeps and whether tol was reached. With zb present, solves the surface
    ! equation instead (u = z_s, F = tau_b/(rho_i g), see godunov_surface).
    subroutine sweep_eikonal(nx, ny, u, F, mask, bc, dfac, dx, dy, max_sweeps, tol, n_sweeps, ok, zb)
        integer,  intent(in)  :: nx, ny
        real(dp), intent(out) :: u(nx, ny)
        real(dp), intent(in)  :: F(nx, ny), bc(nx, ny), dfac(nx, ny)
        logical,  intent(in)  :: mask(nx, ny)
        real(dp), intent(in)  :: dx, dy, tol
        integer,  intent(in)  :: max_sweeps
        integer,  intent(out) :: n_sweeps
        logical,  intent(out) :: ok
        real(dp), optional, intent(in) :: zb(nx, ny)
        integer, parameter :: xd(4) = [1, -1, -1, 1]
        integer, parameter :: yd(4) = [1, 1, -1, -1]
        integer  :: i, j, s, k
        real(dp) :: change

        do j = 1, ny
            do i = 1, nx
                if (mask(i, j)) then
                    u(i, j) = BIG
                else
                    u(i, j) = bc(i, j)
                end if
            end do
        end do

        ok = .false.
        n_sweeps = 0
        do s = 1, max_sweeps
            k = mod(s - 1, 4) + 1
            if (present(zb)) then
                change = sweep_once(KERNEL_SURFACE, nx, ny, u, F, zb, mask, dfac, dx, dy, xd(k), yd(k))
            else
                change = sweep_once(KERNEL_EIKONAL, nx, ny, u, F, F, mask, dfac, dx, dy, xd(k), yd(k))
            end if
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
                                                  dfac, dx, dy)
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
                                                  dfac, dx, dy)
                else
                    r(i, j) = u(i, j) - bc(i, j)
                end if
            end do
        end do
    end subroutine godunov_residual_surface

end module plastic_eikonal
