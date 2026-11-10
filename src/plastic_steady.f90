! Steady SMB-driven reconstruction, direct formulation (docs/fortran-port.md, "direct
! steady formulation"): steady mass conservation div q = smb on the mask with the
! membrane-free DIVA flux q = -H U(tau, H) grad s/|grad s|, tau = rho_i g H |grad s|, and
! the obstacle H >= 0, solved by pseudo-transient continuation with one semismooth
! (active-set) Newton step per pseudo-time step; a cell is active (H = 0) when its own
! Newton step would empty it, H <= F/(1/dt + dF/dH).
!
! Finite volumes with staggered face fluxes: faces between ice cells with the mean
! thickness limited towards the upstream value (face_H; keeps F non-increasing in the
! neighbours' H), and a
! one-sided face flux from an ice cell to the margin (at the margin position of the
! plastic solve, margin_frac). The implied SMB of a geometry, plastic_smb_flux_divergence,
! is the divergence of the same face fluxes, so a geometry forced by its own implied SMB
! is a root of the residual.
module plastic_steady

    use plastic_kinds,    only: dp
    use plastic_eikonal,  only: plastic_flotation_thickness
    use plastic_solver,   only: plastic_t, check_shape, margin_spacing, &
                                PLASTIC_OK, PLASTIC_NOT_CONVERGED, PLASTIC_BAD_INPUT
    use plastic_flow,     only: plastic_glen_t, plastic_sliding_t, LAW_COULOMB, &
                                sliding_law_code, basal_velocity_code, deformational_velocity
    use plastic_linsolve, only: lin_reserve, lin_solve, lin_has_lis

    implicit none

    private
    public :: steady_solve, plastic_smb_flux_divergence, steady_margin_values, steady_tau

contains

    ! --- flux law ---------------------------------------------------------------------

    ! Normal flux -Hf U(tau, Hf) gn/|g| through a face, g = (gn, gt) the surface slope; bf the sliding beta of the face.
    pure function face_q(Hf, gn, gt, rg, taucap, code, slide, rheo, bf) result(q)
        real(dp),                intent(in) :: Hf, gn, gt, rg, taucap, bf
        integer,                 intent(in) :: code
        type(plastic_sliding_t), intent(in) :: slide
        type(plastic_glen_t),    intent(in) :: rheo
        real(dp) :: q, g, tau, u
        q = 0.0_dp
        g = sqrt(gn * gn + gt * gt)
        if (.not. (Hf > 0.0_dp .and. g > 0.0_dp)) return
        tau = min(rg * Hf * g, taucap)
        u = basal_velocity_code(code, slide, tau, bf) + deformational_velocity(rheo, tau, Hf)
        q = -Hf * u * gn / g
    end function face_q

    ! Face thickness between ice cells A and B: the mean (Mahaffy) where the surface drop
    ! across the face is small against the upstream thickness, limited towards the upstream
    ! value where it is not. The downstream weight w = min(1/2, H_up/(4 ds)) keeps the face
    ! flux non-increasing in the downstream thickness (dq/dH_down <= 0 for q ~ H^(p+1) g^p,
    ! p >= 1, i.e. Glen and the sliding laws), so F_P is non-increasing in its neighbours'
    ! H: an M-function, for which the active-set Newton step is consistent with H >= 0.
    pure function face_H(HA, HB, sA, sB) result(Hf)
        real(dp), intent(in) :: HA, HB, sA, sB
        real(dp) :: Hf, Hup, Hdn, ds, w
        if (sA >= sB) then
            Hup = HA; Hdn = HB
        else
            Hup = HB; Hdn = HA
        end if
        ds = abs(sA - sB)
        w = 0.5_dp
        if (4.0_dp * ds * 0.5_dp > Hup) w = Hup / (4.0_dp * ds)
        Hf = (1.0_dp - w) * Hup + w * Hdn
    end function face_H

    ! Outflow from ice cell P to the margin at distance dc (one-sided slope, segment-mean
    ! thickness with the margin thickness hmM; bP the sliding beta of P), never inflow.
    pure function margin_q(HP, hmM, sP, sM, dc, gt, rg, taucap, code, slide, rheo, bP) result(q)
        real(dp),                intent(in) :: HP, hmM, sP, sM, dc, gt, rg, taucap, bP
        integer,                 intent(in) :: code
        type(plastic_sliding_t), intent(in) :: slide
        type(plastic_glen_t),    intent(in) :: rheo
        real(dp) :: q
        q = max(face_q(0.5_dp * (HP + hmM), (sM - sP) / dc, gt, rg, taucap, code, slide, rheo, bP), 0.0_dp)
    end function margin_q

    ! Central slope of s at (i,j) along y (dir = 2) or x (dir = 1), clamp-to-edge.
    pure function cslope(nx, ny, s, i, j, d, dir) result(g)
        integer,  intent(in) :: nx, ny, i, j, dir
        real(dp), intent(in) :: s(nx, ny), d
        real(dp) :: g
        integer  :: a, b
        g = 0.0_dp
        if (dir == 1) then
            a = max(i - 1, 1); b = min(i + 1, nx)
            if (b > a) g = (s(b, j) - s(a, j)) / (real(b - a, dp) * d)
        else
            a = max(j - 1, 1); b = min(j + 1, ny)
            if (b > a) g = (s(i, b) - s(i, a)) / (real(b - a, dp) * d)
        end if
    end function cslope

    ! Residual F = div q - smb on the ice cells (0 elsewhere); s, qx, qy are outputs. bm: the sliding beta per cell
    ! (a face takes the mean of its two cells, a margin face that of its ice cell).
    subroutine steady_residual(nx, ny, H, zb, sm, hm, mask, dfac, smb, dx, dy, rg, taucap, &
                               code, slide, rheo, bm, s, qx, qy, F)
        integer,                 intent(in)  :: nx, ny, code
        real(dp),                intent(in)  :: H(nx, ny), zb(nx, ny), sm(nx, ny), hm(nx, ny), bm(nx, ny)
        real(dp),                intent(in)  :: dfac(nx, ny), smb(nx, ny), dx, dy, rg, taucap
        logical,                 intent(in)  :: mask(nx, ny)
        type(plastic_sliding_t), intent(in)  :: slide
        type(plastic_glen_t),    intent(in)  :: rheo
        real(dp),                intent(out) :: s(nx, ny), qx(0:nx, ny), qy(nx, 0:ny), F(nx, ny)
        integer  :: i, j

        !$omp parallel private(i, j)
        !$omp do schedule(static)
        do j = 1, ny
            do i = 1, nx
                if (mask(i, j)) then
                    s(i, j) = zb(i, j) + H(i, j)
                else
                    s(i, j) = sm(i, j)
                end if
            end do
        end do
        !$omp end do
        !$omp do schedule(static)
        do j = 1, ny
            qx(0, j) = 0.0_dp
            qx(nx, j) = 0.0_dp
            do i = 1, nx - 1
                if (mask(i, j) .and. mask(i + 1, j)) then
                    qx(i, j) = face_q(face_H(H(i, j), H(i + 1, j), s(i, j), s(i + 1, j)), (s(i + 1, j) - s(i, j)) / dx, &
                                      0.5_dp * (cslope(nx, ny, s, i, j, dy, 2) + cslope(nx, ny, s, i + 1, j, dy, 2)), &
                                      rg, taucap, code, slide, rheo, 0.5_dp * (bm(i, j) + bm(i + 1, j)))
                else if (mask(i, j)) then
                    qx(i, j) = margin_q(H(i, j), hm(i + 1, j), s(i, j), s(i + 1, j), dfac(i + 1, j) * dx, &
                                        cslope(nx, ny, s, i, j, dy, 2), rg, taucap, code, slide, rheo, bm(i, j))
                else if (mask(i + 1, j)) then
                    qx(i, j) = -margin_q(H(i + 1, j), hm(i, j), s(i + 1, j), s(i, j), dfac(i, j) * dx, &
                                         cslope(nx, ny, s, i + 1, j, dy, 2), rg, taucap, code, slide, rheo, bm(i + 1, j))
                else
                    qx(i, j) = 0.0_dp
                end if
            end do
        end do
        !$omp end do nowait
        !$omp do schedule(static)
        do j = 0, ny
            do i = 1, nx
                if (j == 0 .or. j == ny) then
                    qy(i, j) = 0.0_dp
                else if (mask(i, j) .and. mask(i, j + 1)) then
                    qy(i, j) = face_q(face_H(H(i, j), H(i, j + 1), s(i, j), s(i, j + 1)), (s(i, j + 1) - s(i, j)) / dy, &
                                      0.5_dp * (cslope(nx, ny, s, i, j, dx, 1) + cslope(nx, ny, s, i, j + 1, dx, 1)), &
                                      rg, taucap, code, slide, rheo, 0.5_dp * (bm(i, j) + bm(i, j + 1)))
                else if (mask(i, j)) then
                    qy(i, j) = margin_q(H(i, j), hm(i, j + 1), s(i, j), s(i, j + 1), dfac(i, j + 1) * dy, &
                                        cslope(nx, ny, s, i, j, dx, 1), rg, taucap, code, slide, rheo, bm(i, j))
                else if (mask(i, j + 1)) then
                    qy(i, j) = -margin_q(H(i, j + 1), hm(i, j), s(i, j + 1), s(i, j), dfac(i, j) * dy, &
                                         cslope(nx, ny, s, i, j + 1, dx, 1), rg, taucap, code, slide, rheo, bm(i, j + 1))
                else
                    qy(i, j) = 0.0_dp
                end if
            end do
        end do
        !$omp end do
        !$omp do schedule(static)
        do j = 1, ny
            do i = 1, nx
                if (mask(i, j)) then
                    F(i, j) = (qx(i, j) - qx(i - 1, j)) / dx + (qy(i, j) - qy(i, j - 1)) / dy - smb(i, j)
                else
                    F(i, j) = 0.0_dp
                end if
            end do
        end do
        !$omp end do
        !$omp end parallel
    end subroutine steady_residual

    ! Complementarity residual max |min(H, F)|: |F| on cells with ice, max(-F, 0) at H = 0.
    function ncp_norm(nx, ny, H, F, mask) result(r)
        integer,  intent(in) :: nx, ny
        real(dp), intent(in) :: H(nx, ny), F(nx, ny)
        logical,  intent(in) :: mask(nx, ny)
        real(dp) :: r
        integer  :: i, j
        r = 0.0_dp
        !$omp parallel do private(i, j) reduction(max:r) schedule(static)
        do j = 1, ny
            do i = 1, nx
                if (.not. mask(i, j)) cycle
                if (H(i, j) > 0.0_dp) then
                    r = max(r, abs(F(i, j)))
                else
                    r = max(r, -F(i, j))
                end if
            end do
        end do
        !$omp end parallel do
    end function ncp_norm

    ! Root-mean-square complementarity residual (steering of dt and step acceptance).
    function ncp_rms(nx, ny, H, F, mask) result(r)
        integer,  intent(in) :: nx, ny
        real(dp), intent(in) :: H(nx, ny), F(nx, ny)
        logical,  intent(in) :: mask(nx, ny)
        real(dp) :: r, e
        integer  :: i, j, m
        r = 0.0_dp
        m = 0
        do j = 1, ny
            do i = 1, nx
                if (.not. mask(i, j)) cycle
                if (H(i, j) > 0.0_dp) then
                    e = F(i, j)
                else
                    e = min(F(i, j), 0.0_dp)
                end if
                r = r + e * e
                m = m + 1
            end do
        end do
        if (m > 0) r = sqrt(r / m)
    end function ncp_rms

    ! Margin surface sm = z_b + H_flt and margin thickness hm = H_flt on all cells.
    subroutine steady_margin_values(ps, z_bed, z_sl, sm, hm)
        type(plastic_t), intent(in)  :: ps
        real(dp),        intent(in)  :: z_bed(:,:), z_sl(:,:)
        real(dp),        intent(out) :: sm(:,:), hm(:,:)
        integer :: i, j
        do j = 1, ps%ny
            do i = 1, ps%nx
                hm(i, j) = plastic_flotation_thickness(z_bed(i, j), z_sl(i, j), ps%par)
                sm(i, j) = z_bed(i, j) + hm(i, j)
            end do
        end do
    end subroutine steady_margin_values

    ! Driving stress rho_i g H |grad s| at the cell centres (central slopes).
    subroutine steady_tau(ps, mask, tau)
        type(plastic_t), intent(in)  :: ps
        logical,         intent(in)  :: mask(:,:)
        real(dp),        intent(out) :: tau(:,:)
        integer  :: i, j
        real(dp) :: rg
        rg = ps%par%rho_i * ps%par%g
        do j = 1, ps%ny
            do i = 1, ps%nx
                tau(i, j) = 0.0_dp
                if (mask(i, j)) tau(i, j) = rg * ps%H(i, j) * &
                    hypot(cslope(ps%nx, ps%ny, ps%wk%sS, i, j, ps%dx, 1), cslope(ps%nx, ps%ny, ps%wk%sS, i, j, ps%dy, 2))
            end do
        end do
    end subroutine steady_tau

    ! Outflow through the margin faces of each ice cell [m ice / yr, >= 0], split by the
    ! ice-free neighbour: marine where its bed lies below sea level (flotation thickness
    ! hm > 0), land elsewhere. Faces on the array edge carry no flux.
    subroutine margin_outflow(nx, ny, mask, hm, qx, qy, dx, dy, ql, qm)
        integer,  intent(in)  :: nx, ny
        logical,  intent(in)  :: mask(nx, ny)
        real(dp), intent(in)  :: hm(nx, ny), qx(0:nx, ny), qy(nx, 0:ny), dx, dy
        real(dp), intent(out) :: ql(nx, ny), qm(nx, ny)
        integer :: i, j
        ql = 0.0_dp
        qm = 0.0_dp
        do j = 1, ny
            do i = 1, nx - 1
                if (mask(i, j) .and. .not. mask(i + 1, j)) call add(i, j, qx(i, j) / dx, hm(i + 1, j))
                if (mask(i + 1, j) .and. .not. mask(i, j)) call add(i + 1, j, -qx(i, j) / dx, hm(i, j))
            end do
        end do
        do j = 1, ny - 1
            do i = 1, nx
                if (mask(i, j) .and. .not. mask(i, j + 1)) call add(i, j, qy(i, j) / dy, hm(i, j + 1))
                if (mask(i, j + 1) .and. .not. mask(i, j)) call add(i, j + 1, -qy(i, j) / dy, hm(i, j))
            end do
        end do
    contains
        subroutine add(ip, jp, r, h)
            integer,  intent(in) :: ip, jp
            real(dp), intent(in) :: r, h
            if (h > 0.0_dp) then
                qm(ip, jp) = qm(ip, jp) + r
            else
                ql(ip, jp) = ql(ip, jp) + r
            end if
        end subroutine add
    end subroutine margin_outflow

    ! Surface velocity at the cell centres [m/yr] of the same closure as the face fluxes:
    ! u_b(tau) + 2A/(n+1) tau^n H along -grad s, tau = rho_i g H |grad s| (central slopes of
    ! s, the margin surface on the ice-free cells), capped as the fluxes; 0 off the mask.
    subroutine surface_velocity(nx, ny, H, s, mask, dx, dy, rg, taucap, code, slide, rheo, bm, ux, uy)
        integer,                 intent(in)  :: nx, ny, code
        real(dp),                intent(in)  :: H(nx, ny), s(nx, ny), dx, dy, rg, taucap, bm(nx, ny)
        logical,                 intent(in)  :: mask(nx, ny)
        type(plastic_sliding_t), intent(in)  :: slide
        type(plastic_glen_t),    intent(in)  :: rheo
        real(dp),                intent(out) :: ux(nx, ny), uy(nx, ny)
        integer  :: i, j
        real(dp) :: gx, gy, g, tau, u
        !$omp parallel do private(i, j, gx, gy, g, tau, u) schedule(static)
        do j = 1, ny
            do i = 1, nx
                ux(i, j) = 0.0_dp
                uy(i, j) = 0.0_dp
                if (.not. (mask(i, j) .and. H(i, j) > 0.0_dp)) cycle
                gx = cslope(nx, ny, s, i, j, dx, 1)
                gy = cslope(nx, ny, s, i, j, dy, 2)
                g = hypot(gx, gy)
                if (.not. g > 0.0_dp) cycle
                tau = min(rg * H(i, j) * g, taucap)
                u = basal_velocity_code(code, slide, tau, bm(i, j)) + &
                    deformational_velocity(rheo, tau, H(i, j)) * (rheo%n + 2.0_dp) / (rheo%n + 1.0_dp)
                ux(i, j) = -u * gx / g
                uy(i, j) = -u * gy / g
            end do
        end do
        !$omp end parallel do
    end subroutine surface_velocity

    ! Implied SMB of the geometry on ps (ps%H): the divergence of the steady solver's face
    ! fluxes [m ice / yr] (closure B of ctrl%method = "steady"). The margin position is
    ! margin_frac or par%margin, as in plastic_solve. Optional: the outflow through the
    ! margin faces of each ice cell, to land (q_land) and to the ocean (q_marine) [m ice / yr,
    ! >= 0; their sum over the mask is that of smb], and the surface velocity of the same
    ! closure at the cell centres (us_x, us_y) [m/yr].
    subroutine plastic_smb_flux_divergence(ps, z_bed, z_sl, mask, rheo, slide, smb, margin_frac, stat, &
                                           q_land, q_marine, us_x, us_y)
        type(plastic_t),         intent(inout) :: ps
        real(dp),                intent(in)    :: z_bed(:,:), z_sl(:,:)
        logical,                 intent(in)    :: mask(:,:)
        type(plastic_glen_t),    intent(in)    :: rheo
        type(plastic_sliding_t), intent(in)    :: slide
        real(dp),                intent(out)   :: smb(:,:)
        real(dp), optional,      intent(in)    :: margin_frac(:,:)
        integer,  optional,      intent(out)   :: stat
        real(dp), optional,      intent(out)   :: q_land(:,:), q_marine(:,:), us_x(:,:), us_y(:,:)
        integer  :: code
        real(dp) :: taucap
        real(dp), allocatable :: ql(:,:), qm(:,:), ux(:,:), uy(:,:), bm(:,:)

        ps%stat = PLASTIC_OK
        ps%msg  = ""
        code = sliding_law_code(slide)
        if (.not. allocated(ps%H) .or. code == 0) then
            ps%stat = PLASTIC_BAD_INPUT
            ps%msg  = "plastic_smb_flux_divergence: not initialized or unknown sliding law"
        else
            call check_shape(ps, size(z_bed, 1), size(z_bed, 2), "z_bed")
            call check_shape(ps, size(z_sl, 1),  size(z_sl, 2),  "z_sl")
            call check_shape(ps, size(mask, 1),  size(mask, 2),  "mask")
            call check_shape(ps, size(smb, 1),   size(smb, 2),   "smb")
            if (present(q_land))   call check_shape(ps, size(q_land, 1),   size(q_land, 2),   "q_land")
            if (present(q_marine)) call check_shape(ps, size(q_marine, 1), size(q_marine, 2), "q_marine")
            if (present(us_x))     call check_shape(ps, size(us_x, 1),     size(us_x, 2),     "us_x")
            if (present(us_y))     call check_shape(ps, size(us_y, 1),     size(us_y, 2),     "us_y")
            if (allocated(slide%beta_map)) call check_shape(ps, size(slide%beta_map, 1), size(slide%beta_map, 2), "beta_map")
            call margin_spacing(ps, mask, margin_frac)
        end if
        if (present(stat)) stat = ps%stat
        if (ps%stat /= PLASTIC_OK) return

        taucap = huge(1.0_dp)
        if (code == LAW_COULOMB) taucap = 0.999_dp * slide%C
        bm = beta_cells(ps, slide)
        call steady_margin_values(ps, z_bed, z_sl, ps%wk%sSm, ps%wk%sHm)
        ps%wk%sFp = 0.0_dp
        call steady_residual(ps%nx, ps%ny, ps%H, z_bed, ps%wk%sSm, ps%wk%sHm, mask, ps%wk%dfac, ps%wk%sFp, &
                             ps%dx, ps%dy, ps%par%rho_i * ps%par%g, taucap, code, slide, rheo, bm, &
                             ps%wk%sS, ps%wk%qx, ps%wk%qy, ps%wk%sF)
        smb = ps%wk%sF
        if (present(q_land) .or. present(q_marine)) then
            allocate(ql(ps%nx, ps%ny), qm(ps%nx, ps%ny))
            call margin_outflow(ps%nx, ps%ny, mask, ps%wk%sHm, ps%wk%qx, ps%wk%qy, ps%dx, ps%dy, ql, qm)
            if (present(q_land))   q_land   = ql
            if (present(q_marine)) q_marine = qm
        end if
        if (present(us_x) .or. present(us_y)) then
            allocate(ux(ps%nx, ps%ny), uy(ps%nx, ps%ny))
            call surface_velocity(ps%nx, ps%ny, ps%H, ps%wk%sS, mask, ps%dx, ps%dy, ps%par%rho_i * ps%par%g, &
                                  taucap, code, slide, rheo, bm, ux, uy)
            if (present(us_x)) us_x = ux
            if (present(us_y)) us_y = uy
        end if
    end subroutine plastic_smb_flux_divergence

    ! The sliding beta per cell: slide%beta_map where allocated (shape checked by the caller), else slide%beta.
    function beta_cells(ps, slide) result(bm)
        type(plastic_t),         intent(in) :: ps
        type(plastic_sliding_t), intent(in) :: slide
        real(dp), allocatable :: bm(:,:)
        allocate(bm(ps%nx, ps%ny))
        if (allocated(slide%beta_map)) then
            bm = slide%beta_map
        else
            bm = slide%beta
        end if
    end function beta_cells

    ! Steady solve from the initial thickness in ps%H (ps%wk%dfac set by the caller).
    ! Returns the number of pseudo-time steps, the final complementarity residual [m/yr]
    ! and whether it is below res_tol; lin_fail counts linear solves that missed lin_tol.
    subroutine steady_solve(ps, z_bed, z_sl, smb, mask, rheo, slide, res_tol, max_steps, dt0, dt_max, &
                            lin_tol, lin_maxit, use_lis, n_steps, res, ok, lin_fail, verbose)
        type(plastic_t),         intent(inout) :: ps
        real(dp),                intent(in)    :: z_bed(:,:), z_sl(:,:), smb(:,:)
        logical,                 intent(in)    :: mask(:,:)
        type(plastic_glen_t),    intent(in)    :: rheo
        type(plastic_sliding_t), intent(in)    :: slide
        real(dp),                intent(in)    :: res_tol, dt0, dt_max, lin_tol
        integer,                 intent(in)    :: max_steps, lin_maxit
        logical,                 intent(in)    :: use_lis
        integer,                 intent(out)   :: n_steps, lin_fail
        real(dp),                intent(out)   :: res
        logical,                 intent(out)   :: ok
        logical,                 intent(in)    :: verbose

        integer  :: nx, ny, code, k, i, j, a, b, c, ci, cj, n, nz, slot, ii, jj, row, kk, loc(2)
        integer  :: nls
        real(dp) :: lam
        real(dp) :: rg, taucap, dt, res_new, dx, dy, d, rms, rms_new
        logical  :: lok
        real(dp), allocatable :: bm(:,:)

        nx = ps%nx; ny = ps%ny; dx = ps%dx; dy = ps%dy
        code = sliding_law_code(slide)
        bm = beta_cells(ps, slide)
        rg = ps%par%rho_i * ps%par%g
        taucap = huge(1.0_dp)
        if (code == LAW_COULOMB) taucap = 0.999_dp * slide%C
        call steady_margin_values(ps, z_bed, z_sl, ps%wk%sSm, ps%wk%sHm)
        call lin_reserve(ps%wk%lin, count(mask), 9 * count(mask))

        where (.not. mask) ps%H = 0.0_dp
        ps%H = max(ps%H, 0.0_dp)
        call resid(ps%H, ps%wk%sF)
        res = ncp_norm(nx, ny, ps%H, ps%wk%sF, mask)
        rms = ncp_rms(nx, ny, ps%H, ps%wk%sF, mask)
        dt = dt0
        n_steps = 0
        lin_fail = 0
        ok = res < res_tol

        do k = 1, max_steps
            if (ok) exit
            n_steps = k


            ! Jacobian dF_P/dH_K over the 3 x 3 block by 9 coloured residual evaluations
            ps%wk%jl = 0.0_dp
            do c = 0, 8
                ci = mod(c, 3); cj = c / 3
                !$omp parallel do private(i, j) schedule(static)
                do j = 1, ny
                    do i = 1, nx
                        ps%wk%sHt(i, j) = ps%H(i, j)
                        if (mask(i, j) .and. mod(i, 3) == ci .and. mod(j, 3) == cj) &
                            ps%wk%sHt(i, j) = ps%H(i, j) + 1.0e-7_dp * max(ps%H(i, j), 10.0_dp)
                    end do
                end do
                !$omp end parallel do
                call resid(ps%wk%sHt, ps%wk%sFp)
                !$omp parallel do private(i, j, a, b, ii, jj, d) schedule(static)
                do j = 1, ny
                    do i = 1, nx
                        if (.not. mask(i, j)) cycle
                        do b = -1, 1
                            jj = j + b
                            if (jj < 1 .or. jj > ny) cycle
                            if (mod(jj, 3) /= cj) cycle
                            do a = -1, 1
                                ii = i + a
                                if (ii < 1 .or. ii > nx) cycle
                                if (mod(ii, 3) /= ci .or. .not. mask(ii, jj)) cycle
                                d = ps%wk%sHt(ii, jj) - ps%H(ii, jj)
                                ps%wk%jl((b + 1) * 3 + a + 2, i, j) = (ps%wk%sFp(i, j) - ps%wk%sF(i, j)) / d
                            end do
                        end do
                    end do
                end do
                !$omp end parallel do
            end do

            ! active set: the cell would empty under its own Newton step,
            ! H <= G/(1/dt + dF_P/dH_P) (G = F at the start of the step); dt where the
            ! diagonal is not positive
            !$omp parallel do private(i, j, d) schedule(static)
            do j = 1, ny
                do i = 1, nx
                    d = 1.0_dp / dt + ps%wk%jl(5, i, j)
                    if (d > 0.0_dp) then
                        ps%wk%sact(i, j) = mask(i, j) .and. ps%H(i, j) * d <= ps%wk%sF(i, j)
                    else
                        ps%wk%sact(i, j) = mask(i, j) .and. ps%H(i, j) <= dt * ps%wk%sF(i, j)
                    end if
                end do
            end do
            !$omp end parallel do

            ! unknowns: inactive ice cells, column-major
            n = 0
            do j = 1, ny
                do i = 1, nx
                    ps%wk%sid(i, j) = 0
                    if (mask(i, j) .and. .not. ps%wk%sact(i, j)) then
                        n = n + 1
                        ps%wk%sid(i, j) = n
                    end if
                end do
            end do

            if (n > 0) then
                ! row pointers
                ps%wk%lin%ia(1) = 1
                do j = 1, ny
                    do i = 1, nx
                        row = ps%wk%sid(i, j)
                        if (row == 0) cycle
                        nz = 0
                        do b = -1, 1
                            do a = -1, 1
                                ii = i + a; jj = j + b
                                if (ii < 1 .or. ii > nx .or. jj < 1 .or. jj > ny) cycle
                                if (ps%wk%sid(ii, jj) > 0) nz = nz + 1
                            end do
                        end do
                        ps%wk%lin%ia(row + 1) = ps%wk%lin%ia(row) + nz
                    end do
                end do
                ps%wk%lin%n = n
                ps%wk%lin%nnz = ps%wk%lin%ia(n + 1) - 1
                ! entries and right-hand side (rows in parallel)
                !$omp parallel do private(i, j, a, b, ii, jj, row, kk, slot) schedule(static)
                do j = 1, ny
                    do i = 1, nx
                        row = ps%wk%sid(i, j)
                        if (row == 0) cycle
                        kk = ps%wk%lin%ia(row)
                        ps%wk%lin%b(row) = -ps%wk%sF(i, j)
                        do b = -1, 1
                            do a = -1, 1
                                ii = i + a; jj = j + b
                                if (ii < 1 .or. ii > nx .or. jj < 1 .or. jj > ny) cycle
                                if (.not. mask(ii, jj)) cycle
                                slot = (b + 1) * 3 + a + 2
                                if (ps%wk%sid(ii, jj) > 0) then
                                    ps%wk%lin%ja(kk) = ps%wk%sid(ii, jj)
                                    ps%wk%lin%a(kk)  = ps%wk%jl(slot, i, j)
                                    if (a == 0 .and. b == 0) ps%wk%lin%a(kk) = ps%wk%lin%a(kk) + 1.0_dp / dt
                                    kk = kk + 1
                                else
                                    ! active neighbour: d_K = -H_K on the right-hand side
                                    ps%wk%lin%b(row) = ps%wk%lin%b(row) + ps%wk%jl(slot, i, j) * ps%H(ii, jj)
                                end if
                            end do
                        end do
                        ps%wk%lin%x(row) = 0.0_dp
                    end do
                end do
                !$omp end parallel do
                call lin_solve(ps%wk%lin, lin_tol, lin_maxit, use_lis .and. lin_has_lis(), lok)
                if (.not. lok) lin_fail = lin_fail + 1
            end if
            !$omp parallel do private(i, j) schedule(static)
            do j = 1, ny
                do i = 1, nx
                    if (.not. mask(i, j)) then
                        ps%wk%sHk(i, j) = 0.0_dp
                    else if (ps%wk%sact(i, j) .or. ps%wk%sid(i, j) == 0) then
                        ps%wk%sHk(i, j) = -ps%H(i, j)
                    else
                        ps%wk%sHk(i, j) = ps%wk%lin%x(ps%wk%sid(i, j))
                    end if
                end do
            end do
            !$omp end parallel do

            ! projected backtracking on the rms residual: H(lam) = max(H + lam d, 0)
            lam = 1.0_dp
            do nls = 1, 6
                !$omp parallel do private(i, j) schedule(static)
                do j = 1, ny
                    do i = 1, nx
                        ps%wk%sHt(i, j) = max(ps%H(i, j) + lam * ps%wk%sHk(i, j), 0.0_dp)
                    end do
                end do
                !$omp end parallel do
                call resid(ps%wk%sHt, ps%wk%sFp)
                rms_new = ncp_rms(nx, ny, ps%wk%sHt, ps%wk%sFp, mask)
                if (rms_new <= (1.0_dp - 1.0e-4_dp * lam) * rms) exit
                lam = 0.5_dp * lam
            end do
            res_new = ncp_norm(nx, ny, ps%wk%sHt, ps%wk%sFp, mask)

            if (verbose) write(*, '(a,i5,a,es9.2,a,es9.2,a,es9.2,a,i7,a,i7,a,i4,a,es9.2)') "  steady ", k, &
                ": dt ", dt, " res ", res, " -> ", res_new, ", unknowns ", n, ", active ", &
                count(ps%wk%sact), ", lin it ", ps%wk%lin%iters, " relres ", ps%wk%lin%relres
            if (verbose) then
                loc = maxloc(merge(abs(ps%wk%sFp), 0.0_dp, mask .and. ps%wk%sHt > 0.0_dp))
                write(*, '(a,2i5,a,f10.2,a,es10.3)') "         max |F| at ", loc, " H ", ps%wk%sHt(loc(1), loc(2)), &
                    " F ", ps%wk%sFp(loc(1), loc(2))
            end if
            if (rms_new > 10.0_dp * rms .or. rms_new /= rms_new) then
                dt = 0.5_dp * dt                     ! reject: shorter pseudo-time step
                if (dt < 1.0e-12_dp * dt0) exit
                cycle
            end if
            ps%H  = ps%wk%sHt
            ps%wk%sF = ps%wk%sFp
            ! SER on the rms residual, growing at least 1.5x after a step that did not
            ! increase it
            if (rms_new <= rms) then
                dt = min(dt * max(rms / max(rms_new, tiny(1.0_dp)), 1.5_dp), dt_max)
            else
                dt = dt * rms / rms_new
            end if
            dt = max(dt, dt0 * 1.0e-6_dp)
            res = res_new
            rms = rms_new
            ok = res < res_tol
        end do
        ! leave s, qx, qy of the final state in the workspace
        call resid(ps%H, ps%wk%sF)

    contains

        subroutine resid(Hx, Fx)
            real(dp), intent(in)  :: Hx(:,:)
            real(dp), intent(out) :: Fx(:,:)
            call steady_residual(nx, ny, Hx, z_bed, ps%wk%sSm, ps%wk%sHm, mask, ps%wk%dfac, smb, dx, dy, &
                                 rg, taucap, code, slide, rheo, bm, ps%wk%sS, ps%wk%qx, ps%wk%qy, Fx)
        end subroutine resid

    end subroutine steady_solve

end module plastic_steady
