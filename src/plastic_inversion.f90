module plastic_inversion
    ! Implied mass balance by joint inversion of the speed and the mass balance on a given plastic geometry (kryos
    ! docs/dev/plastic-smb-inversion-design.md). On the geometry of ps (ps%H on mask, the bed z_bed, sea level z_sl):
    !
    !   unknowns   the depth-averaged speed U of each ice cell along the plastic flow direction -grad s/|grad s|
    !              (fixed by the geometry), x = ln U;
    !   mass       a = dH/dt + div(H U n) with the face geometry of the steady fluxes (face thickness face_H, face
    !              slope; a face's speed that of its upstream cell, as a balance flux: the mean of the two cells
    !              decouples odd and even cells along a flowline; a margin face's outflow at its ice cell's speed);
    !              the outflow to land is ablated in its margin cell, so it is left out of a (as smb - q_land of
    !              plastic_smb_flux_divergence), the outflow to the ocean is in a (and in q_marine);
    !   objective  J = 1/2 sum ((x - x_law)/sigma_u)^2 + 1/2 sum_faces (L/d)^2 ((e_A - e_B)/sigma_u)^2
    !                + 1/2 sum ((a - a_prior)/sigma_a)^2 + 1/2 sum_cold (min(a, 0)/sigma_cold)^2,
    !              e = x - x_law, U_law the closure's speed at the cell (u_b + u_def at tau = min(rho_i g H |grad s|,
    !              the stress cap)), sigma_a = max(sigma_smb_rel |a_prior|, sigma_smb_min), L = ctrl%smooth, d the
    !              face's spacing; cold: the ice cells without a margin face where a_prior >= 0;
    !   bound      x >= ln U_def (the deformation speed at that stress: sliding >= 0).
    !
    ! Gauss-Newton on x: each step the sparse SPD system (13-point) of the normal equations, scaled symmetrically by its
    ! diagonal (the terms span many orders of magnitude) and solved inexactly with the library's linear solver; the bound
    ! by projection (an active set: variables at the bound with a positive gradient are held), backtracking on J. Limits: sigma_u -> 0 gives the flow law's speed; sigma_a -> 0 the balance flux of the prior.
    ! Units: m, yr, Pa; a, a_prior, dH/dt in m ice/yr.

    use plastic_kinds,    only: dp
    use plastic_solver,   only: plastic_t, check_shape, margin_spacing, PLASTIC_OK, PLASTIC_NOT_CONVERGED, &
                                PLASTIC_BAD_INPUT
    use plastic_flow,     only: plastic_glen_t, plastic_sliding_t, sliding_law_code, basal_velocity_code, &
                                deformational_velocity
    use plastic_steady,   only: steady_margin_values, face_H, cslope, cap_cells, fric_cells
    use plastic_linsolve, only: lin_reserve, lin_solve, lin_has_lis

    implicit none

    private
    public :: plastic_inv_ctrl_t, plastic_inv_info_t, plastic_smb_inversion

    type plastic_inv_ctrl_t
        real(dp) :: sigma_u = 0.7_dp              ! sd of ln U about the flow law (0.7: a factor 2)
        real(dp) :: smooth = 0.0_dp               ! [m] the length L of the smoothness of ln U - ln U_law (0: none)
        real(dp) :: sigma_smb_rel = 0.5_dp        ! sd of a relative to |a_prior|
        real(dp) :: sigma_smb_min = 0.1_dp        ! [m ice/yr] its floor
        logical  :: cold = .true.                 ! penalise a < 0 where a_prior >= 0 (ice cells without a margin face)
        real(dp) :: sigma_cold = 0.01_dp          ! [m ice/yr] the scale of that penalty
        real(dp) :: u_min = 1.0e-3_dp             ! [m/yr] floor of U_law and of the bound
        integer  :: n_iter = 100                  ! Gauss-Newton steps
        real(dp) :: tol = 1.0e-6_dp               ! stop at a relative decrease of J below tol with a full step
        character(len=8) :: lin = "auto"         ! "auto" | "lis" | "builtin"
        real(dp) :: lin_tol = 1.0e-6_dp           ! relative residual of the (diagonally scaled) step system
        integer  :: lin_maxit = 2000
        logical  :: verbose = .false.
    end type plastic_inv_ctrl_t

    type plastic_inv_info_t
        integer  :: iters = 0                     ! Gauss-Newton steps taken
        real(dp) :: J = 0.0_dp                    ! the objective at the solution, its four terms
        real(dp) :: J_law = 0.0_dp, J_smooth = 0.0_dp, J_smb = 0.0_dp, J_cold = 0.0_dp
        real(dp) :: J_start = 0.0_dp              ! at U = U_law (the forward closure)
        integer  :: n_bound = 0                   ! cells at the bound (no sliding)
        integer  :: n_cold_neg = 0                ! cold cells with a < 0
    end type plastic_inv_info_t

    ! the 13-point stencil of the normal equations, ordered by (dj, di): increasing cell ids
    integer, parameter :: NS = 13
    integer, parameter :: OFF_I(NS) = [0, -1, 0, 1, -2, -1, 0, 1, 2, -1, 0, 1, 0]
    integer, parameter :: OFF_J(NS) = [-2, -1, -1, -1, 0, 0, 0, 0, 0, 1, 1, 1, 2]
    ! the 5 cells a mass balance row depends on: self, west, east, south, north
    integer, parameter :: D_I(0:4) = [0, -1, 1, 0, 0], D_J(0:4) = [0, 0, 0, -1, 1]

contains

    subroutine plastic_smb_inversion(ps, z_bed, z_sl, mask, rheo, slide, smb_prior, dHdt, smb, speed, ctrl, &
                                     margin_frac, tau_cap, q_land, q_marine, u_law, u_def, info, stat)
        type(plastic_t),          intent(inout) :: ps
        real(dp),                 intent(in)    :: z_bed(:,:), z_sl(:,:), smb_prior(:,:), dHdt(:,:)
        logical,                  intent(in)    :: mask(:,:)
        type(plastic_glen_t),     intent(in)    :: rheo
        type(plastic_sliding_t),  intent(in)    :: slide
        real(dp),                 intent(out)   :: smb(:,:), speed(:,:)
        type(plastic_inv_ctrl_t), intent(in),  optional :: ctrl
        real(dp),                 intent(in),  optional :: margin_frac(:,:), tau_cap(:,:)
        real(dp),                 intent(out), optional :: q_land(:,:), q_marine(:,:), u_law(:,:), u_def(:,:)
        type(plastic_inv_info_t), intent(out), optional :: info
        integer,                  intent(out), optional :: stat

        type(plastic_inv_ctrl_t) :: c
        type(plastic_inv_info_t) :: inf
        integer  :: nx, ny, code, n, i, j, k, p, it, nh, kk, q, s, t, ii, jj
        integer, allocatable :: id(:,:), ci(:), cj(:), sid(:,:)
        real(dp), allocatable :: sm(:,:), hm(:,:), sf(:,:), cm(:,:), bm(:,:)
        real(dp), allocatable :: D(:,:), cmar(:), clan(:), xlaw(:), xlo(:), x(:), xn(:), dx_(:), g(:), w(:), r(:)
        real(dp), allocatable :: ap(:), sa(:), dh(:), u(:), a(:), h(:,:), udf(:), sc(:)
        logical,  allocatable :: cold(:), act(:)
        real(dp) :: rg, tau, ub, ud, gx, gy, sux, suy, Jc, Jn, alpha, lsum(4), lsn(4)
        logical  :: ok, use_lis, conv

        if (present(ctrl)) c = ctrl
        ps%stat = PLASTIC_OK
        ps%msg  = ""
        nx = ps%nx;  ny = ps%ny
        code = sliding_law_code(slide)
        if (.not. allocated(ps%H) .or. code == 0) then
            ps%stat = PLASTIC_BAD_INPUT
            ps%msg  = "plastic_smb_inversion: not initialized or unknown sliding law"
        else if (.not. (c%sigma_u > 0.0_dp .and. c%sigma_smb_min > 0.0_dp .and. c%sigma_smb_rel >= 0.0_dp .and. &
                        c%sigma_cold > 0.0_dp .and. c%smooth >= 0.0_dp .and. c%u_min > 0.0_dp)) then
            ps%stat = PLASTIC_BAD_INPUT
            ps%msg  = "plastic_smb_inversion: sigma_u, sigma_smb_min, sigma_cold, u_min > 0; sigma_smb_rel, smooth >= 0"
        else
            call check_shape(ps, size(z_bed, 1), size(z_bed, 2), "z_bed")
            call check_shape(ps, size(z_sl, 1), size(z_sl, 2), "z_sl")
            call check_shape(ps, size(mask, 1), size(mask, 2), "mask")
            call check_shape(ps, size(smb_prior, 1), size(smb_prior, 2), "smb_prior")
            call check_shape(ps, size(dHdt, 1), size(dHdt, 2), "dHdt")
            call check_shape(ps, size(smb, 1), size(smb, 2), "smb")
            call check_shape(ps, size(speed, 1), size(speed, 2), "speed")
            if (present(tau_cap))  call check_shape(ps, size(tau_cap, 1), size(tau_cap, 2), "tau_cap")
            if (present(q_land))   call check_shape(ps, size(q_land, 1), size(q_land, 2), "q_land")
            if (present(q_marine)) call check_shape(ps, size(q_marine, 1), size(q_marine, 2), "q_marine")
            if (present(u_law))    call check_shape(ps, size(u_law, 1), size(u_law, 2), "u_law")
            if (present(u_def))    call check_shape(ps, size(u_def, 1), size(u_def, 2), "u_def")
            if (allocated(slide%fric_map)) call check_shape(ps, size(slide%fric_map, 1), size(slide%fric_map, 2), &
                                                            "fric_map")
            call margin_spacing(ps, mask, margin_frac)
        end if
        if (present(stat)) stat = ps%stat
        if (ps%stat /= PLASTIC_OK) return

        ! --- geometry and the flow law at the cells ----------------------------------------------------------------
        rg = ps%par%rho_i*ps%par%g
        allocate(sm(nx, ny), hm(nx, ny), sf(nx, ny), id(nx, ny))
        call steady_margin_values(ps, z_bed, z_sl, sm, hm)
        where (mask)
            sf = z_bed + ps%H
        elsewhere
            sf = sm
        end where
        cm = cap_cells(ps, code, slide, tau_cap)
        bm = fric_cells(ps, code, slide)
        n = count(mask)
        allocate(ci(n), cj(n), xlaw(n), xlo(n), ap(n), sa(n), dh(n), cold(n), D(0:4, n), cmar(n), clan(n), udf(n))
        id = 0;  k = 0
        do j = 1, ny
            do i = 1, nx
                if (.not. mask(i, j)) cycle
                k = k + 1;  id(i, j) = k;  ci(k) = i;  cj(k) = j
            end do
        end do
        do k = 1, n
            i = ci(k);  j = cj(k)
            tau = min(rg*ps%H(i, j)*hypot(cslope(nx, ny, sf, i, j, ps%dx, 1), cslope(nx, ny, sf, i, j, ps%dy, 2)), cm(i, j))
            ub = basal_velocity_code(code, slide, tau, bm(i, j))
            if (.not. ub >= 0.0_dp) ub = 0.0_dp
            ud = deformational_velocity(rheo, tau, ps%H(i, j))
            udf(k) = ud
            xlaw(k) = log(max(ub + ud, c%u_min))
            xlo(k)  = min(log(max(ud, c%u_min)), xlaw(k))
            ap(k) = smb_prior(i, j);  dh(k) = dHdt(i, j)
            sa(k) = max(c%sigma_smb_rel*abs(ap(k)), c%sigma_smb_min)
            cold(k) = c%cold .and. ap(k) >= 0.0_dp .and. interior(i, j)
        end do

        ! --- the mass-balance operator a = dh + D U (5 cells per row) ---------------------------------------------
        D = 0.0_dp;  cmar = 0.0_dp;  clan = 0.0_dp
        do j = 1, ny
            do i = 1, nx - 1
                if (mask(i, j) .and. mask(i + 1, j)) then
                    gx = (sf(i + 1, j) - sf(i, j))/ps%dx
                    gy = 0.5_dp*(cslope(nx, ny, sf, i, j, ps%dy, 2) + cslope(nx, ny, sf, i + 1, j, ps%dy, 2))
                    call inner(id(i, j), id(i + 1, j), 2, 1, coef(face_H(ps%H(i, j), ps%H(i + 1, j), sf(i, j), &
                               sf(i + 1, j)), gx, gy)/ps%dx)
                else if (mask(i, j)) then
                    call outer(id(i, j), i + 1, j, coef(0.5_dp*(ps%H(i, j) + hm(i + 1, j)), &
                               (sf(i + 1, j) - sf(i, j))/(ps%wk%dfac(i + 1, j)*ps%dx), cslope(nx, ny, sf, i, j, ps%dy, 2)), ps%dx)
                else if (mask(i + 1, j)) then
                    call outer(id(i + 1, j), i, j, coef(0.5_dp*(ps%H(i + 1, j) + hm(i, j)), &
                               (sf(i, j) - sf(i + 1, j))/(ps%wk%dfac(i, j)*ps%dx), cslope(nx, ny, sf, i + 1, j, ps%dy, 2)), ps%dx)
                end if
            end do
        end do
        do j = 1, ny - 1
            do i = 1, nx
                if (mask(i, j) .and. mask(i, j + 1)) then
                    gy = (sf(i, j + 1) - sf(i, j))/ps%dy
                    gx = 0.5_dp*(cslope(nx, ny, sf, i, j, ps%dx, 1) + cslope(nx, ny, sf, i, j + 1, ps%dx, 1))
                    call inner(id(i, j), id(i, j + 1), 4, 3, coef(face_H(ps%H(i, j), ps%H(i, j + 1), sf(i, j), &
                               sf(i, j + 1)), gy, gx)/ps%dy)
                else if (mask(i, j)) then
                    call outer(id(i, j), i, j + 1, coef(0.5_dp*(ps%H(i, j) + hm(i, j + 1)), &
                               (sf(i, j + 1) - sf(i, j))/(ps%wk%dfac(i, j + 1)*ps%dy), cslope(nx, ny, sf, i, j, ps%dx, 1)), ps%dy)
                else if (mask(i, j + 1)) then
                    call outer(id(i, j + 1), i, j, coef(0.5_dp*(ps%H(i, j + 1) + hm(i, j)), &
                               (sf(i, j) - sf(i, j + 1))/(ps%wk%dfac(i, j)*ps%dy), cslope(nx, ny, sf, i, j + 1, ps%dx, 1)), ps%dy)
                end if
            end do
        end do
        ! the stencil cells of each row (0 where not an ice cell)
        allocate(sid(0:4, n))
        do k = 1, n
            do s = 0, 4
                ii = ci(k) + D_I(s);  jj = cj(k) + D_J(s)
                sid(s, k) = 0
                if (ii >= 1 .and. ii <= nx .and. jj >= 1 .and. jj <= ny) sid(s, k) = id(ii, jj)
            end do
        end do
        sux = 0.0_dp;  suy = 0.0_dp
        if (c%smooth > 0.0_dp) then
            sux = (c%smooth/ps%dx)**2/c%sigma_u**2
            suy = (c%smooth/ps%dy)**2/c%sigma_u**2
        end if

        ! --- Gauss-Newton -------------------------------------------------------------------------------------------
        use_lis = trim(c%lin) == "lis" .or. (trim(c%lin) == "auto" .and. lin_has_lis())
        allocate(x(n), xn(n), dx_(n), g(n), w(n), r(n), u(n), a(n), act(n), h(NS, n), sc(n))
        x = xlaw
        call objective(x, Jc, lsum)
        inf%J_start = Jc
        conv = .false.
        call lin_reserve(ps%wk%lin, n, NS*n)
        do it = 1, c%n_iter
            inf%iters = it
            ! gradient and Hessian at x
            u = exp(x)
            call mass(u, a)
            do k = 1, n
                w(k) = 1.0_dp/sa(k)**2
                r(k) = (a(k) - ap(k))*w(k)
                if (cold(k) .and. a(k) < 0.0_dp) then
                    w(k) = w(k) + 1.0_dp/c%sigma_cold**2
                    r(k) = r(k) + a(k)/c%sigma_cold**2
                end if
            end do
            g = (x - xlaw)/c%sigma_u**2
            h = 0.0_dp
            h(7, :) = 1.0_dp/c%sigma_u**2
            call smooth_terms(x - xlaw, g, h)
            do k = 1, n
                do s = 0, 4
                    p = sid(s, k)
                    if (p == 0 .or. D(s, k) == 0.0_dp) cycle
                    g(p) = g(p) + D(s, k)*u(p)*r(k)
                    do t = 0, 4
                        q = sid(t, k)
                        if (q == 0 .or. D(t, k) == 0.0_dp) cycle
                        kk = slot(D_I(t) - D_I(s), D_J(t) - D_J(s))
                        h(kk, p) = h(kk, p) + D(s, k)*u(p)*w(k)*D(t, k)*u(q)
                    end do
                end do
            end do
            act = x <= xlo + 1.0e-12_dp .and. g > 0.0_dp
            ! the system H dx = -g (held variables: identity rows and columns), scaled: (S H S) y = -S g, dx = S y,
            ! S = diag(H)^(-1/2)
            sc = 1.0_dp/sqrt(h(7, :))
            nh = 0
            ps%wk%lin%ia(1) = 1
            do k = 1, n
                do kk = 1, NS
                    ii = ci(k) + OFF_I(kk);  jj = cj(k) + OFF_J(kk)
                    if (ii < 1 .or. ii > nx .or. jj < 1 .or. jj > ny) cycle
                    q = id(ii, jj)
                    if (q == 0) cycle
                    nh = nh + 1
                    ps%wk%lin%ja(nh) = q
                    if (act(k) .or. act(q)) then
                        ps%wk%lin%a(nh) = merge(1.0_dp, 0.0_dp, q == k)
                    else
                        ps%wk%lin%a(nh) = sc(k)*h(kk, k)*sc(q)
                    end if
                end do
                ps%wk%lin%ia(k + 1) = nh + 1
                ps%wk%lin%b(k) = merge(0.0_dp, -sc(k)*g(k), act(k))
                ps%wk%lin%x(k) = 0.0_dp
            end do
            ps%wk%lin%n = n
            ps%wk%lin%nnz = nh
            call lin_solve(ps%wk%lin, c%lin_tol, c%lin_maxit, use_lis, ok)
            dx_ = sc*ps%wk%lin%x(1:n)
            where (act) dx_ = 0.0_dp
            ! backtracking on J
            alpha = 1.0_dp
            do kk = 1, 30
                xn = max(x + alpha*dx_, xlo)
                call objective(xn, Jn, lsn)
                if (Jn < Jc) exit
                alpha = 0.5_dp*alpha
            end do
            if (c%verbose) write(*, "(a,i4,a,es12.5,a,es10.3,a,i0,a,l1,a,i0,a,es9.2)") "plastic_smb_inversion: it ", it, &
                " J ", Jn, " step ", alpha, " held ", count(act), " lin ok ", ok, " it ", ps%wk%lin%iters, " relres ", &
                ps%wk%lin%relres
            if (.not. Jn < Jc) then
                conv = .true.
                exit
            end if
            conv = alpha == 1.0_dp .and. (Jc - Jn) < c%tol*Jc      ! a full step that no longer gains
            x = xn;  Jc = Jn;  lsum = lsn
            if (conv) exit
        end do

        ! --- outputs -------------------------------------------------------------------------------------------------
        u = exp(x)
        call mass(u, a)
        smb = 0.0_dp;  speed = 0.0_dp
        if (present(q_land))   q_land = 0.0_dp
        if (present(q_marine)) q_marine = 0.0_dp
        if (present(u_law))    u_law = 0.0_dp
        if (present(u_def))    u_def = 0.0_dp
        do k = 1, n
            i = ci(k);  j = cj(k)
            smb(i, j) = a(k);  speed(i, j) = u(k)
            if (present(q_land))   q_land(i, j) = clan(k)*u(k)
            if (present(q_marine)) q_marine(i, j) = cmar(k)*u(k)
            if (present(u_law))    u_law(i, j) = exp(xlaw(k))
            if (present(u_def))    u_def(i, j) = udf(k)
        end do
        inf%J = Jc
        inf%J_law = lsum(1);  inf%J_smooth = lsum(2);  inf%J_smb = lsum(3);  inf%J_cold = lsum(4)
        inf%n_bound = count(x <= xlo + 1.0e-12_dp)
        inf%n_cold_neg = count(cold .and. a < 0.0_dp)
        if (present(info)) info = inf
        if (.not. conv) then
            ps%stat = PLASTIC_NOT_CONVERGED
            ps%msg  = "plastic_smb_inversion: no convergence in n_iter steps"
        end if
        if (present(stat)) stat = ps%stat

    contains

        logical function interior(ip, jp)
            ! an ice cell with four ice neighbours
            integer, intent(in) :: ip, jp
            interior = ip > 1 .and. ip < nx .and. jp > 1 .and. jp < ny
            if (interior) interior = mask(ip - 1, jp) .and. mask(ip + 1, jp) .and. mask(ip, jp - 1) .and. mask(ip, jp + 1)
        end function interior

        pure function coef(Hf, gn, gt) result(cf)
            ! flux per unit speed through a face (along +n): -Hf gn / |g|
            real(dp), intent(in) :: Hf, gn, gt
            real(dp) :: cf, gg
            cf = 0.0_dp
            gg = sqrt(gn*gn + gt*gt)
            if (Hf > 0.0_dp .and. gg > 0.0_dp) cf = -Hf*gn/gg
        end function coef

        subroutine inner(kp, kq, sp, sq, cf)
            ! an inner face from cell kp to kq (kq the neighbour in slot sp of kp, kp in slot sq of kq); cf the flux
            ! per unit face speed divided by the spacing (> 0: from kp to kq); the face speed that of the upstream cell
            integer,  intent(in) :: kp, kq, sp, sq
            real(dp), intent(in) :: cf
            if (cf >= 0.0_dp) then
                D(0, kp) = D(0, kp) + cf;  D(sq, kq) = D(sq, kq) - cf
            else
                D(sp, kp) = D(sp, kp) + cf;  D(0, kq) = D(0, kq) - cf
            end if
        end subroutine inner

        subroutine outer(kp, im, jm, cf, dd)
            ! a margin face of ice cell kp towards the ice-free cell (im, jm): outflow max(cf, 0) U_p / dd, marine
            ! where the neighbour's flotation thickness is positive (in a), else to land (ablated: not in a)
            integer,  intent(in) :: kp, im, jm
            real(dp), intent(in) :: cf, dd
            real(dp) :: co
            co = max(cf, 0.0_dp)/dd
            if (hm(im, jm) > 0.0_dp) then
                D(0, kp) = D(0, kp) + co;  cmar(kp) = cmar(kp) + co
            else
                clan(kp) = clan(kp) + co
            end if
        end subroutine outer

        integer function slot(di, dj)
            ! the position of offset (di, dj) in the 13-point stencil
            integer, intent(in) :: di, dj
            integer :: m
            slot = 0
            do m = 1, NS
                if (OFF_I(m) == di .and. OFF_J(m) == dj) then
                    slot = m
                    return
                end if
            end do
        end function slot

        subroutine mass(uu, aa)
            real(dp), intent(in)  :: uu(:)
            real(dp), intent(out) :: aa(:)
            integer :: m, z, pp
            do m = 1, n
                aa(m) = dh(m)
                do z = 0, 4
                    pp = sid(z, m)
                    if (pp > 0) aa(m) = aa(m) + D(z, m)*uu(pp)
                end do
            end do
        end subroutine mass

        subroutine smooth_terms(e, gg, hh)
            ! the smoothness term: gradient and Hessian (5-point)
            real(dp), intent(in)    :: e(:)
            real(dp), intent(inout) :: gg(:), hh(:,:)
            integer :: m, mq
            if (c%smooth <= 0.0_dp) return
            do m = 1, n
                if (ci(m) < nx) then
                    mq = id(ci(m) + 1, cj(m))
                    if (mq > 0) call smooth_pair(e, gg, hh, m, mq, 8, 6, sux)
                end if
                if (cj(m) < ny) then
                    mq = id(ci(m), cj(m) + 1)
                    if (mq > 0) call smooth_pair(e, gg, hh, m, mq, 11, 3, suy)
                end if
            end do
        end subroutine smooth_terms

        subroutine smooth_pair(e, gg, hh, ma, mb, sab, sba, wt)
            ! one face of the smoothness term: cells ma, mb, mb in slot sab of ma, ma in slot sba of mb
            real(dp), intent(in)    :: e(:), wt
            real(dp), intent(inout) :: gg(:), hh(:,:)
            integer,  intent(in)    :: ma, mb, sab, sba
            gg(ma) = gg(ma) + wt*(e(ma) - e(mb));  gg(mb) = gg(mb) - wt*(e(ma) - e(mb))
            hh(7, ma) = hh(7, ma) + wt;  hh(7, mb) = hh(7, mb) + wt
            hh(sab, ma) = hh(sab, ma) - wt;  hh(sba, mb) = hh(sba, mb) - wt
        end subroutine smooth_pair

        subroutine objective(xx, jj_, parts)
            real(dp), intent(in)  :: xx(:)
            real(dp), intent(out) :: jj_, parts(4)
            real(dp), allocatable :: uu(:), aa(:), e(:)
            integer :: m, mq
            uu = exp(xx)
            allocate(aa(n))
            call mass(uu, aa)
            e = xx - xlaw
            parts = 0.0_dp
            parts(1) = 0.5_dp*sum(e**2)/c%sigma_u**2
            if (c%smooth > 0.0_dp) then
                do m = 1, n
                    if (ci(m) < nx) then
                        mq = id(ci(m) + 1, cj(m))
                        if (mq > 0) parts(2) = parts(2) + 0.5_dp*sux*(e(m) - e(mq))**2
                    end if
                    if (cj(m) < ny) then
                        mq = id(ci(m), cj(m) + 1)
                        if (mq > 0) parts(2) = parts(2) + 0.5_dp*suy*(e(m) - e(mq))**2
                    end if
                end do
            end if
            parts(3) = 0.5_dp*sum(((aa - ap)/sa)**2)
            parts(4) = 0.5_dp*sum(merge(min(aa, 0.0_dp)/c%sigma_cold, 0.0_dp, cold)**2)
            jj_ = sum(parts)
        end subroutine objective

    end subroutine plastic_smb_inversion

end module plastic_inversion
