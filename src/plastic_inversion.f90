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
    !   objective  J = 1/2 sum ((x - x_law)/sigma_u)^2 + 1/2 sum_faces (L_u/d)^2 ((e_A - e_B)/sigma_u)^2
    !                + 1/2 sum eps^2 + 1/2 sum_faces (L_a/d)^2 (eps_A - eps_B)^2 + 1/2 sum_cold (min(a, 0)/sigma_cold)^2,
    !              e = x - x_law, eps = (a - a_prior)/sigma_a, U_law the closure's speed at the cell (u_b + u_def at
    !              tau = min(rho_i g H |grad s|, the stress cap)), sigma_a = max(sigma_smb_rel |a_prior|, sigma_smb_min),
    !              L_u = ctrl%smooth (the speed correction), L_a = ctrl%smooth_smb (the error of a: spatially correlated),
    !              d the face's spacing; cold: the ice cells without a margin face where a_prior >= 0;
    !   soft bound + 1/2 sum (min(x - ln U_def, 0)/sigma_def)^2: the speed at least the deformation speed at that
    !              stress (sliding >= 0) within sigma_def, the uncertainty of the rate factor.
    !
    ! Gauss-Newton on x: each step the SPD system H dx = -g of the normal equations, H = Lambda + S_u + B' W B
    ! (Lambda the flow-law and bound terms, S_u the speed smoothness, B = D diag(U) the linearised mass balance, W =
    ! Sigma^-1 (C + S_a) Sigma^-1 its weights, C the diagonal of the misfit and the cold penalty, S_a the smoothness of
    ! eps), applied matrix-free. The mass-balance part dominates by orders of magnitude, so H is solved by conjugate
    ! gradients preconditioned with the inverse of B' W B: with upwind face speeds D is lower triangular when the cells
    ! are ordered by descending surface, so B^-1 W^-1 B^-T costs two sweeps (the land outflow on the diagonal of the
    ! preconditioner's D, so every row has one) and W^-1 an inner solve of the screened Laplacian C + S_a (Jacobi-
    ! preconditioned CG; with L_a = 0 a division). Backtracking on J. Limits: sigma_u -> 0 gives the flow law's speed; sigma_a -> 0 the balance flux of the prior.
    ! Units: m, yr, Pa; a, a_prior, dH/dt in m ice/yr.

    use plastic_kinds,    only: dp
    use plastic_solver,   only: plastic_t, check_shape, margin_spacing, PLASTIC_OK, PLASTIC_NOT_CONVERGED, &
                                PLASTIC_BAD_INPUT
    use plastic_flow,     only: plastic_glen_t, plastic_sliding_t, sliding_law_code, basal_velocity_code, &
                                deformational_velocity
    use plastic_steady,   only: steady_margin_values, face_H, cslope, cap_cells, fric_cells

    implicit none

    private
    public :: plastic_inv_ctrl_t, plastic_inv_info_t, plastic_smb_inversion

    type plastic_inv_ctrl_t
        real(dp) :: sigma_u = 0.7_dp              ! sd of ln U about the flow law (0.7: a factor 2)
        real(dp) :: smooth = 0.0_dp               ! [m] the length L_u of the smoothness of ln U - ln U_law (0: none)
        real(dp) :: smooth_smb = 0.0_dp           ! [m] the length L_a of the smoothness of (a - a_prior)/sigma_a (0: none)
        real(dp) :: sigma_smb_rel = 0.5_dp        ! sd of a relative to |a_prior|
        real(dp) :: sigma_smb_min = 0.1_dp        ! [m ice/yr] its floor
        logical  :: cold = .true.                 ! penalise a < 0 where a_prior >= 0 (ice cells without a margin face)
        real(dp) :: sigma_cold = 0.01_dp          ! [m ice/yr] the scale of that penalty
        real(dp) :: sigma_def = 0.1_dp            ! sd of ln U below the deformation speed (the soft no-sliding bound)
        real(dp) :: u_min = 1.0e-3_dp             ! [m/yr] floor of U_law and of the bound
        integer  :: n_iter = 100                  ! Gauss-Newton steps
        real(dp) :: tol = 1.0e-6_dp               ! stop at a relative decrease of J below tol with a full step
        real(dp) :: lin_tol = 1.0e-6_dp           ! relative residual of the step system (conjugate gradients)
        integer  :: lin_maxit = 500
        logical  :: verbose = .false.
    end type plastic_inv_ctrl_t

    type plastic_inv_info_t
        integer  :: iters = 0                     ! Gauss-Newton steps taken
        real(dp) :: J = 0.0_dp                    ! the objective at the solution, its four terms
        real(dp) :: J_law = 0.0_dp, J_smooth = 0.0_dp, J_smb = 0.0_dp, J_cold = 0.0_dp, J_def = 0.0_dp
        real(dp) :: J_smooth_smb = 0.0_dp
        real(dp) :: J_start = 0.0_dp              ! at U = U_law (the forward closure)
        integer  :: n_bound = 0                   ! cells below the deformation speed (negative sliding)
        integer  :: n_cold_neg = 0                ! cold cells with a < 0
    end type plastic_inv_info_t

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
        integer  :: nx, ny, code, n, i, j, k, it, kk, s, ii, jj
        integer, allocatable :: id(:,:), ci(:), cj(:), sid(:,:)
        real(dp), allocatable :: sm(:,:), hm(:,:), sf(:,:), cm(:,:), bm(:,:)
        real(dp), allocatable :: D(:,:), cmar(:), clan(:), xlaw(:), xlo(:), x(:), xn(:), dx_(:), g(:), w(:), r(:)
        real(dp), allocatable :: ap(:), sa(:), dh(:), u(:), a(:), udf(:), sc(:), dg(:), eps(:), t1(:), t2(:)
        logical,  allocatable :: cold(:)
        real(dp) :: rg, tau, ub, ud, gx, gy, sux, suy, tax, tay, Jc, Jn, alpha, lsum(6), lsn(6), lres
        integer  :: lit
        integer,  allocatable :: ord(:)
        real(dp), allocatable :: dpre(:), wk(:)
        logical  :: ok, conv

        if (present(ctrl)) c = ctrl
        ps%stat = PLASTIC_OK
        ps%msg  = ""
        nx = ps%nx;  ny = ps%ny
        code = sliding_law_code(slide)
        if (.not. allocated(ps%H) .or. code == 0) then
            ps%stat = PLASTIC_BAD_INPUT
            ps%msg  = "plastic_smb_inversion: not initialized or unknown sliding law"
        else if (.not. (c%sigma_u > 0.0_dp .and. c%sigma_smb_min > 0.0_dp .and. c%sigma_smb_rel >= 0.0_dp .and. &
                        c%sigma_cold > 0.0_dp .and. c%sigma_def > 0.0_dp .and. c%smooth >= 0.0_dp .and. &
                        c%smooth_smb >= 0.0_dp .and. c%u_min > 0.0_dp)) then
            ps%stat = PLASTIC_BAD_INPUT
            ps%msg  = "plastic_smb_inversion: sigma_u, sigma_smb_min, sigma_cold, sigma_def, u_min > 0; sigma_smb_rel, "// &
                      "smooth, smooth_smb >= 0"
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
        sux = (c%smooth/ps%dx)**2/c%sigma_u**2;  suy = (c%smooth/ps%dy)**2/c%sigma_u**2
        tax = (c%smooth_smb/ps%dx)**2;           tay = (c%smooth_smb/ps%dy)**2

        ! --- Gauss-Newton -------------------------------------------------------------------------------------------
        ! the preconditioner: cells by descending surface (the upwind mass balance is lower triangular in that order),
        ! its diagonal with the outflow to land (so every row has one; a cell without any outflow takes a small one)
        ord = order_desc([(sf(ci(k), cj(k)), k = 1, n)])
        dpre = D(0, :) + clan
        where (.not. dpre > 0.0_dp) dpre = 1.0e-3_dp*max(maxval(dpre), tiny(1.0_dp))
        allocate(wk(n))
        allocate(x(n), xn(n), dx_(n), g(n), w(n), r(n), u(n), a(n), sc(n), dg(n), eps(n), t1(n), t2(n))
        x = xlaw
        call objective(x, Jc, lsum)
        inf%J_start = Jc
        conv = .false.
        do it = 1, c%n_iter
            inf%iters = it
            ! the gradient at x and the diagonal parts of H: w = sigma_a^2 C (>= 1; the cold penalty where active),
            ! dg = Lambda
            u = exp(x)
            call mass(u, a)
            eps = (a - ap)/sa
            w = 1.0_dp
            where (cold .and. a < 0.0_dp) w = 1.0_dp + (sa/c%sigma_cold)**2
            call lap(eps, tax, tay, t1)                       ! S_a eps
            r = (eps + t1)/sa                                 ! dJ/da: the misfit and its smoothness,
            where (cold .and. a < 0.0_dp) r = r + a/c%sigma_cold**2   ! the cold penalty
            call bt(r, g)                                     ! B' r
            call lap(x - xlaw, sux, suy, t1)                  ! S_u e
            g = g + (x - xlaw)/c%sigma_u**2 + t1
            dg = 1.0_dp/c%sigma_u**2
            where (x < xlo)
                g = g + (x - xlo)/c%sigma_def**2
                dg = dg + 1.0_dp/c%sigma_def**2
            end where
            ! the system H dx = -g, by conjugate gradients preconditioned with the inverse of its mass-balance part
            sc = -g
            call pcg(sc, dx_, ok, lit, lres)
            ! backtracking on J
            alpha = 1.0_dp
            do kk = 1, 30
                xn = x + alpha*dx_
                call objective(xn, Jn, lsn)
                if (Jn < Jc) exit
                alpha = 0.5_dp*alpha
            end do
            if (c%verbose) write(*, "(a,i4,a,es12.5,a,es10.3,a,i0,a,l1,a,i0,a,es9.2)") "plastic_smb_inversion: it ", it, &
                " J ", Jn, " step ", alpha, " below u_def ", count(xn < xlo), " cg ok ", ok, " it ", lit, " relres ", lres
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
        inf%J_law = lsum(1);  inf%J_smooth = lsum(2);  inf%J_smb = lsum(3);  inf%J_cold = lsum(4);  inf%J_def = lsum(5)
        inf%J_smooth_smb = lsum(6)
        inf%n_bound = count(x < xlo - c%sigma_def)
        inf%n_cold_neg = count(cold .and. a < 0.0_dp)
        if (present(info)) info = inf
        if (.not. conv) then
            ps%stat = PLASTIC_NOT_CONVERGED
            ps%msg  = "plastic_smb_inversion: no convergence in n_iter steps"
        end if
        if (present(stat)) stat = ps%stat

    contains

        subroutine pcg(b, xs, okk, iters, relres)
            ! H xs = b by preconditioned conjugate gradients (H applied by matvec), from xs = 0
            real(dp), intent(in)  :: b(:)
            real(dp), intent(out) :: xs(:)
            logical,  intent(out) :: okk
            integer,  intent(out) :: iters
            real(dp), intent(out) :: relres
            real(dp), allocatable :: rr(:), zz(:), pp(:), qq(:)
            real(dp) :: rz, rzn, bn, al
            integer  :: m
            allocate(rr(n), zz(n), pp(n), qq(n))
            xs = 0.0_dp;  rr = b
            bn = sqrt(sum(b**2));  okk = .true.;  iters = 0;  relres = 0.0_dp
            if (.not. bn > 0.0_dp) return
            call precond(rr, zz)
            pp = zz;  rz = sum(rr*zz)
            do m = 1, c%lin_maxit
                iters = m
                call matvec(pp, qq)
                al = rz/sum(pp*qq)
                xs = xs + al*pp;  rr = rr - al*qq
                relres = sqrt(sum(rr**2))/bn
                if (relres < c%lin_tol) return
                call precond(rr, zz)
                rzn = sum(rr*zz)
                pp = zz + (rzn/rz)*pp;  rz = rzn
            end do
            okk = .false.
        end subroutine pcg

        subroutine matvec(v, hv)
            ! hv = H v = dg v + S_u v + B' Sigma^-1 (w + S_a) Sigma^-1 B v
            real(dp), intent(in)  :: v(:)
            real(dp), intent(out) :: hv(:)
            real(dp), allocatable :: y(:), ly(:)
            allocate(y(n), ly(n))
            call bv(v, y)
            y = y/sa
            call lap(y, tax, tay, ly)
            y = (w*y + ly)/sa
            call bt(y, hv)
            call lap(v, sux, suy, ly)
            hv = hv + dg*v + ly
        end subroutine matvec

        subroutine bv(v, y)
            ! y = B v, B = D diag(U)
            real(dp), intent(in)  :: v(:)
            real(dp), intent(out) :: y(:)
            integer :: m, z, pz
            !$omp parallel do private(m, z, pz) schedule(static)
            do m = 1, n
                y(m) = 0.0_dp
                do z = 0, 4
                    pz = sid(z, m)
                    if (pz > 0) y(m) = y(m) + D(z, m)*u(pz)*v(pz)
                end do
            end do
            !$omp end parallel do
        end subroutine bv

        subroutine bt(y, v)
            ! v = B' y
            real(dp), intent(in)  :: y(:)
            real(dp), intent(out) :: v(:)
            integer :: m, z, pz
            v = 0.0_dp
            do m = 1, n
                do z = 0, 4
                    pz = sid(z, m)
                    if (pz > 0) v(pz) = v(pz) + D(z, m)*u(pz)*y(m)
                end do
            end do
        end subroutine bt

        subroutine lap(v, wx, wy, lv)
            ! lv = L v, L the graph Laplacian over the faces between ice cells, weights wx (x faces), wy (y faces)
            real(dp), intent(in)  :: v(:), wx, wy
            real(dp), intent(out) :: lv(:)
            integer :: m, z, pz
            real(dp) :: wt
            lv = 0.0_dp
            if (wx <= 0.0_dp .and. wy <= 0.0_dp) return
            !$omp parallel do private(m, z, pz, wt) schedule(static)
            do m = 1, n
                do z = 1, 4
                    pz = sid(z, m)
                    if (pz == 0) cycle
                    wt = merge(wx, wy, z <= 2)
                    lv(m) = lv(m) + wt*(v(m) - v(pz))
                end do
            end do
            !$omp end parallel do
        end subroutine lap

        subroutine winv(y, z)
            ! z = W^-1 y = Sigma (w + S_a)^-1 Sigma y: a division, or Jacobi-preconditioned CG on w + S_a
            real(dp), intent(in)  :: y(:)
            real(dp), intent(out) :: z(:)
            real(dp), allocatable :: bb(:), rr(:), pp(:), qq(:), dd(:), zz(:)
            real(dp) :: rz, rzn, bn, al
            integer :: m, z_
            if (c%smooth_smb <= 0.0_dp) then
                z = sa*sa*y/w
                return
            end if
            allocate(bb(n), rr(n), pp(n), qq(n), dd(n), zz(n))
            bb = sa*y
            dd = w
            do m = 1, n
                do z_ = 1, 4
                    if (sid(z_, m) > 0) dd(m) = dd(m) + merge(tax, tay, z_ <= 2)
                end do
            end do
            z = bb/dd;  call lap(z, tax, tay, qq);  rr = bb - (w*z + qq)
            bn = sqrt(sum(bb**2))
            if (bn > 0.0_dp) then
                zz = rr/dd;  pp = zz;  rz = sum(rr*zz)
                do m = 1, 2000
                    call lap(pp, tax, tay, qq);  qq = w*pp + qq
                    al = rz/sum(pp*qq)
                    z = z + al*pp;  rr = rr - al*qq
                    if (sqrt(sum(rr**2)) < 1.0e-10_dp*bn) exit
                    zz = rr/dd;  rzn = sum(rr*zz)
                    pp = zz + (rzn/rz)*pp;  rz = rzn
                end do
            end if
            z = sa*z
        end subroutine winv

        subroutine precond(v, zv)
            ! zv = U^-1 Dp^-1 W^-1 Dp^-T U^-1 v; Dp: D with dpre on
            ! the diagonal, lower triangular in the order ord (row m: its own column and its upstream neighbours')
            real(dp), intent(in)  :: v(:)
            real(dp), intent(out) :: zv(:)
            integer :: o, m, z, pz
            ! y = Dp^-T (U^-1 v): from the lowest cell up; wk accumulates sum over downstream rows of D(z, row) y(row)
            wk = 0.0_dp
            do o = n, 1, -1
                m = ord(o)
                zv(m) = (v(m)/u(m) - wk(m))/dpre(m)
                do z = 1, 4
                    pz = sid(z, m)
                    if (pz > 0 .and. D(z, m) /= 0.0_dp) wk(pz) = wk(pz) + D(z, m)*zv(m)
                end do
            end do
            ! W^-1, then Dp^-1: from the highest cell down
            call winv(zv, wk)
            zv = wk
            do o = 1, n
                m = ord(o)
                do z = 1, 4
                    pz = sid(z, m)
                    if (pz > 0 .and. D(z, m) /= 0.0_dp) wk(m) = wk(m) - D(z, m)*zv(pz)
                end do
                zv(m) = wk(m)/dpre(m)
            end do
            zv = zv/u
        end subroutine precond

        function order_desc(key) result(o)
            ! the indices of key in descending order (merge sort)
            real(dp), intent(in) :: key(:)
            integer, allocatable :: o(:), tmp(:)
            integer :: nn, width, lo, mid, hi, ia, ib, ko
            nn = size(key)
            o = [(ko, ko = 1, nn)];  allocate(tmp(nn))
            width = 1
            do while (width < nn)
                do lo = 1, nn, 2*width
                    mid = min(lo + width - 1, nn);  hi = min(lo + 2*width - 1, nn)
                    ia = lo;  ib = mid + 1;  ko = lo
                    do while (ia <= mid .and. ib <= hi)
                        if (key(o(ia)) >= key(o(ib))) then
                            tmp(ko) = o(ia);  ia = ia + 1
                        else
                            tmp(ko) = o(ib);  ib = ib + 1
                        end if
                        ko = ko + 1
                    end do
                    do while (ia <= mid);  tmp(ko) = o(ia);  ia = ia + 1;  ko = ko + 1;  end do
                    do while (ib <= hi);   tmp(ko) = o(ib);  ib = ib + 1;  ko = ko + 1;  end do
                end do
                o = tmp
                width = 2*width
            end do
        end function order_desc

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

        subroutine objective(xx, jj_, parts)
            real(dp), intent(in)  :: xx(:)
            real(dp), intent(out) :: jj_, parts(6)
            real(dp), allocatable :: uu(:), aa(:), e(:), ee(:), le(:)
            uu = exp(xx)
            allocate(aa(n))
            call mass(uu, aa)
            e = xx - xlaw
            parts = 0.0_dp
            allocate(le(n))
            parts(1) = 0.5_dp*sum(e**2)/c%sigma_u**2
            call lap(e, sux, suy, le)
            parts(2) = 0.5_dp*sum(e*le)                     ! e' L e = sum over faces of w (e_A - e_B)^2
            ee = (aa - ap)/sa
            parts(3) = 0.5_dp*sum(ee**2)
            call lap(ee, tax, tay, le)
            parts(6) = 0.5_dp*sum(ee*le)
            parts(4) = 0.5_dp*sum(merge(min(aa, 0.0_dp)/c%sigma_cold, 0.0_dp, cold)**2)
            parts(5) = 0.5_dp*sum((min(xx - xlo, 0.0_dp)/c%sigma_def)**2)
            jj_ = sum(parts)
        end subroutine objective

    end subroutine plastic_smb_inversion

end module plastic_inversion
