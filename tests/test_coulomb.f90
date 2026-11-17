! The Coulomb yield law in mode "surface" (plastic_solve_coulomb): tau_b = mu N_eff(H),
! 1/N_eff = 1/N + 1/(rho_i g H_N), N = rho_i g H (1 - H_f/H)^p. Kernel properties (the
! update solves its equation, lies at or above flotation, is monotone), the exact flat-land
! profile H^2/2 + H_N H = mu H_N x, the wedge limit (H << H_N: surface slope mu on any bed),
! the Nye limit (H >> H_N: the constant-tau surface solve with tau_b = mu rho_i g H_N), a
! marine bed (the profile thins to flotation at the margin; against an ODE reference), the
! wavefront sweep and the input checks.
program test_coulomb

    use plastic_ice
    use plastic_kinds,      only: BIG
    use plastic_eikonal,    only: godunov_coulomb, coulomb_neff
    use plastic_test_utils

    implicit none

    real(dp), parameter :: rho_i = 917.0_dp, rho_w = 1025.0_dp, grav = 9.80665_dp

    write(*, '(a)') "test_coulomb"
    call test_kernel()
    call test_flat_land()
    call test_wedge()
    call test_nye()
    call test_marine()
    call test_wavefront()
    call test_inputs()
    call finish("test_coulomb")

contains

    ! sqrt(T_a^2 + T_b^2) - 1 of godunov_coulomb's equation at u (with its clamp of the
    ! neighbours on a marine cell), from coulomb_neff.
    function coulomb_eq(u, a, da, b, db, zb, mu, hf, p, H_N) result(g)
        real(dp), intent(in) :: u, a, da, b, db, zb, mu, hf, p, H_N
        real(dp) :: g
        g = sqrt(eq_term(u, a, da, zb, mu, hf, p, H_N) + eq_term(u, b, db, zb, mu, hf, p, H_N)) - 1.0_dp
    end function coulomb_eq

    ! ((u - a)+ / (d s(h)))^2 of one axis (the neighbour clamped to the flotation surface).
    function eq_term(u, aa, dd, zb, mu, hf, p, H_N) result(term)
        real(dp), intent(in) :: u, aa, dd, zb, mu, hf, p, H_N
        real(dp) :: term, ac, h, sl
        term = 0.0_dp
        if (aa >= BIG) return
        ac = aa
        if (hf > 0.0_dp) ac = max(aa, zb + hf)
        if (.not. u > ac) return
        h = 0.5_dp * (u - zb + max(ac - zb, 0.0_dp))
        if (h > 0.0_dp) then
            sl = mu * coulomb_neff(h, hf, p, H_N) / h
        else
            sl = mu
        end if
        term = ((u - ac) / (dd * sl))**2
    end function eq_term

    ! Random land and marine cells (neighbours above and below the bed and the flotation
    ! surface, dx /= dy, p in [0, 1]): the update solves its equation, or brackets its root
    ! within 1e-9 H (the left side of the equation changes sign across u -+ 1e-9 H: near
    ! flotation with p < 1 the left side rises like (H - H_f)^(1-p), steeply; with p = 1 it
    ! jumps at a neighbour at flotation, and u = the lower bound), or is the lower bound where
    ! the left side already exceeds 1 (a land cell above a cliff: H = 0), lies at or above flotation,
    ! and is nondecreasing in a and b.
    subroutine test_kernel()
        integer  :: n, nbad_eq, nbad_mono, nbad_flt, nlo, nseed
        integer, allocatable :: seed(:)
        real(dp) :: r(10), a, b, zb, mu, hf, p, H_N, da, db, u, u2, emax, lo, e, dl
        write(*, '(a)') "Coulomb kernel"
        call check_close(godunov_coulomb(0.0_dp, 1000.0_dp, BIG, 1000.0_dp, 0.0_dp, 0.01_dp, 0.0_dp, 1.0_dp, 1.0e12_dp), &
                         10.0_dp, 1.0e-9_dp, 0.0_dp, "land margin cell, H << H_N: u = a + mu d (wedge)")
        call check(godunov_coulomb(BIG, 1.0_dp, BIG, 1.0_dp, 0.0_dp, 0.01_dp, 0.0_dp, 1.0_dp, 1.0e3_dp) >= BIG, &
                   "no upwind information -> BIG")
        call random_seed(size=nseed)
        allocate(seed(nseed))
        seed = 4242
        call random_seed(put=seed)
        nbad_eq = 0; nbad_mono = 0; nbad_flt = 0; nlo = 0; emax = 0.0_dp
        do n = 1, 200000
            call random_number(r)
            zb  = 2000.0_dp * (r(1) - 0.6_dp)
            hf  = (rho_w / rho_i) * max(-zb, 0.0_dp)
            a   = zb + 3000.0_dp * (r(2) - 0.3_dp)
            b   = zb + 3000.0_dp * (r(3) - 0.3_dp)
            if (r(10) < 0.1_dp) b = BIG
            mu  = 1.0e-3_dp + 0.05_dp * r(4)
            da  = 500.0_dp + 4.0e4_dp * r(5)
            db  = 500.0_dp + 4.0e4_dp * r(6)
            p   = r(7)
            if (r(8) < 0.3_dp) p = 1.0_dp
            H_N = 10.0_dp * 10.0_dp**(3.0_dp * r(9))
            u = godunov_coulomb(a, da, b, db, zb, mu, hf, p, H_N)
            if (u < zb + hf) nbad_flt = nbad_flt + 1
            lo = max(min(a, b), zb)
            if (hf > 0.0_dp) lo = max(min(max(a, zb + hf), max(b, zb + hf)), zb + hf)
            e  = coulomb_eq(u, a, da, b, db, zb, mu, hf, p, H_N)
            dl = 1.0e-9_dp * max(1.0_dp, u - zb)
            if (u == lo) nlo = nlo + 1
            if (abs(e) <= 1.0e-9_dp) then
                emax = max(emax, abs(e))
            else if (u == lo .and. e > 0.0_dp) then        ! the lower bound already steep enough (H = 0 on land)
                continue
            else if (.not. (coulomb_eq(u - dl, a, da, b, db, zb, mu, hf, p, H_N) <= 0.0_dp .and. &
                            coulomb_eq(u + dl, a, da, b, db, zb, mu, hf, p, H_N) >= 0.0_dp)) then
                nbad_eq = nbad_eq + 1
                if (nbad_eq <= 3) write(*, '(a,11es12.4)') "        bad: a b zb mu hf p H_N da db u e ", &
                    a, b, zb, mu, hf, p, H_N, da, db, u, e
            end if
            u2 = godunov_coulomb(a + 1.0_dp, da, b, db, zb, mu, hf, p, H_N)
            if (u2 < u - 1.0e-9_dp * max(1.0_dp, abs(u))) nbad_mono = nbad_mono + 1
            if (b < BIG) then
                u2 = godunov_coulomb(a, da, b + 1.0_dp, db, zb, mu, hf, p, H_N)
                if (u2 < u - 1.0e-9_dp * max(1.0_dp, abs(u))) nbad_mono = nbad_mono + 1
            end if
        end do
        write(*, '(a,es10.3,a,i0,a)') "        max relative equation residual where <= 1e-9: ", emax, " (", nlo, &
            " cases at the lower bound)"
        call check(nbad_eq == 0, "update solves its equation to 1e-9 or brackets its root within 1e-9 H (2e5 cases)")
        call check(nbad_flt == 0, "update at or above flotation (z_b + H_f)")
        call check(nbad_mono == 0, "update nondecreasing in a and b (monotone)")
    end subroutine test_kernel

    subroutine solve_c(nx, ny, dx, z_b, z_sl, mu, mask, p, H_N, H, tau, sweep)
        integer,  intent(in)  :: nx, ny
        real(dp), intent(in)  :: dx, z_b(nx, ny), z_sl(nx, ny), mu(nx, ny), p, H_N
        logical,  intent(in)  :: mask(nx, ny)
        real(dp), intent(out) :: H(nx, ny)
        real(dp), optional, intent(out) :: tau(nx, ny)
        character(len=*), optional, intent(in) :: sweep
        type(plastic_t)     :: ps
        type(plastic_par_t) :: par
        integer :: st
        par%mode = "surface"; par%max_sweeps = 1000; par%tol = 1.0e-10_dp
        if (present(sweep)) par%sweep = sweep
        call plastic_init(ps, nx, ny, dx, dx, par)
        call plastic_solve_coulomb(ps, z_b, z_sl, mu, mask, p, H_N, stat=st)
        if (st /= PLASTIC_OK) write(*, '(a,a)') "        ", trim(ps%msg)
        H = ps%H
        if (present(tau)) tau = ps%tau_b
        call plastic_end(ps)
    end subroutine solve_c

    ! 1-D strip on flat land, margin at the face: dH/dx = mu H_N/(H + H_N) integrates to
    ! H^2/2 + H_N H = mu H_N x, x the distance to the margin. The edge-mean update
    ! reproduces it exactly (the discrete difference of H^2/2 + H_N H), for any dx and H_N.
    subroutine test_flat_land()
        real(dp), parameter :: L = 200000.0_dp, mu0 = 0.01_dp
        real(dp) :: dxs(3), hns(3), dx, x, ex, emax
        real(dp), allocatable :: z_b(:,:), z_sl(:,:), mu(:,:), H(:,:)
        logical,  allocatable :: mask(:,:)
        integer  :: m, k, nx, i
        write(*, '(a)') "flat land: H^2/2 + H_N H = mu H_N x"
        dxs = [8000.0_dp, 2000.0_dp, 500.0_dp]
        hns = [100.0_dp, 1000.0_dp, 5000.0_dp]
        emax = 0.0_dp
        do m = 1, 3
            dx = dxs(m)
            nx = nint(L / dx) + 1
            allocate(z_b(nx, 3), z_sl(nx, 3), mu(nx, 3), H(nx, 3), mask(nx, 3))
            z_b = 0.0_dp; z_sl = -1.0e4_dp; mu = mu0
            mask = .true.; mask(1, :) = .false.; mask(nx, :) = .false.
            do k = 1, 3
                call solve_c(nx, 3, dx, z_b, z_sl, mu, mask, 1.0_dp, hns(k), H)
                do i = 2, nx - 1
                    x  = dx * (min(i - 1, nx - i) - 0.5_dp)
                    ex = -hns(k) + sqrt(hns(k)**2 + 2.0_dp * mu0 * hns(k) * x)
                    emax = max(emax, abs(H(i, 2) - ex) / ex)
                end do
            end do
            deallocate(z_b, z_sl, mu, H, mask)
        end do
        write(*, '(a,es10.3)') "        max rel err ", emax
        call check(emax < 1.0e-10_dp, "flat land: exact for every dx and H_N")
    end subroutine test_flat_land

    ! Wedge limit, H << H_N on land: |grad z_s| = mu whatever the bed. A bed rising inland
    ! (slope 0.004 < mu = 0.01) and one with a 300 m step: z_s(x) = z_s(margin) + mu x
    ! (margin at the face, z_s(margin) the ice-free cell's bed).
    subroutine test_wedge()
        real(dp), parameter :: dx = 2000.0_dp, mu0 = 0.01_dp
        integer,  parameter :: nx = 101
        real(dp) :: z_b(nx, 3), z_sl(nx, 3), mu(nx, 3), H(nx, 3), tau(nx, 3), x, ex, emax, etau
        logical  :: mask(nx, 3)
        integer  :: i, c
        write(*, '(a)') "wedge limit (H_N = 1e15 m): surface slope mu on any bed"
        emax = 0.0_dp; etau = 0.0_dp
        do c = 1, 2
            do i = 1, nx
                z_b(i, :) = 0.004_dp * dx * (i - 1)
                if (c == 2 .and. i > 40) z_b(i, :) = z_b(i, :) + 300.0_dp
            end do
            z_sl = -1.0e4_dp; mu = mu0
            mask = .true.; mask(1, :) = .false.
            call solve_c(nx, 3, dx, z_b, z_sl, mu, mask, 1.0_dp, 1.0e15_dp, H, tau)
            do i = 2, nx
                x  = dx * (i - 1.5_dp)
                ex = z_b(1, 2) + mu0 * x
                if (ex <= z_b(i, 2)) cycle                    ! (not in these cases)
                emax = max(emax, abs(H(i, 2) + z_b(i, 2) - ex))
                etau = max(etau, abs(tau(i, 2) - mu0 * rho_i * grav * H(i, 2)) / (mu0 * rho_i * grav * H(i, 2)))
            end do
        end do
        write(*, '(a,es10.3,a,es10.3)') "        max |z_s - z_m - mu x| ", emax, " m; tau_b vs mu rho_i g H ", etau
        call check(emax < 1.0e-6_dp, "wedge: z_s linear with slope mu (sloped bed, bed step)")
        call check(etau < 1.0e-8_dp, "wedge: tau_b = mu rho_i g H")
    end subroutine test_wedge

    ! Nye limit, H >> H_N: the Coulomb solve equals the constant-tau surface solve with
    ! tau_b = mu rho_i g H_N, on the bed bump of test_surface (2-D, land, sub-cell margin).
    subroutine test_nye()
        integer,  parameter :: n = 80
        real(dp), parameter :: dx = 1000.0_dp, HN = 1.0e-5_dp, tb0 = 1.0e5_dp
        type(plastic_t)     :: ps
        type(plastic_par_t) :: par
        real(dp) :: z_b(n, n), z_sl(n, n), mu(n, n), tb(n, n), Hs(n, n), Hc(n, n), sd(n, n), mf(n, n), x, y, e
        logical  :: mask(n, n)
        integer  :: i, j
        write(*, '(a)') "Nye limit (H_N = 1e-5 m): constant-tau surface solve with tau_b = mu rho_i g H_N"
        do j = 1, n
            do i = 1, n
                x = (i - 0.5_dp) * dx - 40000.0_dp
                y = (j - 0.5_dp) * dx - 40000.0_dp
                z_b(i, j)  = 300.0_dp * exp(-(x * x + y * y) / (2.0_dp * 2.0e4_dp**2))
                tb(i, j)   = tb0 * (1.0_dp + 0.3_dp * cos(j / 9.0_dp))
                sd(i, j)   = 30000.0_dp - hypot(x - 700.0_dp, y + 300.0_dp)
                mask(i, j) = sd(i, j) > 0.0_dp
            end do
        end do
        z_sl = -1.0e4_dp
        mu = tb / (rho_i * grav * HN)
        call plastic_margin_frac_from_sdist(sd, mask, dx, dx, mf)
        par%mode = "surface"; par%max_sweeps = 1000; par%tol = 1.0e-10_dp
        call plastic_init(ps, n, n, dx, dx, par)
        call plastic_solve(ps, z_b, z_sl, tb, mask, margin_frac=mf)
        Hs = ps%H
        call plastic_solve_coulomb(ps, z_b, z_sl, mu, mask, 1.0_dp, HN, margin_frac=mf)
        Hc = ps%H
        e = maxval(abs(Hc - Hs) / max(Hs, 1.0_dp), mask=mask)
        write(*, '(a,es10.3,a,f7.1,a)') "        max |H_c - H_s|/H_s ", e, " (max H ", maxval(Hs), " m)"
        call check(e < 1.0e-5_dp, "Nye limit: matches the constant-tau surface solve (bed bump, sub-cell margin)")
        call check(maxval(abs(ps%tau_b - tb) / tb, mask=mask) < 1.0e-5_dp, "Nye limit: tau_b = mu rho_i g H_N")
        call plastic_end(ps)
    end subroutine test_nye

    ! Marine margin: a 1-D bed rising inland from 800 m below sea level, z_b = -800 + 0.004 x
    ! (x from the margin face), p = 1, H_N = 1000 m, mu = 0.01. Reference: the ODE
    ! dz_s/dx = mu N_eff(H)/(rho_i g H), z_s(0) the flotation surface (RK4, fine step).
    ! The ice is grounded at or above flotation everywhere, thins to flotation at the margin
    ! (H - H_f of the first cell -> 0 with dx), tau_b -> 0 there, and converges to the ODE.
    subroutine test_marine()
        real(dp), parameter :: L = 300000.0_dp, s = 0.004_dp, z0 = -800.0_dp, mu0 = 0.01_dp, HN = 1000.0_dp
        real(dp), parameter :: hstep = 5.0_dp
        real(dp) :: dxs(4), err(4), gap(4), tau1(4), dx, x, zs, k1, k2, k3, k4
        real(dp), allocatable :: z_b(:,:), z_sl(:,:), mu(:,:), H(:,:), tau(:,:), zref(:), hfl(:,:)
        logical,  allocatable :: mask(:,:)
        logical  :: grounded
        integer  :: m, nx, i, nref
        write(*, '(a)') "marine bed (z_b = -800 m + 0.004 x, p = 1, H_N = 1000 m): vs ODE"
        nref = nint(L / hstep)
        allocate(zref(0:nref))
        zref(0) = z0 + (rho_w / rho_i) * max(-z0, 0.0_dp)
        do i = 1, nref
            x  = (i - 1) * hstep
            zs = zref(i - 1)
            k1 = marine_slope(z0 + s * x, zs, mu0, HN)
            k2 = marine_slope(z0 + s * (x + 0.5_dp * hstep), zs + 0.5_dp * hstep * k1, mu0, HN)
            k3 = marine_slope(z0 + s * (x + 0.5_dp * hstep), zs + 0.5_dp * hstep * k2, mu0, HN)
            k4 = marine_slope(z0 + s * (x + hstep), zs + hstep * k3, mu0, HN)
            zref(i) = zs + hstep * (k1 + 2.0_dp * k2 + 2.0_dp * k3 + k4) / 6.0_dp
        end do
        dxs = [8000.0_dp, 4000.0_dp, 2000.0_dp, 1000.0_dp]
        grounded = .true.
        do m = 1, 4
            dx = dxs(m)
            nx = nint(L / dx) + 1
            allocate(z_b(nx, 3), z_sl(nx, 3), mu(nx, 3), H(nx, 3), tau(nx, 3), hfl(nx, 3), mask(nx, 3))
            do i = 1, nx
                z_b(i, :) = z0 + s * dx * (i - 1.5_dp)
            end do
            z_sl = 0.0_dp; mu = mu0
            hfl = (rho_w / rho_i) * max(z_sl - z_b, 0.0_dp)
            mask = .true.; mask(1, :) = .false.
            call solve_c(nx, 3, dx, z_b, z_sl, mu, mask, 1.0_dp, HN, H, tau)
            grounded = grounded .and. all(H(2:, 2) >= hfl(2:, 2) - 1.0e-9_dp)
            err(m) = 0.0_dp
            do i = 2, nx
                x = dx * (i - 1.5_dp)
                err(m) = max(err(m), abs(H(i, 2) + z_b(i, 2) - zref(nint(x / hstep))))
            end do
            gap(m)  = H(2, 2) - hfl(2, 2)
            tau1(m) = tau(2, 2) / maxval(tau(:, 2))
            write(*, '(a,f6.0,a,f8.2,a,f8.3,a,f7.4)') "        dx ", dx, ": max |z_s - ref| ", err(m), &
                " m; first cell H - H_f ", gap(m), " m, tau_b / max ", tau1(m)
            deallocate(z_b, z_sl, mu, H, tau, hfl, mask)
        end do
        call check(grounded, "marine: H >= H_f on every ice cell")
        call check(err(4) < err(3) .and. err(3) < err(2) .and. err(2) < err(1), "marine: converges to the ODE")
        call check(gap(4) < gap(3) .and. gap(3) < gap(2) .and. gap(2) < gap(1), &
                   "marine: thins to flotation at the margin (H - H_f of the first cell -> 0)")
        call check(tau1(4) < 0.05_dp, "marine: tau_b -> 0 at the margin")
    end subroutine test_marine

    ! The marine ODE's slope mu N_eff/(rho_i g H) at bed zb, surface zz (z_sl = 0).
    real(dp) function marine_slope(zb, zz, mu0, HN) result(sl)
        real(dp), intent(in) :: zb, zz, mu0, HN
        real(dp) :: hh
        hh = zz - zb
        sl = 0.0_dp
        if (hh > 0.0_dp) sl = mu0 * coulomb_neff(hh, (rho_w / rho_i) * max(-zb, 0.0_dp), 1.0_dp, HN) / hh
    end function marine_slope

    ! The wavefront sweep: bit-identical to the serial sweep (land and marine, p < 1).
    subroutine test_wavefront()
        integer,  parameter :: n = 90
        real(dp), parameter :: dx = 4000.0_dp
        real(dp) :: z_b(n, n), z_sl(n, n), mu(n, n), H1(n, n), H2(n, n), x, y
        logical  :: mask(n, n)
        integer  :: i, j
        write(*, '(a)') "wavefront sweep == serial sweep"
        do j = 1, n
            do i = 1, n
                x = (i - 0.5_dp) * dx - 0.5_dp * n * dx
                y = (j - 0.5_dp) * dx - 0.5_dp * n * dx
                z_b(i, j)  = 600.0_dp * exp(-(x * x + y * y) / (2.0_dp * 6.0e4_dp**2)) - 400.0_dp + 80.0_dp * sin(i / 7.0_dp)
                mu(i, j)   = 0.01_dp * (1.0_dp + 0.5_dp * cos(j / 9.0_dp))
                mask(i, j) = hypot(x + 3000.0_dp, y - 1700.0_dp) < 0.37_dp * n * dx
            end do
        end do
        z_sl = 0.0_dp
        call solve_c(n, n, dx, z_b, z_sl, mu, mask, 0.5_dp, 800.0_dp, H1, sweep="serial")
        call solve_c(n, n, dx, z_b, z_sl, mu, mask, 0.5_dp, 800.0_dp, H2, sweep="wavefront")
        call check(all(H1 == H2), "wavefront sweep bit-identical to the serial sweep (Coulomb)")
    end subroutine test_wavefront

    ! Mode other than "surface", p outside [0, 1], H_N <= 0, negative mu: PLASTIC_BAD_INPUT.
    subroutine test_inputs()
        integer, parameter :: n = 10
        type(plastic_t)     :: ps
        type(plastic_par_t) :: par
        real(dp) :: z_b(n, n), z_sl(n, n), mu(n, n)
        logical  :: mask(n, n)
        integer  :: st
        write(*, '(a)') "input checks"
        z_b = 0.0_dp; z_sl = 0.0_dp; mu = 0.01_dp
        mask = .false.; mask(3:8, 3:8) = .true.
        par%mode = "full"
        call plastic_init(ps, n, n, 1000.0_dp, 1000.0_dp, par)
        call plastic_solve_coulomb(ps, z_b, z_sl, mu, mask, 1.0_dp, 1000.0_dp, stat=st)
        call check(st == PLASTIC_BAD_INPUT .and. index(ps%msg, "surface") > 0, &
                   'mode "full" -> PLASTIC_BAD_INPUT ("'//trim(ps%msg)//'")')
        ps%par%mode = "surface"
        call plastic_solve_coulomb(ps, z_b, z_sl, mu, mask, 1.5_dp, 1000.0_dp, stat=st)
        call check(st == PLASTIC_BAD_INPUT, "p > 1 -> PLASTIC_BAD_INPUT")
        call plastic_solve_coulomb(ps, z_b, z_sl, mu, mask, 1.0_dp, 0.0_dp, stat=st)
        call check(st == PLASTIC_BAD_INPUT, "H_N = 0 -> PLASTIC_BAD_INPUT")
        mu(5, 5) = -1.0_dp
        call plastic_solve_coulomb(ps, z_b, z_sl, mu, mask, 1.0_dp, 1000.0_dp, stat=st)
        call check(st == PLASTIC_BAD_INPUT, "mu < 0 on the ice -> PLASTIC_BAD_INPUT")
        mu(5, 5) = 0.01_dp
        call plastic_solve_coulomb(ps, z_b, z_sl, mu, mask, 1.0_dp, 1000.0_dp, stat=st)
        call check(st == PLASTIC_OK, "valid input -> PLASTIC_OK")
        call plastic_end(ps)
    end subroutine test_inputs

end program test_coulomb
