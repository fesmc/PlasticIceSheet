! mode "surface": the plastic surface equation solved for z_s directly (no outer
! iteration). Kernel properties, convergence to the analytic Nye profile (against the
! flat w = H^2 solve and a cell-centred variant of the surface update), a sloped bed
! against an ODE reference (flat / full / surface), and the bed bump where full
! converges (surface vs full under refinement).
program test_surface

    use plastic_ice
    use plastic_kinds,      only: BIG
    use plastic_eikonal,    only: godunov_surface
    use plastic_test_utils

    implicit none

    real(dp), parameter :: rho_g = 917.0_dp * 9.80665_dp

    write(*, '(a)') "test_surface"
    call test_kernel()
    call test_nye_strip()
    call test_nye_disc()
    call test_sloped_bed()
    call test_bump()
    call finish("test_surface")

contains

    ! Discrete equation of godunov_surface at u (zero at the update).
    pure function surface_eq(u, a, da, b, db, zb, k) result(g)
        real(dp), intent(in) :: u, a, da, b, db, zb, k
        real(dp) :: g, h, al, be, pa, pb
        h  = u - zb
        al = a - zb
        be = b - zb
        pa = max(h - al, 0.0_dp) * (h + max(al, 0.0_dp)) / 2.0_dp
        pb = max(h - be, 0.0_dp) * (h + max(be, 0.0_dp)) / 2.0_dp
        g  = sqrt((pa / da)**2 + (pb / db)**2) / k - 1.0_dp
    end function surface_eq

    ! Random inputs (neighbours above and below the bed, dx /= dy): the update solves
    ! its equation, lies above the bed, and is nondecreasing in a and b (monotone).
    subroutine test_kernel()
        integer  :: n, nbad_eq, nbad_mono, nbad_bed, nseed
        integer, allocatable :: seed(:)
        real(dp) :: r(6), a, b, zb, k, da, db, u, u2, emax, kk
        write(*, '(a)') "surface kernel"
        kk = 1.0e5_dp / rho_g
        call check_close(godunov_surface(0.0_dp, 1000.0_dp, BIG, 1000.0_dp, 0.0_dp, kk), &
                         sqrt(2.0_dp * kk * 1000.0_dp), 1.0e-14_dp, 0.0_dp, &
                         "flat-bed margin cell: H = sqrt(2 k dx) (Nye)")
        call check(godunov_surface(BIG, 1.0_dp, BIG, 1.0_dp, 0.0_dp, kk) >= BIG, "no upwind information -> BIG")
        call random_seed(size=nseed)
        allocate(seed(nseed))
        seed = 12345
        call random_seed(put=seed)
        nbad_eq = 0; nbad_mono = 0; nbad_bed = 0; emax = 0.0_dp
        do n = 1, 200000
            call random_number(r)
            zb = 2000.0_dp * (r(1) - 0.5_dp)
            a  = zb + 3000.0_dp * (r(2) - 0.3_dp)
            b  = zb + 3000.0_dp * (r(3) - 0.3_dp)
            k  = (5.0e3_dp + 2.0e5_dp * r(4)) / rho_g
            da = 500.0_dp + 4.0e4_dp * r(5)
            db = 500.0_dp + 4.0e4_dp * r(6)
            u = godunov_surface(a, da, b, db, zb, k)
            if (.not. (u > zb)) nbad_bed = nbad_bed + 1
            emax = max(emax, abs(surface_eq(u, a, da, b, db, zb, k)))
            if (abs(surface_eq(u, a, da, b, db, zb, k)) > 1.0e-10_dp) nbad_eq = nbad_eq + 1
            u2 = godunov_surface(a + 1.0_dp, da, b, db, zb, k)
            if (u2 < u - 1.0e-9_dp * abs(u)) nbad_mono = nbad_mono + 1
            u2 = godunov_surface(a, da, b + 1.0_dp, db, zb, k)
            if (u2 < u - 1.0e-9_dp * abs(u)) nbad_mono = nbad_mono + 1
        end do
        write(*, '(a,es10.3)') "        max relative equation residual ", emax
        call check(nbad_eq == 0, "update solves its equation (2e5 random cases)")
        call check(nbad_bed == 0, "update lies above the bed")
        call check(nbad_mono == 0, "update nondecreasing in a and b (monotone)")
    end subroutine test_kernel

    subroutine solve_mode(mode, nx, ny, dx, z_b, tb, mask, H, n_sweeps)
        character(len=*), intent(in)  :: mode
        integer,          intent(in)  :: nx, ny
        real(dp),         intent(in)  :: dx, z_b(nx, ny), tb(nx, ny)
        logical,          intent(in)  :: mask(nx, ny)
        real(dp),         intent(out) :: H(nx, ny)
        integer, optional, intent(out) :: n_sweeps
        type(plastic_t)     :: ps
        type(plastic_par_t) :: par
        real(dp) :: z_sl(nx, ny)
        z_sl = -1.0e4_dp
        par%mode = mode; par%max_sweeps = 1000; par%tol = 1.0e-10_dp
        par%n_outer = 3000; par%outer_tol = 1.0e-7_dp; par%relax = 0.4_dp
        call plastic_init(ps, nx, ny, dx, dx, par)
        call plastic_solve(ps, z_b, z_sl, tb, mask)
        if (.not. ps%converged) write(*, '(a,a,a,a)') "        ", mode, ": ", trim(ps%msg)
        H = ps%H
        if (present(n_sweeps)) n_sweeps = ps%n_sweeps
        call plastic_end(ps)
    end subroutine solve_mode

    ! (i) 1-D strip between land margins: H = sqrt(2 k d). surface and flat reproduce it
    ! to rounding at every dx; the cell-centred update (u - a)(u - z_b) = k dx, solved
    ! here by its 1-D recursion, is sqrt(1/2) low at the first cell.
    subroutine test_nye_strip()
        real(dp), parameter :: L = 64000.0_dp, k = 1.0e5_dp / rho_g, d_fix = 8000.0_dp
        real(dp) :: dxs(5), es(5), ef(5), ec(5), ec_fix(5), dx, d, ex, hc, hprev
        real(dp), allocatable :: z_b(:,:), tb(:,:), Hs(:,:), Hf(:,:)
        logical,  allocatable :: mask(:,:)
        integer  :: m, nx, i, j

        write(*, '(a)') "(i) Nye profile, 1-D strip: max rel err (interior), and at d = 8 km"
        dxs = [4000.0_dp, 2000.0_dp, 1000.0_dp, 500.0_dp, 250.0_dp]
        do m = 1, 5
            dx = dxs(m)
            nx = nint(L / dx) + 1
            allocate(z_b(nx, 4), tb(nx, 4), Hs(nx, 4), Hf(nx, 4), mask(nx, 4))
            z_b = 0.0_dp; tb = 1.0e5_dp
            mask = .true.; mask(1, :) = .false.; mask(nx, :) = .false.
            call solve_mode("surface", nx, 4, dx, z_b, tb, mask, Hs)
            call solve_mode("flat",    nx, 4, dx, z_b, tb, mask, Hf)
            es(m) = 0.0_dp; ef(m) = 0.0_dp; ec(m) = 0.0_dp
            hprev = 0.0_dp
            j = 2
            do i = 2, nx / 2
                d  = dx * (i - 1)
                ex = sqrt(2.0_dp * k * d)
                ! cell-centred recursion: h (h - hprev) = k dx
                hc = (hprev + sqrt(hprev**2 + 4.0_dp * k * dx)) / 2.0_dp
                hprev = hc
                es(m) = max(es(m), abs(Hs(i, j) - ex) / ex)
                ef(m) = max(ef(m), abs(Hf(i, j) - ex) / ex)
                ec(m) = max(ec(m), abs(hc - ex) / ex)
                if (abs(d - d_fix) < 0.5_dp * dx) ec_fix(m) = abs(hc - ex) / ex
            end do
            write(*, '(a,f7.0,a,es9.2,a,es9.2,a,es9.2,a,es9.2)') "        dx ", dx, ": surface ", es(m), &
                ", flat ", ef(m), ", cell-centred ", ec(m), " (", ec_fix(m), " at 8 km)"
            deallocate(z_b, tb, Hs, Hf, mask)
        end do
        write(*, '(a,4f6.2)') "        cell-centred order at 8 km: ", (log(ec_fix(m) / ec_fix(m + 1)) / log(2.0_dp), m = 1, 4)
        call check(maxval(es) < 1.0e-12_dp, "surface reproduces Nye exactly (1-D, all dx)")
        call check(maxval(ef) < 1.0e-12_dp, "flat reproduces Nye exactly (1-D, all dx)")
        call check(abs(ec(5) - (1.0_dp - sqrt(0.5_dp))) < 1.0e-12_dp, &
                   "cell-centred variant: first cell sqrt(1/2) low at every dx")
    end subroutine test_nye_strip

    ! (i) circular cap, land margin: H = sqrt(2 k (R - r)); surface == flat on a flat bed,
    ! and both converge (staircase margin, divide kink).
    subroutine test_nye_disc()
        real(dp), parameter :: R = 30000.0_dp, k = 1.0e5_dp / rho_g
        real(dp) :: dxs(4), es(4), ef(4), dsf(4), dx, rr, ex, x0
        real(dp), allocatable :: z_b(:,:), tb(:,:), Hs(:,:), Hf(:,:)
        logical,  allocatable :: mask(:,:), inner(:,:)
        integer  :: m, n, i, j, ns

        write(*, '(a)') "(i) Nye profile, circular cap: mean abs err over r < R/2 [m]"
        dxs = [2000.0_dp, 1000.0_dp, 500.0_dp, 250.0_dp]
        do m = 1, 4
            dx = dxs(m)
            n = 2 * nint(R / dx) + 6
            x0 = 0.5_dp * (n + 1) * dx
            allocate(z_b(n, n), tb(n, n), Hs(n, n), Hf(n, n), mask(n, n), inner(n, n))
            z_b = 0.0_dp; tb = 1.0e5_dp
            do j = 1, n
                do i = 1, n
                    rr = hypot(i * dx - x0, j * dx - x0)
                    mask(i, j)  = rr < R
                    inner(i, j) = rr < 0.5_dp * R
                end do
            end do
            call solve_mode("surface", n, n, dx, z_b, tb, mask, Hs, ns)
            call solve_mode("flat",    n, n, dx, z_b, tb, mask, Hf)
            es(m) = 0.0_dp; ef(m) = 0.0_dp
            do j = 1, n
                do i = 1, n
                    if (.not. inner(i, j)) cycle
                    ex = sqrt(2.0_dp * k * (R - hypot(i * dx - x0, j * dx - x0)))
                    es(m) = es(m) + abs(Hs(i, j) - ex)
                    ef(m) = ef(m) + abs(Hf(i, j) - ex)
                end do
            end do
            es(m) = es(m) / count(inner); ef(m) = ef(m) / count(inner)
            dsf(m) = maxval(abs(Hs - Hf))
            write(*, '(a,f7.0,a,f8.3,a,f8.3,a,es9.2,a,i0)') "        dx ", dx, ": surface ", es(m), ", flat ", &
                ef(m), ", max|H_surface - H_flat| ", dsf(m), " m, sweeps ", ns
            deallocate(z_b, tb, Hs, Hf, mask, inner)
        end do
        write(*, '(a,3f6.2)') "        order (surface): ", (log(es(m) / es(m + 1)) / log(2.0_dp), m = 1, 3)
        call check(maxval(dsf) < 1.0e-9_dp, "surface == flat on a flat bed (2-D)")
        call check(es(4) < es(3) .and. es(3) < es(2) .and. es(2) < es(1), "error decreases with dx")
    end subroutine test_nye_disc

    ! (ii) 1-D sloped bed z_b = s x rising inland from a land margin at x = 0 (divide at
    ! the far edge): reference from the ODE dw/dx = 2k - 2 s sqrt(w), w = H^2.
    subroutine test_sloped_bed()
        real(dp), parameter :: L = 100000.0_dp, s = 0.005_dp, k = 1.0e5_dp / rho_g, hstep = 0.25_dp
        real(dp) :: dxs(4), es(4), eh(4), ef(4), dx, x, ref, w0, k1, k2, k3, k4
        real(dp), allocatable :: z_b(:,:), tb(:,:), Hs(:,:), Hh(:,:), Hf(:,:), wref(:)
        logical,  allocatable :: mask(:,:)
        integer  :: m, nx, i, nref

        write(*, '(a)') "(ii) sloped bed (s = 5e-3), max abs err for x >= 5 km [m]"
        ! reference w(x) on a fine grid (RK4)
        nref = nint(L / hstep)
        allocate(wref(0:nref))
        wref(0) = 0.0_dp
        do i = 1, nref
            w0 = wref(i - 1)
            k1 = 2.0_dp * k - 2.0_dp * s * sqrt(max(w0, 0.0_dp))
            k2 = 2.0_dp * k - 2.0_dp * s * sqrt(max(w0 + 0.5_dp * hstep * k1, 0.0_dp))
            k3 = 2.0_dp * k - 2.0_dp * s * sqrt(max(w0 + 0.5_dp * hstep * k2, 0.0_dp))
            k4 = 2.0_dp * k - 2.0_dp * s * sqrt(max(w0 + hstep * k3, 0.0_dp))
            wref(i) = w0 + hstep * (k1 + 2.0_dp * k2 + 2.0_dp * k3 + k4) / 6.0_dp
        end do

        dxs = [4000.0_dp, 2000.0_dp, 1000.0_dp, 500.0_dp]
        do m = 1, 4
            dx = dxs(m)
            nx = nint(L / dx) + 1
            allocate(z_b(nx, 3), tb(nx, 3), Hs(nx, 3), Hh(nx, 3), Hf(nx, 3), mask(nx, 3))
            do i = 1, nx
                z_b(i, :) = s * dx * (i - 1)
            end do
            tb = 1.0e5_dp
            mask = .true.; mask(1, :) = .false.
            call solve_mode("surface", nx, 3, dx, z_b, tb, mask, Hs)
            call solve_mode("full",    nx, 3, dx, z_b, tb, mask, Hh)
            call solve_mode("flat",    nx, 3, dx, z_b, tb, mask, Hf)
            es(m) = 0.0_dp; eh(m) = 0.0_dp; ef(m) = 0.0_dp
            do i = 2, nx
                x = dx * (i - 1)
                if (x < 5000.0_dp) cycle
                ref = sqrt(wref(nint(x / hstep)))
                es(m) = max(es(m), abs(Hs(i, 2) - ref))
                eh(m) = max(eh(m), abs(Hh(i, 2) - ref))
                ef(m) = max(ef(m), abs(Hf(i, 2) - ref))
            end do
            write(*, '(a,f7.0,a,f8.2,a,f8.2,a,f8.2)') "        dx ", dx, ": surface ", es(m), ", full ", eh(m), &
                ", flat ", ef(m)
            deallocate(z_b, tb, Hs, Hh, Hf, mask)
        end do
        write(*, '(a,3f6.2)') "        order surface: ", (log(es(m) / es(m + 1)) / log(2.0_dp), m = 1, 3)
        write(*, '(a,3f6.2)') "        order full:    ", (log(eh(m) / eh(m + 1)) / log(2.0_dp), m = 1, 3)
        call check(es(4) < es(3) .and. es(3) < es(2) .and. es(2) < es(1), "surface converges on a sloped bed")
        call check(es(4) < 0.02_dp * ef(4), "surface error << flat (bed relief) error")


    end subroutine test_sloped_bed

    ! (ii) the bed bump of the Julia tests (full converges there): surface and full
    ! approach each other under refinement; flat ignores the bed.
    subroutine test_bump()
        real(dp) :: dxs(3), dsh(3), rsh(3), dfh(3), dx, x, y
        real(dp), allocatable :: z_b(:,:), tb(:,:), Hs(:,:), Hh(:,:), Hf(:,:)
        logical,  allocatable :: mask(:,:), inner(:,:)
        integer  :: m, n, i, j

        write(*, '(a)') "(ii) bed bump (300 m, 80 km domain): surface vs full"
        dxs = [2000.0_dp, 1000.0_dp, 500.0_dp]
        do m = 1, 3
            dx = dxs(m)
            n = nint(80000.0_dp / dx)
            allocate(z_b(n, n), tb(n, n), Hs(n, n), Hh(n, n), Hf(n, n), mask(n, n), inner(n, n))
            do j = 1, n
                do i = 1, n
                    x = (i - 0.5_dp) * dx - 40000.0_dp
                    y = (j - 0.5_dp) * dx - 40000.0_dp
                    z_b(i, j)   = 300.0_dp * exp(-(x * x + y * y) / (2.0_dp * 2.0e4_dp**2))
                    mask(i, j)  = hypot(x, y) <= 30000.0_dp
                    inner(i, j) = hypot(x, y) <= 24000.0_dp
                end do
            end do
            tb = 1.0e5_dp
            call solve_mode("surface", n, n, dx, z_b, tb, mask, Hs)
            call solve_mode("full",    n, n, dx, z_b, tb, mask, Hh)
            call solve_mode("flat",    n, n, dx, z_b, tb, mask, Hf)
            dsh(m) = maxval(abs(Hs - Hh), mask=inner)
            rsh(m) = sqrt(sum((Hs - Hh)**2, mask=mask) / count(mask))
            dfh(m) = maxval(abs(Hf - Hh), mask=inner)
            write(*, '(a,f6.0,a,f7.2,a,f7.2,a,f7.2,a,f7.1)') "        dx ", dx, ": max|surface - full| (r < 24 km) ", &
                dsh(m), " m, rms (all ice) ", rsh(m), " m; max|flat - full| ", dfh(m), " m; max H ", maxval(Hh)
            deallocate(z_b, tb, Hs, Hh, Hf, mask, inner)
        end do
        call check(dsh(3) < dsh(2) .and. dsh(2) < dsh(1), "surface - full shrinks with dx")
        call check(dsh(3) < 0.2_dp * dfh(3), "surface - full << flat - full (bed effect)")
    end subroutine test_bump

end program test_surface
