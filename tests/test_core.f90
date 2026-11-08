! Core solver tests: ports of PlasticIceSheet.jl test/runtests.jl (godunov kernel,
! Nye profile, flat vs full, full over bed relief with its residual), plus the marine
! margin with a 2-D sea level, warm == cold start, and the status reporting.
program test_core

    use plastic_ice
    use plastic_eikonal,    only: godunov_eikonal
    use plastic_test_utils

    implicit none

    real(dp), parameter :: rho_g = 917.0_dp * 9.80665_dp

    write(*, '(a)') "test_core"
    call test_godunov_kernel()
    call test_nye_strip()
    call test_flat_vs_full()
    call test_full_bed_relief()
    call test_marine_margin()
    call test_warm_start()
    call test_status()
    call finish("test_core")

contains

    subroutine circle_mask(mask, cx, cy, R)
        logical,  intent(out) :: mask(:,:)
        real(dp), intent(in)  :: cx, cy, R
        integer :: i, j
        do j = 1, size(mask, 2)
            do i = 1, size(mask, 1)
                mask(i, j) = hypot(real(i, dp) - cx, real(j, dp) - cy) <= R
            end do
        end do
    end subroutine circle_mask

    subroutine test_godunov_kernel()
        real(dp), parameter :: big = 1.0e300_dp
        write(*, '(a)') "godunov kernel"
        call check_close(godunov_eikonal(0.0_dp, 1.0_dp, big, 1.0_dp, 2.0_dp), 2.0_dp, 1.0e-14_dp, 0.0_dp, &
                         "one-sided x")
        call check_close(godunov_eikonal(big, 1.0_dp, 3.0_dp, 2.0_dp, 1.5_dp), 3.0_dp + 1.5_dp * 2.0_dp, &
                         1.0e-14_dp, 0.0_dp, "one-sided y")
        call check_close(godunov_eikonal(0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp), 2.0_dp / sqrt(2.0_dp), &
                         1.0e-14_dp, 0.0_dp, "symmetric two-sided")
        call check_close(godunov_eikonal(1.0_dp, 1.0_dp, 4.0_dp, 1.0_dp, 0.0_dp), 1.0_dp, 1.0e-14_dp, 0.0_dp, &
                         "F = 0 gives the min neighbour")
    end subroutine test_godunov_kernel

    ! 1-D strip between land margins: H(d) = sqrt(2 tau_b/(rho g) d).
    subroutine test_nye_strip()
        integer, parameter :: nx = 60, ny = 8
        real(dp), parameter :: dx = 1000.0_dp, tau = 1.0e5_dp
        type(plastic_t)     :: ps
        type(plastic_par_t) :: par
        real(dp) :: z_b(nx, ny), z_sl(nx, ny), tb(nx, ny), c, d, err
        logical  :: mask(nx, ny)
        integer  :: i, j

        write(*, '(a)') "1-D strip vs analytic Nye profile"
        z_b = 0.0_dp; z_sl = 0.0_dp; tb = tau
        mask = .true.
        mask(1, :) = .false.; mask(nx, :) = .false.
        par%mode = "flat"; par%max_sweeps = 400; par%tol = 1.0e-9_dp
        call plastic_init(ps, nx, ny, dx, dx, par)
        call plastic_solve(ps, z_b, z_sl, tb, mask)
        call check(ps%converged, "flat solve converged")

        c = 2.0_dp * tau / rho_g
        j = ny / 2
        err = 0.0_dp
        do i = 3, nx - 2
            d = dx * min(i - 1, nx - i)
            err = max(err, abs(ps%H(i, j) - sqrt(c * d)) / sqrt(c * d))
        end do
        write(*, '(a,es10.3)') "        max rel err ", err
        call check(err < 1.0e-3_dp, "H = sqrt(c d), rtol 1e-3")
        call check(all(ps%H(1, :) == 0.0_dp) .and. all(ps%H(nx, :) == 0.0_dp), "margins carry H = 0")
        call plastic_end(ps)
    end subroutine test_nye_strip

    ! On a constant bed the flat and full formulations coincide.
    subroutine test_flat_vs_full()
        integer, parameter :: nx = 50, ny = 50
        real(dp), parameter :: dx = 2000.0_dp
        type(plastic_t)     :: pf, ph
        type(plastic_par_t) :: par
        real(dp) :: z_b(nx, ny), z_sl(nx, ny), tb(nx, ny), dmax
        logical  :: mask(nx, ny)

        write(*, '(a)') "flat and full agree on a flat bed"
        z_b = -50.0_dp; z_sl = 0.0_dp; tb = 0.8e5_dp
        call circle_mask(mask, 25.5_dp, 25.5_dp, 18.0_dp)
        par%max_sweeps = 400; par%tol = 1.0e-9_dp
        par%mode = "flat"
        call plastic_init(pf, nx, ny, dx, dx, par)
        call plastic_solve(pf, z_b, z_sl, tb, mask)
        par%mode = "full"; par%n_outer = 40; par%outer_tol = 1.0e-5_dp
        call plastic_init(ph, nx, ny, dx, dx, par)
        call plastic_solve(ph, z_b, z_sl, tb, mask)
        dmax = maxval(abs(pf%H - ph%H))
        write(*, '(a,es10.3,a,i0)') "        max|H_flat - H_full| ", dmax, " m, n_outer ", ph%n_outer_used
        call check(dmax < 2.0_dp, "max|H_flat - H_full| < 2 m")
        call check(ph%converged, "full solve converged")
        call plastic_end(pf); call plastic_end(ph)
    end subroutine test_flat_vs_full

    ! A central bedrock high lifts the surface; the damped fixed point is reached.
    subroutine test_full_bed_relief()
        integer, parameter :: nx = 40, ny = 40
        real(dp), parameter :: dx = 2000.0_dp
        type(plastic_t)     :: p0, pb
        type(plastic_par_t) :: par
        real(dp) :: z_flat(nx, ny), z_bump(nx, ny), z_sl(nx, ny), tb(nx, ny), r(nx, ny)
        logical  :: mask(nx, ny)
        integer  :: i, j

        write(*, '(a)') "full responds to bed relief (damped fixed point converges)"
        call circle_mask(mask, 20.5_dp, 20.5_dp, 15.0_dp)
        z_flat = 0.0_dp; z_sl = 0.0_dp; tb = 1.0e5_dp
        do j = 1, ny
            do i = 1, nx
                z_bump(i, j) = 300.0_dp * exp(-((i - 20.5_dp)**2 + (j - 20.5_dp)**2) / 200.0_dp)
            end do
        end do
        par%mode = "full"; par%n_outer = 600; par%outer_tol = 1.0e-7_dp; par%relax = 0.4_dp
        call plastic_init(p0, nx, ny, dx, dx, par)
        call plastic_init(pb, nx, ny, dx, dx, par)
        call plastic_solve(p0, z_flat, z_sl, tb, mask)
        call plastic_solve(pb, z_bump, z_sl, tb, mask)
        write(*, '(a,i0)') "        n_outer ", pb%n_outer_used
        call check(pb%converged, "full solve over the bump converged")
        call check(pb%z_srf(20, 20) > p0%H(20, 20), "bedrock high lifts the surface")
        call plastic_residual(pb, z_bump, z_sl, tb, mask, r)
        write(*, '(a,es10.3)') "        max|residual| ", maxval(abs(r))
        call check(maxval(abs(r)) < 1.0e-2_dp, "Godunov residual of the :full solution < 1e-2")
        call plastic_end(p0); call plastic_end(pb)
    end subroutine test_full_bed_relief

    ! Marine margin: bed below sea level, margin value H_flt; with a 2-D sea level the
    ! two margins carry different flotation thickness. Flat bed: w = H_flt^2 + c d.
    subroutine test_marine_margin()
        integer, parameter :: nx = 60, ny = 8
        real(dp), parameter :: dx = 1000.0_dp, tau = 1.0e5_dp
        type(plastic_t)     :: pf, ph
        type(plastic_par_t) :: par
        real(dp) :: z_b(nx, ny), z_sl(nx, ny), tb(nx, ny), c, d, hl, hr, err, ex
        logical  :: mask(nx, ny)
        integer  :: i, j

        write(*, '(a)') "marine margin, 2-D sea level"
        z_b = -500.0_dp; tb = tau
        z_sl(1:nx/2, :) = 0.0_dp
        z_sl(nx/2+1:nx, :) = -200.0_dp
        mask = .true.
        mask(1, :) = .false.; mask(nx, :) = .false.
        hl = plastic_flotation_thickness(-500.0_dp, 0.0_dp, par)
        hr = plastic_flotation_thickness(-500.0_dp, -200.0_dp, par)
        call check_close(hl, 1025.0_dp / 917.0_dp * 500.0_dp, 1.0e-14_dp, 0.0_dp, "flotation thickness")
        call check(plastic_flotation_thickness(10.0_dp, 0.0_dp, par) == 0.0_dp, "no flotation on land")

        par%mode = "flat"; par%max_sweeps = 400; par%tol = 1.0e-9_dp
        call plastic_init(pf, nx, ny, dx, dx, par)
        call plastic_solve(pf, z_b, z_sl, tb, mask)
        c = 2.0_dp * tau / rho_g
        j = ny / 2
        err = 0.0_dp
        do i = 2, 6
            d = dx * (i - 1)
            ex = sqrt(hl**2 + c * d)
            err = max(err, abs(pf%H(i, j) - ex) / ex)
            ex = sqrt(hr**2 + c * d)
            err = max(err, abs(pf%H(nx - i + 1, j) - ex) / ex)
        end do
        write(*, '(a,es10.3)') "        max rel err near the margins ", err
        call check(err < 1.0e-6_dp, "H = sqrt(H_flt^2 + c d) at both marine margins")
        call check(abs(pf%H(1, j)) == 0.0_dp, "ocean cell carries H = 0")

        ! Constant bed: full == flat (bed gradient vanishes).
        par%mode = "full"; par%outer_tol = 1.0e-6_dp
        call plastic_init(ph, nx, ny, dx, dx, par)
        call plastic_solve(ph, z_b, z_sl, tb, mask)
        write(*, '(a,es10.3)') "        max|H_full - H_flat| ", maxval(abs(ph%H - pf%H))
        call check(maxval(abs(ph%H - pf%H)) < 1.0e-6_dp, "full == flat on a constant marine bed")
        call plastic_end(pf); call plastic_end(ph)
    end subroutine test_marine_margin

    ! Warm start. Both warm and cold solves stop on an outer step < outer_tol, so each
    ! sits ~ rho/(1-rho) outer_tol from the fixed point (rho ~ 0.7 is the linear rate of
    ! the damped map): at the default tolerance they agree to a few outer_tol, and with
    ! a tight tolerance they converge to the same fixed point. Mask growth takes the
    ! flat solve on the new cells; a large mask change falls back to a cold start.
    subroutine test_warm_start()
        integer, parameter :: nx = 80, ny = 80
        real(dp), parameter :: dx = 2000.0_dp
        type(plastic_t)     :: pw, pc
        type(plastic_par_t) :: par
        real(dp) :: z_b(nx, ny), z_b1(nx, ny), z_b2(nx, ny), z_sl(nx, ny), tb(nx, ny), tb2(nx, ny)
        real(dp) :: H_prev(nx, ny)
        logical  :: mask(nx, ny), mask2(nx, ny), mask3(nx, ny)
        integer  :: i, j

        write(*, '(a)') "warm start == cold start"
        call circle_mask(mask, 40.5_dp, 40.5_dp, 30.0_dp)
        z_sl = 0.0_dp
        do j = 1, ny
            do i = 1, nx
                z_b(i, j)  = 300.0_dp * exp(-((i - 35.0_dp)**2 + (j - 45.0_dp)**2) / 300.0_dp) &
                             + 50.0_dp * sin(i / 7.0_dp) * cos(j / 9.0_dp) - 100.0_dp
                tb(i, j)   = 1.0e5_dp * (1.0_dp + 0.3_dp * sin(i / 11.0_dp) * cos(j / 13.0_dp))
                ! one coupling step: a 1 m GIA-like change of the bed
                z_b1(i, j) = z_b(i, j) - 1.0_dp * exp(-((i - 40.0_dp)**2 + (j - 40.0_dp)**2) / 400.0_dp)
                ! a large change: 20 m of bed and 2 % of tau_b
                z_b2(i, j) = z_b(i, j) - 20.0_dp * exp(-((i - 40.0_dp)**2 + (j - 40.0_dp)**2) / 400.0_dp)
                tb2(i, j)  = 1.02_dp * tb(i, j)
            end do
        end do
        par%mode = "full"; par%outer_tol = 1.0e-4_dp; par%n_outer = 500; par%relax = 0.4_dp
        call plastic_init(pw, nx, ny, dx, dx, par)
        call plastic_init(pc, nx, ny, dx, dx, par)

        call plastic_solve(pw, z_b, z_sl, tb, mask)
        call check(pw%converged, "initial cold solve converged")

        ! unchanged inputs: the stored solution is already the fixed point
        H_prev = pw%H
        call plastic_solve(pw, z_b, z_sl, tb, mask, warm=.true.)
        call check(pw%warm_used .and. pw%n_outer_used == 1, "unchanged inputs: one outer iteration")
        call check(maxval(abs(pw%H - H_prev)) < par%outer_tol, "unchanged inputs: H unchanged")

        call warm_vs_cold(pw, pc, z_sl, "1 m bed change", z_b1, tb, mask, 5.0_dp * par%outer_tol)
        call warm_vs_cold(pw, pc, z_sl, "20 m bed + 2 % tau_b change", z_b2, tb2, mask, 5.0_dp * par%outer_tol)

        ! tight tolerance: the same fixed point
        pw%par%outer_tol = 1.0e-7_dp; pc%par%outer_tol = 1.0e-7_dp
        call plastic_solve(pw, z_b, z_sl, tb, mask, warm=.true.)
        call warm_vs_cold(pw, pc, z_sl, "tight outer_tol = 1e-7", z_b2, tb2, mask, 1.0e-6_dp)
        pw%par%outer_tol = par%outer_tol; pc%par%outer_tol = par%outer_tol

        ! grown mask (margin moved by one cell): the new ring takes the flat solve
        call circle_mask(mask2, 40.5_dp, 40.5_dp, 31.0_dp)
        call warm_vs_cold(pw, pc, z_sl, "grown mask", z_b2, tb2, mask2, 5.0_dp * par%outer_tol, fewer=.false.)
        call check(pw%warm_used .and. pw%n_new_cells > 0, "warm start used with new cells")

        ! a large mask change falls back to a cold start
        call circle_mask(mask3, 40.5_dp, 40.5_dp, 38.0_dp)
        call plastic_solve(pw, z_b2, z_sl, tb2, mask3, warm=.true.)
        call plastic_solve(pc, z_b2, z_sl, tb2, mask3, warm=.false.)
        call check(.not. pw%warm_used, "large mask change falls back to cold")
        call check(maxval(abs(pw%H - pc%H)) == 0.0_dp, "fallback is bit-identical to cold")

        ! reset drops the warm start
        call plastic_reset(pw)
        call plastic_solve(pw, z_b2, z_sl, tb2, mask3, warm=.true.)
        call check(.not. pw%warm_used, "plastic_reset drops the warm start")

        call plastic_end(pw); call plastic_end(pc)


    end subroutine test_warm_start

    ! Warm solve on pw (holding the previous solution) and cold solve on pc; compare.
    ! fewer = .false.: only require the warm solve to be no slower (a moved margin
    ! changes H everywhere, so the stored solution is no better than the flat one).
    subroutine warm_vs_cold(pw, pc, z_sl, label, zb, tbx, msk, tol_eq, fewer)
        type(plastic_t),   intent(inout) :: pw, pc
        real(dp),          intent(in)    :: z_sl(:,:), zb(:,:), tbx(:,:), tol_eq
        character(len=*),  intent(in)    :: label
        logical,           intent(in)    :: msk(:,:)
        logical, optional, intent(in)    :: fewer
        real(dp) :: t0, tw, tc, dmax
        logical  :: need_fewer

        t0 = wall_time()
        call plastic_solve(pw, zb, z_sl, tbx, msk, warm=.true.)
        tw = wall_time() - t0
        t0 = wall_time()
        call plastic_solve(pc, zb, z_sl, tbx, msk, warm=.false.)
        tc = wall_time() - t0
        dmax = maxval(abs(pw%H - pc%H))
        write(*, '(a,a,a,i0,a,i0,a,i0,a,f5.1,a,es9.2,a,es9.2)') "        ", label, &
            ": new cells ", pw%n_new_cells, ", n_outer warm ", pw%n_outer_used, &
            ", cold ", pc%n_outer_used, ", speedup ", tc / max(tw, 1.0e-9_dp), &
            "x, max|H_warm - H_cold| ", dmax, " m, outer_tol ", pw%par%outer_tol
        call check(pw%warm_used .and. pw%converged .and. pc%converged, label//": warm and cold converged")
        need_fewer = .true.
        if (present(fewer)) need_fewer = fewer
        if (need_fewer) then
            call check(pw%n_outer_used < pc%n_outer_used, label//": warm needs fewer outer iterations")
        else
            call check(pw%n_outer_used <= pc%n_outer_used, label//": warm needs no more outer iterations")
        end if
        call check(dmax < tol_eq, label//": warm == cold")
    end subroutine warm_vs_cold

    ! Status reporting: bad input, unknown mode, unreachable ice; no abort.
    subroutine test_status()
        type(plastic_t)     :: ps
        type(plastic_par_t) :: par
        real(dp) :: a(10, 10), b(9, 10)
        logical  :: mask(10, 10)
        integer  :: st

        write(*, '(a)') "status reporting"
        call plastic_init(ps, 10, 10, 1000.0_dp, 1000.0_dp, par)
        a = 0.0_dp; b = 0.0_dp
        mask = .true.
        call plastic_solve(ps, b, a, a, mask, stat=st)
        call check(st == PLASTIC_BAD_INPUT .and. .not. ps%converged, "shape mismatch -> PLASTIC_BAD_INPUT")
        ps%par%mode = "steep"
        call plastic_solve(ps, a, a, a + 1.0e5_dp, mask, stat=st)
        call check(st == PLASTIC_BAD_INPUT, "unknown mode -> PLASTIC_BAD_INPUT")
        ps%par%mode = "full"
        call plastic_solve(ps, a, a, a + 1.0e5_dp, mask, stat=st)
        write(*, '(a,a)') "        msg: ", trim(ps%msg)
        call check(st == PLASTIC_NOT_CONVERGED .and. all(ps%H == 0.0_dp), &
                   "ice with no margin -> PLASTIC_NOT_CONVERGED, H = 0")
        call plastic_init(ps, 0, 10, 1000.0_dp, 1000.0_dp, par, stat=st)
        call check(st == PLASTIC_BAD_INPUT, "plastic_init with nx = 0 -> PLASTIC_BAD_INPUT")
        call plastic_end(ps)
    end subroutine test_status

end program test_core
