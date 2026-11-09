! Closure tests: ports of PlasticIceSheet.jl test/runtests.jl (constitutive laws, DIVA
! velocity and implied SMB, balance flux/velocity, SMB-driven reconstruction), plus
! the warm start of the SMB-driven reconstruction.
program test_closures

    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use plastic_ice
    use plastic_test_utils

    implicit none

    write(*, '(a)') "test_closures"
    call test_laws()
    call test_velocity_closure_b()
    call test_balance_closure_a()
    call test_reconstruct()
    call finish("test_closures")

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

    ! Flat-mode plastic solve of a uniform tau_b on a flat bed (geometry for the tests),
    ! margin at the ice-free cell centre as in PlasticIceSheet.jl.
    subroutine flat_cap(ps, nx, ny, dx, tau, mask, max_sweeps)
        type(plastic_t), intent(inout) :: ps
        integer,         intent(in)    :: nx, ny, max_sweeps
        real(dp),        intent(in)    :: dx, tau
        logical,         intent(in)    :: mask(nx, ny)
        type(plastic_par_t) :: par
        real(dp) :: zero(nx, ny), tb(nx, ny)
        par%mode = "flat"; par%max_sweeps = max_sweeps; par%tol = 1.0e-9_dp
        par%margin = "centre"           ! the PlasticIceSheet.jl tests (exact symmetry at the divide)
        zero = 0.0_dp; tb = tau
        call plastic_init(ps, nx, ny, dx, dx, par)
        call plastic_solve(ps, zero, zero, tb, mask)
    end subroutine flat_cap

    subroutine test_laws()
        type(plastic_glen_t)    :: r
        type(plastic_sliding_t) :: lin, wm, wm1, rc
        real(dp) :: tau, ub, small

        write(*, '(a)') "constitutive laws: deformation and sliding inversions"
        r%A = 1.0e-16_dp; r%n = 3.0_dp
        call check_close(plastic_deformational_velocity(r, 8.0e4_dp, 2000.0_dp), &
                         (2.0_dp * 1.0e-16_dp / 5.0_dp) * (8.0e4_dp)**3 * 2000.0_dp, 1.0e-14_dp, 0.0_dp, &
                         "Glen: 2A/(n+2) tau^n H")
        call check(plastic_deformational_velocity(r, 0.0_dp, 2000.0_dp) == 0.0_dp, "Glen: no stress, no flow")

        tau = 8.0e4_dp
        lin%law = "linear"; lin%beta = 2.0e3_dp
        call check_close(lin%beta * plastic_basal_velocity(lin, tau), tau, 1.0e-14_dp, 0.0_dp, "linear inverts")
        wm%law = "weertman"; wm%beta = 1.0e4_dp; wm%m = 3.0_dp
        call check_close(wm%beta * plastic_basal_velocity(wm, tau)**(1.0_dp / wm%m), tau, 1.0e-12_dp, 0.0_dp, &
                         "weertman inverts")
        wm1%law = "weertman"; wm1%beta = 2.0e3_dp; wm1%m = 1.0_dp
        call check_close(plastic_basal_velocity(wm1, tau), plastic_basal_velocity(lin, tau), 1.0e-14_dp, 0.0_dp, &
                         "weertman m = 1 == linear")
        rc%law = "coulomb"; rc%C = 1.0e5_dp; rc%u_0 = 100.0_dp; rc%m = 3.0_dp
        ub = plastic_basal_velocity(rc, tau)
        call check_close(rc%C * (ub / (ub + rc%u_0))**(1.0_dp / rc%m), tau, 1.0e-12_dp, 0.0_dp, &
                         "coulomb forward relation holds")
        call check(ieee_is_nan(plastic_basal_velocity(rc, rc%C)), "coulomb: NaN at the cap")
        call check(ieee_is_nan(plastic_basal_velocity(rc, 1.2e5_dp)), "coulomb: NaN above the cap")
        small = 1.0e3_dp
        call check_close(plastic_basal_velocity(rc, small), rc%u_0 * (small / rc%C)**rc%m, 1.0e-2_dp, 0.0_dp, &
                         "coulomb low-speed limit is the power law")
    end subroutine test_laws

    subroutine test_velocity_closure_b()
        integer, parameter :: nx = 61, ny = 61
        real(dp), parameter :: dx = 2000.0_dp, tau = 8.0e4_dp
        type(plastic_t)         :: ps
        type(plastic_glen_t)    :: r
        type(plastic_sliding_t) :: sl
        real(dp) :: tb(nx, ny), ux(nx, ny), uy(nx, ny), speed(nx, ny), smb(nx, ny), Q(nx, ny), out(nx, ny)
        real(dp) :: err
        logical  :: mask(nx, ny)

        write(*, '(a)') "DIVA velocity and implied SMB (closure B)"
        call circle_mask(mask, 31.0_dp, 31.0_dp, 24.0_dp)
        call flat_cap(ps, nx, ny, dx, tau, mask, 600)
        tb = tau
        r%A = 1.0e-16_dp; r%n = 3.0_dp
        sl%law = "weertman"; sl%beta = 1.0e4_dp; sl%m = 3.0_dp
        call plastic_velocity(ps, tb, mask, r, sl, ux, uy, speed)

        call check_close(speed(31, 31), plastic_basal_velocity(sl, tau) + &
                         plastic_deformational_velocity(r, tau, ps%H(31, 31)), 1.0e-14_dp, 0.0_dp, &
                         "speed = u_b + u_def at tau = tau_b")
        write(*, '(a,es10.3)') "        |u| at the divide ", hypot(ux(31, 31), uy(31, 31))
        call check(hypot(ux(31, 31), uy(31, 31)) == 0.0_dp, "zero horizontal velocity at the divide")
        call check(ux(41, 31) > 0.0_dp .and. uy(31, 41) > 0.0_dp, "flow is outward")
        call check(all(speed == 0.0_dp .or. mask), "ice-free cells carry no velocity")

        ! MFD divergence is the exact inverse of balance_flux.
        call plastic_smb_implied(ps, ux, uy, mask, smb, speed=speed)
        call plastic_balance_flux(ps%z_srf, mask, smb, dx, dx, Q)
        out = speed * ps%H * sqrt(dx * dx)
        err = maxval(abs(Q - out), mask=mask) / maxval(out)
        write(*, '(a,es10.3)') "        round trip rel err (speed) ", err
        call check(err < 1.0e-10_dp, "balance_flux(smb_implied(speed)) == speed H w")

        ! Without speed, |u| = hypot(ux, uy) and the round trip holds for that.
        call plastic_smb_implied(ps, ux, uy, mask, smb)
        call plastic_balance_flux(ps%z_srf, mask, smb, dx, dx, Q)
        out = hypot(ux, uy) * ps%H * sqrt(dx * dx)
        err = maxval(abs(Q - out), mask=mask) / maxval(out)
        write(*, '(a,es10.3)') "        round trip rel err (hypot) ", err
        call check(err < 1.0e-10_dp, "balance_flux(smb_implied(ux, uy)) == |u| H w")
        call plastic_end(ps)
    end subroutine test_velocity_closure_b

    subroutine test_balance_closure_a()
        integer, parameter :: nx = 41, ny = 6, nx2 = 61, ny2 = 61
        real(dp), parameter :: dx = 1000.0_dp, adot = 0.3_dp, dx2 = 2000.0_dp, tau2 = 8.0e4_dp
        type(plastic_t)         :: ps, p2
        type(plastic_glen_t)    :: r
        type(plastic_sliding_t) :: sl
        real(dp) :: smb(nx, ny), Q(nx, ny), cellarea, err, rr
        real(dp) :: tb2(nx2, ny2), ux(nx2, ny2), uy(nx2, ny2), speed(nx2, ny2), smbB(nx2, ny2)
        real(dp) :: bsp(nx2, ny2), bux(nx2, ny2), buy(nx2, ny2), u_basal(nx2, ny2), ub
        logical  :: mask(nx, ny), mask2(nx2, ny2), sel(nx2, ny2)
        integer  :: i, j

        write(*, '(a)') "balance flux/velocity and consistency (closure A)"
        ! 1-D strip, uniform accumulation: exact catchment integration.
        mask = .true.
        mask(1, :) = .false.; mask(nx, :) = .false.
        call flat_cap(ps, nx, ny, dx, 1.0e5_dp, mask, 600)
        smb = adot
        call plastic_balance_flux(ps%z_srf, mask, smb, dx, dx, Q)
        cellarea = dx * dx
        j = 3
        err = 0.0_dp
        do i = 22, nx - 1
            err = max(err, abs(Q(i, j) / (adot * cellarea) - (1.5_dp + (i - 22))))
        end do
        call check(err < 1.0e-10_dp, "right flank reads 1.5, 2.5, ..., 19.5")
        err = 0.0_dp
        do i = 23, nx - 1
            err = max(err, abs((Q(i, j) - Q(i - 1, j)) - adot * cellarea) / (adot * cellarea))
        end do
        call check(err < 1.0e-10_dp, "unit increment per cell")
        call check_close(sum(Q(2, :)) + sum(Q(nx - 1, :)), adot * cellarea * count(mask), 1.0e-12_dp, 0.0_dp, &
                         "mass conservation: all accumulation exits at the margins")
        call plastic_end(ps)

        ! A <-> B on a radial cap.
        call circle_mask(mask2, 31.0_dp, 31.0_dp, 24.0_dp)
        call flat_cap(p2, nx2, ny2, dx2, tau2, mask2, 800)
        tb2 = tau2
        r%A = 1.0e-16_dp; r%n = 3.0_dp
        sl%law = "weertman"; sl%beta = 1.0e4_dp; sl%m = 3.0_dp
        call plastic_velocity(p2, tb2, mask2, r, sl, ux, uy, speed)
        call plastic_smb_implied(p2, ux, uy, mask2, smbB, speed=speed)
        call plastic_balance_velocity(p2%z_srf, p2%H, mask2, smbB, dx2, dx2, bsp, bux, buy)
        do j = 1, ny2
            do i = 1, nx2
                rr = hypot(i - 31.0_dp, j - 31.0_dp)
                sel(i, j) = mask2(i, j) .and. rr > 6.0_dp .and. rr < 20.0_dp
            end do
        end do
        err = maxval(abs(bsp - speed), mask=sel) / maxval(speed, mask=sel)
        write(*, '(a,es10.3)') "        balance vs DIVA speed rel err ", err
        call check(err < 1.0e-9_dp, "balance velocity of B's SMB recovers B's speed")

        call plastic_implied_basal_velocity(p2%z_srf, p2%H, tb2, mask2, smbB, dx2, dx2, r, u_basal)
        ub = plastic_basal_velocity(sl, tau2)
        call check(all(u_basal > 0.0_dp .or. .not. sel), "implied basal velocity > 0")
        call check(maxval(abs(u_basal - ub), mask=sel) < 1.0e-6_dp, "implied basal velocity == u_b")
        call plastic_end(p2)
    end subroutine test_balance_closure_a

    subroutine test_reconstruct()
        integer, parameter :: nx = 41, ny = 41
        real(dp), parameter :: dx = 2500.0_dp, tau = 8.0e4_dp
        type(plastic_t)          :: p0, pr, pw, pc
        type(plastic_glen_t)     :: r
        type(plastic_sliding_t)  :: sl
        type(plastic_smb_ctrl_t) :: ctrl
        type(plastic_par_t)      :: pcen
        real(dp) :: z_b(nx, ny), z_sl(nx, ny), tb(nx, ny), ux(nx, ny), uy(nx, ny), speed(nx, ny)
        real(dp) :: smb(nx, ny), H0(nx, ny), Hg(nx, ny), smb2(nx, ny), t0, tw, tc
        real(dp), allocatable :: rel(:), srt(:)
        logical  :: mask(nx, ny), sel(nx, ny)
        integer  :: i, j, st

        write(*, '(a)') "SMB-driven reconstruction"
        call circle_mask(mask, 21.0_dp, 21.0_dp, 16.0_dp)
        call flat_cap(p0, nx, ny, dx, tau, mask, 600)
        pcen%margin = "centre"          ! as the reference geometry (flat_cap)
        H0 = p0%H
        z_b = 0.0_dp; z_sl = 0.0_dp; tb = tau
        r%A = 1.0e-16_dp; r%n = 3.0_dp
        sl%law = "weertman"; sl%beta = 1.0e4_dp; sl%m = 3.0_dp
        call plastic_velocity(p0, tb, mask, r, sl, ux, uy, speed)
        call plastic_smb_implied(p0, ux, uy, mask, smb, speed=speed)
        do j = 1, ny
            do i = 1, nx
                sel(i, j) = mask(i, j) .and. hypot(i - 21.0_dp, j - 21.0_dp) < 12.0_dp
            end do
        end do

        ! The plastic solution is a fixed point of the SMB-driven map.
        call plastic_init(pr, nx, ny, dx, dx, pcen)
        call plastic_solve_smb(pr, z_b, z_sl, smb, mask, r, sl, H_init=H0, stat=st)
        write(*, '(a,es10.3,a,i0)') "        fixed point max|H - H0| ", maxval(abs(pr%H - H0), mask=mask), &
            " m, n_outer ", pr%n_outer_used
        call check(st == PLASTIC_OK, "fixed point: converged")
        call check(maxval(abs(pr%H - H0), mask=mask) < 1.0e-6_dp, "plastic solution is a fixed point (inner surface)")
        ctrl%mode = "flat"
        call plastic_solve_smb(pr, z_b, z_sl, smb, mask, r, sl, ctrl=ctrl, H_init=H0, stat=st)
        call check(st == PLASTIC_OK .and. maxval(abs(pr%H - H0), mask=mask) < 1.0e-6_dp, &
                   "plastic solution is a fixed point (inner flat, Julia parity)")
        ctrl%mode = "steep"
        call plastic_solve_smb(pr, z_b, z_sl, smb, mask, r, sl, ctrl=ctrl, H_init=H0, stat=st)
        call check(st == PLASTIC_BAD_INPUT, "unknown ctrl%mode -> PLASTIC_BAD_INPUT")
        ctrl = plastic_smb_ctrl_t()

        ! Recovered from a flat 400 m start by the under-relaxed iteration.
        Hg = 0.0_dp
        where (mask) Hg = 400.0_dp
        ctrl%relax = 0.05_dp; ctrl%n_outer = 800
        call plastic_solve_smb(pr, z_b, z_sl, smb, mask, r, sl, ctrl=ctrl, H_init=Hg)
        rel = pack(abs(pr%H - H0) / max(H0, 1.0e-30_dp), sel)
        srt = sort_vec(rel)
        write(*, '(a,es10.3,a,i0)') "        cold start: lower-median rel err ", srt(size(srt) / 2), &
            ", n_outer ", pr%n_outer_used
        call check(srt(size(srt) / 2) < 0.02_dp, "recovered from a flat start (median rel err < 2 %)")
        call check_close(maxval(pr%H), maxval(H0), 0.05_dp, 0.0_dp, "max H recovered within 5 %")

        ! Warm start: after a 5 % change of smb, warm and cold reach the same geometry.
        ctrl = plastic_smb_ctrl_t()
        ctrl%tol = 1.0e-2_dp
        call plastic_init(pw, nx, ny, dx, dx, pcen)
        call plastic_init(pc, nx, ny, dx, dx, pcen)
        call plastic_solve_smb(pw, z_b, z_sl, smb, mask, r, sl, ctrl=ctrl)
        smb2 = 1.05_dp * smb
        t0 = wall_time()
        call plastic_solve_smb(pw, z_b, z_sl, smb2, mask, r, sl, warm=.true., ctrl=ctrl)
        tw = wall_time() - t0
        t0 = wall_time()
        call plastic_solve_smb(pc, z_b, z_sl, smb2, mask, r, sl, warm=.false., ctrl=ctrl)
        tc = wall_time() - t0
        write(*, '(a,i0,a,i0,a,f6.1,a,es10.3,a)') "        n_outer warm ", pw%n_outer_used, ", cold ", &
            pc%n_outer_used, ", speedup ", tc / max(tw, 1.0e-9_dp), "x, max|H_warm - H_cold| ", &
            maxval(abs(pw%H - pc%H)), " m"
        call check(pw%warm_used .and. pw%converged .and. pc%converged, "warm and cold converged")
        call check(pw%n_outer_used < pc%n_outer_used, "warm start needs fewer iterations")
        ! Both stop at max|dH| < tol with relax 0.05: they agree to ~ tol/relax.
        call check(maxval(abs(pw%H - pc%H)) < ctrl%tol / ctrl%relax, "warm == cold within tol/relax")
        call plastic_end(p0); call plastic_end(pr); call plastic_end(pw); call plastic_end(pc)
    end subroutine test_reconstruct

end program test_closures
