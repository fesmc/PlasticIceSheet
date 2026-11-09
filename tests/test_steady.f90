! Steady SMB-driven solver (ctrl%method = "steady"): the plastic solution forced by its
! own implied SMB is recovered from a cold start; the Vialov steady profile (uniform
! accumulation, no sliding) against its analytic solution; built-in vs Lis; the same
! result for 1 and 4 OpenMP threads; a 150 x 150 cap at 2 km with timing.
program test_steady

    !$ use omp_lib
    use plastic_ice
    use plastic_test_utils

    implicit none

    real(dp), parameter :: rho_g = 917.0_dp * 9.80665_dp

    write(*, '(a)') "test_steady"
    call test_fixed_point(40, 2000.0_dp, "bump 40 x 40, dx 2 km")
    call test_fixed_point(150, 2000.0_dp, "cap 150 x 150, dx 2 km")
    call test_vialov()
    call test_margin_outflow()
    call finish("test_steady")

contains

    ! Synthetic cap of n x n cells: circular mask, 300 m bedrock high, tau_b = 1e5 Pa;
    ! H_p from mode "surface"; smb = its implied SMB (plastic_smb_flux_divergence).
    subroutine test_fixed_point(n, dx, label)
        integer,          intent(in) :: n
        real(dp),         intent(in) :: dx
        character(len=*), intent(in) :: label
        type(plastic_t)          :: p0, ps, pl
        type(plastic_par_t)      :: par
        type(plastic_glen_t)     :: r
        type(plastic_sliding_t)  :: sl
        type(plastic_smb_ctrl_t) :: ctrl
        real(dp), allocatable :: z_b(:,:), z_sl(:,:), tb(:,:), smb(:,:), Hp(:,:), H1(:,:)
        logical,  allocatable :: mask(:,:)
        real(dp) :: x, y, L, t0, t1, dmax
        integer  :: i, j, st, nt

        write(*, '(a,a)') "plastic fixed point under its implied SMB: ", label
        allocate(z_b(n, n), z_sl(n, n), tb(n, n), smb(n, n), Hp(n, n), H1(n, n), mask(n, n))
        L = n * dx
        do j = 1, n
            do i = 1, n
                x = (i - 0.5_dp) * dx - 0.5_dp * L
                y = (j - 0.5_dp) * dx - 0.5_dp * L
                z_b(i, j)  = 300.0_dp * exp(-(x * x + y * y) / (2.0_dp * (0.25_dp * L)**2))
                mask(i, j) = hypot(x, y) <= 0.375_dp * L
            end do
        end do
        z_sl = -1.0e4_dp; tb = 1.0e5_dp
        r%A = 1.0e-16_dp; r%n = 3.0_dp
        sl%law = "weertman"; sl%beta = 1.0e4_dp; sl%m = 3.0_dp
        par%mode = "surface"
        call plastic_init(p0, n, n, dx, dx, par)
        call plastic_solve(p0, z_b, z_sl, tb, mask)
        Hp = p0%H
        call plastic_smb_flux_divergence(p0, z_b, z_sl, mask, r, sl, smb, stat=st)
        call check(st == PLASTIC_OK, "implied SMB computed")

        ! started at H_p: already converged
        call plastic_init(ps, n, n, dx, dx, par)
        ctrl%lin = "builtin"
        call plastic_solve_smb(ps, z_b, z_sl, smb, mask, r, sl, ctrl=ctrl, H_init=Hp, stat=st)
        call check(st == PLASTIC_OK .and. ps%n_outer_used == 0 .and. maxval(abs(ps%H - Hp)) == 0.0_dp, &
                   "H_p is a root (no step from H_p)")

        ! cold start (surface solve at tau_init)
        t0 = wall_time()
        call plastic_solve_smb(ps, z_b, z_sl, smb, mask, r, sl, ctrl=ctrl, stat=st)
        t1 = wall_time() - t0
        dmax = maxval(abs(ps%H - Hp), mask=mask)
        nt = 1
        !$ nt = omp_get_max_threads()
        write(*, '(a,i0,a,es9.2,a,f8.3,a,i0,a,es9.2,a)') "        cold start: ", ps%n_outer_used, &
            " steps, residual ", ps%outer_change, " m/yr, ", t1, " s (", nt, " threads), max|H - H_p| ", dmax, " m"
        call check(st == PLASTIC_OK, "cold start converged")
        call check(dmax < 1.0e-3_dp * maxval(Hp), "plastic solution recovered (max|dH| < 1e-3 max H)")
        H1 = ps%H

        if (plastic_has_lis()) then
            call plastic_init(pl, n, n, dx, dx, par)
            ctrl%lin = "lis"
            t0 = wall_time()
            call plastic_solve_smb(pl, z_b, z_sl, smb, mask, r, sl, ctrl=ctrl, stat=st)
            t1 = wall_time() - t0
            write(*, '(a,i0,a,f8.3,a,es9.2,a)') "        Lis: ", pl%n_outer_used, " steps, ", t1, &
                " s, max|H_lis - H_builtin| ", maxval(abs(pl%H - H1)), " m"
            call check(st == PLASTIC_OK .and. maxval(abs(pl%H - H1)) < 1.0e-3_dp * maxval(Hp), &
                       "Lis and built-in agree")
            call plastic_end(pl)
            ctrl%lin = "builtin"
        end if

        ! the built-in solver gives the same bits for 1 and 4 threads
        !$ if (n <= 40) then
        !$     call omp_set_num_threads(1)
        !$     call plastic_solve_smb(ps, z_b, z_sl, smb, mask, r, sl, ctrl=ctrl, stat=st)
        !$     H1 = ps%H
        !$     call omp_set_num_threads(4)
        !$     call plastic_solve_smb(ps, z_b, z_sl, smb, mask, r, sl, ctrl=ctrl, stat=st)
        !$     call check(maxval(abs(ps%H - H1)) == 0.0_dp, "built-in: bit-identical for 1 and 4 threads")
        !$     call omp_set_num_threads(nt)
        !$ end if
        call plastic_end(p0); call plastic_end(ps)
    end subroutine test_fixed_point

    ! Margin outflow and surface velocity of plastic_smb_flux_divergence on a 40 x 40 cap: the
    ! outflow through the margin faces sums to the divergence over the mask (interior faces
    ! cancel); with the sea above the bed west of the cap the outflow there is marine, east of it
    ! land; the surface speed is positive on the ice and points down the surface.
    subroutine test_margin_outflow()
        integer,  parameter :: n = 40
        real(dp), parameter :: dx = 2000.0_dp
        type(plastic_t)         :: ps
        type(plastic_par_t)     :: par
        type(plastic_glen_t)    :: r
        type(plastic_sliding_t) :: sl
        real(dp), allocatable :: z_b(:,:), z_sl(:,:), tb(:,:), smb(:,:), ql(:,:), qm(:,:), ux(:,:), uy(:,:)
        logical,  allocatable :: mask(:,:)
        real(dp) :: x, y, L, sdiv, sout
        integer  :: i, j, st
        logical  :: west_marine, east_land, down

        write(*, '(a)') "margin outflow and surface velocity of the implied SMB"
        allocate(z_b(n, n), z_sl(n, n), tb(n, n), smb(n, n), ql(n, n), qm(n, n), ux(n, n), uy(n, n), mask(n, n))
        L = n * dx
        do j = 1, n
            do i = 1, n
                x = (i - 0.5_dp) * dx - 0.5_dp * L
                y = (j - 0.5_dp) * dx - 0.5_dp * L
                z_b(i, j)  = 300.0_dp * exp(-(x * x + y * y) / (2.0_dp * (0.25_dp * L)**2))
                mask(i, j) = hypot(x, y) <= 0.375_dp * L
                z_sl(i, j) = -1.0e4_dp
                if (x < 0.0_dp) z_sl(i, j) = z_b(i, j) + 50.0_dp     ! the sea west of the divide
            end do
        end do
        tb = 1.0e5_dp
        r%A = 1.0e-16_dp; r%n = 3.0_dp
        sl%law = "weertman"; sl%beta = 1.0e4_dp; sl%m = 3.0_dp
        par%mode = "surface"
        call plastic_init(ps, n, n, dx, dx, par)
        call plastic_solve(ps, z_b, z_sl, tb, mask)
        call plastic_smb_flux_divergence(ps, z_b, z_sl, mask, r, sl, smb, stat=st, q_land=ql, q_marine=qm, us_x=ux, us_y=uy)
        call check(st == PLASTIC_OK, "implied SMB with the margin outflow computed")
        sdiv = sum(smb, mask = mask); sout = sum(ql + qm, mask = mask)
        call check_close(sout, sdiv, 1.0e-10_dp, 0.0_dp, "margin outflow sums to the divergence")
        call check(all(ql >= 0.0_dp) .and. all(qm >= 0.0_dp), "outflow non-negative")
        west_marine = .true.; east_land = .true.; down = .true.
        do j = 1, n
            do i = 1, n
                x = (i - 0.5_dp) * dx - 0.5_dp * L
                if (x < -dx .and. ql(i, j) > 0.0_dp) west_marine = .false.
                if (x > dx .and. qm(i, j) > 0.0_dp) east_land = .false.
                if (mask(i, j) .and. ps%H(i, j) > 0.0_dp) then
                    if (.not. hypot(ux(i, j), uy(i, j)) > 0.0_dp) down = .false.
                end if
            end do
        end do
        call check(west_marine .and. sum(qm) > 0.0_dp, "west of the divide the outflow is marine")
        call check(east_land .and. sum(ql) > 0.0_dp, "east of it the outflow is to land")
        call check(down .and. ux(n / 2 + 8, n / 2) > 0.0_dp .and. ux(n / 2 - 8, n / 2) < 0.0_dp, &
                   "surface speed positive on the ice, away from the divide")
        call plastic_end(ps)
    end subroutine test_margin_outflow

    ! Vialov profile: divide at x = 0 (closed array edge), margin at x = L (cell face),
    ! smb = a, no sliding (linear law, beta = 1e30): H = [2 (a/G)^(1/n) (L^p - x^p)]^(n/(2n+2)),
    ! p = (n+1)/n, G = 2A (rho g)^n/(n+2).
    subroutine test_vialov()
        real(dp), parameter :: L = 300000.0_dp, a = 0.3_dp
        type(plastic_t)          :: ps
        type(plastic_glen_t)     :: r
        type(plastic_sliding_t)  :: sl
        type(plastic_smb_ctrl_t) :: ctrl
        real(dp) :: dxs(4), ei(4), ev(4), dx, x, ex, G, pw, vex, vh, t0
        real(dp), allocatable :: z_b(:,:), z_sl(:,:), smb(:,:)
        logical,  allocatable :: mask(:,:)
        integer  :: m, nx, i, st

        write(*, '(a)') "Vialov steady profile (n = 3, a = 0.3 m/yr, L = 300 km)"
        r%A = 1.0e-16_dp; r%n = 3.0_dp
        sl%law = "linear"; sl%beta = 1.0e30_dp
        G  = 2.0_dp * r%A * rho_g**r%n / (r%n + 2.0_dp)
        pw = (r%n + 1.0_dp) / r%n
        ctrl%lin = "builtin"
        dxs = [20000.0_dp, 10000.0_dp, 5000.0_dp, 2500.0_dp]
        do m = 1, 4
            dx = dxs(m)
            nx = nint(L / dx) + 3
            allocate(z_b(nx, 3), z_sl(nx, 3), smb(nx, 3), mask(nx, 3))
            z_b = 0.0_dp; z_sl = -1.0e4_dp; smb = a
            do i = 1, nx
                mask(i, :) = (i - 0.5_dp) * dx < L
            end do
            call plastic_init(ps, nx, 3, dx, dx)
            t0 = wall_time()
            call plastic_solve_smb(ps, z_b, z_sl, smb, mask, r, sl, ctrl=ctrl, stat=st)
            ei(m) = 0.0_dp; vex = 0.0_dp; vh = 0.0_dp
            do i = 1, nx
                if (.not. mask(i, 2)) cycle
                x  = (i - 0.5_dp) * dx
                ex = (2.0_dp * (a / G)**(1.0_dp / r%n) * (L**pw - x**pw))**(r%n / (2.0_dp * r%n + 2.0_dp))
                if (x <= 0.9_dp * L) ei(m) = max(ei(m), abs(ps%H(i, 2) - ex) / ex)
                vex = vex + ex
                vh  = vh + ps%H(i, 2)
            end do
            ev(m) = vh / vex - 1.0_dp
            write(*, '(a,f7.0,a,i0,a,es9.2,a,es10.3,a,f7.4,a,f6.2,a)') "        dx ", dx, ": ", ps%n_outer_used, &
                " steps, residual ", ps%outer_change, ", max rel err (x <= 0.9 L) ", ei(m), &
                ", volume rel err ", ev(m), ", ", wall_time() - t0, " s"
            call check(st == PLASTIC_OK, "Vialov solve converged")
            deallocate(z_b, z_sl, smb, mask)
            call plastic_end(ps)
        end do
        write(*, '(a,3f6.2)') "        order (max rel err, x <= 0.9 L): ", &
            (log(ei(m) / ei(m + 1)) / log(2.0_dp), m = 1, 3)
        call check(ei(4) < ei(3) .and. ei(3) < ei(2) .and. ei(2) < ei(1), "Vialov error decreases with dx")
        call check(ei(4) < 0.02_dp, "Vialov within 2 % at dx = 2.5 km")
    end subroutine test_vialov

end program test_steady
