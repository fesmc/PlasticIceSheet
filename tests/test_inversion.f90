! Joint inversion of speed and mass balance (plastic_smb_inversion) on a synthetic dome, marine west of the divide and
! land to the east: its limits (sigma_u -> 0: the flow law's speed; sigma_u large: the prior's mass balance on the
! marine half), mass closure, the cold penalty and the smoothness of the speed correction. A uniform positive prior
! cannot hold on the land half (no ablation zone: the inflow of a land margin cell is ablated there), so the prior
! limit is checked on the marine half.
program test_inversion

    use plastic_ice
    use plastic_test_utils
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    implicit none

    integer,  parameter :: n = 48
    real(dp), parameter :: dx = 1.0e4_dp
    type(plastic_t)          :: ps
    type(plastic_par_t)      :: par
    type(plastic_glen_t)     :: r
    type(plastic_sliding_t)  :: sl
    type(plastic_inv_ctrl_t) :: c
    type(plastic_inv_info_t) :: inf
    real(dp), allocatable :: z_b(:,:), z_sl(:,:), tb(:,:), ap(:,:), dh(:,:), smb(:,:), sp(:,:), ul(:,:), ql(:,:), qm(:,:)
    real(dp), allocatable :: smb0(:,:)
    logical,  allocatable :: mask(:,:), inner(:,:), west(:,:)
    real(dp) :: x, y, L, rough0, rough1, f0, f1
    integer  :: i, j, st

    write(*, '(a)') "test_inversion"
    allocate(z_b(n, n), z_sl(n, n), tb(n, n), ap(n, n), dh(n, n), smb(n, n), sp(n, n), ul(n, n), ql(n, n), qm(n, n), &
             smb0(n, n), mask(n, n), inner(n, n), west(n, n))
    L = n * dx
    do j = 1, n
        do i = 1, n
            x = (i - 0.5_dp) * dx - 0.5_dp * L
            y = (j - 0.5_dp) * dx - 0.5_dp * L
            z_b(i, j)  = 300.0_dp * exp(-(x * x + y * y) / (2.0_dp * (0.25_dp * L)**2))
            mask(i, j) = hypot(x, y) <= 0.375_dp * L
            z_sl(i, j) = -1.0e4_dp
            if (x < 0.0_dp) z_sl(i, j) = z_b(i, j) + 50.0_dp
            west(i, j) = x < 0.0_dp
        end do
    end do
    inner = .false.
    do j = 2, n - 1
        do i = 2, n - 1
            inner(i, j) = mask(i, j) .and. mask(i - 1, j) .and. mask(i + 1, j) .and. mask(i, j - 1) .and. mask(i, j + 1)
        end do
    end do
    tb = 1.0e5_dp
    r%A = 1.0e-18_dp; r%n = 3.0_dp                         ! deformation ~1 m/yr: the prior is reachable
    sl%law = "weertman"; sl%beta = 3.0e4_dp; sl%m = 3.0_dp ! sliding ~40 m/yr
    par%mode = "surface"
    call plastic_init(ps, n, n, dx, dx, par)
    call plastic_solve(ps, z_b, z_sl, tb, mask)
    ap = 0.3_dp;  dh = 0.0_dp

    ! the flow-law limit
    c%sigma_u = 1.0e-6_dp
    call plastic_smb_inversion(ps, z_b, z_sl, mask, r, sl, ap, dh, smb0, sp, c, tau_cap=tb, u_law=ul, info=inf, stat=st)
    call check(st == PLASTIC_OK .and. all(ieee_is_finite(smb0)), "sigma_u small: converged, finite")
    call check(maxval(abs(log(sp / ul)), mask = mask) < 1.0e-3_dp, "sigma_u small: the speed of the flow law")

    ! the prior limit, with mass closure
    c%sigma_u = 30.0_dp;  c%cold = .false.
    call plastic_smb_inversion(ps, z_b, z_sl, mask, r, sl, ap, dh, smb, sp, c, tau_cap=tb, q_land=ql, q_marine=qm, &
                               info=inf, stat=st)
    write(*, '(a,i0,a,es10.3,a,es10.3,a,f7.4)') "        prior limit: ", inf%iters, " steps, J ", inf%J, " from ", &
        inf%J_start, ", rms(a - prior) inner marine half ", sqrt(sum((smb - ap)**2, mask = inner .and. west) / &
        count(inner .and. west))
    call check(st == PLASTIC_OK, "sigma_u large: converged")
    call check(sqrt(sum((smb - ap)**2, mask = inner .and. west) / count(inner .and. west)) < 0.1_dp * 0.3_dp, &
               "sigma_u large: the prior's mass balance on the inner cells of the marine half")
    call check_close(sum(smb, mask = mask), sum(qm), 1.0e-8_dp, 1.0e-10_dp, "sum(a) = the marine outflow (dH/dt = 0)")
    call check(all(sp > 0.0_dp .or. .not. mask), "speed positive on the ice")

    ! the cold penalty: fewer inner cells with a < 0 than the flow law alone
    c%sigma_u = 0.7_dp;  c%cold = .true.
    call plastic_smb_inversion(ps, z_b, z_sl, mask, r, sl, ap, dh, smb, sp, c, tau_cap=tb, info=inf, stat=st)
    f0 = real(count(inner .and. smb0 < 0.0_dp), dp) / count(inner)
    f1 = real(count(inner .and. smb < 0.0_dp), dp) / count(inner)
    write(*, '(a,f6.3,a,f6.3,a,i0)') "        inner cells with a < 0: flow law ", f0, ", inversion ", f1, &
        "; cold cells < 0: ", inf%n_cold_neg
    call check(st == PLASTIC_OK .and. f1 <= f0 .and. inf%J < inf%J_start, "cold: fewer negative cells, lower J")

    ! the smoothness of the speed correction
    c%smooth = 0.0_dp
    call plastic_smb_inversion(ps, z_b, z_sl, mask, r, sl, ap, dh, smb, sp, c, tau_cap=tb, u_law=ul, stat=st)
    rough0 = roughness(log(sp / ul))
    c%smooth = 5.0_dp * dx
    call plastic_smb_inversion(ps, z_b, z_sl, mask, r, sl, ap, dh, smb, sp, c, tau_cap=tb, u_law=ul, stat=st)
    rough1 = roughness(log(sp / ul))
    write(*, '(a,es10.3,a,es10.3)') "        roughness of ln U - ln U_law: smooth 0 ", rough0, ", 5 dx ", rough1
    call check(st == PLASTIC_OK .and. rough1 < 0.5_dp * rough0, "smooth = 5 dx: a smoother speed correction")

    call plastic_end(ps)
    call finish("test_inversion")

contains

    real(dp) function roughness(e)
        ! mean squared difference between inner-cell x-neighbours
        real(dp), intent(in) :: e(:,:)
        integer :: ii, jj, m
        roughness = 0.0_dp;  m = 0
        do jj = 1, n
            do ii = 1, n - 1
                if (inner(ii, jj) .and. inner(ii + 1, jj)) then
                    roughness = roughness + (e(ii + 1, jj) - e(ii, jj))**2;  m = m + 1
                end if
            end do
        end do
        roughness = roughness / max(m, 1)
    end function roughness

end program test_inversion
