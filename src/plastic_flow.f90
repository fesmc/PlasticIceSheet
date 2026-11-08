! Closures on a fixed plastic geometry: Glen deformation, basal sliding laws, the
! local DIVA velocity, and the multiple-flow-direction (MFD) flux operator shared by
! the implied SMB (closure B) and the balance flux / velocity (closure A).
!
! The plastic state is the membrane-free limit of DIVA (tau_d = tau_b), so the
! velocity is local: it follows surface descent with speed u_b(tau_b) + u_def(tau_b,H).
! Units follow the rate factor and sliding coefficients (per year -> m/yr).
!
! balance_flux (source -> throughput) and flux_divergence (throughput -> source) use
! the same MFD weights and are exact inverses, so routing the implied SMB back down
! the surface recovers u H to machine precision.
module plastic_flow

    use plastic_kinds,   only: dp, NAN_DP
    use plastic_eikonal, only: grad_central
    use plastic_solver,  only: plastic_t, PLASTIC_OK, PLASTIC_BAD_INPUT

    implicit none

    private
    public :: plastic_glen_t, plastic_sliding_t
    public :: plastic_deformational_velocity, plastic_basal_velocity
    public :: plastic_flow_direction, plastic_velocity, plastic_smb_implied
    public :: plastic_balance_flux, plastic_balance_velocity, plastic_implied_basal_velocity
    ! for plastic_reconstruct
    public :: LAW_LINEAR, LAW_WEERTMAN, LAW_COULOMB
    public :: sliding_law_code, basal_velocity_code, deformational_velocity, balance_flux_core

    ! Isothermal Glen flow law: depth-averaged deformational speed 2A/(n+2) tau^n H.
    type plastic_glen_t
        real(dp) :: A = 1.0e-16_dp      ! rate factor [Pa^-n yr^-1]
        real(dp) :: n = 3.0_dp          ! flow-law exponent
    end type plastic_glen_t

    ! Basal sliding law, inverted for the sliding speed at the drag tau_b:
    !   "linear"   tau_b = beta u_b                 -> u_b = tau_b / beta
    !   "weertman" tau_b = beta u_b^(1/m)           -> u_b = (tau_b/beta)^m
    !   "coulomb"  tau_b = C (u_b/(u_b+u_0))^(1/m)  -> u_b = u_0 x/(1-x), x = (tau_b/C)^m,
    !              NaN where tau_b >= C (bed at its Coulomb cap, speed undetermined).
    type plastic_sliding_t
        character(len=8) :: law = "weertman"
        real(dp) :: beta = 1.0e4_dp     ! linear [Pa yr m-1], weertman [Pa (yr m-1)^(1/m)]
        real(dp) :: m    = 3.0_dp       ! weertman / coulomb stress exponent
        real(dp) :: C    = 1.0e5_dp     ! coulomb cap [Pa]
        real(dp) :: u_0  = 100.0_dp     ! coulomb threshold speed [m yr-1]
    end type plastic_sliding_t

    integer, parameter :: LAW_LINEAR = 1, LAW_WEERTMAN = 2, LAW_COULOMB = 3

contains

    ! --- constitutive laws ---------------------------------------------------------

    ! x^p, with an integer power where p is integral (exact, and cheap in loops).
    pure function powr(x, p) result(y)
        real(dp), intent(in) :: x, p
        real(dp) :: y
        if (p == aint(p) .and. abs(p) <= 16.0_dp) then
            y = x**nint(p)
        else
            y = x**p
        end if
    end function powr

    pure function deformational_velocity(rheo, tau, H) result(u)
        type(plastic_glen_t), intent(in) :: rheo
        real(dp),             intent(in) :: tau, H
        real(dp) :: u
        u = (2.0_dp * rheo%A / (rheo%n + 2.0_dp)) * powr(tau, rheo%n) * H
    end function deformational_velocity

    ! Depth-averaged internal-deformation speed 2A/(n+2) tau^n H.
    elemental function plastic_deformational_velocity(rheo, tau, H) result(u)
        type(plastic_glen_t), intent(in) :: rheo
        real(dp),             intent(in) :: tau, H
        real(dp) :: u
        u = deformational_velocity(rheo, tau, H)
    end function plastic_deformational_velocity

    ! LAW_* code of a sliding law, 0 if unknown.
    pure function sliding_law_code(slide) result(code)
        type(plastic_sliding_t), intent(in) :: slide
        integer :: code
        select case (trim(slide%law))
        case ("linear")
            code = LAW_LINEAR
        case ("weertman")
            code = LAW_WEERTMAN
        case ("coulomb")
            code = LAW_COULOMB
        case default
            code = 0
        end select
    end function sliding_law_code

    elemental function basal_velocity_code(code, slide, tau_b) result(u)
        integer,                 intent(in) :: code
        type(plastic_sliding_t), intent(in) :: slide
        real(dp),                intent(in) :: tau_b
        real(dp) :: u, x
        select case (code)
        case (LAW_LINEAR)
            u = tau_b / slide%beta
        case (LAW_WEERTMAN)
            u = powr(tau_b / slide%beta, slide%m)
        case (LAW_COULOMB)
            if (tau_b < slide%C) then
                x = powr(tau_b / slide%C, slide%m)
                u = slide%u_0 * x / (1.0_dp - x)
            else
                u = NAN_DP
            end if
        case default
            u = NAN_DP
        end select
    end function basal_velocity_code

    ! Sliding speed u_b >= 0 that produces basal drag tau_b (NaN for "coulomb" at or
    ! above its cap, and for an unknown law).
    elemental function plastic_basal_velocity(slide, tau_b) result(u)
        type(plastic_sliding_t), intent(in) :: slide
        real(dp),                intent(in) :: tau_b
        real(dp) :: u
        u = basal_velocity_code(sliding_law_code(slide), slide, tau_b)
    end function plastic_basal_velocity

    ! --- flow direction and DIVA velocity ------------------------------------------

    ! Unit flow direction (ex, ey) = -grad z_srf / sqrt(|grad z_srf|^2 + slope_eps^2):
    ! tapers to zero at divides.
    subroutine flow_direction(nx, ny, z_srf, dx, dy, slope_eps, ex, ey)
        integer,  intent(in)  :: nx, ny
        real(dp), intent(in)  :: z_srf(nx, ny), dx, dy, slope_eps
        real(dp), intent(out) :: ex(nx, ny), ey(nx, ny)
        integer  :: i, j
        real(dp) :: gm

        call grad_central(nx, ny, z_srf, dx, dy, ex, ey)
        do j = 1, ny
            do i = 1, nx
                gm = sqrt(ex(i, j)**2 + ey(i, j)**2 + slope_eps**2)
                ex(i, j) = -ex(i, j) / gm
                ey(i, j) = -ey(i, j) / gm
            end do
        end do
    end subroutine flow_direction

    subroutine plastic_flow_direction(z_srf, dx, dy, ex, ey, slope_eps)
        real(dp),           intent(in)  :: z_srf(:,:), dx, dy
        real(dp),           intent(out) :: ex(:,:), ey(:,:)
        real(dp), optional, intent(in)  :: slope_eps
        real(dp) :: eps
        eps = 1.0e-6_dp
        if (present(slope_eps)) eps = slope_eps
        call flow_direction(size(z_srf, 1), size(z_srf, 2), z_srf, dx, dy, eps, ex, ey)
    end subroutine plastic_flow_direction

    ! Local DIVA depth-averaged velocity on the geometry of ps (ps%z_srf, ps%H): speed
    ! u_b(tau_b) + u_def(tau_b, H) along the surface descent; zero off the mask.
    ! u_b is NaN where a "coulomb" bed is at its cap (take the balance velocity there).
    subroutine plastic_velocity(ps, tau_b, mask, rheo, slide, ux, uy, speed, slope_eps, &
                                u_b, u_def, stat)
        type(plastic_t),         intent(inout) :: ps
        real(dp),                intent(in)    :: tau_b(:,:)
        logical,                 intent(in)    :: mask(:,:)
        type(plastic_glen_t),    intent(in)    :: rheo
        type(plastic_sliding_t), intent(in)    :: slide
        real(dp),                intent(out)   :: ux(:,:), uy(:,:)
        real(dp), optional,      intent(out)   :: speed(:,:), u_b(:,:), u_def(:,:)
        real(dp), optional,      intent(in)    :: slope_eps
        integer,  optional,      intent(out)   :: stat
        integer  :: i, j, code, st
        real(dp) :: eps, ud, ub, s

        st = PLASTIC_OK
        code = sliding_law_code(slide)
        if (.not. allocated(ps%H) .or. code == 0) st = PLASTIC_BAD_INPUT
        if (st == PLASTIC_OK) then
            if (.not. (same_shape(ps, tau_b) .and. same_shape_l(ps, mask) .and. &
                       same_shape(ps, ux) .and. same_shape(ps, uy))) st = PLASTIC_BAD_INPUT
            if (present(speed)) then
                if (.not. same_shape(ps, speed)) st = PLASTIC_BAD_INPUT
            end if
            if (present(u_b)) then
                if (.not. same_shape(ps, u_b)) st = PLASTIC_BAD_INPUT
            end if
            if (present(u_def)) then
                if (.not. same_shape(ps, u_def)) st = PLASTIC_BAD_INPUT
            end if
        end if
        if (present(stat)) stat = st
        if (st /= PLASTIC_OK) return

        eps = 1.0e-6_dp
        if (present(slope_eps)) eps = slope_eps
        call flow_direction(ps%nx, ps%ny, ps%z_srf, ps%dx, ps%dy, eps, ps%wk%ex, ps%wk%ey)

        do j = 1, ps%ny
            do i = 1, ps%nx
                if (mask(i, j)) then
                    ud = deformational_velocity(rheo, tau_b(i, j), ps%H(i, j))
                    ub = basal_velocity_code(code, slide, tau_b(i, j))
                    s  = ub + ud
                    ux(i, j) = s * ps%wk%ex(i, j)
                    uy(i, j) = s * ps%wk%ey(i, j)
                else
                    ud = 0.0_dp
                    ub = 0.0_dp
                    s  = 0.0_dp
                    ux(i, j) = 0.0_dp
                    uy(i, j) = 0.0_dp
                end if
                if (present(speed)) speed(i, j) = s
                if (present(u_b))   u_b(i, j)   = ub
                if (present(u_def)) u_def(i, j) = ud
            end do
        end do
    end subroutine plastic_velocity

    ! --- shared MFD operator -------------------------------------------------------

    ! Downhill surface drops from (i,j) to its four neighbours (i-1,j), (i+1,j),
    ! (i,j-1), (i,j+1); off-grid, ice-free and uphill neighbours give zero.
    pure subroutine downhill_drops(nx, ny, z_s, mask, i, j, ni, nj, d, wsum)
        integer,  intent(in)  :: nx, ny, i, j
        real(dp), intent(in)  :: z_s(nx, ny)
        logical,  intent(in)  :: mask(nx, ny)
        integer,  intent(out) :: ni(4), nj(4)
        real(dp), intent(out) :: d(4), wsum
        integer  :: n
        real(dp) :: dd

        ni = [i - 1, i + 1, i, i]
        nj = [j, j, j - 1, j + 1]
        do n = 1, 4
            d(n) = 0.0_dp
            if (ni(n) >= 1 .and. ni(n) <= nx .and. nj(n) >= 1 .and. nj(n) <= ny) then
                if (mask(ni(n), nj(n))) then
                    dd = z_s(i, j) - z_s(ni(n), nj(n))
                    if (dd > 0.0_dp) d(n) = dd
                end if
            end if
        end do
        wsum = ((d(1) + d(2)) + d(3)) + d(4)
    end subroutine downhill_drops

    ! Net volume-flux divergence S = out - inflow of the per-cell throughput out
    ! (the adjoint of balance_flux). S = smb * cell area.
    subroutine flux_divergence(nx, ny, out, z_s, mask, inflow, S)
        integer,  intent(in)  :: nx, ny
        real(dp), intent(in)  :: out(nx, ny), z_s(nx, ny)
        logical,  intent(in)  :: mask(nx, ny)
        real(dp), intent(out) :: inflow(nx, ny), S(nx, ny)
        integer  :: i, j, n, ni(4), nj(4)
        real(dp) :: d(4), wsum, oc

        inflow = 0.0_dp
        do j = 1, ny
            do i = 1, nx
                if (.not. mask(i, j)) cycle
                call downhill_drops(nx, ny, z_s, mask, i, j, ni, nj, d, wsum)
                if (.not. (wsum > 0.0_dp)) cycle
                oc = out(i, j)
                do n = 1, 4
                    if (d(n) > 0.0_dp) inflow(ni(n), nj(n)) = inflow(ni(n), nj(n)) + oc * d(n) / wsum
                end do
            end do
        end do
        do j = 1, ny
            do i = 1, nx
                if (mask(i, j)) then
                    S(i, j) = out(i, j) - inflow(i, j)
                else
                    S(i, j) = 0.0_dp
                end if
            end do
        end do
    end subroutine flux_divergence

    ! Sort key(1:n) descending, permuting idx alongside; stable bottom-up merge sort
    ! with the work arrays tk, ti.
    subroutine sort_desc(n, key, idx, tk, ti)
        integer,  intent(in)    :: n
        real(dp), intent(inout) :: key(n), tk(n)
        integer,  intent(inout) :: idx(n), ti(n)
        integer :: width, lo, mid, hi, a, b, k

        width = 1
        do while (width < n)
            lo = 1
            do while (lo <= n)
                mid = min(lo + width - 1, n)
                hi  = min(lo + 2*width - 1, n)
                a = lo
                b = mid + 1
                k = lo
                do while (a <= mid .and. b <= hi)
                    if (key(a) >= key(b)) then
                        tk(k) = key(a); ti(k) = idx(a); a = a + 1
                    else
                        tk(k) = key(b); ti(k) = idx(b); b = b + 1
                    end if
                    k = k + 1
                end do
                do while (a <= mid)
                    tk(k) = key(a); ti(k) = idx(a); a = a + 1; k = k + 1
                end do
                do while (b <= hi)
                    tk(k) = key(b); ti(k) = idx(b); b = b + 1; k = k + 1
                end do
                lo = lo + 2*width
            end do
            key(1:n) = tk(1:n)
            idx(1:n) = ti(1:n)
            width = 2*width
        end do
    end subroutine sort_desc

    ! Volume flux Q [m3/time] through each ice cell of the steady balance div q = smb:
    ! smb is rolled downhill on z_s with the MFD weights, highest cell first (Budd &
    ! Warner 1996); flux reaching the margin terminates there. idx, key, tk, ti are
    ! work arrays of length >= count(mask).
    subroutine balance_flux_core(nx, ny, z_s, mask, smb, dx, dy, Q, idx, key, ti, tk)
        integer,  intent(in)    :: nx, ny
        real(dp), intent(in)    :: z_s(nx, ny), smb(nx, ny), dx, dy
        logical,  intent(in)    :: mask(nx, ny)
        real(dp), intent(out)   :: Q(nx, ny)
        integer,  intent(inout) :: idx(:), ti(:)
        real(dp), intent(inout) :: key(:), tk(:)
        integer  :: i, j, k, n, nn, lin, ni(4), nj(4)
        real(dp) :: cellarea, d(4), wsum, Qk

        cellarea = dx * dy
        n = 0
        do j = 1, ny
            do i = 1, nx
                if (mask(i, j)) then
                    Q(i, j) = smb(i, j) * cellarea
                    n = n + 1
                    idx(n) = (j - 1) * nx + i
                    key(n) = z_s(i, j)
                else
                    Q(i, j) = 0.0_dp
                end if
            end do
        end do

        call sort_desc(n, key, idx, tk, ti)

        do k = 1, n
            lin = idx(k)
            i = mod(lin - 1, nx) + 1
            j = (lin - 1) / nx + 1
            call downhill_drops(nx, ny, z_s, mask, i, j, ni, nj, d, wsum)
            if (.not. (wsum > 0.0_dp)) cycle      ! pit or margin: throughput terminates
            Qk = Q(i, j)
            do nn = 1, 4
                if (d(nn) > 0.0_dp) Q(ni(nn), nj(nn)) = Q(ni(nn), nj(nn)) + Qk * d(nn) / wsum
            end do
        end do
    end subroutine balance_flux_core

    ! --- closure B: implied SMB -------------------------------------------------------

    ! Steady-state SMB implied by a velocity on the geometry of ps, smb = div(u H)
    ! [m ice / yr], as the MFD divergence of the throughput |u| H sqrt(dx dy) down
    ! ps%z_srf. |u| = speed if given (the signed speed of plastic_velocity: this makes
    ! plastic_balance_flux recover speed*H*sqrt(dx dy) exactly), else hypot(ux, uy).
    subroutine plastic_smb_implied(ps, ux, uy, mask, smb, speed, stat)
        type(plastic_t),    intent(inout) :: ps
        real(dp),           intent(in)    :: ux(:,:), uy(:,:)
        logical,            intent(in)    :: mask(:,:)
        real(dp),           intent(out)   :: smb(:,:)
        real(dp), optional, intent(in)    :: speed(:,:)
        integer,  optional, intent(out)   :: stat
        integer  :: i, j, st
        real(dp) :: wf

        st = PLASTIC_OK
        if (.not. allocated(ps%H)) then
            st = PLASTIC_BAD_INPUT
        else if (.not. (same_shape(ps, ux) .and. same_shape(ps, uy) .and. &
                        same_shape_l(ps, mask) .and. same_shape(ps, smb))) then
            st = PLASTIC_BAD_INPUT
        else if (present(speed)) then
            if (.not. same_shape(ps, speed)) st = PLASTIC_BAD_INPUT
        end if
        if (present(stat)) stat = st
        if (st /= PLASTIC_OK) return

        wf = sqrt(ps%dx * ps%dy)
        do j = 1, ps%ny
            do i = 1, ps%nx
                if (present(speed)) then
                    ps%wk%f1(i, j) = speed(i, j) * ps%H(i, j) * wf
                else
                    ps%wk%f1(i, j) = hypot(ux(i, j), uy(i, j)) * ps%H(i, j) * wf
                end if
            end do
        end do
        call flux_divergence(ps%nx, ps%ny, ps%wk%f1, ps%z_srf, mask, ps%wk%f2, smb)
        smb = smb / (ps%dx * ps%dy)
    end subroutine plastic_smb_implied

    ! --- closure A: balance flux and velocity -------------------------------------

    ! Balance volume flux Q [m3/yr] of smb [m/yr] routed down z_srf (see
    ! balance_flux_core). Allocates its sort work arrays (not for a hot loop).
    subroutine plastic_balance_flux(z_srf, mask, smb, dx, dy, Q, stat)
        real(dp),          intent(in)  :: z_srf(:,:), smb(:,:), dx, dy
        logical,           intent(in)  :: mask(:,:)
        real(dp),          intent(out) :: Q(:,:)
        integer, optional, intent(out) :: stat
        integer,  allocatable :: idx(:), ti(:)
        real(dp), allocatable :: key(:), tk(:)
        integer :: nx, ny, n

        nx = size(z_srf, 1)
        ny = size(z_srf, 2)
        if (any(shape(mask) /= [nx, ny]) .or. any(shape(smb) /= [nx, ny]) .or. &
            any(shape(Q) /= [nx, ny])) then
            if (present(stat)) stat = PLASTIC_BAD_INPUT
            return
        end if
        n = max(count(mask), 1)
        allocate(idx(n), ti(n), key(n), tk(n))
        call balance_flux_core(nx, ny, z_srf, mask, smb, dx, dy, Q, idx, key, ti, tk)
        if (present(stat)) stat = PLASTIC_OK
    end subroutine plastic_balance_flux

    ! Balance velocity: speed = Q / (sqrt(dx dy) max(H, H_min)) along the flow
    ! direction. Net-ablation catchments can give a negative speed (flux reversal).
    subroutine plastic_balance_velocity(z_srf, H, mask, smb, dx, dy, speed, ux, uy, Q, &
                                        slope_eps, H_min, stat)
        real(dp),           intent(in)  :: z_srf(:,:), H(:,:), smb(:,:), dx, dy
        logical,            intent(in)  :: mask(:,:)
        real(dp),           intent(out) :: speed(:,:), ux(:,:), uy(:,:)
        real(dp), optional, intent(out) :: Q(:,:)
        real(dp), optional, intent(in)  :: slope_eps, H_min
        integer,  optional, intent(out) :: stat
        real(dp), allocatable :: Qw(:,:), ex(:,:), ey(:,:)
        integer  :: nx, ny, i, j, st
        real(dp) :: eps, hmin, wf

        nx = size(z_srf, 1)
        ny = size(z_srf, 2)
        st = PLASTIC_OK
        if (any(shape(H) /= [nx, ny]) .or. any(shape(speed) /= [nx, ny]) .or. &
            any(shape(ux) /= [nx, ny]) .or. any(shape(uy) /= [nx, ny])) st = PLASTIC_BAD_INPUT
        if (present(Q)) then
            if (any(shape(Q) /= [nx, ny])) st = PLASTIC_BAD_INPUT
        end if
        allocate(Qw(nx, ny), ex(nx, ny), ey(nx, ny))
        if (st == PLASTIC_OK) call plastic_balance_flux(z_srf, mask, smb, dx, dy, Qw, st)
        if (present(stat)) stat = st
        if (st /= PLASTIC_OK) return

        eps  = 1.0e-6_dp
        hmin = 1.0_dp
        if (present(slope_eps)) eps  = slope_eps
        if (present(H_min))     hmin = H_min
        call flow_direction(nx, ny, z_srf, dx, dy, eps, ex, ey)
        wf = sqrt(dx * dy)
        do j = 1, ny
            do i = 1, nx
                if (mask(i, j)) then
                    speed(i, j) = Qw(i, j) / (wf * max(H(i, j), hmin))
                else
                    speed(i, j) = 0.0_dp
                end if
                ux(i, j) = speed(i, j) * ex(i, j)
                uy(i, j) = speed(i, j) * ey(i, j)
            end do
        end do
        if (present(Q)) Q = Qw
    end subroutine plastic_balance_velocity

    ! Basal velocity the geometry and smb demand: u_basal = u_balance - u_def(tau_b, H).
    ! Negative where the balance flux is below what deformation alone carries.
    subroutine plastic_implied_basal_velocity(z_srf, H, tau_b, mask, smb, dx, dy, rheo, &
                                              u_basal, u_balance, u_def, Q, slope_eps, H_min, stat)
        real(dp),             intent(in)  :: z_srf(:,:), H(:,:), tau_b(:,:), smb(:,:), dx, dy
        logical,              intent(in)  :: mask(:,:)
        type(plastic_glen_t), intent(in)  :: rheo
        real(dp),             intent(out) :: u_basal(:,:)
        real(dp), optional,   intent(out) :: u_balance(:,:), u_def(:,:), Q(:,:)
        real(dp), optional,   intent(in)  :: slope_eps, H_min
        integer,  optional,   intent(out) :: stat
        real(dp), allocatable :: spd(:,:), ux(:,:), uy(:,:)
        integer  :: nx, ny, i, j, st
        real(dp) :: ud

        nx = size(z_srf, 1)
        ny = size(z_srf, 2)
        allocate(spd(nx, ny), ux(nx, ny), uy(nx, ny))
        st = PLASTIC_OK
        if (any(shape(tau_b) /= [nx, ny]) .or. any(shape(u_basal) /= [nx, ny])) st = PLASTIC_BAD_INPUT
        if (st == PLASTIC_OK) call plastic_balance_velocity(z_srf, H, mask, smb, dx, dy, spd, ux, uy, &
                                                            Q, slope_eps, H_min, st)
        if (present(stat)) stat = st
        if (st /= PLASTIC_OK) return

        do j = 1, ny
            do i = 1, nx
                if (mask(i, j)) then
                    ud = deformational_velocity(rheo, tau_b(i, j), H(i, j))
                    u_basal(i, j) = spd(i, j) - ud
                else
                    ud = 0.0_dp
                    u_basal(i, j) = 0.0_dp
                end if
                if (present(u_def)) u_def(i, j) = ud
            end do
        end do
        if (present(u_balance)) u_balance = spd
    end subroutine plastic_implied_basal_velocity

    ! --- shape helpers -------------------------------------------------------------

    pure logical function same_shape(ps, a)
        type(plastic_t), intent(in) :: ps
        real(dp),        intent(in) :: a(:,:)
        same_shape = (size(a, 1) == ps%nx .and. size(a, 2) == ps%ny)
    end function same_shape

    pure logical function same_shape_l(ps, a)
        type(plastic_t), intent(in) :: ps
        logical,         intent(in) :: a(:,:)
        same_shape_l = (size(a, 1) == ps%nx .and. size(a, 2) == ps%ny)
    end function same_shape_l

end module plastic_flow
