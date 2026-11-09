! The plastic ice-sheet solver object: geometry from (z_bed, z_sl, tau_b, mask) by the
! damped fixed point on w = H^2, with an optional warm start from the previous solve.
!
! All work arrays live in the plastic_t instance and are allocated in plastic_init,
! so a solve allocates nothing and several instances can coexist (no module state).
module plastic_solver

    use plastic_kinds,   only: dp, BIG, MSG_LEN
    use plastic_eikonal
    use plastic_linsolve, only: plastic_lin_t

    implicit none

    private
    public :: plastic_t, plastic_work_t
    public :: plastic_init, plastic_solve, plastic_reset, plastic_end
    public :: plastic_volume, plastic_area, plastic_surface_misfit, plastic_residual
    public :: plastic_stress_residual, plastic_margin_frac_from_sdist
    public :: margin_spacing
    public :: PLASTIC_OK, PLASTIC_NOT_CONVERGED, PLASTIC_BAD_INPUT
    public :: check_shape, cold_flat_solve, finalize_geometry, warm_start_cells

    ! Status codes (ps%stat and the optional stat arguments).
    integer, parameter :: PLASTIC_OK            = 0   ! converged
    integer, parameter :: PLASTIC_NOT_CONVERGED = 1   ! result usable, tolerance not met
    integer, parameter :: PLASTIC_BAD_INPUT     = 2   ! nothing computed, state unchanged

    ! Work arrays (internal; allocated in plastic_init).
    type plastic_work_t
        real(dp), allocatable :: G(:,:)     ! flat RHS 2 tau_b / (rho_i g)
        real(dp), allocatable :: F(:,:)     ! full RHS
        real(dp), allocatable :: wc(:,:)    ! inner sweep result
        real(dp), allocatable :: bc(:,:)    ! margin values H_flt^2
        real(dp), allocatable :: zbx(:,:), zby(:,:)   ! bed gradient
        logical,  allocatable :: mask_prev(:,:)       ! mask of the stored solution
        real(dp), allocatable :: dfac(:,:)            ! spacing factor to each cell as a neighbour
        ! closures (plastic_flow, plastic_reconstruct)
        real(dp), allocatable :: ex(:,:), ey(:,:)     ! unit flow direction
        real(dp), allocatable :: f1(:,:), f2(:,:)     ! scratch: flux, inflow / driving stress
        integer,  allocatable :: idx(:), ti(:)        ! balance-flux ordering (nx*ny)
        real(dp), allocatable :: key(:), tk(:)
        ! steady SMB-driven solver (plastic_steady)
        real(dp), allocatable :: sF(:,:), sFp(:,:), sS(:,:), sHk(:,:), sHt(:,:), sHm(:,:), sSm(:,:)
        real(dp), allocatable :: qx(:,:), qy(:,:)     ! face fluxes (0:nx,ny), (nx,0:ny)
        real(dp), allocatable :: jl(:,:,:)            ! dF_P/dH over the 3 x 3 block, (9,nx,ny)
        integer,  allocatable :: sid(:,:)             ! unknown number of an ice cell (0: none)
        logical,  allocatable :: sact(:,:)            ! active set (H = 0)
        type(plastic_lin_t)   :: lin
    end type plastic_work_t

    type plastic_t
        type(plastic_par_t) :: par
        integer  :: nx = 0, ny = 0
        real(dp) :: dx = 0.0_dp, dy = 0.0_dp
        real(dp), allocatable :: w(:,:)      ! H^2 of the last solve (warm-start state)
        real(dp), allocatable :: H(:,:)      ! ice thickness [m]
        real(dp), allocatable :: z_srf(:,:)  ! ice surface elevation [m] (z_bed off the mask)
        integer  :: n_outer_used = 0         ! outer iterations of the last solve
        logical  :: converged    = .false.
        ! Diagnostics of the last solve
        integer  :: stat = PLASTIC_OK        ! PLASTIC_OK | PLASTIC_NOT_CONVERGED | PLASTIC_BAD_INPUT
        character(len=MSG_LEN) :: msg = ""
        integer  :: n_sweeps     = 0         ! inner sweeps in total
        real(dp) :: outer_change = 0.0_dp    ! last outer max|dH| [m]
        logical  :: warm_used    = .false.   ! the last solve started from the stored w
        integer  :: n_new_cells  = 0         ! ice cells with no valid stored w (warm start)
        logical  :: has_warm     = .false.   ! w, mask_prev hold a usable solution
        type(plastic_work_t) :: wk
    end type plastic_t

contains

    subroutine plastic_init(ps, nx, ny, dx, dy, par, stat)
        type(plastic_t),               intent(inout) :: ps
        integer,                       intent(in)    :: nx, ny
        real(dp),                      intent(in)    :: dx, dy
        type(plastic_par_t), optional, intent(in)    :: par
        integer,             optional, intent(out)   :: stat

        call plastic_end(ps)

        if (present(par)) ps%par = par
        if (nx < 1 .or. ny < 1 .or. .not. (dx > 0.0_dp) .or. .not. (dy > 0.0_dp)) then
            ps%stat = PLASTIC_BAD_INPUT
            ps%msg  = "plastic_init: need nx, ny >= 1 and dx, dy > 0"
            if (present(stat)) stat = ps%stat
            return
        end if

        ps%nx = nx
        ps%ny = ny
        ps%dx = dx
        ps%dy = dy

        allocate(ps%w(nx, ny), ps%H(nx, ny), ps%z_srf(nx, ny))
        ps%w     = 0.0_dp
        ps%H     = 0.0_dp
        ps%z_srf = 0.0_dp

        allocate(ps%wk%G(nx, ny), ps%wk%F(nx, ny), ps%wk%wc(nx, ny), ps%wk%bc(nx, ny))
        allocate(ps%wk%zbx(nx, ny), ps%wk%zby(nx, ny), ps%wk%mask_prev(nx, ny), ps%wk%dfac(nx, ny))
        ps%wk%mask_prev = .false.
        allocate(ps%wk%ex(nx, ny), ps%wk%ey(nx, ny), ps%wk%f1(nx, ny), ps%wk%f2(nx, ny))
        allocate(ps%wk%idx(nx*ny), ps%wk%ti(nx*ny), ps%wk%key(nx*ny), ps%wk%tk(nx*ny))
        allocate(ps%wk%sF(nx, ny), ps%wk%sFp(nx, ny), ps%wk%sS(nx, ny), ps%wk%sHk(nx, ny), &
                 ps%wk%sHt(nx, ny), ps%wk%sHm(nx, ny), ps%wk%sSm(nx, ny), ps%wk%qx(0:nx, ny), &
                 ps%wk%qy(nx, 0:ny), ps%wk%jl(9, nx, ny), ps%wk%sid(nx, ny), ps%wk%sact(nx, ny))

        ps%has_warm = .false.
        ps%stat = PLASTIC_OK
        ps%msg  = ""
        if (present(stat)) stat = ps%stat
    end subroutine plastic_init

    ! Drop the warm start: the next solve starts cold.
    subroutine plastic_reset(ps)
        type(plastic_t), intent(inout) :: ps
        ps%has_warm = .false.
        if (allocated(ps%wk%mask_prev)) ps%wk%mask_prev = .false.
    end subroutine plastic_reset

    subroutine plastic_end(ps)
        type(plastic_t), intent(inout) :: ps
        if (allocated(ps%w))            deallocate(ps%w)
        if (allocated(ps%H))            deallocate(ps%H)
        if (allocated(ps%z_srf))        deallocate(ps%z_srf)
        if (allocated(ps%wk%G))         deallocate(ps%wk%G)
        if (allocated(ps%wk%F))         deallocate(ps%wk%F)
        if (allocated(ps%wk%wc))        deallocate(ps%wk%wc)
        if (allocated(ps%wk%bc))        deallocate(ps%wk%bc)
        if (allocated(ps%wk%zbx))       deallocate(ps%wk%zbx)
        if (allocated(ps%wk%zby))       deallocate(ps%wk%zby)
        if (allocated(ps%wk%mask_prev)) deallocate(ps%wk%mask_prev)
        if (allocated(ps%wk%dfac))      deallocate(ps%wk%dfac)
        if (allocated(ps%wk%ex))        deallocate(ps%wk%ex)
        if (allocated(ps%wk%ey))        deallocate(ps%wk%ey)
        if (allocated(ps%wk%f1))        deallocate(ps%wk%f1)
        if (allocated(ps%wk%f2))        deallocate(ps%wk%f2)
        if (allocated(ps%wk%idx))       deallocate(ps%wk%idx)
        if (allocated(ps%wk%ti))        deallocate(ps%wk%ti)
        if (allocated(ps%wk%key))       deallocate(ps%wk%key)
        if (allocated(ps%wk%tk))        deallocate(ps%wk%tk)
        if (allocated(ps%wk%sF)) deallocate(ps%wk%sF, ps%wk%sFp, ps%wk%sS, ps%wk%sHk, ps%wk%sHt, &
                                            ps%wk%sHm, ps%wk%sSm, ps%wk%qx, ps%wk%qy, ps%wk%jl, &
                                            ps%wk%sid, ps%wk%sact)
        ps%wk%lin = plastic_lin_t()
        ps%nx = 0
        ps%ny = 0
        ps%has_warm = .false.
        ps%n_outer_used = 0
        ps%converged = .false.
    end subroutine plastic_end

    ! Plastic reconstruction: H, z_srf (and w = H^2) on ps from bed z_bed [m], sea level
    ! z_sl [m], basal shear stress tau_b [Pa] and the grounded-ice mask.
    !
    ! mode "flat": one eikonal solve |grad w| = 2 tau_b/(rho_i g).
    ! mode "full": damped fixed point on the bed-corrected RHS, started from the flat
    ! solution (cold) or from the stored w of the last solve (warm=.true.).
    ! mode "surface": one fast-sweeping solve of the surface equation
    ! |grad z_s| = tau_b / (rho_i g (z_s - z_bed)) for z_s itself (godunov_surface):
    ! full bed relief, no outer iteration; warm has no effect. Not in PlasticIceSheet.jl.
    !
    ! Warm start (full mode only; a flat solve is a single sweep from +BIG either way):
    ! ice cells that were also ice in the stored solution keep their stored w; ice-free
    ! cells take the new margin value H_flt^2. Cells new to the mask take the flat
    ! solution, which costs one extra eikonal solve. If more than par%warm_frac_max of
    ! the ice cells are new, the solve starts cold. Each inner sweep still starts from
    ! +BIG on the ice cells, as the monotone scheme requires.
    !
    ! Non-convergence does not abort: ps%converged, ps%stat and ps%msg report it.
    subroutine plastic_solve(ps, z_bed, z_sl, tau_b, mask, warm, stat, margin_frac)
        type(plastic_t),    intent(inout) :: ps
        real(dp),           intent(in)    :: z_bed(:,:), z_sl(:,:), tau_b(:,:)
        logical,            intent(in)    :: mask(:,:)
        logical,  optional, intent(in)    :: warm
        integer,  optional, intent(out)   :: stat
        real(dp), optional, intent(in)    :: margin_frac(:,:)

        integer  :: nx, ny, it, i, j, ns, n_new, n_ice
        logical  :: full, want_warm, inner_ok, ok, outer_ok
        real(dp) :: r, conv, wn

        ps%stat = PLASTIC_OK
        ps%msg  = ""
        call check_inputs(ps, z_bed, z_sl, tau_b, mask)
        call margin_spacing(ps, mask, margin_frac)
        if (ps%stat /= PLASTIC_OK) then
            ps%converged = .false.
            if (present(stat)) stat = ps%stat
            return
        end if

        if (trim(ps%par%mode) == "surface") then
            call solve_surface(ps, z_bed, z_sl, tau_b, mask, stat)
            return
        end if

        nx = ps%nx
        ny = ps%ny
        full = (trim(ps%par%mode) == "full")
        want_warm = .false.
        if (present(warm)) want_warm = warm

        ps%n_sweeps     = 0
        ps%n_outer_used = 0
        ps%outer_change = 0.0_dp
        ps%warm_used    = .false.
        ps%n_new_cells  = 0
        inner_ok = .true.

        call flotation_bc(nx, ny, z_bed, z_sl, ps%par, ps%wk%bc)
        call rhs_flat(nx, ny, tau_b, ps%par, ps%wk%G)

        if (.not. full) then
            call sweep_eikonal(nx, ny, ps%w, ps%wk%G, mask, ps%wk%bc, ps%wk%dfac, ps%dx, ps%dy, &
                               ps%par%max_sweeps, ps%par%tol, ns, ok)
            ps%n_sweeps = ns
            ps%converged = ok
            if (.not. ok) ps%msg = "plastic_solve: inner sweep reached max_sweeps"
        else
            call grad_central(nx, ny, z_bed, ps%dx, ps%dy, ps%wk%zbx, ps%wk%zby)

            ! Initial w: warm (stored w, flat solve on the new cells) or cold (flat).
            if (want_warm .and. ps%has_warm) then
                call warm_start_cells(ps, mask, n_new, n_ice)
                ps%n_new_cells = n_new
                if (real(n_new, dp) <= ps%par%warm_frac_max * real(n_ice, dp)) then
                    ps%warm_used = .true.
                    if (n_new > 0) then
                        call cold_flat_solve(ps, mask, ps%wk%wc, inner_ok)
                    end if
                    do j = 1, ny
                        do i = 1, nx
                            if (.not. mask(i, j)) then
                                ps%w(i, j) = ps%wk%bc(i, j)
                            else if (.not. (ps%wk%mask_prev(i, j) .and. ps%w(i, j) < BIG)) then
                                ps%w(i, j) = ps%wk%wc(i, j)
                            end if
                        end do
                    end do
                end if
            end if
            if (.not. ps%warm_used) then
                call cold_flat_solve(ps, mask, ps%w, inner_ok)
            end if

            ! Damped outer fixed point: freeze F(w), re-solve, under-relax.
            r = ps%par%relax
            outer_ok = .false.
            do it = 1, ps%par%n_outer
                call rhs_full(nx, ny, ps%w, ps%wk%G, ps%wk%zbx, ps%wk%zby, ps%dx, ps%dy, &
                              ps%par%H_min, ps%wk%F)
                call sweep_eikonal(nx, ny, ps%wk%wc, ps%wk%F, mask, ps%wk%bc, ps%wk%dfac, ps%dx, ps%dy, &
                                   ps%par%max_sweeps, ps%par%tol, ns, ok)
                ps%n_sweeps = ps%n_sweeps + ns
                inner_ok = inner_ok .and. ok
                conv = 0.0_dp
                do j = 1, ny
                    do i = 1, nx
                        wn = (1.0_dp - r) * ps%w(i, j) + r * ps%wk%wc(i, j)
                        conv = max(conv, abs(sqrt(max(wn, 0.0_dp)) - sqrt(max(ps%w(i, j), 0.0_dp))))
                        ps%w(i, j) = wn
                    end do
                end do
                ps%n_outer_used = it
                ps%outer_change = conv
                if (conv < ps%par%outer_tol) then
                    outer_ok = .true.
                    exit
                end if
            end do

            ps%converged = outer_ok .and. inner_ok
            if (.not. outer_ok) then
                ps%msg = "plastic_solve: outer fixed point reached n_outer"
            else if (.not. inner_ok) then
                ps%msg = "plastic_solve: an inner sweep reached max_sweeps"
            end if
        end if

        call finalize_geometry(ps, z_bed, mask)
        if (present(stat)) stat = ps%stat
    end subroutine plastic_solve

    ! mode "surface": z_s from one fast-sweeping solve of the surface equation.
    ! Ice-free cells hold z_s = z_bed + H_flt; H = z_s - z_bed on the ice.
    subroutine solve_surface(ps, z_bed, z_sl, tau_b, mask, stat)
        type(plastic_t),   intent(inout) :: ps
        real(dp),          intent(in)    :: z_bed(:,:), z_sl(:,:), tau_b(:,:)
        logical,           intent(in)    :: mask(:,:)
        integer, optional, intent(out)   :: stat
        integer :: nx, ny, i, j, ns, n_unreached
        logical :: ok

        nx = ps%nx
        ny = ps%ny
        ps%n_outer_used = 0
        ps%outer_change = 0.0_dp
        ps%warm_used    = .false.
        ps%n_new_cells  = 0

        call surface_bc(ps, z_bed, z_sl, tau_b)
        call sweep_eikonal(nx, ny, ps%wk%wc, ps%wk%G, mask, ps%wk%bc, ps%wk%dfac, ps%dx, ps%dy, &
                           ps%par%max_sweeps, ps%par%tol, ns, ok, zb=z_bed)
        ps%n_sweeps  = ns
        ps%converged = ok
        if (.not. ok) ps%msg = "plastic_solve: surface sweep reached max_sweeps"

        n_unreached = 0
        do j = 1, ny
            do i = 1, nx
                if (mask(i, j) .and. ps%wk%wc(i, j) < BIG) then
                    ps%z_srf(i, j) = ps%wk%wc(i, j)
                    ps%H(i, j)     = ps%wk%wc(i, j) - z_bed(i, j)
                    ps%w(i, j)     = ps%H(i, j) * ps%H(i, j)
                else if (mask(i, j)) then
                    n_unreached    = n_unreached + 1
                    ps%z_srf(i, j) = z_bed(i, j)
                    ps%H(i, j)     = 0.0_dp
                    ps%w(i, j)     = BIG
                else
                    ps%z_srf(i, j) = z_bed(i, j)
                    ps%H(i, j)     = 0.0_dp
                    ps%w(i, j)     = (ps%wk%bc(i, j) - z_bed(i, j))**2
                end if
            end do
        end do
        call finalize_status(ps, mask, n_unreached)
        if (present(stat)) stat = ps%stat
    end subroutine solve_surface

    ! Surface-mode boundary values bc = z_bed + H_flt (all cells) and K = tau_b/(rho_i g)
    ! in ps%wk%bc and ps%wk%G.
    subroutine surface_bc(ps, z_bed, z_sl, tau_b)
        type(plastic_t), intent(inout) :: ps
        real(dp),        intent(in)    :: z_bed(:,:), z_sl(:,:), tau_b(:,:)
        integer  :: i, j
        real(dp) :: rg

        rg = ps%par%rho_i * ps%par%g
        do j = 1, ps%ny
            do i = 1, ps%nx
                ps%wk%bc(i, j) = z_bed(i, j) + plastic_flotation_thickness(z_bed(i, j), z_sl(i, j), ps%par)
                ps%wk%G(i, j)  = tau_b(i, j) / rg
            end do
        end do
    end subroutine surface_bc

    ! Count the ice cells without a valid stored w (new to the mask) and all ice cells.
    subroutine warm_start_cells(ps, mask, n_new, n_ice)
        type(plastic_t), intent(in)  :: ps
        logical,         intent(in)  :: mask(:,:)
        integer,         intent(out) :: n_new, n_ice
        integer :: i, j

        n_new = 0
        n_ice = 0
        do j = 1, ps%ny
            do i = 1, ps%nx
                if (.not. mask(i, j)) cycle
                n_ice = n_ice + 1
                if (.not. (ps%wk%mask_prev(i, j) .and. ps%w(i, j) < BIG)) n_new = n_new + 1
            end do
        end do
    end subroutine warm_start_cells

    ! Flat eikonal solve with the RHS in ps%wk%G and margins in ps%wk%bc.
    subroutine cold_flat_solve(ps, mask, u, ok)
        type(plastic_t), intent(inout) :: ps
        logical,         intent(in)    :: mask(:,:)
        real(dp),        intent(out)   :: u(:,:)
        logical,         intent(inout) :: ok
        integer :: ns
        logical :: ok1

        call sweep_eikonal(ps%nx, ps%ny, u, ps%wk%G, mask, ps%wk%bc, ps%wk%dfac, ps%dx, ps%dy, &
                           ps%par%max_sweeps, ps%par%tol, ns, ok1)
        ps%n_sweeps = ps%n_sweeps + ns
        ok = ok .and. ok1
    end subroutine cold_flat_solve

    ! H and z_srf from ps%w; flag ice cells the sweep never reached (an ice region
    ! with no ice-free cell); store the mask for the next warm start.
    subroutine finalize_geometry(ps, z_bed, mask)
        type(plastic_t), intent(inout) :: ps
        real(dp),        intent(in)    :: z_bed(:,:)
        logical,         intent(in)    :: mask(:,:)
        integer :: i, j, n_unreached

        n_unreached = 0
        do j = 1, ps%ny
            do i = 1, ps%nx
                if (mask(i, j) .and. ps%w(i, j) < BIG) then
                    ps%H(i, j)     = sqrt(max(ps%w(i, j), 0.0_dp))
                    ps%z_srf(i, j) = ps%H(i, j) + z_bed(i, j)
                else
                    if (mask(i, j)) n_unreached = n_unreached + 1
                    ps%H(i, j)     = 0.0_dp
                    ps%z_srf(i, j) = z_bed(i, j)
                end if
            end do
        end do
        call finalize_status(ps, mask, n_unreached)
    end subroutine finalize_geometry

    ! Status from ps%converged and the unreached ice cells; store the mask for the
    ! next warm start.
    subroutine finalize_status(ps, mask, n_unreached)
        type(plastic_t), intent(inout) :: ps
        logical,         intent(in)    :: mask(:,:)
        integer,         intent(in)    :: n_unreached
        character(len=32) :: s

        if (n_unreached > 0) then
            ps%converged = .false.
            write(s, '(i0)') n_unreached
            ps%msg = "plastic: "//trim(s)//" ice cells not connected to an ice-free cell (H set to 0)"
        end if

        if (ps%converged) then
            ps%stat = PLASTIC_OK
        else
            ps%stat = PLASTIC_NOT_CONVERGED
        end if

        ps%wk%mask_prev = mask
        ps%has_warm = .true.
    end subroutine finalize_status

    ! Shape of an (nx,ny) input; sets ps%stat/msg on mismatch.
    subroutine check_shape(ps, nx, ny, name)
        type(plastic_t),  intent(inout) :: ps
        integer,          intent(in)    :: nx, ny
        character(len=*), intent(in)    :: name
        if (ps%stat /= PLASTIC_OK) return
        if (nx /= ps%nx .or. ny /= ps%ny) then
            ps%stat = PLASTIC_BAD_INPUT
            ps%msg  = "plastic: "//name//" does not have shape (nx,ny)"
        end if
    end subroutine check_shape

    subroutine check_inputs(ps, z_bed, z_sl, tau_b, mask)
        type(plastic_t), intent(inout) :: ps
        real(dp),        intent(in)    :: z_bed(:,:), z_sl(:,:), tau_b(:,:)
        logical,         intent(in)    :: mask(:,:)

        if (.not. allocated(ps%w)) then
            ps%stat = PLASTIC_BAD_INPUT
            ps%msg  = "plastic: not initialized (call plastic_init)"
            return
        end if
        if (trim(ps%par%mode) /= "full" .and. trim(ps%par%mode) /= "flat" .and. &
            trim(ps%par%mode) /= "surface") then
            ps%stat = PLASTIC_BAD_INPUT
            ps%msg  = "plastic: par%mode must be ""full"", ""flat"" or ""surface"""
            return
        end if
        call check_shape(ps, size(z_bed, 1), size(z_bed, 2), "z_bed")
        call check_shape(ps, size(z_sl, 1),  size(z_sl, 2),  "z_sl")
        call check_shape(ps, size(tau_b, 1), size(tau_b, 2), "tau_b")
        call check_shape(ps, size(mask, 1),  size(mask, 2),  "mask")
    end subroutine check_inputs

    ! --- diagnostics --------------------------------------------------------------

    ! Total ice volume [m3].
    function plastic_volume(ps) result(V)
        type(plastic_t), intent(in) :: ps
        real(dp) :: V
        V = sum(ps%H) * ps%dx * ps%dy
    end function plastic_volume

    ! Grounded-ice area [m2].
    function plastic_area(ps, mask) result(A)
        type(plastic_t), intent(in) :: ps
        logical,         intent(in) :: mask(:,:)
        real(dp) :: A
        A = real(count(mask), dp) * ps%dx * ps%dy
    end function plastic_area

    ! Mean-squared misfit of ps%z_srf to z_srf_obs over the mask [m2].
    function plastic_surface_misfit(ps, z_srf_obs, mask) result(mis)
        type(plastic_t), intent(in) :: ps
        real(dp),        intent(in) :: z_srf_obs(:,:)
        logical,         intent(in) :: mask(:,:)
        real(dp) :: mis
        integer  :: i, j, n
        real(dp) :: d

        mis = 0.0_dp
        n = 0
        do j = 1, ps%ny
            do i = 1, ps%nx
                if (mask(i, j)) then
                    d = ps%z_srf(i, j) - z_srf_obs(i, j)
                    mis = mis + d * d
                    n = n + 1
                end if
            end do
        end do
        if (n > 0) mis = mis / real(n, dp)
    end function plastic_surface_misfit

    ! Godunov residual of the stored solution for the current par%mode (convergence
    ! check): r = w - (Godunov update) on ice cells, w - H_flt^2 elsewhere (in mode
    ! "surface" the same for z_s, with z_bed + H_flt); ~0 when converged.
    subroutine plastic_residual(ps, z_bed, z_sl, tau_b, mask, r, stat, margin_frac)
        type(plastic_t),    intent(inout) :: ps
        real(dp),           intent(in)    :: z_bed(:,:), z_sl(:,:), tau_b(:,:)
        logical,            intent(in)    :: mask(:,:)
        real(dp),           intent(out)   :: r(:,:)
        integer,  optional, intent(out)   :: stat
        real(dp), optional, intent(in)    :: margin_frac(:,:)
        integer :: nx, ny, stat_solve
        character(len=MSG_LEN) :: msg_solve

        ! Keep the status of the last solve unless the inputs are bad.
        stat_solve = ps%stat
        msg_solve  = ps%msg
        ps%stat = PLASTIC_OK
        call check_inputs(ps, z_bed, z_sl, tau_b, mask)
        call check_shape(ps, size(r, 1), size(r, 2), "r")
        call margin_spacing(ps, mask, margin_frac)
        if (present(stat)) stat = ps%stat
        if (ps%stat /= PLASTIC_OK) return
        ps%stat = stat_solve
        ps%msg  = msg_solve

        nx = ps%nx
        ny = ps%ny
        if (trim(ps%par%mode) == "surface") then
            call surface_bc(ps, z_bed, z_sl, tau_b)
            call surface_field(ps, mask, ps%wk%wc)
            call godunov_residual_surface(nx, ny, ps%wk%wc, ps%wk%G, z_bed, mask, ps%wk%bc, &
                                          ps%wk%dfac, ps%dx, ps%dy, r)
            return
        end if
        call flotation_bc(nx, ny, z_bed, z_sl, ps%par, ps%wk%bc)
        call rhs_flat(nx, ny, tau_b, ps%par, ps%wk%G)
        if (trim(ps%par%mode) == "full") then
            call grad_central(nx, ny, z_bed, ps%dx, ps%dy, ps%wk%zbx, ps%wk%zby)
            call rhs_full(nx, ny, ps%w, ps%wk%G, ps%wk%zbx, ps%wk%zby, ps%dx, ps%dy, &
                          ps%par%H_min, ps%wk%F)
            call godunov_residual(nx, ny, ps%w, ps%wk%F, mask, ps%wk%bc, ps%wk%dfac, ps%dx, ps%dy, r)
        else
            call godunov_residual(nx, ny, ps%w, ps%wk%G, mask, ps%wk%bc, ps%wk%dfac, ps%dx, ps%dy, r)
        end if
    end subroutine plastic_residual

    ! Surface field of the stored solution with the margin values: ps%z_srf on the
    ! ice, ps%wk%bc (= z_bed + H_flt, from surface_bc) elsewhere.
    subroutine surface_field(ps, mask, u)
        type(plastic_t), intent(in)  :: ps
        logical,         intent(in)  :: mask(:,:)
        real(dp),        intent(out) :: u(:,:)
        integer :: i, j
        do j = 1, ps%ny
            do i = 1, ps%nx
                if (mask(i, j)) then
                    u(i, j) = ps%z_srf(i, j)
                else
                    u(i, j) = ps%wk%bc(i, j)
                end if
            end do
        end do
    end subroutine surface_field

    ! Stress residual of the stored geometry (any mode): on ice cells with tau_b > 0
    ! and H > 0, r = | |grad z_s| rho_i g H - tau_b | / tau_b, with the Godunov upwind
    ! gradient |grad z_s|^2 = (max_s (z_s - a_s)+/d_s)^2 + (same in y) over the two sides
    ! (z_bed + H_flt on ice-free neighbours, at the margin spacing). rmax = max r; r = 0
    ! elsewhere. A cell-centred measure: on the exact Nye profile it is O(1) at the first
    ! ice cell (the one-sided slope of a sqrt profile), so read it away from the margin.
    subroutine plastic_stress_residual(ps, z_bed, z_sl, tau_b, mask, rmax, r, stat, margin_frac)
        type(plastic_t),    intent(inout) :: ps
        real(dp),           intent(in)    :: z_bed(:,:), z_sl(:,:), tau_b(:,:)
        logical,            intent(in)    :: mask(:,:)
        real(dp),           intent(out)   :: rmax
        real(dp), optional, intent(out)   :: r(:,:)
        integer,  optional, intent(out)   :: stat
        real(dp), optional, intent(in)    :: margin_frac(:,:)
        integer  :: i, j, im, ip, jm, jp, stat_solve
        real(dp) :: u, gx, gy, rr, rg
        character(len=MSG_LEN) :: msg_solve

        stat_solve = ps%stat
        msg_solve  = ps%msg
        ps%stat = PLASTIC_OK
        call check_inputs(ps, z_bed, z_sl, tau_b, mask)
        if (present(r)) call check_shape(ps, size(r, 1), size(r, 2), "r")
        call margin_spacing(ps, mask, margin_frac)
        if (present(stat)) stat = ps%stat
        rmax = 0.0_dp
        if (ps%stat /= PLASTIC_OK) return
        ps%stat = stat_solve
        ps%msg  = msg_solve

        call surface_bc(ps, z_bed, z_sl, tau_b)
        call surface_field(ps, mask, ps%wk%wc)
        rg = ps%par%rho_i * ps%par%g
        do j = 1, ps%ny
            jm = max(j - 1, 1)
            jp = min(j + 1, ps%ny)
            do i = 1, ps%nx
                rr = 0.0_dp
                if (mask(i, j) .and. tau_b(i, j) > 0.0_dp .and. ps%H(i, j) > 0.0_dp) then
                    im = max(i - 1, 1)
                    ip = min(i + 1, ps%nx)
                    u  = ps%wk%wc(i, j)
                    gx = max(u - ps%wk%wc(im, j), 0.0_dp) / (ps%dx * ps%wk%dfac(im, j))
                    gx = max(gx, max(u - ps%wk%wc(ip, j), 0.0_dp) / (ps%dx * ps%wk%dfac(ip, j)))
                    gy = max(u - ps%wk%wc(i, jm), 0.0_dp) / (ps%dy * ps%wk%dfac(i, jm))
                    gy = max(gy, max(u - ps%wk%wc(i, jp), 0.0_dp) / (ps%dy * ps%wk%dfac(i, jp)))
                    rr = abs(sqrt(gx * gx + gy * gy) * rg * ps%H(i, j) - tau_b(i, j)) / tau_b(i, j)
                end if
                rmax = max(rmax, rr)
                if (present(r)) r(i, j) = rr
            end do
        end do
    end subroutine plastic_stress_residual

    ! Spacing factor ps%wk%dfac to each cell as a neighbour: 1 on ice cells; on ice-free
    ! cells max(1 - f, MARGIN_DMIN), f = margin_frac if given, else 0.5 (par%margin =
    ! "face") or 0 ("centre"). PLASTIC_BAD_INPUT for an unknown par%margin, a wrong shape,
    ! or margin_frac not in [0, 1) on an ice-free cell. Does nothing if ps%stat is
    ! already an error.
    subroutine margin_spacing(ps, mask, margin_frac)
        type(plastic_t),    intent(inout) :: ps
        logical,            intent(in)    :: mask(:,:)
        real(dp), optional, intent(in)    :: margin_frac(:,:)
        integer  :: i, j
        real(dp) :: f0

        if (ps%stat /= PLASTIC_OK) return
        if (trim(ps%par%margin) /= "face" .and. trim(ps%par%margin) /= "centre") then
            ps%stat = PLASTIC_BAD_INPUT
            ps%msg  = 'plastic: par%margin must be "face" or "centre"'
            return
        end if
        if (present(margin_frac)) then
            call check_shape(ps, size(margin_frac, 1), size(margin_frac, 2), "margin_frac")
            if (ps%stat /= PLASTIC_OK) return
            do j = 1, ps%ny
                do i = 1, ps%nx
                    if (mask(i, j)) cycle
                    if (.not. (margin_frac(i, j) >= 0.0_dp .and. margin_frac(i, j) < 1.0_dp)) then
                        ps%stat = PLASTIC_BAD_INPUT
                        ps%msg  = "plastic: margin_frac must be in [0,1) on ice-free cells"
                        return
                    end if
                end do
            end do
        end if
        f0 = 0.0_dp
        if (trim(ps%par%margin) == "face") f0 = 0.5_dp
        do j = 1, ps%ny
            do i = 1, ps%nx
                if (mask(i, j)) then
                    ps%wk%dfac(i, j) = 1.0_dp
                else if (present(margin_frac)) then
                    ps%wk%dfac(i, j) = max(1.0_dp - margin_frac(i, j), MARGIN_DMIN)
                else
                    ps%wk%dfac(i, j) = max(1.0_dp - f0, MARGIN_DMIN)
                end if
            end do
        end do
    end subroutine margin_spacing

    ! margin_frac from the signed distance to the margin, sdist [m], positive in the ice,
    ! negative outside, at the cell centres: on an ice-free cell with an ice neighbour,
    ! frac = -sdist / h clamped to [0, 1 - MARGIN_DMIN], h = sqrt(dx dy) (= dx on a square
    ! grid); 0 elsewhere. The margin is then at (1 - frac) d from the centre of each
    ! neighbouring ice cell (isotropic: exact for margins along a grid axis).
    subroutine plastic_margin_frac_from_sdist(sdist, mask, dx, dy, frac)
        real(dp), intent(in)  :: sdist(:,:), dx, dy
        logical,  intent(in)  :: mask(:,:)
        real(dp), intent(out) :: frac(:,:)
        integer  :: nx, ny, i, j
        real(dp) :: h
        logical  :: next_to_ice

        nx = size(mask, 1)
        ny = size(mask, 2)
        h = sqrt(dx * dy)
        do j = 1, ny
            do i = 1, nx
                frac(i, j) = 0.0_dp
                if (mask(i, j)) cycle
                next_to_ice = mask(max(i-1,1), j) .or. mask(min(i+1,nx), j) .or. &
                              mask(i, max(j-1,1)) .or. mask(i, min(j+1,ny))
                if (next_to_ice) frac(i, j) = min(max(-sdist(i, j) / h, 0.0_dp), 1.0_dp - MARGIN_DMIN)
            end do
        end do
    end subroutine plastic_margin_frac_from_sdist

end module plastic_solver
