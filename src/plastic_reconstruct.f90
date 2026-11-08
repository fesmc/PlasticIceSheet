! SMB-driven reconstruction (PlasticIceSheet.jl diva_reconstruct): the geometry on a
! fixed mask whose DIVA flux balances a surface mass balance. tau_b is derived, not
! prescribed. Each outer iteration
!
!   1. routes smb down the current surface        -> required throughput Q
!   2. inverts the DIVA flux law per cell for tau_d with Q/w = H (u_b(tau_d) + u_def(tau_d,H))
!   3. sweeps the plastic eikonal (flat RHS) with tau_d -> new thickness
!   4. under-relaxes H and repeats.
!
! The plastic solution is a fixed point when smb is its own implied SMB.
module plastic_reconstruct

    use plastic_kinds,   only: dp, BIG
    use plastic_eikonal, only: flotation_bc, sweep_eikonal
    use plastic_solver,  only: plastic_t, check_shape, warm_start_cells, &
                               PLASTIC_OK, PLASTIC_NOT_CONVERGED, PLASTIC_BAD_INPUT
    use plastic_flow,    only: plastic_glen_t, plastic_sliding_t, LAW_COULOMB, &
                               sliding_law_code, basal_velocity_code, deformational_velocity, &
                               balance_flux_core

    implicit none

    private
    public :: plastic_smb_ctrl_t, plastic_solve_smb, plastic_invert_driving_stress

    ! Controls of plastic_solve_smb (defaults as diva_reconstruct).
    type plastic_smb_ctrl_t
        integer  :: n_outer    = 1000       ! outer iterations cap
        real(dp) :: relax      = 0.05_dp    ! under-relaxation in (0,1]
        real(dp) :: tol        = 1.0_dp     ! exit on max|dH| [m] (before relaxation)
        real(dp) :: tau_hi     = 5.0e5_dp   ! upper bracket of the tau_d inversion [Pa]
        integer  :: max_sweeps = 400        ! inner eikonal controls
        real(dp) :: sweep_tol  = 1.0e-9_dp
        real(dp) :: tau_init   = 5.0e4_dp   ! uniform tau_b of the cold-start flat solve [Pa]
    end type plastic_smb_ctrl_t

contains

    ! Driving stress tau >= 0 with q = H (u_b(tau) + u_def(tau, H)) (monotone in tau),
    ! by 60 bisections on [0, tau_hi] (capped at 0.999 C for "coulomb"); saturates at
    ! the bracket top where even that cannot carry q.
    function invert_driving_stress(q, H, code, slide, rheo, tau_hi) result(tau)
        real(dp),                intent(in) :: q, H, tau_hi
        integer,                 intent(in) :: code
        type(plastic_sliding_t), intent(in) :: slide
        type(plastic_glen_t),    intent(in) :: rheo
        real(dp) :: tau, lo, hi, mid
        integer  :: k

        if (q <= 0.0_dp .or. H <= 0.0_dp) then
            tau = 0.0_dp
            return
        end if
        hi = tau_hi
        if (code == LAW_COULOMB) hi = min(tau_hi, 0.999_dp * slide%C)
        if (flux(hi) <= q) then
            tau = hi
            return
        end if
        lo = 0.0_dp
        do k = 1, 60
            mid = (lo + hi) / 2.0_dp
            if (flux(mid) < q) then
                lo = mid
            else
                hi = mid
            end if
        end do
        tau = (lo + hi) / 2.0_dp

    contains

        real(dp) function flux(t)
            real(dp), intent(in) :: t
            flux = H * (basal_velocity_code(code, slide, t) + deformational_velocity(rheo, t, H))
        end function flux

    end function invert_driving_stress

    ! Public scalar form of the inversion (per-cell flux q = Q / sqrt(dx dy) [m2/yr]).
    function plastic_invert_driving_stress(q, H, slide, rheo, tau_hi) result(tau)
        real(dp),                intent(in) :: q, H, tau_hi
        type(plastic_sliding_t), intent(in) :: slide
        type(plastic_glen_t),    intent(in) :: rheo
        real(dp) :: tau
        tau = invert_driving_stress(q, H, sliding_law_code(slide), slide, rheo, tau_hi)
    end function plastic_invert_driving_stress

    ! SMB-driven reconstruction of H, z_srf (and w = H^2) on ps, from bed z_bed, sea
    ! level z_sl, smb [m ice / yr] and the grounded-ice mask, for the DIVA flux law
    ! (rheo, slide). The inner sweep uses the flat RHS (bed relief to first order).
    !
    ! Start: H_init if given; else, with warm=.true. and a stored solution, the stored
    ! H on cells that were ice before and a cold-start value on the new cells (cold
    ! start if more than par%warm_frac_max of the ice cells are new); else cold, a flat
    ! solve at the uniform tau_b = ctrl%tau_init.
    !
    ! tau_d (optional) returns the derived driving stress [Pa] of the last iteration.
    ! Non-convergence does not abort: ps%converged, ps%stat and ps%msg report it.
    subroutine plastic_solve_smb(ps, z_bed, z_sl, smb, mask, rheo, slide, warm, ctrl, &
                                 H_init, tau_d, stat)
        type(plastic_t),                    intent(inout) :: ps
        real(dp),                           intent(in)    :: z_bed(:,:), z_sl(:,:), smb(:,:)
        logical,                            intent(in)    :: mask(:,:)
        type(plastic_glen_t),               intent(in)    :: rheo
        type(plastic_sliding_t),            intent(in)    :: slide
        logical,                  optional, intent(in)    :: warm
        type(plastic_smb_ctrl_t), optional, intent(in)    :: ctrl
        real(dp),                 optional, intent(in)    :: H_init(:,:)
        real(dp),                 optional, intent(out)   :: tau_d(:,:)
        integer,                  optional, intent(out)   :: stat

        type(plastic_smb_ctrl_t) :: c
        integer  :: nx, ny, i, j, it, code, ns, n_new, n_ice, n_unreached
        logical  :: ok, inner_ok, outer_ok, start_warm
        real(dp) :: rg, wf, delta, Hn, r
        character(len=32) :: s

        if (present(ctrl)) c = ctrl
        code = sliding_law_code(slide)

        ! --- input checks -------------------------------------------------------------
        ps%stat = PLASTIC_OK
        ps%msg  = ""
        if (.not. allocated(ps%w)) then
            ps%stat = PLASTIC_BAD_INPUT
            ps%msg  = "plastic: not initialized (call plastic_init)"
        else if (code == 0) then
            ps%stat = PLASTIC_BAD_INPUT
            ps%msg  = "plastic_solve_smb: slide%law must be linear, weertman or coulomb"
        else if (.not. (c%relax > 0.0_dp .and. c%relax <= 1.0_dp)) then
            ps%stat = PLASTIC_BAD_INPUT
            ps%msg  = "plastic_solve_smb: ctrl%relax must be in (0,1]"
        else
            call check_shape(ps, size(z_bed, 1), size(z_bed, 2), "z_bed")
            call check_shape(ps, size(z_sl, 1),  size(z_sl, 2),  "z_sl")
            call check_shape(ps, size(smb, 1),   size(smb, 2),   "smb")
            call check_shape(ps, size(mask, 1),  size(mask, 2),  "mask")
            if (present(H_init)) call check_shape(ps, size(H_init, 1), size(H_init, 2), "H_init")
            if (present(tau_d))  call check_shape(ps, size(tau_d, 1),  size(tau_d, 2),  "tau_d")
        end if
        if (ps%stat /= PLASTIC_OK) then
            ps%converged = .false.
            if (present(stat)) stat = ps%stat
            return
        end if

        nx = ps%nx
        ny = ps%ny
        rg = ps%par%rho_i * ps%par%g
        wf = sqrt(ps%dx * ps%dy)
        ps%n_sweeps     = 0
        ps%n_outer_used = 0
        ps%outer_change = 0.0_dp
        ps%warm_used    = .false.
        ps%n_new_cells  = 0
        inner_ok = .true.

        call flotation_bc(nx, ny, z_bed, z_sl, ps%par, ps%wk%bc)

        ! --- initial thickness ----------------------------------------------------------
        start_warm = .false.
        if (present(H_init)) then
            ps%H = H_init
        else
            if (present(warm)) start_warm = warm .and. ps%has_warm
            if (start_warm) then
                call warm_start_cells(ps, mask, n_new, n_ice)
                ps%n_new_cells = n_new
                start_warm = real(n_new, dp) <= ps%par%warm_frac_max * real(n_ice, dp)
            end if
            if (start_warm) then
                ps%warm_used = .true.
                if (n_new > 0) call cold_start(ps%wk%wc)
                do j = 1, ny
                    do i = 1, nx
                        if (.not. mask(i, j)) then
                            ps%H(i, j) = 0.0_dp
                        else if (.not. (ps%wk%mask_prev(i, j) .and. ps%w(i, j) < BIG)) then
                            ps%H(i, j) = sqrt(max(ps%wk%wc(i, j), 0.0_dp))
                        end if
                    end do
                end do
            else
                call cold_start(ps%wk%wc)
                do j = 1, ny
                    do i = 1, nx
                        if (mask(i, j)) then
                            ps%H(i, j) = sqrt(max(ps%wk%wc(i, j), 0.0_dp))
                        else
                            ps%H(i, j) = 0.0_dp
                        end if
                    end do
                end do
            end if
        end if

        ! --- outer iteration: route, invert, sweep, relax --------------------------------
        r = c%relax
        outer_ok = .false.
        do it = 1, c%n_outer
            do j = 1, ny
                do i = 1, nx
                    if (mask(i, j)) then
                        ps%z_srf(i, j) = ps%H(i, j) + z_bed(i, j)
                    else
                        ps%z_srf(i, j) = z_bed(i, j)
                    end if
                end do
            end do
            call balance_flux_core(nx, ny, ps%z_srf, mask, smb, ps%dx, ps%dy, ps%wk%f1, &
                                   ps%wk%idx, ps%wk%key, ps%wk%ti, ps%wk%tk)
            do j = 1, ny
                do i = 1, nx
                    if (mask(i, j)) then
                        ps%wk%f2(i, j) = invert_driving_stress(ps%wk%f1(i, j) / wf, ps%H(i, j), &
                                                               code, slide, rheo, c%tau_hi)
                    else
                        ps%wk%f2(i, j) = 0.0_dp
                    end if
                    ps%wk%G(i, j) = 2.0_dp * ps%wk%f2(i, j) / rg
                end do
            end do
            call sweep_eikonal(nx, ny, ps%wk%wc, ps%wk%G, mask, ps%wk%bc, ps%dx, ps%dy, &
                               c%max_sweeps, c%sweep_tol, ns, ok)
            ps%n_sweeps = ps%n_sweeps + ns
            inner_ok = inner_ok .and. ok

            delta = 0.0_dp
            do j = 1, ny
                do i = 1, nx
                    if (mask(i, j)) then
                        Hn = sqrt(max(ps%wk%wc(i, j), 0.0_dp))
                    else
                        Hn = 0.0_dp
                    end if
                    delta = max(delta, abs(Hn - ps%H(i, j)))
                    ps%H(i, j) = (1.0_dp - r) * ps%H(i, j) + r * Hn
                end do
            end do
            ps%n_outer_used = it
            ps%outer_change = delta
            if (delta < c%tol) then
                outer_ok = .true.
                exit
            end if
        end do

        ! --- finalize ------------------------------------------------------------------
        ps%converged = outer_ok .and. inner_ok
        if (.not. outer_ok) then
            ps%msg = "plastic_solve_smb: outer iteration reached n_outer"
        else if (.not. inner_ok) then
            ps%msg = "plastic_solve_smb: an inner sweep reached max_sweeps"
        end if

        n_unreached = 0
        do j = 1, ny
            do i = 1, nx
                if (mask(i, j) .and. ps%wk%wc(i, j) >= BIG) then
                    n_unreached = n_unreached + 1
                    ps%H(i, j) = 0.0_dp
                end if
                if (mask(i, j)) then
                    ps%z_srf(i, j) = ps%H(i, j) + z_bed(i, j)
                    ps%w(i, j)     = ps%H(i, j) * ps%H(i, j)
                else
                    ps%z_srf(i, j) = z_bed(i, j)
                    ps%w(i, j)     = ps%wk%bc(i, j)
                end if
            end do
        end do
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
        if (present(tau_d)) tau_d = ps%wk%f2
        ps%wk%mask_prev = mask
        ps%has_warm = .true.
        if (present(stat)) stat = ps%stat

    contains

        ! Flat plastic solve at the uniform tau_b = c%tau_init, into u (= w).
        subroutine cold_start(u)
            real(dp), intent(out) :: u(:,:)
            integer :: ns0
            logical :: ok0
            ps%wk%G = 2.0_dp * c%tau_init / rg
            call sweep_eikonal(nx, ny, u, ps%wk%G, mask, ps%wk%bc, ps%dx, ps%dy, &
                               c%max_sweeps, c%sweep_tol, ns0, ok0)
            ps%n_sweeps = ps%n_sweeps + ns0
            inner_ok = inner_ok .and. ok0
        end subroutine cold_start

    end subroutine plastic_solve_smb

end module plastic_reconstruct
