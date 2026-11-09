! PlasticIceSheet: perfectly-plastic, steady-state ice-sheet reconstruction.
!
! Public API. Fortran port of PlasticIceSheet.jl: the plastic surface equation
! |grad z_s| = tau_b / (rho_i g H) is solved for w = H^2 as a static eikonal problem
! by fast sweeping, on a grid with a prescribed grounded-ice mask.
module plastic_ice

    use plastic_kinds,   only: plastic_dp => dp
    use plastic_eikonal, only: plastic_par_t, plastic_flotation_thickness
    use plastic_solver,  only: plastic_t, plastic_init, plastic_solve, plastic_reset, &
                               plastic_end, plastic_volume, plastic_area, &
                               plastic_surface_misfit, plastic_residual, &
                               plastic_stress_residual, plastic_margin_frac_from_sdist, &
                               PLASTIC_OK, PLASTIC_NOT_CONVERGED, PLASTIC_BAD_INPUT
    use plastic_flow,    only: plastic_glen_t, plastic_sliding_t, &
                               plastic_deformational_velocity, plastic_basal_velocity, &
                               plastic_flow_direction, plastic_velocity, plastic_smb_implied, &
                               plastic_balance_flux, plastic_balance_velocity, &
                               plastic_implied_basal_velocity
    use plastic_reconstruct, only: plastic_smb_ctrl_t, plastic_solve_smb, &
                                   plastic_invert_driving_stress

    implicit none

    private
    public :: plastic_dp            ! = real64; not exported as dp to avoid clashing with host kinds
    public :: plastic_par_t, plastic_t
    public :: plastic_init, plastic_solve, plastic_reset, plastic_end
    public :: plastic_flotation_thickness
    public :: plastic_volume, plastic_area, plastic_surface_misfit, plastic_residual
    public :: plastic_stress_residual, plastic_margin_frac_from_sdist
    public :: PLASTIC_OK, PLASTIC_NOT_CONVERGED, PLASTIC_BAD_INPUT
    ! closures
    public :: plastic_glen_t, plastic_sliding_t
    public :: plastic_deformational_velocity, plastic_basal_velocity
    public :: plastic_flow_direction, plastic_velocity, plastic_smb_implied
    public :: plastic_balance_flux, plastic_balance_velocity, plastic_implied_basal_velocity
    public :: plastic_smb_ctrl_t, plastic_solve_smb, plastic_invert_driving_stress

end module plastic_ice
