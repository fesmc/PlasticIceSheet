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
                               PLASTIC_OK, PLASTIC_NOT_CONVERGED, PLASTIC_BAD_INPUT

    implicit none

    private
    public :: plastic_dp            ! = real64; not exported as dp to avoid clashing with host kinds
    public :: plastic_par_t, plastic_t
    public :: plastic_init, plastic_solve, plastic_reset, plastic_end
    public :: plastic_flotation_thickness
    public :: plastic_volume, plastic_area, plastic_surface_misfit, plastic_residual
    public :: PLASTIC_OK, PLASTIC_NOT_CONVERGED, PLASTIC_BAD_INPUT

end module plastic_ice
