! Kinds and numeric constants shared by all PlasticIceSheet modules.
module plastic_kinds

    use, intrinsic :: iso_fortran_env, only: real64, int64

    implicit none

    private
    public :: dp, BIG, NAN_DP, MSG_LEN

    integer,  parameter :: dp = real64

    ! "No upwind information yet" in the eikonal sweep (the Julia code uses Inf).
    ! A finite sentinel keeps the scheme safe under fast-math compiler flags.
    real(dp), parameter :: BIG = 1.0e300_dp

    ! Quiet NaN from its bit pattern (ieee_value can trap under -ffpe-trap=invalid).
    real(dp), parameter :: NAN_DP = transfer(9221120237041090560_int64, 1.0_dp)

    integer,  parameter :: MSG_LEN = 256

end module plastic_kinds
