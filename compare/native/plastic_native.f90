! Native-grid check against ICESHEET (PaleoMIST 1.0, 5 km North America; see native.jl):
! plastic_native.x <dir> <mode> <margin> <rho_i> <g> reads <dir>/hdr.txt (nx ny dx) and the
! raw Float64 (nx,ny) fields bed.bin, tau.bin, mask.bin (ice where > 0.5), optionally zsl.bin
! and, for margin = "sdist", sdist.bin (signed distance to the margin [m], positive in ice).
! margin: "centre" | "face" (par%margin) | "sdist" (margin_frac from
! plastic_margin_frac_from_sdist). Writes <dir>/H_<mode>_<margin>.bin.
! Adapted from /work/ba1442/robinson/plastic-diag/drv/plastic_raw.f90.
program plastic_native
    use plastic_ice
    implicit none
    integer, parameter :: dp = plastic_dp
    character(len=512) :: dir, arg
    character(len=16) :: mode, margin
    integer :: nx, ny, stat, u, c0, c1, cr
    real(dp) :: dx, rho, g, rmax
    real(dp), allocatable :: zb(:,:), tau(:,:), m(:,:), zsl(:,:), sd(:,:), mf(:,:)
    logical, allocatable :: mask(:,:)
    logical :: ex
    type(plastic_par_t) :: par
    type(plastic_t) :: ps
    call get_command_argument(1, dir); call get_command_argument(2, mode); call get_command_argument(3, margin)
    call get_command_argument(4, arg); read(arg, *) rho
    call get_command_argument(5, arg); read(arg, *) g
    open(newunit=u, file=trim(dir)//"/hdr.txt", status="old"); read(u, *) nx, ny, dx; close(u)
    allocate(zb(nx,ny), tau(nx,ny), m(nx,ny), zsl(nx,ny), mask(nx,ny), sd(nx,ny), mf(nx,ny))
    call rd("bed.bin", zb); call rd("tau.bin", tau); call rd("mask.bin", m)
    zsl = 0.0_dp
    inquire(file=trim(dir)//"/zsl.bin", exist=ex); if (ex) call rd("zsl.bin", zsl)
    mask = m > 0.5_dp
    par%rho_i = rho; par%g = g; par%mode = trim(mode); par%max_sweeps = 2000
    if (trim(margin) == "centre" .or. trim(margin) == "face") par%margin = trim(margin)
    call plastic_init(ps, nx, ny, dx, dx, par, stat=stat)
    call system_clock(c0, cr)
    if (trim(margin) == "sdist") then
        call rd("sdist.bin", sd)
        call plastic_margin_frac_from_sdist(sd, mask, dx, dx, mf)
        call plastic_solve(ps, zb, zsl, tau, mask, stat=stat, margin_frac=mf)
        call system_clock(c1)
        call plastic_stress_residual(ps, zb, zsl, tau, mask, rmax, margin_frac=mf)
    else
        call plastic_solve(ps, zb, zsl, tau, mask, stat=stat)
        call system_clock(c1)
        call plastic_stress_residual(ps, zb, zsl, tau, mask, rmax)
    end if
    print "(a,a,a,i3,a,i6,a,l2,a,f8.3,a,a)", trim(margin), ": stat ", "", stat, " sweeps ", ps%n_sweeps, &
        " conv ", ps%converged, " t ", real(c1 - c0, dp) / cr, " s ", trim(ps%msg)
    open(newunit=u, file=trim(dir)//"/H_"//trim(mode)//"_"//trim(margin)//".bin", access="stream", &
         form="unformatted", status="replace")
    write(u) ps%H; close(u)
contains
    subroutine rd(f, a)
        character(len=*), intent(in) :: f
        real(dp), intent(out) :: a(:,:)
        integer :: v
        open(newunit=v, file=trim(dir)//"/"//f, access="stream", form="unformatted", status="old"); read(v) a; close(v)
    end subroutine rd
end program plastic_native
