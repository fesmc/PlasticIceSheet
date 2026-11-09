! Thread scaling (run with OMP_NUM_THREADS = 1, 4, 8, 16; compare/bench_threads.sh):
! the surface-mode eikonal, serial vs wavefront sweep (best of 5, and bit-identity), on
! synth300 and nh32_20ka; the steady SMB-driven solve (built-in and Lis) on grl32 (MAR),
! synth150 (its input SMB) and nh32_0ka (the implied SMB of its plastic solution).
! Usage: bench_threads.x <outdir>
program bench_threads

    !$ use omp_lib
    use ncio
    use plastic_ice

    implicit none

    integer, parameter :: dp = plastic_dp
    character(len=512) :: outdir
    integer :: nt

    call get_command_argument(1, outdir)
    nt = 1
    !$ nt = omp_get_max_threads()
    write(*, '(a,i0)') "threads ", nt
    call bench_eikonal("synth300")
    call bench_eikonal("nh32_20ka")
    call bench_steady("grl32", "file")
    call bench_steady("synth150", "file")
    call bench_steady("nh32_0ka", "implied")

contains

    function wall_time() result(t)
        real(dp) :: t
        integer(8) :: c, r
        call system_clock(c, r)
        t = real(c, dp) / real(r, dp)
    end function wall_time

    subroutine read_case(name, nx, ny, dx, z_bed, z_sl, tau_b, smb, mask, rheo, slide)
        character(len=*),        intent(in)  :: name
        integer,                 intent(out) :: nx, ny
        real(dp),                intent(out) :: dx
        real(dp), allocatable,   intent(out) :: z_bed(:,:), z_sl(:,:), tau_b(:,:), smb(:,:)
        logical,  allocatable,   intent(out) :: mask(:,:)
        type(plastic_glen_t),    intent(out) :: rheo
        type(plastic_sliding_t), intent(out) :: slide
        character(len=1024) :: fin
        character(len=32)   :: law
        integer, allocatable :: im(:,:)
        fin = trim(outdir)//"/in_"//name//".nc"
        nx = nc_size(fin, "x"); ny = nc_size(fin, "y")
        allocate(z_bed(nx, ny), z_sl(nx, ny), tau_b(nx, ny), smb(nx, ny), mask(nx, ny), im(nx, ny))
        call nc_read(fin, "z_bed", z_bed); call nc_read(fin, "z_sl", z_sl); call nc_read(fin, "tau_b", tau_b)
        call nc_read(fin, "smb", smb); call nc_read(fin, "mask", im); mask = im == 1
        call nc_read_attr(fin, "dx", dx)
        call nc_read_attr(fin, "glen_A", rheo%A); call nc_read_attr(fin, "glen_n", rheo%n)
        call nc_read_attr(fin, "slide_law", law); slide%law = trim(law)
        call nc_read_attr(fin, "slide_beta", slide%beta); call nc_read_attr(fin, "slide_m", slide%m)
        call nc_read_attr(fin, "slide_C", slide%C); call nc_read_attr(fin, "slide_u0", slide%u_0)
    end subroutine read_case

    subroutine bench_eikonal(name)
        character(len=*), intent(in) :: name
        type(plastic_t)         :: ps
        type(plastic_par_t)     :: par
        type(plastic_glen_t)    :: rheo
        type(plastic_sliding_t) :: slide
        real(dp), allocatable :: z_bed(:,:), z_sl(:,:), tau_b(:,:), smb(:,:), Hs(:,:)
        logical,  allocatable :: mask(:,:)
        real(dp) :: dx, t0, ts, tw
        integer  :: nx, ny, k

        call read_case(name, nx, ny, dx, z_bed, z_sl, tau_b, smb, mask, rheo, slide)
        par%mode = "surface"
        call plastic_init(ps, nx, ny, dx, dx, par)
        ts = huge(1.0_dp)
        do k = 1, 5
            t0 = wall_time(); call plastic_solve(ps, z_bed, z_sl, tau_b, mask); ts = min(ts, wall_time() - t0)
        end do
        Hs = ps%H
        ps%par%sweep = "wavefront"
        tw = huge(1.0_dp)
        do k = 1, 5
            t0 = wall_time(); call plastic_solve(ps, z_bed, z_sl, tau_b, mask); tw = min(tw, wall_time() - t0)
        end do
        write(*, '(a,a,a,i0,a,i0,a,f9.5,a,f9.5,a,l1)') "  eikonal ", name, " (", nx, "x", ny, "): serial ", ts, &
            " s, wavefront ", tw, " s, bit-identical ", all(ps%H == Hs)
        call plastic_end(ps)
    end subroutine bench_eikonal

    subroutine bench_steady(name, forcing)
        character(len=*), intent(in) :: name, forcing
        type(plastic_t)          :: ps, p0
        type(plastic_par_t)      :: par
        type(plastic_glen_t)     :: rheo
        type(plastic_sliding_t)  :: slide
        type(plastic_smb_ctrl_t) :: ctrl
        real(dp), allocatable :: z_bed(:,:), z_sl(:,:), tau_b(:,:), smb(:,:), Hb(:,:)
        logical,  allocatable :: mask(:,:)
        real(dp) :: dx, t0, t
        integer  :: nx, ny, st, k
        character(len=8) :: lins(2) = [character(len=8) :: "builtin", "lis"]

        call read_case(name, nx, ny, dx, z_bed, z_sl, tau_b, smb, mask, rheo, slide)
        if (forcing == "implied") then
            par%mode = "surface"
            call plastic_init(p0, nx, ny, dx, dx, par)
            call plastic_solve(p0, z_bed, z_sl, tau_b, mask)
            call plastic_smb_flux_divergence(p0, z_bed, z_sl, mask, rheo, slide, smb)
        end if
        do k = 1, 2
            if (k == 2 .and. .not. plastic_has_lis()) cycle
            ctrl%lin = lins(k)
            call plastic_init(ps, nx, ny, dx, dx)
            t0 = wall_time()
            call plastic_solve_smb(ps, z_bed, z_sl, smb, mask, rheo, slide, ctrl=ctrl, stat=st)
            t = wall_time() - t0
            write(*, '(a,a,a,a,a,a,a,i0,a,es9.2,a,f8.3,a,f8.1,a,i0,a)', advance="no") "  steady ", name, " (", forcing, &
                ", ", trim(lins(k)), "): stat ", st, ", residual ", ps%outer_change, " m/yr, ", t, " s, mean H ", &
                sum(ps%H, mask=mask) / count(mask), " m, steps ", ps%n_outer_used, ""
            if (k == 1) then
                Hb = ps%H
                if (forcing == "implied") write(*, '(a,es9.2,a)', advance="no") ", max|H - H_plastic| ", &
                    maxval(abs(ps%H - p0%H), mask=mask), " m"
            else
                write(*, '(a,es9.2,a)', advance="no") ", max|H_lis - H_builtin| ", maxval(abs(ps%H - Hb)), " m"
            end if
            write(*, '(a,i0,a,i0)') ", ice-free cells in mask ", count(mask .and. ps%H <= 0.0_dp), " of ", count(mask)
            call plastic_end(ps)
        end do
    end subroutine bench_steady

end program bench_threads
