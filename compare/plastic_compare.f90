! Fortran side of the comparison with PlasticIceSheet.jl (see compare/julia_ref.jl).
!
! For each case named on the command line, reads out/in_<case>.nc, runs the flat and
! full solves and the SMB-driven reconstruction with the same settings as the Julia
! run, times them, measures the warm-start speedup on a coupling-step-like change
! (1 m GIA-like bed change for the full solve, 1 % SMB change for the SMB-driven
! one), and writes out/fortran_<case>.nc.
!
! Usage: plastic_compare.x <outdir> <case> [<case> ...]
program plastic_compare

    use ncio
    use plastic_ice

    implicit none

    integer, parameter :: dp = plastic_dp
    character(len=512) :: outdir, case_name
    integer :: k

    call get_command_argument(1, outdir)
    do k = 2, command_argument_count()
        call get_command_argument(k, case_name)
        call run_case(trim(outdir), trim(case_name))
    end do

contains

    function wall_time() result(t)
        real(dp) :: t
        integer(8) :: c, r
        call system_clock(c, r)
        t = real(c, dp) / real(r, dp)
    end function wall_time

    subroutine run_case(outdir, name)
        character(len=*), intent(in) :: outdir, name
        character(len=1024) :: fin, fout
        character(len=32)   :: law
        type(plastic_t)          :: ps, pw, pc
        type(plastic_par_t)      :: par
        type(plastic_glen_t)     :: rheo
        type(plastic_sliding_t)  :: slide
        type(plastic_smb_ctrl_t) :: ctrl
        real(dp), allocatable :: x(:), y(:), z_bed(:,:), z_sl(:,:), tau_b(:,:), smb(:,:)
        real(dp), allocatable :: H_flat(:,:), H_full(:,:), H_smb(:,:), z_bed1(:,:), smb1(:,:)
        integer,  allocatable :: imask(:,:)
        logical,  allocatable :: mask(:,:)
        real(dp) :: dx, dy, relax_full, t_flat, t_full, t_smb, t0, tw, tc, tsw, tsc
        integer  :: nx, ny, n_outer_full, do_smb, rep, i, j
        integer  :: no_full, no_smb, nw, nc, nsw, nsc
        logical  :: cv_flat, cv_full, cv_smb

        fin  = outdir//"/in_"//name//".nc"
        fout = outdir//"/fortran_"//name//".nc"
        nx = nc_size(fin, "x")
        ny = nc_size(fin, "y")
        allocate(x(nx), y(ny), z_bed(nx, ny), z_sl(nx, ny), tau_b(nx, ny), smb(nx, ny), imask(nx, ny))
        allocate(mask(nx, ny), H_flat(nx, ny), H_full(nx, ny), H_smb(nx, ny), z_bed1(nx, ny), smb1(nx, ny))
        call nc_read(fin, "x", x)
        call nc_read(fin, "y", y)
        call nc_read(fin, "z_bed", z_bed)
        call nc_read(fin, "z_sl", z_sl)
        call nc_read(fin, "tau_b", tau_b)
        call nc_read(fin, "smb", smb)
        call nc_read(fin, "mask", imask)
        mask = imask == 1
        call nc_read_attr(fin, "dx", dx)
        call nc_read_attr(fin, "dy", dy)
        call nc_read_attr(fin, "relax_full", relax_full)
        call nc_read_attr(fin, "n_outer_full", n_outer_full)
        call nc_read_attr(fin, "do_smb", do_smb)
        call nc_read_attr(fin, "glen_A", rheo%A)
        call nc_read_attr(fin, "glen_n", rheo%n)
        call nc_read_attr(fin, "slide_law", law)
        call nc_read_attr(fin, "slide_beta", slide%beta)
        call nc_read_attr(fin, "slide_m", slide%m)
        call nc_read_attr(fin, "slide_C", slide%C)
        call nc_read_attr(fin, "slide_u0", slide%u_0)
        slide%law = trim(law)
        write(*, '(a,": ",i0," x ",i0,", dx = ",f0.0," m, ice cells ",i0)') name, nx, ny, dx, count(mask)

        ! flat
        par%mode = "flat"
        call plastic_init(ps, nx, ny, dx, dy, par)
        t_flat = huge(1.0_dp)
        do rep = 1, 3
            t0 = wall_time()
            call plastic_solve(ps, z_bed, z_sl, tau_b, mask)
            t_flat = min(t_flat, wall_time() - t0)
        end do
        H_flat = ps%H
        cv_flat = ps%converged

        ! full (cold every repetition)
        ps%par%mode = "full"
        ps%par%relax = relax_full
        ps%par%n_outer = n_outer_full
        t_full = huge(1.0_dp)
        do rep = 1, 3
            t0 = wall_time()
            call plastic_solve(ps, z_bed, z_sl, tau_b, mask, warm=.false.)
            t_full = min(t_full, wall_time() - t0)
        end do
        H_full = ps%H
        cv_full = ps%converged
        no_full = ps%n_outer_used
        write(*, '(a,f9.4,a,l1,a,f9.4,a,l1,a,i0)') "  flat ", t_flat, " s (conv ", cv_flat, &
            "), full ", t_full, " s (conv ", cv_full, ", n_outer ", no_full
        if (.not. cv_full) write(*, '(a,a,a,es9.2)') "  full: ", trim(ps%msg), ", last max|dH| ", ps%outer_change

        ! warm start on a coupling-step-like change: 1 m GIA-like bed lowering
        do j = 1, ny
            do i = 1, nx
                z_bed1(i, j) = z_bed(i, j) - exp(-((i - nx / 2.0_dp)**2 + (j - ny / 2.0_dp)**2) &
                                                 / (2.0_dp * (0.2_dp * max(nx, ny))**2))
            end do
        end do
        call plastic_init(pw, nx, ny, dx, dy, ps%par)
        call plastic_init(pc, nx, ny, dx, dy, ps%par)
        call plastic_solve(pw, z_bed, z_sl, tau_b, mask)
        t0 = wall_time()
        call plastic_solve(pw, z_bed1, z_sl, tau_b, mask, warm=.true.)
        tw = wall_time() - t0
        nw = pw%n_outer_used
        t0 = wall_time()
        call plastic_solve(pc, z_bed1, z_sl, tau_b, mask, warm=.false.)
        tc = wall_time() - t0
        nc = pc%n_outer_used
        write(*, '(a,i0,a,f8.4,a,i0,a,f8.4,a,f5.2,a,es9.2,a)') "  full warm (1 m bed change): n_outer ", nw, &
            " (", tw, " s) vs cold ", nc, " (", tc, " s), speedup ", tc / tw, "x, max|H_warm - H_cold| ", &
            maxval(abs(pw%H - pc%H)), " m"

        ! SMB-driven reconstruction (defaults as diva_reconstruct)
        t_smb = -1.0_dp; no_smb = 0; cv_smb = .false.
        tsw = 0.0_dp; tsc = 0.0_dp; nsw = 0; nsc = 0
        H_smb = -1.0_dp
        if (do_smb == 1) then
            call plastic_init(ps, nx, ny, dx, dy)
            t0 = wall_time()
            call plastic_solve_smb(ps, z_bed, z_sl, smb, mask, rheo, slide, ctrl=ctrl)
            t_smb = wall_time() - t0
            H_smb = ps%H
            no_smb = ps%n_outer_used
            cv_smb = ps%converged
            write(*, '(a,f9.4,a,l1,a,i0)') "  smb  ", t_smb, " s (conv ", cv_smb, ", n_outer ", no_smb

            smb1 = 1.01_dp * smb
            t0 = wall_time()
            call plastic_solve_smb(ps, z_bed, z_sl, smb1, mask, rheo, slide, warm=.true., ctrl=ctrl)
            tsw = wall_time() - t0
            nsw = ps%n_outer_used
            call plastic_init(pc, nx, ny, dx, dy)
            t0 = wall_time()
            call plastic_solve_smb(pc, z_bed, z_sl, smb1, mask, rheo, slide, ctrl=ctrl)
            tsc = wall_time() - t0
            nsc = pc%n_outer_used
            write(*, '(a,i0,a,f8.4,a,i0,a,f8.4,a,f5.2,a,es9.2,a)') "  smb warm (1 % smb change): n_outer ", nsw, &
                " (", tsw, " s) vs cold ", nsc, " (", tsc, " s), speedup ", tsc / tsw, &
                "x, max|H_warm - H_cold| ", maxval(abs(ps%H - pc%H)), " m"
        end if

        ! output
        call nc_create(fout, overwrite=.true.)
        call nc_write_dim(fout, "x", x)
        call nc_write_dim(fout, "y", y)
        call nc_write(fout, "H_flat", H_flat, dim1="x", dim2="y")
        call nc_write(fout, "H_full", H_full, dim1="x", dim2="y")
        call nc_write(fout, "H_smb",  H_smb,  dim1="x", dim2="y")
        call nc_write_attr(fout, "t_flat", t_flat)
        call nc_write_attr(fout, "t_full", t_full)
        call nc_write_attr(fout, "t_smb", t_smb)
        call nc_write_attr(fout, "n_outer_full", no_full)
        call nc_write_attr(fout, "n_outer_smb", no_smb)
        call nc_write_attr(fout, "converged_full", merge(1, 0, cv_full))
        call nc_write_attr(fout, "converged_smb", merge(1, 0, cv_smb))
        call nc_write_attr(fout, "warm_full_n_outer", nw)
        call nc_write_attr(fout, "cold_full_n_outer", nc)
        call nc_write_attr(fout, "warm_full_t", tw)
        call nc_write_attr(fout, "cold_full_t", tc)
        call nc_write_attr(fout, "warm_smb_n_outer", nsw)
        call nc_write_attr(fout, "cold_smb_n_outer", nsc)
        call nc_write_attr(fout, "warm_smb_t", tsw)
        call nc_write_attr(fout, "cold_smb_t", tsc)
        call plastic_end(ps); call plastic_end(pw); call plastic_end(pc)
    end subroutine run_case

end program plastic_compare
