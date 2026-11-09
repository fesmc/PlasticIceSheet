! Prototype of a direct formulation of the SMB-driven reconstruction: steady mass
! conservation div q = smb on the fixed mask, with the membrane-free DIVA flux of the
! closures, q = -H (u_b(tau) + u_def(tau, H)) grad z_s / |grad z_s|, tau = rho_i g H |grad z_s|
! (steady SIA with sliding), solved by explicit pseudo-time stepping to steady state
! (staggered Mahaffy fluxes, H >= 0, ice-free cells fixed at z_b + H_flt).
!
! Usage: smb_sia_proto.x <outdir> <case> <forcing: file|implied> <max_steps> <tol [m/yr]>
program smb_sia_proto

    use ncio
    use plastic_ice

    implicit none

    integer, parameter :: dp = plastic_dp
    character(len=512) :: outdir, cname, forcing, arg
    character(len=1024) :: fin
    character(len=32) :: law
    type(plastic_t)         :: p0
    type(plastic_par_t)     :: par
    type(plastic_glen_t)    :: rheo
    type(plastic_sliding_t) :: slide
    real(dp), allocatable :: z_bed(:,:), z_sl(:,:), tau_b(:,:), smb(:,:), H(:,:), H0(:,:), s(:,:)
    real(dp), allocatable :: qx(:,:), qy(:,:), dHdt(:,:), ux(:,:), uy(:,:), speed(:,:)
    integer,  allocatable :: imask(:,:)
    logical,  allocatable :: mask(:,:)
    real(dp) :: dx, rg, tau_cap, dt, t, dmax, res, tol, Hf, sx, sy, Dmx
    integer  :: nx, ny, i, j, it, max_steps, c0, c1, cr

    call get_command_argument(1, outdir)
    call get_command_argument(2, cname)
    call get_command_argument(3, forcing)
    call get_command_argument(4, arg); read(arg, *) max_steps
    call get_command_argument(5, arg); read(arg, *) tol

    fin = trim(outdir)//"/in_"//trim(cname)//".nc"
    nx = nc_size(fin, "x"); ny = nc_size(fin, "y")
    allocate(z_bed(nx,ny), z_sl(nx,ny), tau_b(nx,ny), smb(nx,ny), imask(nx,ny), mask(nx,ny), H(nx,ny), &
             H0(nx,ny), s(nx,ny), qx(0:nx,ny), qy(nx,0:ny), dHdt(nx,ny), ux(nx,ny), uy(nx,ny), &
             speed(nx,ny))
    call nc_read(fin, "z_bed", z_bed); call nc_read(fin, "z_sl", z_sl); call nc_read(fin, "tau_b", tau_b)
    call nc_read(fin, "smb", smb); call nc_read(fin, "mask", imask); mask = imask == 1
    call nc_read_attr(fin, "dx", dx)
    call nc_read_attr(fin, "glen_A", rheo%A); call nc_read_attr(fin, "glen_n", rheo%n)
    call nc_read_attr(fin, "slide_law", law); slide%law = trim(law)
    call nc_read_attr(fin, "slide_beta", slide%beta); call nc_read_attr(fin, "slide_m", slide%m)
    call nc_read_attr(fin, "slide_C", slide%C); call nc_read_attr(fin, "slide_u0", slide%u_0)

    par%mode = "surface"
    call plastic_init(p0, nx, ny, dx, dx, par)
    call plastic_solve(p0, z_bed, z_sl, tau_b, mask)
    H0 = p0%H
    if (trim(forcing) == "implied") then
        call plastic_velocity(p0, tau_b, mask, rheo, slide, ux, uy, speed)
        where (speed /= speed) speed = 0.0_dp
        call plastic_smb_implied(p0, ux, uy, mask, smb, speed=speed)
    end if

    rg = par%rho_i * par%g
    tau_cap = huge(1.0_dp)
    if (trim(law) == "coulomb") tau_cap = 0.999_dp * slide%C
    ! start: the plastic solution at 5e4 Pa (as the cold start of plastic_solve_smb)
    call plastic_solve(p0, z_bed, z_sl, 0.0_dp * tau_b + 5.0e4_dp, mask)
    H = p0%H
    do j = 1, ny
        do i = 1, nx
            if (.not. mask(i, j)) H(i, j) = plastic_flotation_thickness(z_bed(i, j), z_sl(i, j), par)
        end do
    end do

    call system_clock(c0, cr)
    t = 0.0_dp
    do it = 1, max_steps
        s = z_bed + H
        ! x faces (i+1/2, j), i = 1..nx-1; Dmx = max diffusivity for the time step
        Dmx = 0.0_dp
        qx = 0.0_dp; qy = 0.0_dp
        do j = 1, ny
            do i = 1, nx - 1
                if (.not. (mask(i, j) .or. mask(i + 1, j))) cycle
                Hf = 0.5_dp * (H(i, j) + H(i + 1, j))
                sx = (s(i + 1, j) - s(i, j)) / dx
                sy = 0.25_dp * (s(i, min(j+1,ny)) - s(i, max(j-1,1)) + s(i+1, min(j+1,ny)) - s(i+1, max(j-1,1))) / dx
                call face_flux(Hf, sx, sy, qx(i, j), Dmx)
            end do
        end do
        do j = 1, ny - 1
            do i = 1, nx
                if (.not. (mask(i, j) .or. mask(i, j + 1))) cycle
                Hf = 0.5_dp * (H(i, j) + H(i, j + 1))
                sy = (s(i, j + 1) - s(i, j)) / dx
                sx = 0.25_dp * (s(min(i+1,nx), j) - s(max(i-1,1), j) + s(min(i+1,nx), j+1) - s(max(i-1,1), j+1)) / dx
                call face_flux(Hf, sy, sx, qy(i, j), Dmx)
            end do
        end do
        dt = 0.2_dp * dx * dx / max(Dmx, 1.0e-30_dp)
        dt = min(dt, 10.0_dp)
        res = 0.0_dp
        do j = 1, ny
            do i = 1, nx
                if (.not. mask(i, j)) cycle
                dHdt(i, j) = smb(i, j) - (qx(i, j) - qx(i - 1, j)) / dx - (qy(i, j) - qy(i, j - 1)) / dx
                if (H(i, j) <= 0.0_dp .and. dHdt(i, j) < 0.0_dp) dHdt(i, j) = 0.0_dp
                res = max(res, abs(dHdt(i, j)))
                H(i, j) = max(H(i, j) + dt * dHdt(i, j), 0.0_dp)
            end do
        end do
        t = t + dt
        if (it == 1 .or. mod(it, max(max_steps / 10, 1)) == 0) &
            write(*, '(a,i8,a,f10.1,a,es9.2,a,es10.3,a)') "  step ", it, ", t ", t, " yr, dt ", dt, &
                " yr, max|dH/dt| ", res, " m/yr"
        if (res < tol) exit
    end do
    call system_clock(c1)
    dmax = maxval(abs(H - H0), mask=mask)
    write(*, '(a,a,a,a,a,i0,a,f10.1,a,es10.3,a,f8.2,a)') trim(cname), " (", trim(forcing), &
        "): steps ", "", it, ", t ", t, " yr, max|dH/dt| ", res, " m/yr, ", real(c1 - c0, dp) / cr, " s"
    write(*, '(a,f8.1,a,f8.1,a,f8.1,a,i0,a,i0)') "  vs plastic (surface, tau_b): max|dH| ", dmax, " m, rms ", &
        sqrt(sum((H - H0)**2, mask=mask) / count(mask)), " m, mean ", sum(H - H0, mask=mask) / count(mask), &
        " m; ice-free (H = 0) cells in the mask ", count(mask .and. H <= 0.0_dp), " of ", count(mask)

contains

    ! Flux normal to a face with normal slope sn and tangential slope st.
    subroutine face_flux(Hf, sn, st, q, Dmx)
        real(dp), intent(in)    :: Hf, sn, st
        real(dp), intent(out)   :: q
        real(dp), intent(inout) :: Dmx
        real(dp) :: gs, tau, u, D
        gs = sqrt(sn * sn + st * st)
        q = 0.0_dp
        if (Hf <= 0.0_dp .or. gs <= 0.0_dp) return
        tau = min(rg * Hf * gs, tau_cap)
        u = plastic_basal_velocity(slide, tau) + plastic_deformational_velocity(rheo, tau, Hf)
        D = Hf * u / gs
        q = -D * sn
        Dmx = max(Dmx, D)
    end subroutine face_flux

end program smb_sia_proto
