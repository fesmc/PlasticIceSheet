! Diagnostics of the SMB-driven reconstruction (plastic_solve_smb) outer iteration
! "route, invert, sweep, relax": step sequence max|T(H) - H|, contraction factor, and
! where the un-relaxed map T moves H (ablation cells with Q <= 0, saturated tau_d,
! margin cells, interior). With forcing "implied" the plastic solution H0 (inner mode)
! is an exact fixed point; its linear stability is measured by power iteration on the
! finite-difference Jacobian of T at H0, and by the damped iteration from H0 + 1 m.
!
! Usage: smb_diag.x <outdir> <case> <forcing> <mode> <n_iter> [relax]
!   forcing: "file" (smb of in_<case>.nc) | "implied" (closure-B SMB of the surface-mode
!            plastic solution at the case's tau_b, with its rheology and sliding law)
!   mode:    inner plastic solve, "surface" | "flat"
program smb_diag

    use ncio
    use plastic_ice

    implicit none

    integer, parameter :: dp = plastic_dp
    character(len=512) :: outdir, cname, forcing, mode, arg
    character(len=1024) :: fin
    character(len=32) :: law
    type(plastic_t)          :: ps, p0
    type(plastic_par_t)      :: par
    type(plastic_glen_t)     :: rheo
    type(plastic_sliding_t)  :: slide
    type(plastic_smb_ctrl_t) :: ctrl
    real(dp), allocatable :: z_bed(:,:), z_sl(:,:), tau_b(:,:), smb(:,:), H(:,:), H0(:,:), tau_d(:,:)
    real(dp), allocatable :: ux(:,:), uy(:,:), speed(:,:), Q(:,:), d(:,:), hist(:)
    integer,  allocatable :: imask(:,:)
    logical,  allocatable :: mask(:,:), margin(:,:), abl(:,:), sat(:,:)
    real(dp) :: dx, t0, t1, rho, hi, dmax, relax
    integer  :: nx, ny, n_iter, it, i, j, c0, c1, cr, loc(2)

    call get_command_argument(1, outdir)
    call get_command_argument(2, cname)
    call get_command_argument(3, forcing)
    call get_command_argument(4, mode)
    call get_command_argument(5, arg)
    read(arg, *) n_iter
    relax = 0.05_dp
    if (command_argument_count() >= 6) then
        call get_command_argument(6, arg)
        read(arg, *) relax
    end if

    fin = trim(outdir)//"/in_"//trim(cname)//".nc"
    nx = nc_size(fin, "x"); ny = nc_size(fin, "y")
    allocate(z_bed(nx,ny), z_sl(nx,ny), tau_b(nx,ny), smb(nx,ny), imask(nx,ny), mask(nx,ny), H(nx,ny), &
             H0(nx,ny), tau_d(nx,ny), ux(nx,ny), uy(nx,ny), speed(nx,ny), Q(nx,ny), d(nx,ny), &
             margin(nx,ny), abl(nx,ny), sat(nx,ny), hist(n_iter))
    call nc_read(fin, "z_bed", z_bed); call nc_read(fin, "z_sl", z_sl); call nc_read(fin, "tau_b", tau_b)
    call nc_read(fin, "smb", smb); call nc_read(fin, "mask", imask); mask = imask == 1
    call nc_read_attr(fin, "dx", dx)
    call nc_read_attr(fin, "glen_A", rheo%A); call nc_read_attr(fin, "glen_n", rheo%n)
    call nc_read_attr(fin, "slide_law", law); slide%law = trim(law)
    call nc_read_attr(fin, "slide_beta", slide%beta); call nc_read_attr(fin, "slide_m", slide%m)
    call nc_read_attr(fin, "slide_C", slide%C); call nc_read_attr(fin, "slide_u0", slide%u_0)

    ! reference plastic solution (inner mode) and, for "implied", its own SMB
    par%mode = mode
    call plastic_init(p0, nx, ny, dx, dx, par)
    call plastic_solve(p0, z_bed, z_sl, tau_b, mask)
    H0 = p0%H
    if (trim(forcing) == "implied") then
        call plastic_velocity(p0, tau_b, mask, rheo, slide, ux, uy, speed)
        where (speed /= speed) speed = 0.0_dp          ! coulomb at the cap
        call plastic_smb_implied(p0, ux, uy, mask, smb, speed=speed)
    end if

    margin = .false.
    do j = 1, ny
        do i = 1, nx
            if (.not. mask(i, j)) cycle
            margin(i, j) = .not. (mask(max(i-1,1), j) .and. mask(min(i+1,nx), j) .and. &
                                  mask(i, max(j-1,1)) .and. mask(i, min(j+1,ny)))
        end do
    end do

    write(*, '(a,a,a,a,a,a,a,i0,a)') trim(cname), ", forcing ", trim(forcing), ", inner ", trim(mode), &
        ", law ", trim(law)//", ", count(mask), " ice cells"

    ! one outer iteration per call (ctrl%n_outer = 1, H_init = previous H)
    ctrl%mode = mode
    ctrl%relax = relax
    ctrl%n_outer = 1
    ctrl%tol = 0.0_dp
    call plastic_init(ps, nx, ny, dx, dx)
    call system_clock(c0, cr)
    call plastic_solve_smb(ps, z_bed, z_sl, smb, mask, rheo, slide, ctrl=ctrl, tau_d=tau_d)
    hist(1) = ps%outer_change
    do it = 2, n_iter
        H = ps%H
        call plastic_solve_smb(ps, z_bed, z_sl, smb, mask, rheo, slide, ctrl=ctrl, H_init=H, tau_d=tau_d)
        hist(it) = ps%outer_change
    end do
    call system_clock(c1)
    t1 = real(c1 - c0, dp) / cr
    do it = 1, n_iter
        if (it == 1 .or. mod(it, max(n_iter / 10, 1)) == 0) &
            write(*, '(a,i6,a,es10.3,a)') "  it ", it, "  max|T(H) - H| ", hist(it), " m"
    end do
    rho = (hist(n_iter) / hist(n_iter / 2)) ** (1.0_dp / real(n_iter - n_iter / 2, dp))
    write(*, '(a,f8.5,a,f8.2,a)') "  observed rate (second half) ", rho, ", ", t1, " s"
    write(*, '(a,es10.3,a,es10.3,a)') "  final max|H - H_plastic| ", maxval(abs(ps%H - H0), mask=mask), &
        " m, rms ", sqrt(sum((ps%H - H0)**2, mask=mask) / count(mask)), " m"

    ! where the un-relaxed map moves H: one more step with relax = 1 from H
    H = ps%H
    ctrl%relax = 1.0_dp
    call plastic_solve_smb(ps, z_bed, z_sl, smb, mask, rheo, slide, ctrl=ctrl, H_init=H, tau_d=tau_d)
    d = ps%H - H
    call plastic_balance_flux(merge(H + z_bed, z_bed, mask), mask, smb, dx, dx, Q)
    abl = mask .and. Q <= 0.0_dp
    hi = ctrl%tau_hi
    if (trim(law) == "coulomb") hi = min(hi, 0.999_dp * slide%C)
    sat = mask .and. tau_d >= hi * (1.0_dp - 1.0e-9_dp)
    loc = maxloc(abs(d), mask=mask)
    write(*, '(a,es10.3,a,i0,a,i0,a,f7.0,a,es10.3,a,f9.0,a,l1,a,l1,a,l1)') "  T(H) - H: max ", maxval(abs(d), mask=mask), &
        " m at (", loc(1), ",", loc(2), "), H ", H(loc(1), loc(2)), " m, Q ", Q(loc(1), loc(2)), &
        " m3/yr, tau_d ", tau_d(loc(1), loc(2)), " Pa; ablation ", abl(loc(1), loc(2)), ", saturated ", &
        sat(loc(1), loc(2)), ", margin ", margin(loc(1), loc(2))
    call cls("ablation (Q <= 0, tau_d = 0)", abl)
    call cls("saturated tau_d", sat .and. .not. abl)
    call cls("margin (not ablation/saturated)", margin .and. .not. abl .and. .not. sat)
    call cls("interior (none of the above)", mask .and. .not. margin .and. .not. abl .and. .not. sat)

    if (trim(forcing) == "implied") call fixed_point_stability()

contains

    ! T(Hx): one un-relaxed outer iteration from Hx.
    subroutine apply_T(Hx, THx)
        real(dp), intent(in)  :: Hx(:,:)
        real(dp), intent(out) :: THx(:,:)
        type(plastic_smb_ctrl_t) :: c1
        c1%mode = mode; c1%n_outer = 1; c1%tol = 0.0_dp; c1%relax = 1.0_dp
        call plastic_solve_smb(ps, z_bed, z_sl, smb, mask, rheo, slide, ctrl=c1, H_init=Hx)
        THx = ps%H
    end subroutine apply_T

    subroutine fixed_point_stability()
        real(dp), allocatable :: T0(:,:), v(:,:), Jv(:,:), Hk(:,:), r(:,:)
        real(dp) :: eps, nv, lam, ray
        integer  :: k
        type(plastic_smb_ctrl_t) :: c1
        allocate(T0(nx,ny), v(nx,ny), Jv(nx,ny), Hk(nx,ny), r(nx,ny))
        call apply_T(H0, T0)
        write(*, '(a,es10.3,a)') "  fixed point: max|T(H0) - H0| ", maxval(abs(T0 - H0), mask=mask), " m"
        ! power iteration on J = dT/dH at H0
        eps = 0.01_dp
        call random_number(r)
        v = merge(r - 0.5_dp, 0.0_dp, mask)
        v = v / sqrt(sum(v**2))
        do k = 1, 40
            call apply_T(H0 + eps * v, Jv)
            Jv = merge((Jv - T0) / eps, 0.0_dp, mask)
            nv  = sqrt(sum(Jv**2))
            ray = sum(v * Jv)
            v = Jv / max(nv, tiny(1.0_dp))
            if (k == 1 .or. mod(k, 10) == 0) write(*, '(a,i3,a,es10.3,a,es10.3)') "    power it ", k, &
                ": |J v| ", nv, ", v.Jv ", ray
        end do
        lam = nv
        loc = maxloc(abs(v), mask=mask)
        write(*, '(a,es10.3,a,f7.3,a,i0,a,i0,a,l1,a)') "  dominant |lambda(J)| ", lam, ", sign(v.Jv)/|J v| ", &
            ray / max(nv, tiny(1.0_dp)), ", eigenvector peak at (", loc(1), ",", loc(2), "), margin ", &
            margin(loc(1), loc(2)), ""
        write(*, '(a,f6.3,a,f8.3,a,f8.4)') "  damped map eigenvalue 1 - r + r lambda, r = ", relax, ": ", &
            1.0_dp - relax + relax * ray, "; stable for r < 2/(1 - lambda) = ", 2.0_dp / (1.0_dp - ray)
        ! damped iteration from H0 + 1 m
        c1%mode = mode; c1%n_outer = 1; c1%tol = 0.0_dp; c1%relax = relax
        Hk = merge(H0 + 1.0_dp, H0, mask)
        do k = 1, 400
            call plastic_solve_smb(ps, z_bed, z_sl, smb, mask, rheo, slide, ctrl=c1, H_init=Hk)
            Hk = ps%H
            if (k == 1 .or. mod(k, 50) == 0) write(*, '(a,i4,a,es10.3,a)') "    from H0 + 1 m, it ", k, &
                ": max|H - H0| ", maxval(abs(Hk - H0), mask=mask), " m"
        end do
    end subroutine fixed_point_stability

    subroutine cls(label, m)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: m(:,:)
        if (count(m) == 0) then
            write(*, '(a,a,a)') "    ", label, ": 0 cells"
        else
            write(*, '(a,a,a,i6,a,es10.3,a,es10.3,a)') "    ", label, ": ", count(m), " cells, max|T(H)-H| ", &
                maxval(abs(d), mask=m), " m, rms ", sqrt(sum(d**2, mask=m) / count(m)), " m"
        end if
    end subroutine cls

end program smb_diag
