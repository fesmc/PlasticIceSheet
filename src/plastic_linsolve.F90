! Sparse linear solves for the steady SMB-driven solver: a CSR matrix (1-based, sorted
! columns, diagonal present) solved by BiCGStab with ILU(0) preconditioning, either
! built in (no dependency) or through Lis (build with lis=1, -DPLASTIC_LIS).
!
! The built-in solver is deterministic for any number of OpenMP threads: row-parallel
! mat-vec, dot products summed in fixed blocks, sequential triangular solves.
module plastic_linsolve

    use plastic_kinds, only: dp

    implicit none

    private
    public :: plastic_lin_t, lin_reserve, lin_solve, lin_has_lis

    integer, parameter :: DOT_BLOCK = 4096

    ! Matrix and work vectors (capacity nmax rows, nzmax entries; reused across solves).
    type plastic_lin_t
        integer :: n = 0, nnz = 0, nmax = 0, nzmax = 0
        integer,  allocatable :: ia(:), ja(:), idiag(:)
        real(dp), allocatable :: a(:), lu(:), b(:), x(:)
        real(dp), allocatable :: r(:), r0(:), p(:), v(:), s(:), t(:), y(:), z(:), blk(:)
        integer,  allocatable :: iw(:)
        integer  :: iters = 0
        real(dp) :: relres = 0.0_dp
    end type plastic_lin_t

contains

    logical function lin_has_lis()
#ifdef PLASTIC_LIS
        lin_has_lis = .true.
#else
        lin_has_lis = .false.
#endif
    end function lin_has_lis

    ! Make room for n rows and nz entries (grows only).
    subroutine lin_reserve(L, n, nz)
        type(plastic_lin_t), intent(inout) :: L
        integer,             intent(in)    :: n, nz
        integer :: nb
        if (n > L%nmax) then
            if (allocated(L%ia)) deallocate(L%ia, L%idiag, L%b, L%x, L%r, L%r0, L%p, L%v, L%s, &
                                            L%t, L%y, L%z, L%iw, L%blk)
            L%nmax = n
            nb = (n + DOT_BLOCK - 1) / DOT_BLOCK
            allocate(L%ia(n + 1), L%idiag(n), L%b(n), L%x(n), L%r(n), L%r0(n), L%p(n), L%v(n), &
                     L%s(n), L%t(n), L%y(n), L%z(n), L%iw(n), L%blk(nb))
        end if
        if (nz > L%nzmax) then
            if (allocated(L%ja)) deallocate(L%ja, L%a, L%lu)
            L%nzmax = nz
            allocate(L%ja(nz), L%a(nz), L%lu(nz))
        end if
    end subroutine lin_reserve

    ! Solve A x = b (L%ia, L%ja, L%a, L%b of size L%n, initial guess L%x) to the relative
    ! residual tol. use_lis selects Lis when compiled in. ok: tolerance reached.
    subroutine lin_solve(L, tol, maxit, use_lis, ok)
        type(plastic_lin_t), intent(inout) :: L
        real(dp),            intent(in)    :: tol
        integer,             intent(in)    :: maxit
        logical,             intent(in)    :: use_lis
        logical,             intent(out)   :: ok
#ifdef PLASTIC_LIS
        if (use_lis) then
            call solve_lis(L, tol, maxit, ok)
            return
        end if
#else
        if (use_lis) continue                 ! Lis not compiled in: built-in solver
#endif
        call ilu0(L)
        call bicgstab(L, tol, maxit, ok)
    end subroutine lin_solve

    ! --- built-in BiCGStab with ILU(0) ------------------------------------------------

    ! ILU(0) factors of A in L%lu (unit lower L, upper U on the same pattern).
    subroutine ilu0(L)
        type(plastic_lin_t), intent(inout) :: L
        integer  :: i, k, kk, jw, jrow
        real(dp) :: piv

        L%lu(1:L%nnz) = L%a(1:L%nnz)
        L%iw(1:L%n) = 0
        do i = 1, L%n
            do k = L%ia(i), L%ia(i + 1) - 1
                L%iw(L%ja(k)) = k
                if (L%ja(k) == i) L%idiag(i) = k
            end do
            do k = L%ia(i), L%ia(i + 1) - 1
                jrow = L%ja(k)
                if (jrow >= i) exit
                piv = L%lu(k) / L%lu(L%idiag(jrow))
                L%lu(k) = piv
                do kk = L%idiag(jrow) + 1, L%ia(jrow + 1) - 1
                    jw = L%iw(L%ja(kk))
                    if (jw /= 0) L%lu(jw) = L%lu(jw) - piv * L%lu(kk)
                end do
            end do
            do k = L%ia(i), L%ia(i + 1) - 1
                L%iw(L%ja(k)) = 0
            end do
            if (L%lu(L%idiag(i)) == 0.0_dp) L%lu(L%idiag(i)) = 1.0e-30_dp
        end do
    end subroutine ilu0

    ! z = (LU)^-1 y
    subroutine ilu_apply(L, y, z)
        type(plastic_lin_t), intent(in)  :: L
        real(dp),            intent(in)  :: y(:)
        real(dp),            intent(out) :: z(:)
        integer  :: i, k
        real(dp) :: acc
        do i = 1, L%n
            acc = y(i)
            do k = L%ia(i), L%idiag(i) - 1
                acc = acc - L%lu(k) * z(L%ja(k))
            end do
            z(i) = acc
        end do
        do i = L%n, 1, -1
            acc = z(i)
            do k = L%idiag(i) + 1, L%ia(i + 1) - 1
                acc = acc - L%lu(k) * z(L%ja(k))
            end do
            z(i) = acc / L%lu(L%idiag(i))
        end do
    end subroutine ilu_apply

    subroutine matvec(L, x, y)
        type(plastic_lin_t), intent(in)  :: L
        real(dp),            intent(in)  :: x(:)
        real(dp),            intent(out) :: y(:)
        integer  :: i, k
        real(dp) :: acc
        !$omp parallel do private(i, k, acc) schedule(static)
        do i = 1, L%n
            acc = 0.0_dp
            do k = L%ia(i), L%ia(i + 1) - 1
                acc = acc + L%a(k) * x(L%ja(k))
            end do
            y(i) = acc
        end do
        !$omp end parallel do
    end subroutine matvec

    ! Dot product summed in fixed blocks (independent of the thread count).
    function dot(L, x, y) result(d)
        type(plastic_lin_t), intent(inout) :: L
        real(dp),            intent(in)    :: x(:), y(:)
        real(dp) :: d, acc
        integer  :: ib, nb, i, i0, i1
        nb = (L%n + DOT_BLOCK - 1) / DOT_BLOCK
        !$omp parallel do private(ib, i, i0, i1, acc) schedule(static)
        do ib = 1, nb
            i0 = (ib - 1) * DOT_BLOCK + 1
            i1 = min(ib * DOT_BLOCK, L%n)
            acc = 0.0_dp
            do i = i0, i1
                acc = acc + x(i) * y(i)
            end do
            L%blk(ib) = acc
        end do
        !$omp end parallel do
        d = 0.0_dp
        do ib = 1, nb
            d = d + L%blk(ib)
        end do
    end function dot

    ! Right-preconditioned BiCGStab (M = ILU(0)) from the initial guess in L%x.
    subroutine bicgstab(L, tol, maxit, ok)
        type(plastic_lin_t), intent(inout) :: L
        real(dp),            intent(in)    :: tol
        integer,             intent(in)    :: maxit
        logical,             intent(out)   :: ok
        integer  :: it, n
        real(dp) :: rho, rho1, alpha, omega, beta, bnorm, rnorm, tt

        n = L%n
        ok = .false.
        L%iters = 0
        call matvec(L, L%x, L%r)
        L%r(1:n) = L%b(1:n) - L%r(1:n)
        bnorm = sqrt(dot(L, L%b, L%b))
        if (bnorm == 0.0_dp) then
            L%x(1:n) = 0.0_dp
            L%relres = 0.0_dp
            ok = .true.
            return
        end if
        rnorm = sqrt(dot(L, L%r, L%r))
        L%relres = rnorm / bnorm
        if (L%relres <= tol) then
            ok = .true.
            return
        end if
        L%r0(1:n) = L%r(1:n)
        L%p(1:n)  = 0.0_dp
        L%v(1:n)  = 0.0_dp
        rho = 1.0_dp; alpha = 1.0_dp; omega = 1.0_dp
        do it = 1, maxit
            rho1 = dot(L, L%r0, L%r)
            if (rho1 == 0.0_dp) exit
            beta = (rho1 / rho) * (alpha / omega)
            L%p(1:n) = L%r(1:n) + beta * (L%p(1:n) - omega * L%v(1:n))
            call ilu_apply(L, L%p, L%y)
            call matvec(L, L%y, L%v)
            alpha = rho1 / dot(L, L%r0, L%v)
            L%s(1:n) = L%r(1:n) - alpha * L%v(1:n)
            call ilu_apply(L, L%s, L%z)
            call matvec(L, L%z, L%t)
            tt = dot(L, L%t, L%t)
            if (tt > 0.0_dp) then
                omega = dot(L, L%t, L%s) / tt
            else
                omega = 0.0_dp
            end if
            L%x(1:n) = L%x(1:n) + alpha * L%y(1:n) + omega * L%z(1:n)
            L%r(1:n) = L%s(1:n) - omega * L%t(1:n)
            rho = rho1
            rnorm = sqrt(dot(L, L%r, L%r))
            L%iters = it
            L%relres = rnorm / bnorm
            if (L%relres <= tol) then
                ok = .true.
                exit
            end if
            if (omega == 0.0_dp) exit
        end do
    end subroutine bicgstab

#ifdef PLASTIC_LIS
    ! Lis: BiCGStab + ILU(0), as yelmo's solver_linear (0-based CSR).
    subroutine solve_lis(L, tol, maxit, ok)
        type(plastic_lin_t), intent(inout) :: L
        real(dp),            intent(in)    :: tol
        integer,             intent(in)    :: maxit
        logical,             intent(out)   :: ok
#include "lisf.h"
        LIS_INTEGER :: ierr, nr, nnz, iter, status
        LIS_MATRIX  :: A
        LIS_VECTOR  :: b, x
        LIS_SOLVER  :: solver
        LIS_REAL    :: resid
        LIS_INTEGER, allocatable :: ptr(:), idx(:)
        LIS_SCALAR,  allocatable :: val(:), bv(:), xv(:)
        character(len=128) :: opt
        integer :: k

        nr  = L%n
        nnz = L%nnz
        allocate(ptr(nr + 1), idx(nnz), val(nnz), bv(nr), xv(nr))
        !$omp parallel private(k)
        !$omp do
        do k = 1, nr + 1
            ptr(k) = L%ia(k) - 1
        end do
        !$omp end do nowait
        !$omp do
        do k = 1, nnz
            idx(k) = L%ja(k) - 1
            val(k) = L%a(k)
        end do
        !$omp end do nowait
        !$omp do
        do k = 1, nr
            bv(k) = L%b(k)
            xv(k) = L%x(k)
        end do
        !$omp end do
        !$omp end parallel

        call lis_initialize(ierr)
        call lis_matrix_create(LIS_COMM_WORLD, A, ierr)
        call lis_vector_create(LIS_COMM_WORLD, b, ierr)
        call lis_vector_create(LIS_COMM_WORLD, x, ierr)
        call lis_matrix_set_size(A, 0, nr, ierr)
        call lis_vector_set_size(b, 0, nr, ierr)
        call lis_vector_set_size(x, 0, nr, ierr)
        call lis_matrix_set_csr(nnz, ptr, idx, val, A, ierr)
        call lis_matrix_assemble(A, ierr)
        call lis_vector_scatter(bv, b, ierr)
        call lis_vector_scatter(xv, x, ierr)
        call lis_solver_create(solver, ierr)
        write(opt, '(a,es10.3,a,i0,a)') "-i bicgstab -p ilu -tol ", tol, " -maxiter ", maxit, &
            " -initx_zeros false"
        call lis_solver_set_option(trim(opt), solver, ierr)
        call lis_solve(A, b, x, solver, ierr)
        call lis_solver_get_iter(solver, iter, ierr)
        call lis_solver_get_status(solver, status, ierr)
        call lis_solver_get_residualnorm(solver, resid, ierr)
        call lis_vector_gather(x, xv, ierr)
        L%x(1:nr) = xv(1:nr)
        L%iters  = iter
        L%relres = resid
        ok = (status == 0)
        call lis_matrix_destroy(A, ierr)
        call lis_vector_destroy(b, ierr)
        call lis_vector_destroy(x, ierr)
        call lis_solver_destroy(solver, ierr)
        call lis_finalize(ierr)
    end subroutine solve_lis
#endif

end module plastic_linsolve
