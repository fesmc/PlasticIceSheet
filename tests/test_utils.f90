! Minimal check/report helpers for the test programs (test-only module state).
module plastic_test_utils

    use plastic_ice, only: dp => plastic_dp

    implicit none

    private
    public :: dp, check, check_close, finish, median, sort_vec, wall_time

    integer :: n_pass = 0, n_fail = 0

contains

    subroutine check(cond, name)
        logical,          intent(in) :: cond
        character(len=*), intent(in) :: name
        if (cond) then
            n_pass = n_pass + 1
            write(*, '(a,a)') "  PASS  ", name
        else
            n_fail = n_fail + 1
            write(*, '(a,a)') "  FAIL  ", name
        end if
    end subroutine check

    ! |a - b| <= max(atol, rtol |b|)
    subroutine check_close(a, b, rtol, atol, name)
        real(dp),         intent(in) :: a, b, rtol, atol
        character(len=*), intent(in) :: name
        logical :: ok
        ok = abs(a - b) <= max(atol, rtol * abs(b))
        if (.not. ok) write(*, '(a,es24.16,a,es24.16)') "        got ", a, "  expected ", b
        call check(ok, name)
    end subroutine check_close

    subroutine finish(prog)
        character(len=*), intent(in) :: prog
        write(*, '(a,a,i0,a,i0,a)') prog, ": ", n_pass, " passed, ", n_fail, " failed"
        if (n_fail > 0) error stop 1
    end subroutine finish

    ! Ascending insertion sort (small test vectors).
    function sort_vec(v) result(s)
        real(dp), intent(in) :: v(:)
        real(dp), allocatable :: s(:)
        real(dp) :: t
        integer :: i, j
        s = v
        do i = 2, size(s)
            t = s(i)
            j = i - 1
            do while (j >= 1)
                if (s(j) <= t) exit
                s(j + 1) = s(j)
                j = j - 1
            end do
            s(j + 1) = t
        end do
    end function sort_vec

    function median(v) result(m)
        real(dp), intent(in) :: v(:)
        real(dp) :: m
        real(dp), allocatable :: s(:)
        integer :: n
        s = sort_vec(v)
        n = size(s)
        if (mod(n, 2) == 1) then
            m = s((n + 1) / 2)
        else
            m = 0.5_dp * (s(n / 2) + s(n / 2 + 1))
        end if
    end function median

    function wall_time() result(t)
        real(dp) :: t
        integer(8) :: c, r
        call system_clock(c, r)
        t = real(c, dp) / real(r, dp)
    end function wall_time

end module plastic_test_utils
