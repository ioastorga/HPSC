
module functions

    use omp_lib
    implicit none
    integer :: fevals(0:7), gevals(0:7)
    real(kind=8) :: k
    save

contains

    real(kind=8) function f(x)
        implicit none
        real(kind=8), intent(in) :: x
!    	real(kind=8), external :: g

        integer thread_num,ny, j
	    real(kind=8) :: ay, by, hy, trap_sum, yj

        ! keep track of number of function evaluations by
        ! each thread:
        thread_num = 0   ! serial mode
        !$ thread_num = omp_get_thread_num()
        fevals(thread_num) = fevals(thread_num) + 1

        ay = 4.d0
        by = 1.d0
        ny = 1000

        hy = (by-ay)/(ny-1)
        trap_sum = 0.5d0*(g(x,ay) + g(x,by))  ! endpoint contributions

        do j=2,ny-1
            yj = ay + (j-1)*hy
            trap_sum = trap_sum + g(x,yj)
            enddo

        f = hy * trap_sum
    end function f

    real(kind=8) function g(x,y)
        implicit none
        real(kind=8), intent(in) :: x, y
        integer thread_num

        ! keep track of number of function evaluations by
        ! each thread:
        thread_num = 0   ! serial mode
        !$ thread_num = omp_get_thread_num()
        gevals(thread_num) = gevals(thread_num) + 1

        g = sin(x+y)

    end function g

end module functions
