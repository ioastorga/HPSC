    
! intersection  of two curves, fortran verison to find the 4 points within
    ! -5 to 5

program Intersections

    use newton, only: solve
    use functions, only: gval, gpval 
    implicit none
    real(kind=8),dimension(4) :: x0
    real(kind=8),dimension(4) :: x
    integer, dimension(4) :: iter
    integer :: i, iters
    real(kind=8) :: fx, xtry, xout
    logical :: debug = .false.

    x0 = (/-2.2d0, -1.6d0, -0.8d0, 1.45d0/) 

    i=1
    do while (i<= 4)
        xtry = x0(i)
        call solve(gval, gpval, xtry,xout,iters,debug)
        fx = 1.d0 - 0.6d0*xout**2
        print *,  "The value of intersection is: ",fx
        print *,  "Point of intersection: ", xout
        x(i) = xout
        iter(i) =iters 
        i = i+1
        enddo
    print *, "Intersections: ", x
    print *, "Iterations: ", iter
    
end program Intersections


