    ! intersection  of two curves, fortran verison to find the 4 points within
    ! -5 to 5

program Intersections

    use newton, only: solve
    use functions, only: gval, gpval 
    implicit none
    real(kind=8),dimension(:), intent(in) :: x0
    real(kind=8),dimension(size(x0)), intent(out):: x
    integer, dimension(size(x0)),intent(out) :: iter
    integer :: i, iters
    real(kind=8) :: fx, xtry, xout

    i=0
    do while (i< size(x0))
        xtry = x0(i)
        solve(gval, gpval, xtry,xout,iters,false)
        fx = 1.d0 - 0.6d0*xout**2.d0
        print *,  "The value of intersection is: ",fx
        print *,  "Point of intersection: ", xout
        x(i) = xout
        iter(i) =iters 
        enddo
    print *, "Intersections: ", x
e
    print *, "Iterations: ", iter
    
end program Intersections


