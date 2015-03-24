    ! intersection  of two curves, fortran verison to find the 4 points within
    ! -5 to 5

program Intersections(x0, x, iter)
    use newton, only: solve
    implicit none
    real(kind=8),dimension(:), intent(in) :: x0
    real(kind=8),dimension(size(x0)), intent(out):: x
    integer, dimension(sixze(x0)),intent(out) :: iter
    integer :: i

    i=0
    do while (i< size(x0))
        xtry = x0(i)
        call f_val(xtry, fx, fxp)
        solve(fx, fxp, xtry,xout,iters,false)
        call g1vla(xout,xg1,g1p)
        print *,  "The value of intersection is: ",fx
        print *,  "Point of intersection: ", xg1
        enddo
    
end program Intersections



subroutine f_val(x,fval,fpval)
    implicit none
    real(kind=8), intent(in) :: x
    real(kind=8), external :: g1,g1p,g2,g2p
    real(kind=8), intent(out):: fval, fpval
    call g1val(x,g1,g1p)
    call g2val(x,g2,g2p)
    fval = g1-g2
    fpval = g1p-g2p
end subroutine f_val
