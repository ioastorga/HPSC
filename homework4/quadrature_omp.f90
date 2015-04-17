! quadrature integration form quadrature.py

module quadrature_omp
    use omp_lib
    ! module parameters:
    implicit none
contains

real(kind=8) function trapezoid(f,a,b,n)

    ! Estimate the  integral of f(x) using quadrature  method. 
    ! Input:
    !   f:  the function to find an integral of
    !   a: lower limit
    !   b: upper limit
    !   n: number of subdivisions to integrate
    ! Returns:
    !   the estimate integral of f: int_f

    implicit none
    real(kind=8), intent(in) :: a, b
    integer, intent(in) :: n
    real(kind=8), external :: f
    !real(kind=8), intent(out) :: int_f

    ! Declare any local variables:
    real(kind=8) :: h, int_tmp, int_ini, int_f
    integer :: k

    h = (b-a)/(n-1)
    int_tmp = 0
    int_ini = -0.5*(f(a)+f(b))
    int_f = 0
    !$omp parallel do reduction(+ : int_f)
    do k=1,n-1 
        int_tmp = f(a+k*h)
        int_f = int_f + int_tmp
        enddo
    !$omp end parallel do

    int_f = h*(int_f+int_ini)
end function trapezoid

subroutine  error_table(f,a,b,nvals, int_true)
    implicit none
    real(kind=8), intent(in) :: a, b, int_true
    integer, dimension(:), intent(in) :: nvals
    real(kind=8), external :: f

    real(kind=8) ::last_error, error, ratio, int_trap 
    integer :: j,n
    
    last_error = 0
    print *, "      n       trapezoid       error       ratio"
    do j=1,size(nvals)
        n = nvals(j)
        int_trap =  trapezoid(f,a,b,n)
        error =abs(int_trap-int_true)
        ratio = last_error/error
        last_error = error
        print 11,  n , int_trap, error, ratio
    11  format(i8, es22.14, es13.3, es13.3)
        enddo
end subroutine error_table

end module quadrature_omp
