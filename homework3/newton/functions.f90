! $UWHPSC/codes/fortran/newton/functions.f90

module functions

contains

real(kind=8) function f_sqrt(x)
    implicit none
    real(kind=8), intent(in) :: x

    f_sqrt = x**2 - 4.d0

end function f_sqrt


real(kind=8) function fprime_sqrt(x)
    implicit none
    real(kind=8), intent(in) :: x
    
    fprime_sqrt = 2.d0 * x

end function fprime_sqrt


real(kind=8) function g1val(x)
    implicit none
    real(kind=8), parameter :: pi = 3.141592653589793d0
    real(kind=8), intent(in)::x

    g1 =  x+cos(pi*x)

end function g1val


real(kind=8) function g1pval(x)
    implicit none
    real(kind=8), parameter :: pi = 3.141592653589793d0
    real(kind=8), intent(in)::x

    g1p = cos(pi*x) -x*pi*sin(pi*x)

end function g1pval


real(kind=8) function g2val(x)
    implicit none
    real(kind=8), intent(in)::x

    g2 = 1.d0-0.6d0*x**2.d0

end function g2val


real(kind=8) function g2pval(x)
    implicit none
    real(kind=8), intent(in)::x

    g2p = -1.2d0*x

end function g2pval



end module functions
