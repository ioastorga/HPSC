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


real(kind=8) function gval(x)
    implicit none
    real(kind=8), parameter :: pi = 3.141592653589793d0
    real(kind=8), intent(in) ::x

    gval =  x*cos(pi*x) - 1.d0 + 0.6d0*x**2

end function gval


real(kind=8) function gpval(x)
    implicit none
    real(kind=8), parameter :: pi = 3.141592653589793d0
    real(kind=8), intent(in) ::x

    gpval = cos(pi*x) - x*pi*sin(pi*x) + 1.2d0*x

end function gpval

end module functions
