
module fun_u

    implicit none
    real(kind=8), parameter :: ax = 0.0d0
    real(kind=8), parameter :: bx = 1.0d0
    real(kind=8), parameter :: ay = 0.4d0
    real(kind=8), parameter :: by = 1.0d0
    integer, parameter :: nx = 19
    integer, parameter :: ny = 11
    real(kind=8), parameter :: dx = (bx - ax) / (nx+1)
    real(kind=8), parameter :: dy = (by - ay) / (ny+1)

contains

real(kind=8) function utrue(x,y)

    ! True solution for comparison, if known.
    implicit none
    real(kind=8), intent(in) :: x,y  

    utrue = x**2 - y**2

end function utrue

real(kind=8) function uboundary(x,y)

    ! Return u(x,y) assuming (x,y) is a boundary point.
    implicit none
    real(kind=8), intent(in) :: x,y
    real(kind=8) :: try_boundary

    try_boundary = (x-ax)*(x-bx)*(y-ay)*(y-by)
    if (try_boundary .ne. 0.d0) then
        print *, "*** Error -- called uboundary at non-boundary point"
        stop
        endif

    uboundary = utrue(x,y)   ! assuming we know this

end function uboundary
    

end module fun_u
