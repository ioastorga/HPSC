
module quadrature_mc

    implicit none

contains
    ! g       is the function defining the integrand. g takes two arguments x and ndim,
    !         where x is an array of length ndim, the number of dimensions we are integrating
    !         over. 
    ! a, b    are arrays of length ndim that have the lower and upper limits of
    !         integration in each dimension.
    ! ndim    is the number of dimensions to integrate over.
    ! npoints is the number of Monte Carlo samples to use.
    function quad_mc(g, a, b, ndim, npoints)

    implicit none
    real(kind=8), external :: g
    integer, intent(in) :: ndim, npoints
    real(kind=8), dimension(ndim), intent(in) :: a,b

    real(kind=8) :: V, sum_npoints, ai,bi
    real(kind=8), dimension(ndim) :: r_ndim
    real(kind=8), dimension(ndim) :: x
    real(kind=8) :: quad_mc
    integer :: i,j,c
    real(kind=8), allocatable :: r(:)
    allocate(r(ndim*npoints))

    call random_number(r)

    c = 1
    V = 1
    sum_npoints = 0

    do j = 1,npoints
        r_ndim = r(c:c+ndim-1)
        do i = 1,ndim
            bi = b(i)
            ai = a(i)
            if (j == 1) then
                V = V*(bi-ai)
            endif
            x(i) = ai + r_ndim(i)*(bi-ai)

            !sum_npoints =  sum_npoints + g(x,ndim)
            enddo
        sum_npoints = sum_npoints + g(x,ndim)
        c = c+ndim
        enddo

      quad_mc = V/npoints*sum_npoints

     end function quad_mc

end module quadrature_mc
