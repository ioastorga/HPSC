
! intersection  of two curves, fortran verison to find the 4 points within
    ! -5 to 5

program testquartic

    use newton, only: solve
    use functions, only: f_quartic, fprime_quartic
    implicit none
    real(kind=8),dimension(3) :: toler, epiN
    integer :: i,j, iters
    real(kind=8) :: fx, xtry, xout, errorX
    logical :: debug = .false.

    xtry = 4.0
    toler = (/1.d-5, 1.d-10, 1.d-14/)
    epiN = (/1.d-4, 1.d-8, 1.d-12/)

    print *, '    epsilon        tol    iters          x                 f(x)
    x-xstar'

    i=1
    do while (i<= 3)
        epsilon = epiN(i)
        j=1
        do while (j <= 3)
            tol = toler(j)
            call solve(f_quartic, fprime_quartic, xtry,xout,iters,debug)
            fx = f_quartic(xout)
            errorX = xout - xtry
            print 11, epsilon, tol, iters, xout, fx, errorX
         11 format(2es13.3, i4, es24.15, 2es13.3)
        i = i+1
        enddo
    enddo
end program testquartic


