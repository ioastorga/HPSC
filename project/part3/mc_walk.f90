
module mc_walk
    implicit none
    integer :: nwalks
    save

contains

    subroutine random_walk(i0, j0, max_steps,ub,iabort)
    use fun_u, only: nx,ny,ax,ay,dx,dy,uboundary
    
    implicit none
    integer, intent(in) :: i0, j0, max_steps
    real(kind=8), intent(out) :: ub
    integer, intent(out) :: iabort
    !use fun_u, only: nx,ny,ax,ay,dx,dy,uboundary
    integer :: istep, i, j
    real(kind=8) :: xb,yb
    real(kind=4) :: rn 
    real(kind=8), allocatable :: r(:)
    allocate(r(max_steps))

    call random_number(r)

    i = i0
    j = j0

    !move
    do istep =1, max_steps
        rn = r(istep)
        if (rn < 0.25) then
            i = i - 1
        else if (rn < 0.5)  then
            i = i + 1
        else if (rn < 0.75) then
            j = j - 1
        else
            j = j + 1
        endif

        !boundary?
        if (i*j*(nx+1-i)*(ny+1-j) == 0) then
            xb = ax + i*dx
            yb = ay + j*dy
            ub =uboundary(xb,yb)
            exit
        endif
    enddo

        !did not hit
    if (istep == (max_steps -1)) then
        iabort = 1
    else
        iabort =0
    endif
    
    
    nwalks = nwalks + 1
    end subroutine random_walk

    subroutine many_walks(i0, j0, max_steps, n_mc, u_mc, n_success)

    implicit none
    integer, intent(in) :: i0, j0, max_steps, n_mc
    real(kind=8), intent(out) :: u_mc
    integer, intent(out) :: n_success

    real(kind=8) :: ub_sum, ub
    integer :: i, j,k, iabort

    ub_sum = 0
    n_success = 0

    do  k=1,n_mc
        i = i0
        j = j0
        call random_walk(i0, j0, max_steps, ub, iabort)
        if (iabort==0) then
            ub_sum = ub_sum + ub
            n_success = n_success +1
        endif
    enddo

    u_mc = ub_sum/n_success

    end subroutine many_walks

end module mc_walk

