
module mc_walk
 !   use mpi
    implicit none
    integer :: nwalks
    save
    !call MPI_INIT(ierr)
!    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
  !  call MPI_COMM_RANK(MPI_COMM_WORLD, proc_num, ierr)

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
    use mpi
    implicit none
    integer, intent(in) :: i0, j0, max_steps, n_mc
    real(kind=8), intent(out) :: u_mc
    integer, intent(out) :: n_success

    real(kind=8) :: ub_sum, ub
    integer :: i, j,k,m, iabort
    integer :: ierr, numprocs, proc_num, numsent, sender
    integer, dimension(MPI_STATUS_SIZE) :: status

    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, proc_num, ierr)
 
    ub_sum = 0.d0
    n_success = 0

    !do  k=1,n_mc
    

    if (proc_num==0) then
        !call MPI_BCAST(i0, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
        !call MPI_BCAST(j0, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
        !call MPI_BCAST(max_steps, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
        numsent =0    
        
        do k=1,min(numprocs-1,n_mc)
            call MPI_SEND(MPI_BOTTOM, 0, MPI_INTEGER, k, 1, &
            MPI_COMM_WORlD, ierr)
            numsent = numsent +1
            enddo

        do k=1, n_mc
            call MPI_RECV(ub,1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
            MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
            sender = status(MPI_SOURCE)
            iabort = status(MPI_TAG)

            if (iabort==0) then
                ub_sum = ub_sum + ub
                n_success = n_success +1
            endif
            
            !u_mc = ub_sum/n_success

            if (numsent < n_mc) then
                call MPI_SEND(MPI_BOTTOM, 0, MPI_INTEGER, sender, 1, &
                            MPI_COMM_WORLD,ierr)
                numsent =numsent+1
            else
                call MPI_SEND(MPI_BOTTOM, 0, MPI_INTEGER, sender, &
                                0, MPI_COMM_WORLD, ierr)
            endif
        !    u_mc = ub_sum/n_success
        enddo
        u_mc = ub_sum/n_success
    endif !proc_num == 0 to send processes

    if (proc_num /= 0) then
        if (proc_num > n_mc) go to 99 ! no work expected

        do while (.true.)
            call MPI_RECV(MPI_BOTTOM,0,MPI_INTEGER,0, &
                         MPI_ANY_TAG, MPI_COMM_WORLD,status, ierr )

            m = status(MPI_TAG)

            if (m==0) go to 99

            call random_walk(i0, j0, max_steps, ub, iabort)
            call MPI_SEND(ub, 1, MPI_DOUBLE_PRECISION, 0, iabort, &
                        MPI_COMM_WORLD, ierr)

            enddo
        endif

    99 continue

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    !do k=1, numprocs-1
      !  call MPI_RECV(nwalks, 1, MPI_INTEGER, 0, MPI_ANY_TAG, &
      !          MPI_COMM_WORLD, status, ierr)
     !   print *, "walks by Process :  ", nwalks
    !enddo

    !call MPI_BARRIER(MPI_COMM_WORLD,ierr) 

    !if (proc_num == 0) then

      !  call MPI_REDUCE(nwalks, nwalksTot, 1, MPI_INTEGER, MPI_SUM, &
     !               0,MPI_COMM_WORLD, ierr)

    !endif

    
    end subroutine many_walks
    !call MPI_FINALIZE(ierr)
end module mc_walk

