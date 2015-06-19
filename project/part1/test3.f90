
program test2

    use mpi

    use quadrature, only: trapezoid
    use functions, only: f, fevals_proc, k

    implicit none
    real(kind=8) :: ab, dx_sub, a,b,int_true, int_approx
    real(kind=8) :: int_sub_each, ab_sub(2), int_sub
    integer :: j,nerr, nsub,numsent, nextinterv,sender
    integer :: proc_num, num_procs, ierr, n, fevals_total
    integer, dimension(MPI_STATUS_SIZE) :: status

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, proc_num, ierr)

    ! All processes set these values so we don't have to broadcast:
    k = 1.d3   ! functions module variable
    a = 0.d0
    b = 2.d0
    int_true = (b-a) + (b**4 - a**4) / 4.d0 - (1.d0/k) * (cos(k*b) - cos(k*a))
    n = 1000
!    nsub = num_procs -1
   
!int_approx = 0
    ! Each process keeps track of number of fevals:
    fevals_proc = 0

    nerr = 0
    if (proc_num==0) then
        int_approx =0.d0 
        print '("Using ",i3," processes")', num_procs
        print '("true integral: ", es22.14)', int_true
        print *, " "  ! blank line
        print *, "How many subintervals? "
        read *, nsub
        !Check for  number of processors larger than 1
        if (num_procs == 1) then
            print *, "*** Error, this version requires at least two processors ***"
            nerr = 1
            endif
        endif

    call MPI_BCAST(nerr, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    if (nerr == 1) then
        call MPI_FINALIZE(ierr)
        stop
        endif
    call MPI_BCAST(nsub, 1, MPI_INT, 0,  MPI_COMM_WORLD, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr) ! wait for process 0 to print

    ! Note: In this version all processes call trap and repeat the
    !       same work, so each should get the same answer.

    if (proc_num==0) then
        numsent = 0
        dx_sub = (b-a) / nsub
        do j=1,min(num_procs-1, nsub)
            ab_sub(1) = a+ (j-1)*dx_sub
            ab_sub(2) = a+ j*dx_sub
            call MPI_SEND(ab_sub, 2, MPI_DOUBLE_PRECISION, j, j, &
                          MPI_COMM_WORLD,ierr)
            numsent = numsent+1
            enddo
        do j=1,nsub
            call MPI_RECV(int_sub, 1, MPI_DOUBLE_PRECISION, &
                          MPI_ANY_SOURCE, MPI_ANY_TAG, &
                          MPI_COMM_WORLD,status,ierr)
            sender = status(MPI_SOURCE)
            int_approx =  int_approx+int_sub

            if (numsent < nsub) then
                nextinterv = numsent + 1
                ab_sub(1) = a +(nextinterv-1)*dx_sub
                ab_sub(2) = a + nextinterv*dx_sub
                call MPI_SEND(ab_sub, 2, MPI_DOUBLE_PRECISION, sender, &
                            nextinterv, MPI_COMM_WORLD, ierr)
                numsent = numsent +1
                else
                !reached last
                call MPI_SEND(MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION, sender, &
                            0, MPI_COMM_WORLD,ierr)
                endif
            enddo
        endif
    !CHANGED:
    !int_approx = trapezoid(f,a,b,n)
    !print '("Process ",i3," with n = ",i8," computes int_approx = ",es22.14)', &
    !        proc_num,n, int_approx
    if (proc_num /= 0) then
        if (proc_num >  nsub) go to 99
            do while (.true.)
            j = status(MPI_TAG)
            if (j==0) go to 99
            call MPI_RECV(ab_sub, 2, MPI_DOUBLE_PRECISION, &
                         0, MPI_ANY_TAG, &
                         MPI_COMM_WORLD,status,ierr)
            int_sub = trapezoid(f,ab_sub(1), ab_sub(2), n)
            !j = status(MPI_TAG)
            !if (j==0) go to 99
            call MPI_SEND(int_sub, 1,  MPI_DOUBLE_PRECISION, &
                          0, j, MPI_COMM_WORLD, ierr)
            enddo
        endif

99  continue

    call MPI_BARRIER(MPI_COMM_WORLD,ierr) ! wait for all process to print
    ! print the number of function evaluations by each thread:
    print '("fevals by Process ",i2,": ",i13)',  proc_num, fevals_proc

    call MPI_BARRIER(MPI_COMM_WORLD,ierr) ! wait for all process to print

    call MPI_REDUCE(fevals_proc,fevals_total,1,MPI_INT,MPI_SUM,0,&
                    MPI_COMM_WORLD,ierr)
    if (proc_num==0) then
        ! This is wrong -- part of homework is to fix this:
        !call MPI_REDUCE(fevals_proc,fevals_total,1, MPI_INT,MPI_SUM,0, &
         !           MPI_COMM_WORLD,ierr) 
        !fevals_total = 0   !! need to fix
        print '("Trapezoid approximation with ",i8, " total points: &
        ",es22.14)', nsub*n, int_approx
        print '("Total number of fevals: ",i10)', fevals_total
        endif

    call MPI_FINALIZE(ierr)

end program test2
