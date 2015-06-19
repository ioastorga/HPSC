
program Laplace
    
    use mpi
!    use fun_u, only: utrue, uboundary
    use mc_walk, only: random_walk, many_walks, nwalks
    use fun_u, only: utrue, uboundary 
    use random_util, only: init_random_seed

    implicit none
    real(kind=8), parameter :: ax = 0.0d0
    real(kind=8), parameter :: bx = 1.0d0
    real(kind=8), parameter :: ay = 0.4d0
    real(kind=8), parameter :: by = 1.0d0
    integer, parameter :: nx = 19
    integer, parameter :: ny = 11
    real(kind=8), parameter :: dx = (bx - ax) / (nx+1)
    real(kind=8), parameter :: dy = (by - ay) / (ny+1)

    real(kind=8) :: x0, y0, u_mc, error, ut
    real(kind=8) :: u_mc_total, u_sum_old, u_sum_new
    integer :: seed1, max_steps, n_mc, n_success, n_total
    integer :: i,i0,j0

    call MPI_INIT(ierr)
    nwalks = 0
    open(unit=25, file='mc_Laplace_error.txt', status='unknown')

    x0 = 0.9d0
    y0 = 0.6d0

    i0 = nint((x0-ax)/dx)
    j0 = nint((y0-ay)/dy)

    x0 = ax + i0*dx
    y0 = ay + j0*dy
    
    ut= utrue(x0,y0)
    seed1 = 12345
    seed1 =seed1 + 97*proc_num
    call init_random_seed(seed1)

    max_steps = 100*max(nx, ny)
    n_mc = 10

    call many_walks(i0, j0, max_steps, n_mc, u_mc, n_success)
    error = abs((u_mc - ut)/ut)
    u_mc_total = u_mc
    n_total = n_success

    write(25,'(i10,e23.15,e15.6)') n_total, u_mc_total, error
    print '(i10,e23.15,e15.6)', n_total, u_mc_total, error

    do i=1, 12
        u_sum_old = u_mc_total * n_total
        call many_walks(i0, j0, max_steps, n_mc, u_mc, n_success)
        u_sum_new = u_mc * n_success
        n_total = n_total + n_success
        u_mc_total = (u_sum_old + u_sum_new)/n_total
        error = abs((u_mc_total -ut)/ut)

        write(25,'(i10,e23.15, e15.6)') n_total, u_mc_total, error
        print '(i10,e23.15,e15.6)', n_total, u_mc_total, error
        
        n_mc = 2*n_mc !double N for next iter 
    enddo

    print '("Final approximation to u(x0, y0):  ", es22.14)', u_mc_total
    print *, "Total number random walks:  ", nwalks

    call MPI_FINALIZE(ierr)

end program Laplace
