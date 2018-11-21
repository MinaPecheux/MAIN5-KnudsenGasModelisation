!--------------------------------------------------------------------
! HEAT APPROXIMATION MODULE
!--------------------------------------------------------------------

module heat_approximation

! avoid automatic variables declarations
implicit none
    integer                             :: P, Q, N
    real*8                              :: R, k, dr, dtheta, dt, eps, dr1
    real*8, parameter                   :: PI  = 4 * atan (1.0_8)
    real*8, dimension(:), allocatable   :: u, u_exact
    real*8, dimension(:,:), allocatable :: M
    
    ! (specific to Gaussian)
    real :: mu = 0.0, sigma = 0.1

contains

    !--------------------------------------------------------------------
    ! UTIL FUNCTIONS
    !--------------------------------------------------------------------
    integer function idx(lp, lq)
        integer, intent(in) :: lp, lq
        idx = lp*Q + mod(lq, Q)
    end function idx

    real function g()
        implicit none
        g = 1.0
    end function g

    subroutine initialize_arrays()
        implicit none
        integer :: lp, lq
        real*8  :: c, tmp
        
        ! allocate arrays
        allocate(u(N))
        allocate(u_exact(N))
        
        c          = 1.0 / (sqrt(2*PI) * sigma)
        u(1)       = c
        u_exact(1) = c
        do lp = 1, P-1
            do lq = 1, Q
                tmp                  = c * exp(-(lp*dr - mu) ** 2 / (2*sigma ** 2))
                u(idx(lp, lq))       = tmp
                u_exact(idx(lp, lq)) = tmp
            end do
        end do
    end subroutine initialize_arrays
    
    real*8 function Ai(i)
        integer, intent(in) :: i
        Ai = 1.0 + 2*k*dt / (dr**2) + 2*k*dt/((i*dr)**2 * (dtheta**2))
    end function Ai
    
    real*8 function Bi(i)
        integer, intent(in) :: i
        Bi = -k*dt / ((i*dr)** 2 * (dtheta**2))
    end function Bi
    
    real*8 function Ci_p(i)
        integer, intent(in) :: i
        Ci_p = -k*dt/(dr**2) * (1.0 + 1/(2*i))
    end function Ci_p
    
    real*8 function Ci_m(i)
        integer, intent(in) :: i
        Ci_m = -k*dt/(dr**2) * (1.0 - 1/(2*i))
    end function Ci_m
    
    real*8 function Ai1(i)
        integer, intent(in) :: i
        Ai1 = 1.0 + 2*k*dt / (dr*dr1) + 2*k*dt/(i*i*dr*dr1 * (dtheta**2))
    end function Ai1
    
    real*8 function Bi1(i)
        integer, intent(in) :: i
        Bi1 = -k*dt / (i*i*dr*dr1 * (dtheta**2))
    end function Bi1
    
    real*8 function Ci_p1(i)
        integer, intent(in) :: i
        Ci_p1 = -k*dt/(dr*dr1) * (1.0 + 1/(2*i))
    end function Ci_p1
    
    real*8 function Ci_m1(i)
        integer, intent(in) :: i
        Ci_m1 = -k*dt/(dr*dr1) * (1.0 - 1/(2*i))
    end function Ci_m1
    
    subroutine initialize_matrix()
        integer :: lp, lq, tmp_q
        
        ! allocate matrix + initialize to zero
        allocate(M(N, N))
        do lp = 1, N
            do lq = 1, N
                M(lp,lq) = 0.0
            end do
        end do
        
        ! central points
        do lq = 1, Q
            M(lq,lq) = 1.0 + 4.0*k*dt/(eps**2)
            M(lq,lq+Q) = -2.0*k*dt/(eps**2)
            if (lq <= Q/2) then
                M(lq,lq+Q/2) = -2.0*k*dt/(eps**2)
            else
                M(lq,lq-Q/2) = -2.0*k*dt/(eps**2)
            end if
        end do
        
        ! specific treatment for first ring (p = 1)
        do tmp_q = 1, Q
            lq         = tmp_q + Q
            M(lq,lq-Q) = Ci_m1(1)
            M(lq,lq)   = Ai1(1)
            if (tmp_q == 1) then
                M(lq, lq-1+Q) = Bi1(1)
            else
                M(lq, lq-1)   = Bi1(1)
            end if
            if (tmp_q == Q) then
                M(lq, Q+1)    = Bi1(1)
            else
                M(lq, lq+1)   = Bi1(1)
            end if
            M(lq, lq+Q) = Ci_p1(1)
        end do
        
        ! middle points (including last ring)
        do lp = 2, P
            do tmp_q = 1, Q
                lq           = tmp_q + lp * Q
                M(lq,lq-Q) = Ci_m(lp)
                M(lq,lq)   = Ai(lp)
                if (tmp_q == 1) then
                    M(lq, lq-1+Q) = Bi(lp)
                else
                    M(lq, lq-1)   = Bi(lp)
                end if
                if (tmp_q == Q) then
                    M(lq, 1+lp*Q) = Bi(lp)
                else
                    M(lq, lq+1)   = Bi(lp)
                end if
            end do
            if (lp < P) then
                M(lq, lq+Q) = Ci_p(lp)
            end if
        end do
        
    end subroutine initialize_matrix

    subroutine initialize_variables(lp, lq, lr, lk, leps)
        implicit none
        integer, intent(in) :: lp, lq
        real*8,  intent(in) :: lr, lk, leps
        integer :: i,j

        ! initialize variables
        P = lp; Q = lq; R = lr; k = lk
        dr = R / P; dtheta = 2*PI / Q
        eps = leps
        dr1 = dr - eps
        dt = 0.1 / (k/(dr**2) + k/(dtheta**2))
        N = (P+1)*Q ! add ring for p = 0
        
        ! initialize arrays
        call initialize_arrays()        
        ! initialize scheme matrix M
        call initialize_matrix()
        
        !do i = 1,N
        !    do j = 1,N
        !        print *, "M[", i, ", ", j, "] =", M(i,j)
        !    end do
        !end do
    end subroutine initialize_variables
    
    subroutine free_memory()
        deallocate(u)
        deallocate(u_exact)
        deallocate(M)
    end subroutine free_memory
    
    !--------------------------------------------------------------------
    ! SCHEME APPROXIMATION FUNCTIONS
    !--------------------------------------------------------------------
    subroutine border_func()
        integer :: lq
        do lq = 1, Q
            u(idx(P-1, lq)) = 0.0    ! with Dirichlet condition
        end do
    end subroutine border_func
    
    subroutine update_exact(t)
        implicit none
        integer, intent(in) :: t
        integer             :: lp, lq
        real*8              :: c
        
        ! update values
        c = 1.0 / (sqrt(2*PI * (sigma ** 2 + 2*k*t*dt)))
        u_exact(1) = c
        do lp = 1, P-1
            do lq = 1, Q
                u_exact(idx(lp, lq)) = c * exp(-(lp*dr - mu) ** 2 / (2*(sigma ** 2 + 2*k*t*dt)))
            end do
        end do
    end subroutine update_exact

    subroutine update()
        implicit none
        external dgesv
        integer                :: info, i, j
        real*8, dimension(N)   :: ipiv
        real*8, dimension(N,N) :: A
        
        ! prepare temporary matrix
        do i = 1,N
            do j = 1,N
                A(i,j) = M(i,j)
            end do
        end do

        ! solve linear system M*next_u = u (result is directly stored in u)
        call dgesv(N, 1, A, N, ipiv, u, N, info)
        
        ! apply border condition
        call border_func()
    end subroutine update

end module heat_approximation

!--------------------------------------------------------------------
! TEST PROGRAM
!--------------------------------------------------------------------

program simulate
    ! import heat_approximation module
    use heat_approximation
    
    integer :: T = 20, lt
    real*8  :: lr = 1.0, lk = 1.0, leps = 0.0001
    
    ! initialize variables
    !call initialize_variables(20, 30, lr, lk, leps)
    call initialize_variables(30, 30, lr, lk, leps)
    
    open(1, file = 'output.dat')
    print *, "t = 0", "     u_exact(1) =", u_exact(1), "u(1) =", u(1)
    do lt = 1, T
        call update()
        call update_exact(lt)
        write(1,*) 0, u(1), u_exact(1)
        do lp = 1, P-1
            write(1,*) lp, u(idx(lp, 1)), u_exact(idx(lp, 1))
        end do
        print *, "t =", lt, "     diff =", u_exact(1)-u(1), "     u_ex(1) =", u_exact(1), "u(1) =", u(1)
    end do
    close(1)
    
    ! free memory
    call free_memory()
end program simulate
