!--------------------------------------------------------------------
! HEAT APPROXIMATION MODULE
!--------------------------------------------------------------------

module heat_approximation

! avoid automatic variables declarations
implicit none
    integer                             :: P, Q, N
    real*8                              :: R, k, dr, dtheta, dt, TP = 1.0
    real*8, parameter                   :: PI  = 4 * atan (1.0_8)
    real*8, dimension(:), allocatable   :: u, u_exact
    real*8, dimension(:,:), allocatable :: M
    
    ! (specific to Gaussian)
    real*8 :: mu = 0.0, sigma = 0.1

contains

    !--------------------------------------------------------------------
    ! UTIL FUNCTIONS
    !--------------------------------------------------------------------
    integer function idx(lp, lq)
        integer, intent(in) :: lp, lq
        if (lp > 0) then
            idx = 2 + (lp-1)*Q + mod(lq-1, Q) ! Fortran is 1-indexed...
        else
            idx = 1
        end if
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
    
    subroutine initialize_matrix()
        integer :: lp, lq, tmp_q
        
        ! allocate matrix + initialize to zero
        allocate(M(N, N))
        do lp = 1, N
            do lq = 1, N
                M(lp,lq) = 0.0
            end do
        end do
        
        ! central point
        M(1,1)     = 1.0 + 4.0*k*dt/(dr**2)
        M(1,2)     = -2.0*k*dt/(dr**2)
        M(1,Q/2+2) = -2.0*k*dt/(dr**2)
        
        ! specific treatment for first ring (p = 1)
        do tmp_q = 1, Q
            lq = 1 + tmp_q
            M(lq,1)  = Ci_m(1)
            M(lq,lq) = Ai(1)
            if (lq == 2) then
                M(lq, Q+1)  = Bi(1)
            else
                M(lq, lq-1) = Bi(1)
            end if
            if (lq == Q+1) then
                M(lq, 2)    = Bi(1)
            else
                M(lq, lq+1) = Bi(1)
            end if
            M(lq, lq+Q) = Ci_p(1)
        end do
        
        ! other points (mid points + last ring)
        do lp = 2, P-1
            do tmp_q = 1, Q
                lq         = 1 + tmp_q + (lp-1) * Q
                M(lq,lq-Q) = Ci_m(lp)
                M(lq,lq)   = Ai(lp)
                if (tmp_q == 1) then
                    M(lq,lq+Q-1)     = Bi(lp)
                else
                    M(lq,lq-1)       = Bi(lp)
                end if
                if (tmp_q == Q) then
                    M(lq,2+(lp-1)*Q) = Bi(lp)
                else
                    M(lq,lq+1)       = Bi(lp)
                end if
                
                if (lp < P-1) then
                    M(lq,lq+Q) = Ci_p(lp)
                end if
            end do
        end do
        
    end subroutine initialize_matrix

    subroutine initialize_variables(lp, lq, lr, lk, ldt)
        implicit none
        integer, intent(in) :: lp, lq
        real*8,  intent(in) :: lr, lk, ldt

        ! initialize variables
        P = lp; Q = lq; R = lr; k = lk
        dr = R / P; dtheta = 2*PI / Q
        dt = ldt
        N = 1 + (P-1) * Q
        
        ! initialize arrays
        call initialize_arrays()        
        ! initialize scheme matrix M
        call initialize_matrix()
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
            !u(idx(P-1, lq)) = 0.0            ! with Dirichlet condition
            u(idx(P-1, lq)) = TP - g()*dr    ! with Neumann condition
        end do
    end subroutine border_func
    
    subroutine update_exact(t)
        implicit none
        real*8, intent(in) :: t
        integer            :: lp, lq
        real*8             :: c
        
        ! update values
        c = 1.0 / (sqrt(2*PI * (sigma ** 2 + 2*k*t)))
        u_exact(1) = c
        do lp = 1, P-1
            do lq = 1, Q
                u_exact(idx(lp, lq)) = c * exp(-(lp*dr - mu) ** 2 / (2*(sigma ** 2 + 2*k*t)))
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
    
    integer :: lt
    real*8  :: lr = 1.0, lk = 1.0, ldt = 0.001, T = 0.02
    
    ! initialize variables
    call initialize_variables(40, 30, lr, lk, ldt)
    
    open(1, file = 'output.dat')
    print *, "t = 0", "     u_exact(1) =", u_exact(1), "u(1) =", u(1)
    do lt = 1, int(T/ldt)
        call update()
        call update_exact(lt*ldt)
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
