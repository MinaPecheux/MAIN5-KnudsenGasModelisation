!--------------------------------------------------------------------
! HEAT APPROXIMATION MODULE
!--------------------------------------------------------------------

module heat_approximation

! avoid automatic variables declarations
implicit none
    integer                           :: P, Q, N
    real                              :: R, k, dr, dtheta, dt
    real(8), parameter                :: PI  = 4 * atan (1.0_8)
    real, dimension(:), allocatable   :: u, u_exact
    real, dimension(:,:), allocatable :: M
    
    ! (specific to Gaussian)
    real :: mu = 0.0, sigma = 0.1

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
        real    :: c, tmp
        
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
    
    real function Ai(i)
        integer, intent(in) :: i
        Ai = 1.0 + 2*k*dt / (dr**2) + 2*k*dt/((i*dr)**2 * (dtheta**2))
    end function Ai
    
    real function Bi(i)
        integer, intent(in) :: i
        Bi = -k*dt / ((i*dr)** 2 * (dtheta**2))
    end function Bi
    
    real function Ci_p(i)
        integer, intent(in) :: i
        Ci_p = -k*dt/(dr**2) * (1.0 + 1/(2*i))
    end function Ci_p
    
    real function Ci_m(i)
        integer, intent(in) :: i
        Ci_m = -k*dt/(dr**2) * (1.0 - 1/(2*i))
    end function Ci_m
    
    subroutine initialize_matrix()
        integer :: lp, lq, tmp_q
        
        print *, "....",Ci_p(1)
        print *, Ci_m(1)
        
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
        do lq = 1, Q
            M(lq+1,1)    = Ci_m(1)
            M(lq+1,lq+1) = Ai(1)
            if (lq == 1) then
                M(lq+1, Q+1)  = Bi(1)
            else
                M(lq+1, lq)   = Bi(1)
            end if
            if (lq == Q) then
                M(lq+1, 1)    = Bi(1)
            else
                M(lq+1, lq+2) = Bi(1)
            end if
            M(lq+1, lq+Q+1) = Ci_p(1)
        end do
        
        ! middle points
        do lp = 2, P-2
            do tmp_q = 1, Q
                lq           = tmp_q + (lp-1) * Q
                M(lq+1,lq+1-Q) = Ci_m(lp)
                M(lq+1,lq+1)   = Ai(lp)
                if (tmp_q == 1) then
                    M(lq+1, lq+Q)       = Bi(lp)
                else
                    M(lq+1, lq)         = Bi(lp)
                end if
                if (tmp_q == Q) then
                    M(lq+1, 2+(lp-1)*Q) = Bi(lp)
                else
                    M(lq+1, lq+2)       = Bi(lp)
                end if
            end do
            M(lq+1, lq+Q+1) = Ci_p(lp)
        end do
            
        ! specific treatment for last ring (p = P-1)
        do tmp_q = 1, Q
            lq             = tmp_q + (P-2) * Q
            M(lq+1,lq+1-Q) = Ci_m(P-1)
            M(lq+1,lq+1)   = Ai(P-1)
            if (tmp_q == 1) then
                M(lq+1, lq+Q)      = Bi(P-1)
            else
                M(lq+1, lq)        = Bi(P-1)
            end if
            if (tmp_q == Q) then
                M(lq+1, 2+(P-2)*Q) = Bi(P-1)
            else
                M(lq+1, lq+2)      = Bi(P-1)
            end if
        end do
        
    end subroutine initialize_matrix

    subroutine initialize_variables(lP, lQ, lR, lk, ldt)
        implicit none
        integer, intent(in) :: lP, lQ
        real,    intent(in) :: lR, lk, ldt
        integer :: i, j

        ! initialize variables
        P = lP; Q = lQ; R = lR; k = lk
        dr = R / P; dtheta = 2*PI / Q
        !dt = ldt
        dt = 0.5 / (k/(dr**2) + k/(dtheta**2))
        N = 1 + (P-1) * Q
        
        ! initialize arrays
        call initialize_arrays()        
        ! initialize scheme matrix M
        call initialize_matrix()
        
        do i = 1,N
            do j = 1,N
                print *, "M[", i, ",", j, "] =", M(i,j)
            end do
        end do
    end subroutine initialize_variables
    
    subroutine free_memory()
        deallocate(u)
        deallocate(u_exact)
        deallocate(M)
    end subroutine free_memory
    
    !--------------------------------------------------------------------
    ! SCHEME APPROXIMATION FUNCTIONS
    !--------------------------------------------------------------------
    subroutine update_exact(t)
        implicit none
        integer, intent(in) :: t
        integer             :: lp, lq
        real                :: c
        
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
        external sgesv
        integer            :: info
        real, dimension(N) :: ipiv

        ! solve linear system M*next_u = u (result is directly stored in u)
        call sgesv(N, 1, M, N, ipiv, u, N, info)
    end subroutine update
    
    real function diff_solutions()
        diff_solutions = u_exact(1) - u(1)
    end function diff_solutions

end module heat_approximation

!--------------------------------------------------------------------
! TEST PROGRAM
!--------------------------------------------------------------------

program simulate
    ! import heat_approximation module
    use heat_approximation
    
    integer :: T = 20, lt
    real    :: d
    
    ! initialize variables
    call initialize_variables(5, 6, 1.0, 1.0, 0.001)
    !call initialize_variables(20, 30, 1.0, 1.0, 0.001)
    
    do lt = 1, T
        call update()
        call update_exact(lt)
        d = diff_solutions()
        print *, "t =", lt, "     diff =", d, "     u_exact(1) =", u_exact(1), "u(1) =", u(1)
    end do
    
    ! free memory
    call free_memory()
end program simulate
