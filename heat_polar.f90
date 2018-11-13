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
    integer function idx(p, q)
        integer, intent(in) :: p, q
        if (p > 0) then
            idx = 2 + (p-1)*Q + mod(q-1, Q) ! Fortran is 1-indexed...
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
        integer :: p, q
        real    :: c
        
        ! allocate arrays
        allocate(u(N))
        allocate(u_exact(N))
        
        c          = 1.0 / (sqrt(2*PI) * sigma)
        u(1)       = c
        u_exact(1) = c
        do p = 1, P
            do q = 1, Q
                u(idx(p, q))       = c * exp(-(p*dr - mu) ** 2 / (2*sigma ** 2))
                u_exact(idx(p, q)) = c * exp(-(p*dr - mu) ** 2 / (2*sigma ** 2))
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
        Ci_p = -k*dt/(dr**2) * (1.0 + 1.0/(2*i))
    end function Ci_p
    
    real function Ci_m(i)
        integer, intent(in) :: i
        Ci_m = -k*dt/(dr**2) * (1.0 - 1.0/(2*i))
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
        M(1,1)     = 1.0 + 4*k*dt/(dr**2)
        M(1,2)     = -2*k*dt/(dr**2)
        M(1,Q/2+2) = -2*k*dt/(dr**2)
        
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
        do p = 2, P-2
            do tmp_q = 1, Q
                lq           = tmp_q + (p-1) * Q
                M(lq+1,lq+1-Q) = Ci_m(p)
                M(lq+1,lq+1)   = Ai(p)
                if (tmp_q == 1) then
                    M(lq+1, lq+Q)      = Bi(p)
                else
                    M(lq+1, lq)        = Bi(p)
                end if
                if (tmp_q == Q) then
                    M(lq+1, 2+(p-1)*Q) = Bi(p)
                else
                    M(lq+1, lq+2)      = Bi(p)
                end if
            end do
            M(lq+1, lq+Q+1) = Ci_p(p)
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

        ! initialize variables
        P = lP; Q = lQ; R = lR
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
    subroutine update_exact(t)
        implicit none
        integer, intent(in) :: t
        integer             :: p, q
        real                :: c
        
        ! update values
        c = 1.0 / (sqrt(2*PI * (sigma ** 2 + 2*k*t*dt)))
        do p = 0, P-1
            do q = 1, Q
                u(idx(p, q)) = c * exp(-(p*dr - mu) ** 2 / (2*(sigma ** 2 + 2*k*t*dt)))
            end do
        end do
    end subroutine update_exact

    subroutine update()
        implicit none
        external dgesv
        integer                         :: info
        integer, dimension(N)           :: ipiv
        real, dimension(:), allocatable :: next_u
        
        ! temporary holder
        allocate(next_u(N))
        
        ! solve linear system M*next_u = u
        call dgesv(N, 1, M, N, ipiv, u, N, info)
        ! update values
        u = next_u
        
        ! free memory
        deallocate(next_u)
    end subroutine update
    
    real function diff_solutions()
        diff_solutions = u_exact(1) - u(1)
    end function diff_solutions

end module heat_approximation

program simulate
    ! import linear algebra packages
    !use la_precision, only: wp => sp
    !use f95_lapack
    ! import heat_approximation module
    use heat_approximation
    
    integer :: T = 10, lt
    real    :: d
    
    ! initialize variables
    call initialize_variables(20, 30, 1.0, 1.0, 0.001)
    
    do lt = 1, T
        call update()
        call update_exact(lt)
        d = diff_solutions()
        print *, lt, d
    end do
    
    ! free memory
    call free_memory()
end program simulate
