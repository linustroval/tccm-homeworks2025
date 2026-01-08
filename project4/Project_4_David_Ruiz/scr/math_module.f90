module math_utils
    use basis_types   
    implicit none

    ! Helper structure for math only (polynomials x^i y^j z^k)
    type :: cart_poly
        integer :: i, j, k
        real(dp) :: prefactor
    end type cart_poly

contains

    ! --- 1. PURE MATH ---

    function fact2(n) result(res)
        integer, intent(in) :: n
        integer :: res, i
        res = 1
        if (n <= 0) return
        do i = n, 1, -2
            res = res * i
        end do
    end function fact2

    function gaussian_int_1d(n, a) result(val)
        integer, intent(in) :: n
        real(dp), intent(in) :: a
        real(dp) :: val, pi
        
        if (mod(n, 2) /= 0) then
            val = 0.0_dp; return
        end if
        pi = acos(-1.0_dp)
        val = real(fact2(n-1), dp) / ((4.0_dp * a)**(n/2)) * sqrt(pi / (2.0_dp * a))
    end function gaussian_int_1d

    ! --- 2. GEOMETRY (TABLE 1) ---

    subroutine get_solid_harmonic(l, m, terms, n_terms)
        integer, intent(in) :: l, m
        type(cart_poly), intent(out), allocatable :: terms(:)
        integer, intent(out) :: n_terms
        real(dp) :: sqrt3

        sqrt3 = sqrt(3.0_dp)

        select case (l)
        case (0) ! S
            n_terms = 1; allocate(terms(1))
            terms(1) = cart_poly(0,0,0, 1.0_dp)
        case (1) ! P
            n_terms = 1; allocate(terms(1))
            select case (m)
                case (1);  terms(1) = cart_poly(1,0,0, 1.0_dp) ! x
                case (-1); terms(1) = cart_poly(0,1,0, 1.0_dp) ! y
                case (0);  terms(1) = cart_poly(0,0,1, 1.0_dp) ! z
            end select
        case (2) ! D
            select case (m)
                case (-2); n_terms=1; allocate(terms(1)); terms(1)=cart_poly(1,1,0, sqrt3)
                case (-1); n_terms=1; allocate(terms(1)); terms(1)=cart_poly(0,1,1, sqrt3)
                case (0) 
                    n_terms=3; allocate(terms(3))
                    terms(1)=cart_poly(0,0,2, 1.0_dp); terms(2)=cart_poly(2,0,0, -0.5_dp); terms(3)=cart_poly(0,2,0, -0.5_dp)
                case (1); n_terms=1; allocate(terms(1)); terms(1)=cart_poly(1,0,1, sqrt3)
                case (2)
                    n_terms=2; allocate(terms(2))
                    terms(1)=cart_poly(2,0,0, 0.5_dp*sqrt3); terms(2)=cart_poly(0,2,0, -0.5_dp*sqrt3)
            end select
        end select
    end subroutine get_solid_harmonic

    ! --- 3. ELEMENTAL OVERLAP CALCULATION ---
    
    function overlap_spherical_component(l, m, a) result(ovlp)
        integer, intent(in) :: l, m
        real(dp), intent(in) :: a
        real(dp) :: ovlp
        type(cart_poly), allocatable :: terms(:)
        integer :: n_terms, i, j, p_x, p_y, p_z

        call get_solid_harmonic(l, m, terms, n_terms)
        ovlp = 0.0_dp

        do i = 1, n_terms
            do j = 1, n_terms
                p_x = terms(i)%i + terms(j)%i
                p_y = terms(i)%j + terms(j)%j
                p_z = terms(i)%k + terms(j)%k
                
                ovlp = ovlp + (terms(i)%prefactor * terms(j)%prefactor) * &
                       gaussian_int_1d(p_x, a) * gaussian_int_1d(p_y, a) * gaussian_int_1d(p_z, a)
            end do
        end do
        if (allocated(terms)) deallocate(terms)
    end function overlap_spherical_component

end module math_utils
