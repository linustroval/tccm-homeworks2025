module basis_types
    implicit none
    integer, parameter :: dp = kind(1.0d0)

    ! Structure for a Primitive Gaussian
    type :: primitive
        real(dp) :: a_nu      ! Exponent
        real(dp) :: d_nu_mu   ! Contraction Coefficient
    end type primitive

    ! Structure for a Contracted Atomic Orbital (AO)
    type :: contracted
        integer :: l, m       ! Angular momentum and magnetic component
        integer :: idx        ! Unique identifier index
        integer :: n_prims    ! Number of primitives
        real(dp) :: center(3) ! Coordinates (X, Y, Z) of the nucleus
        type(primitive), allocatable :: prims(:)
    end type contracted

contains

    ! Subroutine to print orbital information
    subroutine print_ao(ao, prefix)
        type(contracted), intent(in) :: ao
        character(len=*), intent(in) :: prefix
        
        print *, trim(prefix), " Orbital", ao%idx, &
                 " (l=", ao%l, "m=", ao%m, ")", &
                 " at [", ao%center(1), ao%center(2), ao%center(3), "]"
    end subroutine print_ao

    ! Helper function to convert L number to character (s, p, d)
    function get_l_char(l) result(ch)
        integer, intent(in) :: l
        character(len=1) :: ch
        select case (l)
            case (0); ch='s'
            case (1); ch='p'
            case (2); ch='d'
            case (3); ch='f'
            case default; ch='?'
        end select
    end function get_l_char

end module basis_types
