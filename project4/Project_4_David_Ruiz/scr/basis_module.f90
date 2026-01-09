module basis_types
    implicit none
    integer, parameter :: dp = kind(1.0d0)

    ! Primitive Gaussian Structure
    type :: primitive
        real(dp) :: a_nu      ! Exponent
        real(dp) :: d_nu_mu   ! Contraction Coefficient
    end type primitive

    ! Contracted Gaussian Structure
    type :: contracted
        integer :: l, m       ! Angular momentum and magnetic component
        integer :: idx        ! Unique identifier index
        integer :: n_prims    ! Number of primitives
        real(dp) :: center(3) ! Coordinates (X, Y, Z) - Required for Ex 3.3
        type(primitive), allocatable :: prims(:)
    end type contracted

contains

    ! Subroutine to print orbital information in formula format
    subroutine print_ao(ao, prefix)
        type(contracted), intent(in) :: ao
        character(len=*), intent(in) :: prefix
        integer :: i
        character(len=20) :: angular_part, m_str
        character(len=20) :: label
        logical :: first_term 

        ! 1. Build label (e.g., 1s, 2p-1...)
        write(label, '(A,I0,A)') trim(prefix), ao%idx, get_l_char(ao%l)

        if (ao%l > 0) then
            write(m_str, '(I0)') ao%m
            label = trim(label) // trim(m_str)
        end if

        print *, "chi " // trim(label) // " ="

        ! 2. Angular Part Construction
        if (ao%l == 0) then
            angular_part = ""
        else
            write(m_str, '(I0)') ao%m
            write(angular_part, '(" * S", I0, A, "(r)")') ao%l, trim(m_str)
        end if

        ! 3. Print primitives (Formula Style)
        first_term = .true.
        do i = 1, ao%n_prims
            ! Only print if coefficient is significant
            if (abs(ao%prims(i)%d_nu_mu) > 1.0e-10) then
                
                if (.not. first_term) then
                    write(*, '(" +")', advance='no')
                    print *, ""
                end if

                if (ao%l == 0) then
                    write(*, '(4X, F16.8, " * exp(", F16.8, " * r^2)")', advance='no') &
                          ao%prims(i)%d_nu_mu, -ao%prims(i)%a_nu
                else
                    write(*, '(4X, F16.8, A, " * exp(", F16.8, " * r^2)")', advance='no') &
                          ao%prims(i)%d_nu_mu, trim(angular_part), -ao%prims(i)%a_nu
                end if
                first_term = .false.
            end if
        end do
        
        print *, ""
        print *, ""

    end subroutine print_ao

    ! Helper function to convert L quantum number to character
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
