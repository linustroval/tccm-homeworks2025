program exercise_3_2
    use basis_types
    use math_utils
    implicit none

    integer :: i, k, m, error, unit_num
    integer :: target_Z, current_Z
    character(len=100) :: line
    integer :: l_curr, n_prim, n_cont
    
    real(dp), allocatable :: buffer(:,:)
    integer :: s_count, p_count, d_count 
    
    unit_num = 10
    
    print *, "=========================================="
    print *, "      EXERCISE 3.2: OVERLAP NORMS         "
    print *, "=========================================="
    print *, "Enter the desired Atomic Number (Z):"
    read(*,*) target_Z
    
    open(unit=unit_num, file='6-31g.1.dalton', status='old', action='read')

    ! 1. FIND ATOM
    do
        read(unit_num, '(A)', iostat=error) line
        if (error /= 0) stop "Error: Atom not found"
        if (line(1:1) == 'a') then
            read(line(3:), *) current_Z
            if (current_Z == target_Z) exit
        end if
    end do

    ! 2. READ DATA
    s_count = 0; p_count = 0; d_count = 0
    
    do
        read(unit_num, '(A)', iostat=error) line
        if (error /= 0) exit 
        if (line(1:1) == 'a') exit 
        
        if (index(line, 's functions') > 0) l_curr = 0
        if (index(line, 'p functions') > 0) l_curr = 1
        if (index(line, 'd functions') > 0) l_curr = 2
        
        if (line(1:1) == 'H') then
            read(line(2:), *) n_prim, n_cont
            
            if (allocated(buffer)) deallocate(buffer)
            allocate(buffer(n_prim, n_cont + 1))
            
            do i = 1, n_prim
                read(unit_num, *) buffer(i, :)
            end do
            
            do k = 1, n_cont
                if (l_curr == 0) s_count = s_count + 1
                if (l_curr == 1) p_count = p_count + 1
                if (l_curr == 2) d_count = d_count + 1
                
                do m = -l_curr, l_curr
                    call process_orbital(l_curr, m, n_prim, buffer(:,1), buffer(:, k+1))
                end do
            end do
        end if
    end do

    close(unit_num)
    if (allocated(buffer)) deallocate(buffer)

contains

    subroutine process_orbital(l, m, n, exps, coefs)
        integer, intent(in) :: l, m, n
        real(dp), intent(in) :: exps(:), coefs(:)
        type(contracted) :: temp_ao
        integer :: i, j
        real(dp) :: norm, exponent_sum, term_ovlp
        
        ! Arrays for Normalized coefficients
        real(dp), allocatable :: norm_coefs(:)
        real(dp) :: self_overlap, N_const

        ! 1. PREPARE DATA
        temp_ao%l = l; temp_ao%m = m; temp_ao%n_prims = n
        if (l==0) temp_ao%idx = s_count
        if (l==1) temp_ao%idx = p_count+1 
        if (l==2) temp_ao%idx = d_count+2 
        
        allocate(temp_ao%prims(n))
        allocate(norm_coefs(n)) 

        ! ---------------------------------------------------------
        ! NORMALIZATION STEP:
        ! 1. Calculate the integral of the primitive alone <Gi|Gi>
        ! 2. N = 1 / sqrt(Integral)
        ! 3. New Coefficient = Old Coefficient * N
        ! ---------------------------------------------------------
        do i = 1, n
            ! Calculate self-overlap of the primitive
            self_overlap = overlap_spherical_component(l, m, exps(i))
            
            ! Normalization Constant
            N_const = 1.0_dp / sqrt(self_overlap)
            
            ! Store normalized coefficient
            norm_coefs(i) = coefs(i) * N_const
            
            ! Update temp_ao for printing
            temp_ao%prims(i)%a_nu    = exps(i)
            temp_ao%prims(i)%d_nu_mu = coefs(i)
        end do
        
        ! Print the orbital formula
        call print_ao(temp_ao, "")

        norm = 0.0_dp
        
        do i = 1, n
            do j = 1, n
                exponent_sum = (exps(i) + exps(j)) * 0.5_dp
                term_ovlp = overlap_spherical_component(l, m, exponent_sum)
                
                norm = norm + (norm_coefs(i) * norm_coefs(j) * term_ovlp)
            end do
        end do
        
        print *, "   -> Calculated Norm: ", norm
        print *, "----------------------------------------"
        
        deallocate(temp_ao%prims)
        deallocate(norm_coefs)
    end subroutine process_orbital

end program exercise_3_2
