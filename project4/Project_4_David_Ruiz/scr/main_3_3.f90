program exercise_3_3_final
    use basis_types
    use mmd_utils
    implicit none

    integer :: i, j, m, error, unit_num, out_file
    integer :: target_Z, current_Z, n_atom_orbitals, n_total
    integer :: l_curr, n_prim, n_cont, ao_idx
    character(len=100) :: line
    real(dp), allocatable :: buffer(:,:)
    
    ! Molecule Database
    type(contracted), allocatable :: mol_basis(:) 
    real(dp), allocatable :: S_matrix(:,:)
    real(dp) :: norm_factor
    
    ! Molecule Configuration
    real(dp) :: dist_R = 1.4_dp 

    ! File Configuration
    unit_num = 10  ! Input 
    out_file = 20  ! Output 

    print *, "==============================================="
    print *, "   EXERCISE 3.3: MMD DIATOMIC OVERLAP"
    print *, "==============================================="
    print *, "Enter Atomic Number (Z) to build dimer (e.g., 1=H2, 6=C2):"
    read(*,*) target_Z

    open(unit=unit_num, file='6-31g.1.dalton', status='old', action='read')
    
    ! ---------------------------------------------------------
    ! 1. FIND THE CORRECT ATOM
    ! ---------------------------------------------------------
    do
        read(unit_num, '(A)', iostat=error) line
        if (error/=0) stop "Error: Atom not found in file."
        if (line(1:1) == 'a') then
             read(line(3:),*) current_Z
             if (current_Z == target_Z) exit
        end if
    end do
    
    ! ---------------------------------------------------------
    ! 2. COUNT ORBITALS 
    ! ---------------------------------------------------------
    n_atom_orbitals = 0
    l_curr = -1
    
    do
        read(unit_num, '(A)', iostat=error) line
        if (error/=0 .or. line(1:1)=='a') exit
        
        ! Detect L type
        if (index(line, 's functions') > 0) l_curr = 0
        if (index(line, 'p functions') > 0) l_curr = 1
        if (index(line, 'd functions') > 0) l_curr = 2
        
        ! Process H block
        if (line(1:1)=='H') then
            read(line(2:),*) n_prim, n_cont
            do i=1,n_prim; read(unit_num,*); end do ! Skip data
            
            if (l_curr >= 0) then
                n_atom_orbitals = n_atom_orbitals + (n_cont * (2*l_curr + 1))
            end if
        end if
    end do
    
    if (n_atom_orbitals == 0) then
        print *, "ERROR: No orbitals detected for this atom."
        stop
    end if

    ! ---------------------------------------------------------
    ! 3. ALLOCATE MEMORY
    ! ---------------------------------------------------------
    n_total = n_atom_orbitals * 2
    allocate(mol_basis(n_total))
    allocate(S_matrix(n_total, n_total))
    
    print *, "-> Orbitals per atom:", n_atom_orbitals
    print *, "-> Building molecule with", n_total, "total orbitals."
    print *, "-> Atom A at (0,0,0) | Atom B at (0,0,1.4)"
    
    ! ---------------------------------------------------------
    ! 4. READ REAL DATA AND BUILD MOLECULE
    ! ---------------------------------------------------------
    rewind(unit_num)
    ! Find atom again
    do; read(unit_num,'(A)') line; if(line(1:1)=='a') then; read(line(3:),*) current_Z; if(current_Z==target_Z) exit; endif; end do
    
    ao_idx = 0
    l_curr = -1
    
    do
        read(unit_num, '(A)', iostat=error) line
        if (error/=0 .or. line(1:1)=='a') exit
        
        if (index(line,'s functions')>0) l_curr=0
        if (index(line,'p functions')>0) l_curr=1
        if (index(line,'d functions')>0) l_curr=2
        
        if (line(1:1) == 'H') then
            read(line(2:),*) n_prim, n_cont
            
            if (allocated(buffer)) deallocate(buffer)
            allocate(buffer(n_prim, n_cont+1))
            
            do i=1, n_prim; read(unit_num, *) buffer(i,:); end do
            
            do m=1, n_cont
               do i = -l_curr, l_curr
                  
                  ! --- ATOM A (Origin) ---
                  ao_idx = ao_idx + 1
                  if (ao_idx > n_atom_orbitals) stop "Error: AO Index out of bounds"
                  
                  mol_basis(ao_idx)%l = l_curr
                  mol_basis(ao_idx)%m = i
                  mol_basis(ao_idx)%idx = ao_idx
                  mol_basis(ao_idx)%n_prims = n_prim
                  mol_basis(ao_idx)%center = [0.0_dp, 0.0_dp, 0.0_dp]
                  
                  allocate(mol_basis(ao_idx)%prims(n_prim))
                  do j=1, n_prim
                      mol_basis(ao_idx)%prims(j)%a_nu = buffer(j,1)
                      mol_basis(ao_idx)%prims(j)%d_nu_mu = buffer(j,m+1)
                  end do
                  
                  ! Normalize Atom A 
                  norm_factor = 1.0_dp/sqrt(overlap_mmd(mol_basis(ao_idx), mol_basis(ao_idx)))
                  mol_basis(ao_idx)%prims(:)%d_nu_mu = mol_basis(ao_idx)%prims(:)%d_nu_mu * norm_factor
                  
                  ! --- ATOM B (Displaced) ---
                  mol_basis(ao_idx + n_atom_orbitals) = mol_basis(ao_idx) 
                  
                  ! Update properties for Atom B
                  mol_basis(ao_idx + n_atom_orbitals)%idx = ao_idx + n_atom_orbitals
                  mol_basis(ao_idx + n_atom_orbitals)%center = [0.0_dp, 0.0_dp, dist_R] 
               end do
            end do
        end if
    end do
    
    if (allocated(buffer)) deallocate(buffer)
    close(unit_num)

    ! ---------------------------------------------------------
    ! 5. CALCULATE FULL OVERLAP MATRIX
    ! ---------------------------------------------------------
    print *, "Calculating Overlap Matrix..."
    do i = 1, n_total
        do j = 1, n_total
            S_matrix(i,j) = overlap_mmd(mol_basis(i), mol_basis(j))
        end do
    end do

    ! ---------------------------------------------------------
    ! 6. SAVE RESULTS TO FILE
    ! ---------------------------------------------------------
    open(unit=out_file, file='overlap_matrix.txt', status='replace', action='write')
    
    print *, "Writing 'overlap_matrix.txt'..."
    
    write(out_file, *) "OVERLAP MATRIX (MMD SCHEME)"
    write(out_file, *) "Atom A (0,0,0) - Atom B (0,0,", dist_R, ")"
    write(out_file, *) "Dimension: ", n_total, "x", n_total
    write(out_file, *) "------------------------------------------------"
    
    do i = 1, n_total
        write(out_file, '(100F9.4)') (S_matrix(i,j), j=1, n_total)
    end do
    
    close(out_file)
    print *, "Done!"

end program exercise_3_3_final
