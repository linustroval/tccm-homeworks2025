program main
        use basis_types
        implicit none
        integer :: i,k,m
        integer :: error
        integer :: unit_num
        integer :: target_Z, current_Z
        character(len=100) :: line
        integer :: l_curr,n_prim,n_cont
        real(dp), allocatable :: buffer(:,:)
        integer :: s_count,p_count,d_count

        unit_num=10

        print *, "=========================================="
        print *, "      DALTON BASIS READER (6-31G)         "
        print *, "=========================================="
        print *, "Enter the desired Atomic Number (Z):"
        print *, "(e.g., 1=H, 6=C, 8=O, 16=S...)"

        read(*,*) target_Z

        print *, "------------------------------------------"
        print *, "Searching data for Z =", target_Z, "..."
        print *, "------------------------------------------"

        open(unit=unit_num,file='6-31g.1.dalton',status='old',action='read')


        do
                read(unit_num,'(A)',iostat=error) line
                if (error/=0) then
                        print *, "CRITICAL ERROR: Atom Z=", target_Z, "not found"
                        print *, "Please ensure the number is between 1 and 36 in the file"
                        stop
                end if

                if (line(1:1) == 'a') then
                        read(line(3:),*) current_Z
                        if (current_Z == target_Z) exit
                end if
        end do


        s_count=0
        p_count=0
        d_count=0

        do 
                read(unit_num,'(A)', iostat=error) line
                if (error /= 0) exit
                if (line(1:1) == 'a') exit

                if (index(line, 's functions') > 0) then
                        l_curr=0
                else if (index(line, 'p functions') > 0) then
                        l_curr=1
                else if (index(line, 'd functions') > 0) then
                        l_curr=2
                else if (line(1:1) == 'H') then
                        read(line(2:),*) n_prim, n_cont

                        if (allocated(buffer)) deallocate(buffer)
                        allocate(buffer(n_prim, n_cont + 1))

                        do i=1,n_prim
                                read(unit_num,*) buffer(i,:)
                        end do


                        do k=1,n_cont
                                if (l_curr == 0) s_count = s_count + 1
                                if (l_curr == 1) p_count = p_count + 1
                                if (l_curr == 2) d_count = d_count + 1

                                do m=-l_curr,l_curr
                                        call create_print(l_curr, m, n_prim, buffer(:,1), buffer(:, k+1))
                                end do
                        end do
                end if  
        end do
        close(unit_num)
        if (allocated(buffer)) deallocate(buffer)
        
        print *, ""
        print *, "--- End of program ---"


 contains
        
        subroutine create_print(l,m,n,exps,coefs)
               integer, intent(in) :: l, m, n
               real(dp), intent(in) :: exps(:), coefs(:)
               type(contracted) :: temp_ao
               integer :: i
               temp_ao%l=l
               temp_ao%m=m
               temp_ao%n_prims=n

               if (l == 0) temp_ao%idx = s_count
               if (l == 1) temp_ao%idx = p_count + 1
               if (l == 2) temp_ao%idx = d_count + 2

               allocate(temp_ao%prims(n))
               do i=1,n
                        temp_ao%prims(i)%a_nu=exps(i)
                        temp_ao%prims(i)%d_nu_mu=coefs(i)
               end do

               call print_ao(temp_ao, "")

               deallocate(temp_ao%prims)

        end subroutine create_print

end program main



