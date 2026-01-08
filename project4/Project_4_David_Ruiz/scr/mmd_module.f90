module mmd_utils
    use basis_types
    implicit none

contains

    ! ----------------------------------------------------------------
    ! MAIN FUNCTION: OVERLAP BETWEEN TWO ORBITALS (MMD)
    ! Calculates <Chi_A | Chi_B> by summing over all primitives
    ! ----------------------------------------------------------------
    function overlap_mmd(ao1, ao2) result(S_val)
        type(contracted), intent(in) :: ao1, ao2
        real(dp) :: S_val
        
        integer :: i, j
        real(dp) :: term, Sx, Sy, Sz, prefactor
        real(dp) :: alpha, beta, p
        real(dp) :: A(3), B(3)
        integer :: ix, iy, iz, jx, jy, jz
        
        S_val = 0.0_dp
        A = ao1%center
        B = ao2%center
        
        ! Get Cartesian powers (nx, ny, nz) from quantum numbers (l, m)
        call get_powers(ao1, ix, iy, iz)
        call get_powers(ao2, jx, jy, jz)

        ! Double loop over primitives
        do i = 1, ao1%n_prims
            do j = 1, ao2%n_prims
                
                alpha = ao1%prims(i)%a_nu
                beta  = ao2%prims(j)%a_nu
                p = alpha + beta
                
                ! 1. Calculate 1D Recursive Overlap for each axis
                Sx = get_1d_overlap(ix, jx, alpha, beta, A(1), B(1))
                Sy = get_1d_overlap(iy, jy, alpha, beta, A(2), B(2))
                Sz = get_1d_overlap(iz, jz, alpha, beta, A(3), B(3))
                
                ! 2. Normalization factor for the Gaussian integral (pi/p)^1.5
                prefactor = (acos(-1.0_dp) / p)**1.5_dp
                
                ! 3. Sum term: Coef1 * Coef2 * Prefactor * Sx * Sy * Sz
                term = ao1%prims(i)%d_nu_mu * ao2%prims(j)%d_nu_mu * prefactor * Sx * Sy * Sz
                
                S_val = S_val + term
            end do
        end do
    end function overlap_mmd

    ! ----------------------------------------------------------------
    ! 1D RECURSION (Obara-Saika / MacMurchie-Davidson Scheme)
    ! Formula: <i|j> = (P-A)*<i-1|j> + 1/2p * ( i*<i-2|j> + j*<i-1|j-1> )
    ! ----------------------------------------------------------------
    recursive function get_1d_overlap(i, j, alpha, beta, Ax, Bx) result(val)
        integer, intent(in) :: i, j
        real(dp), intent(in) :: alpha, beta, Ax, Bx
        real(dp) :: val, p, Px, mu

        p = alpha + beta
        Px = (alpha * Ax + beta * Bx) / p
        
        ! --- Base Case: <0|0> (s-s orbitals) ---
        if (i == 0 .and. j == 0) then
            mu = (alpha * beta) / p
            val = exp(-mu * (Ax - Bx)**2)
            return
        end if
        
        val = 0.0_dp
        
        ! --- Recursive Step: Reduce index i ---
        if (i > 0) then
            val = (Px - Ax) * get_1d_overlap(i-1, j, alpha, beta, Ax, Bx)
            if (i > 1) val = val + (real(i-1,dp)/(2.0_dp*p)) * get_1d_overlap(i-2, j, alpha, beta, Ax, Bx)
            if (j > 0) val = val + (real(j,dp)/(2.0_dp*p))   * get_1d_overlap(i-1, j-1, alpha, beta, Ax, Bx)
            return
        end if
        
        ! --- Recursive Step: Reduce index j (if i is already 0) ---
        if (j > 0) then
            val = (Px - Bx) * get_1d_overlap(i, j-1, alpha, beta, Ax, Bx)
            if (j > 1) val = val + (real(j-1,dp)/(2.0_dp*p)) * get_1d_overlap(i, j-2, alpha, beta, Ax, Bx)
            return
        end if
        
    end function get_1d_overlap

    ! Helper: Converts (l, m) to Cartesian powers (nx, ny, nz)
    subroutine get_powers(ao, nx, ny, nz)
        type(contracted), intent(in) :: ao
        integer, intent(out) :: nx, ny, nz
        
        nx = 0; ny = 0; nz = 0
        
        ! Simple mapping for the exercise
        if (ao%l == 0) return ! s -> (0,0,0)
        
        if (ao%l == 1) then
            if (ao%m == 1)  nx = 1 ! x
            if (ao%m == -1) ny = 1 ! y
            if (ao%m == 0)  nz = 1 ! z
        else if (ao%l == 2) then
            ! Simplified d-shell mapping:
            if (ao%m == 0) nz = 2        ! z^2
            if (ao%m == 1) then; nx=1; nz=1; end if ! xz
            if (ao%m == -1) then; ny=1; nz=1; end if ! yz
            if (ao%m == 2) then; nx=2; end if ! x^2 - y^2 (approx x^2)
            if (ao%m == -2) then; nx=1; ny=1; end if ! xy
        end if
    end subroutine get_powers

end module mmd_utils
