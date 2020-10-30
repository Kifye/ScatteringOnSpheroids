module matrix
	use regime
	implicit none
	private
	
	public :: solve_system, inverse_matrix, multiply_by_diag_left, multiply_by_diag_right, check_symmetric
contains
	
!   Swaps rows i and j of the matrix A nxn
!   Used in LU_gaussian	
	subroutine swap_rows(A, n, i, j)
		complex(knd) A(n,n), v
		integer :: n, i, j, k

		do k = 1, n
			v = A(i, k)
			A(i, k) = A(j, k)
			A(j, k) = v
		enddo
		
	end subroutine swap_rows
	
!   Performs LU decomposition of the matrix A nxn to lower L and upper U
!   with P as a row permutation matrix
	subroutine LU_gaussian(A, n, L, U, P)
		integer :: n, i, j, k
		complex(knd) :: A(n, n), L(n, n), U(n, n), P(n, n), v
		real(knd), parameter :: EPS = 1d-6
		
		U = A
		L = 0
		P = 0
		do k = 1, n
			P(k, K) = qcmplx(1q0, 0q0)
		enddo
		do k = 1, n
			L(k, k) = qcmplx(1q0, 0q0)
			
			if (abs(u(k, k)) < EPS) then
				do i = k + 1, n
					if (abs(u(i, k)) > EPS) then
!						write(*,*) 'we found 1!'
						call swap_rows(U, n, i, k)
						call swap_rows(P, n, i, k)
						exit
					endif
				enddo
				if (i > n) then
					write(*,*) 'matrix is degenerate'
				endif
			endif
			
			do i = k + 1, n
!				write(*,*) 'i = ', i, 'k = ', k, 'u = ', u(i, k)
				L(i, k) = u(i, k) / u(k, k)
				do j = 1, n
					u(i, j) = u(i, j) - L(i, k) * u(k, j)
				enddo
			enddo
		enddo
		
	end subroutine LU_gaussian

!   Solves the linear system Ax = b, A:nxn, b:n
!   result in the last parameter res	
	subroutine solve_system(A, n, b, res)
		integer :: n, i, j, k
		complex(knd) :: A(n, n), L(n, n), U(n, n), P(n, n), b(n), res(n), res_low(n)
		
		call LU_gaussian(A, n, L, U, P)
		b = matmul(P, b)
		write(*,*) 'check LU'
		! write(*,*) 'L = '
		! write(*,*) L
		! write(*,*) 'U = '
		! write(*,*) U
		! write(*,*) 'LU = '
		! write(*,*) matmul(L, U)
		! write(*,*) 'A = '
		! write(*,*) A
		! write(*,*)
		call solve_lower_system(L, n, b, res_low)
		call solve_upper_system(U, n, res_low, res)
		! write(*,*) 'check LUres'
		! write(*,*) 'LUres = '
		! write(*,*) matmul(a, res) - b
		! write(*,*) 'b = '
		! write(*,*) b
		! write(*,*)

		! write(*,*) 'endres = '
		! write(*,*) res
		! write(*,*)
		
	end subroutine solve_system

	subroutine solve_lower_system(A, n, b, res)
		integer :: n, i, j
		complex(knd) :: A(n, n), b(n), res(n)
		
		res = b
		
		do i = 1, n
			do j = 1, i - 1
				res(i) = res(i) - res(j) * A(i, j)
			enddo
			res(i) = res(i) / A(i, i)
		enddo
		
		! write(*,*) 'check lower_solver'
		! write(*,*) 'Lres = '
		! write(*,*) matmul(A, res)
		! write(*,*) 'b = '
		! write(*,*) b
		! write(*,*)
		
		
	end subroutine solve_lower_system

	subroutine solve_upper_system(A, n, b, res)
		integer :: n, i, j
		complex(knd) :: A(n, n), b(n), res(n)
		
		res = b
		
		do i = n, 1, -1
			do j = i + 1, n
				res(i) = res(i) - res(j) * A(i, j)
			enddo
			res(i) = res(i) / A(i, i)
		enddo
		! write(*,*) 'check upper_solver'
		! write(*,*) 'Ures = '
		! write(*,*) matmul(A, res)
		! write(*,*) 'b = '
		! write(*,*) b
		! write(*,*)
		
	end subroutine solve_upper_system

	subroutine inverse_matrix(A, n, res)
		integer :: n, i, j, k
		complex(knd) :: A(n, n), L(n, n), U(n, n), P(n, n), res(n, n), res_low(n)
		
		call LU_gaussian(A, n, L, U, P)
				
		do k = 1, n
			
			call solve_lower_system(L, n, P(k,:), res_low)
			call solve_upper_system(U, n, res_low, res(:, k))

		enddo
		
	end subroutine inverse_matrix
	
	subroutine multiply_by_diag_right(A, n, diag_matr)
		complex(knd) A(n, n), diag_matr(n)
		integer :: n, i, j
		do i = 1, n
			do j = 1, n
!				write(*,*) 'Aj before'
!				write(*,*) A(1:4,j)
				A(i, j) = A(i, j) * diag_matr(j)
!				write(*,*) 'Aj after'
!				write(*,*) A(1:4,j)
			enddo
		enddo
	end subroutine multiply_by_diag_right
	
	subroutine multiply_by_diag_left(A, n, diag_matr)
		complex(knd) A(n, n), diag_matr(n)
		integer :: n, i, j
		do i = 1, n
			do j = 1, n
				A(i,j) = diag_matr(i) * A(i,j)
			enddo
		enddo
	end subroutine multiply_by_diag_left

	subroutine check_symmetric(A, n)
		integer :: n, i, j
		complex(knd) :: A(n,n)
		logical symm
		
		symm = .true.
		do i = 1, n
			do j = i + 1, n
				if (A(i,j) /= A(j,i)) then
					symm = .false.
					write(*,*) 'i =', i, 'j =', j, A(i,j), '!=', A(j,i)
				endif
			enddo
		enddo
		if (symm) then
			write(*,*) 'symmetrical'
		else
			write(*,*) 'not symmetrical'
		endif
	end subroutine check_symmetric
end module matrix