module scattering

	use regime
	use constants
	use matrix
	use integrals
	use spheroidal
	
	implicit none
	private

	type,public :: Scatterer
	private
		
		type(SpheroidalCalculation) :: layer(0:NUMBER_OF_LAYERS)
		
		real(knd) :: ksi(NUMBER_OF_LAYERS), alpha, &
		mu(0:NUMBER_OF_LAYERS), eps(0:NUMBER_OF_LAYERS), ext(MODES), sca(MODES)
		
		complex(knd), allocatable, dimension(:,:,:) :: A11, A31, A31inv, Tmatr
		complex(knd), allocatable, dimension(:,:) :: solution, Delta
		integer :: matrix_size, f, accuracy = 100
		
		logical functions_calculated, matrices_calculated(MODES), tmatrix_calculated(MODES), &
		solution_calculated(MODES), ext_calculated(MODES), sca_calculated(MODES)

	contains
	
		procedure, public :: set, get_tmatrix, get_solution, &
		get_extinction_factor, get_scattering_factor
	
		procedure, private :: calculate_functions, calculate_matrices, calculate_tmatrix, &
		calculate_solution, get_initial, &
		common_factor, calculate_extinction_factor, calculate_scattering_factor
		
		final :: delete_scatterer
	end type
					
contains
	subroutine set(this, ni, lambda, rv, ab, alpha, matrix_size, eps, mu, f)
	
		class(Scatterer) :: this
		real(knd) :: lambda, rv, ab, alpha, eps(MODES), mu(MODES), arg(1)
		complex(knd) :: ni, c0, c
		integer :: matrix_size, f
	
		this%f = f
		this%mu = mu
		this%eps = eps
		this%alpha = alpha
		if (f == 1) then
			this%ksi(1) = ab / qsqrt(ab * ab - 1) 
			c0 = (rv / ab**2q0)**(1q0 / 3q0)
		else
			this%ksi(1) = 1q0 / qsqrt(ab * ab - 1)
			c0 = (rv * ab)**(1q0 / 3q0)
		endif
		c0 = 2q0 * PI / lambda * qsqrt(ab**2q0 - 1q0) * c0
		c = ni * c0
		
		if (allocated(this%A11) .and. matrix_size /= this%matrix_size) then 
			deallocate(this%A11, this%A31, this%A31inv, this%Tmatr, this%Delta, &
			this%solution)
		endif
		write(*,*) 'ready to allocate scatterer'
		if (.not. allocated(this%A11)) then
			allocate(this%A11(matrix_size, matrix_size, MODES))
			allocate(this%A31(matrix_size, matrix_size, MODES))
			allocate(this%A31inv(matrix_size, matrix_size, MODES))
			
			allocate(this%Tmatr(matrix_size, matrix_size, MODES))
			allocate(this%Delta(matrix_size, matrix_size))
			
			allocate(this%solution(matrix_size, MODES))
		endif
		write(*,*) 'allocated scatterer'

		this%matrix_size = matrix_size
		
		arg(1) = qcos(alpha)
		call this%layer(0)%set(1, matrix_size, c0, this%ksi(1), 1, arg, f)
		call this%layer(1)%set(1, matrix_size, c, this%ksi(1), 1, arg, f)
		
		this%functions_calculated = .false.
		this%matrices_calculated = .false.
		this%tmatrix_calculated = .false.
		this%solution_calculated = .false.
		this%ext_calculated = .false.
		this%sca_calculated = .false.
		
	end subroutine set

	subroutine calculate_functions(this)
		class(Scatterer) :: this
		integer :: i
		do i = 0, NUMBER_OF_LAYERS
			call this%layer(i)%calculate()
		enddo
		
		this%functions_calculated = .true.
		
	end subroutine calculate_functions
	
	subroutine calculate_matrices(this, mode)
	
		class(Scatterer) :: this
		
		integer :: n, mode
		complex(knd), allocatable, dimension(:,:) :: adder
		complex(knd), allocatable, dimension(:) :: R11, R31, R12, W1
		complex(knd) :: mu(0:1)
		
		allocate(adder(this%matrix_size, this%matrix_size), R11(this%matrix_size), &
		R31(this%matrix_size), R12(this%matrix_size), W1(this%matrix_size))
		
		if (.not.(this%functions_calculated)) then 
			call calculate_functions(this)
		endif
		
		mu = this%mu
		if (mode == 2) then
			mu = this%eps
		endif
		
		R11 = this%layer(0)%r1d / this%layer(0)%r1
		R31 = this%layer(0)%r3d / this%layer(0)%r3
		R12 = this%layer(1)%r1d / this%layer(1)%r1
		W1 =  -1q0 / (R31 - R11)
		
		call calculate_delta(this%layer(0), this%layer(1), this%Delta, this%accuracy)
		
		this%A11(:,:,mode) = (1q0 - mu(0) / mu(1)) * this%ksi(1) / (this%ksi(1) ** 2 - this%f) * &
		this%Delta
		adder = -(mu(0) / mu(1)) * this%Delta	
		call multiply_by_diag_right(adder, this%matrix_size, R12)
		this%A11(:,:,mode) = this%A11(:,:,mode) + adder
		this%A31(:,:,mode) = this%A11(:,:,mode)
		
		adder = this%Delta
		call multiply_by_diag_left(adder, this%matrix_size, R11)
		this%A11(:,:,mode) = this%A11(:,:,mode) + adder

		adder = this%Delta
		call multiply_by_diag_left(adder, this%matrix_size, R31)
		this%A31(:,:,mode) = this%A31(:,:,mode) + adder
		
		call multiply_by_diag_left(this%A11(:,:,mode), this%matrix_size, W1)
		call multiply_by_diag_left(this%A31(:,:,mode), this%matrix_size, W1)

		deallocate(adder, R11, R31, R12, W1)

!		write(*,*) 'A11 = '
!		write(*,*) this%A11(1:4, 1:4)
!		write(*,*) 'A31 = '
!		write(*,*) this%A31(1:4, 1:4)
		this%matrices_calculated(mode) = .true.

	end subroutine calculate_matrices

	subroutine calculate_tmatrix(this, mode)
		class(Scatterer) :: this
		integer :: mode
		
		if (.not. this%matrices_calculated(mode)) then
			call this%calculate_matrices(mode)
		endif
		
		call inverse_matrix(this%A31(:,:,mode), this%matrix_size, this%A31inv(:,:,mode))

		!  write(*,*) 'check inverse'
		!  write(*,*) matmul(this%A31, this%A31inv)
		!  write(*,*)
		this%Tmatr(:,:,mode) = -matmul(this%A11(:,:,mode), this%A31inv(:,:,mode))
		
		this%tmatrix_calculated(mode) = .true.
	end subroutine calculate_tmatrix

	function get_tmatrix(this, mode) result(tmatr)
		class(Scatterer) :: this
		integer :: mode
		complex(knd) :: tmatr(this%matrix_size, this%matrix_size)
		
		if (.not. this%tmatrix_calculated(mode)) then
			call this%calculate_tmatrix(mode)
		endif
		
		tmatr = this%Tmatr(:,:,mode)
		
	end function get_tmatrix

		
!   mode == 1 for TE
!   mode == 2 for TM
	subroutine get_initial(this, res, mode)
		class(Scatterer) :: this
		complex(knd) :: res(this%matrix_size)
		integer :: mode, i
		
		do i = 1, this%matrix_size
			if (mode == 1) then
				res(i) = -2q0 * (qcmplx(0q0, 1q0) ** i) * &
				this%layer(0)%s1(i, 1) * this%layer(0)%r1(i)
!				write(*,*) 'assigned res(i) = ', res(i)
!				write(*,*) 's = ', this%layer(0)%s1(i, 1), 'r = ', this%layer(0)%r1(i)
			else
				res(i) = 2q0 * qsqrt(this%eps(0) / this%mu(0)) * (qcmplx(0q0, 1q0) ** i) * &
				this%layer(0)%s1(i, 1) * this%layer(0)%r1(i)
			endif
		enddo
	end subroutine get_initial
	
	
	subroutine calculate_solution(this, mode)
		class(Scatterer) :: this
		complex(knd) init(this%matrix_size)
		integer :: mode

		if (.not. this%tmatrix_calculated(mode)) then 
			call calculate_tmatrix(this, mode)
		endif
		
		call get_initial(this, init, mode)
		
		call solve_system(this%Tmatr(:,:,mode), this%matrix_size, init, this%solution(:, mode))
		
		! write(*,*) 'check solution:'
		! write(*,*) 'init:'
		! write(*,*) init
		! write(*,*) 'result:'
		! write(*,*) matmul(this%Tmatr, this%solution(:, mode))
		! write(*,*) 'check solution'
		
		this%solution(:,mode) = matmul(this%Tmatr(:,:,mode), init)
		
		
		this%solution_calculated(mode) = .true.
		
	end subroutine calculate_solution
	
	function get_solution(this, mode) result(solution)
		
		class(Scatterer) :: this
		integer :: mode
		complex(knd) :: solution(this%matrix_size)
		
		if (.not. this%solution_calculated(mode)) then
			call this%calculate_solution(mode)
		endif
		
		solution = this%solution(:, mode)
	end function get_solution
		
	real(knd) function common_factor(this)
		class(Scatterer), intent(in) :: this
		
		common_factor = 1q0 / (cqabs(this%layer(0)%c) * &
		qsqrt((this%ksi(1)**2 - this%f) * (this%ksi(1)**2 - this%f * qcos(this%alpha)**2)))
	
	end function common_factor
	
	subroutine calculate_extinction_factor(this, mode)
		class(Scatterer) :: this
		integer :: mode, i
		complex(knd) :: ideg
		
		if (.not. this%solution_calculated(mode)) then
			call this%calculate_solution(mode)
		endif
		
		this%ext(mode) = 0q0
		
		ideg = qcmplx(0q0, -1q0)
		
		do i = 1, this%matrix_size
			this%ext(mode) = this%ext(mode) + &
			real(this%solution(i, mode) * ideg**i * this%layer(0)%s1(i, 1) &
			/ this%layer(0)%r3(i), 16)
!			ideg = ideg * qcmplx(0q0, -1q0)
		enddo
		
		this%ext(mode) = -this%ext(mode) * 4q0 * this%common_factor()
		
	end subroutine calculate_extinction_factor

	function get_extinction_factor(this, mode) result(ext)
		class(Scatterer) :: this
		integer :: mode
		real(knd) :: ext

		if (.not. this%ext_calculated(mode)) then
			call this%calculate_extinction_factor(mode)
		endif
		ext = this%ext(mode)
	end function get_extinction_factor

	subroutine calculate_scattering_factor(this, mode)
		class(Scatterer) :: this
		integer :: mode, i
		complex(knd) :: ideg
		
		if (.not. this%solution_calculated(mode)) then
			call this%calculate_solution(mode)
		endif
		
		this%sca(mode) = 0q0
		
		do i = 1, this%matrix_size
			this%sca(mode) = this%sca(mode) + &
			cqabs(this%solution(i, mode))**2 / cqabs(this%layer(0)%r3(i))**2
		enddo
		
		this%sca(mode) = this%sca(mode) * 2q0 * this%common_factor()
		
	end subroutine calculate_scattering_factor

	function get_scattering_factor(this, mode) result(sca)
		class(Scatterer) :: this
		real(knd) :: sca
		integer :: mode

		if (.not. this%sca_calculated(mode)) then
			call this%calculate_scattering_factor(mode)
		endif
		sca = this%sca(mode)
	end function get_scattering_factor
	
	subroutine delete_scatterer(this)

		type(Scatterer), intent(inout) :: this

		if (allocated(this%A11)) then
			deallocate(this%A11)
		endif
		if (allocated(this%A31)) then
			deallocate(this%A31)
		endif
		if (allocated(this%A31inv)) then
			deallocate(this%A31inv)
		endif
		if (allocated(this%Delta)) then
			deallocate(this%Delta)
		endif
		if (allocated(this%Tmatr)) then
			deallocate(this%Tmatr)
		endif
		if	(allocated(this%solution)) then
			deallocate(this%solution)
		endif

	end subroutine delete_scatterer
end module scattering