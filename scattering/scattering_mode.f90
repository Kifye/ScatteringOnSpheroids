module scattering_mode
    use matrix
    use regime
    use scattering
    use spheroidal
    use matrix
    use elemfunc
    use constants
    use integrals
    implicit none
    private

    type, public, abstract :: ScatteringMode
        type(Scatterer), pointer :: scatterer
        complex(knd), allocatable, dimension(:, :) :: initial, solution  !  [1..(1 or maxm)][1..matrix_size]
        complex(knd), allocatable, dimension(:, :, :) :: tmatr, A11, A31, A31inv   !  [1..(1 or maxm)][1..matrix_size][1..matrix_size]
        integer :: matrix_size, maxm

        logical tmatrix_calculated, initial_calculated, solution_calculated, matrices_calculated

    contains
        procedure, public, non_overridable :: get_tmatrix, get_solution
        procedure, public :: set => set_scattering_mode
        procedure, private, non_overridable :: calculate_solution, calculate_tmatrix
        procedure(calculate_parameters), private, deferred :: calculate_initial, calculate_matrices
        procedure(get_parameters), private, deferred :: get_extinction, get_scattering
        procedure, public :: delete_mode => delete_scattering_mode

    end type ScatteringMode

    abstract interface
        subroutine calculate_parameters(this)
            import ScatteringMode
            class(ScatteringMode) :: this
        end subroutine
    end interface

    abstract interface
        function get_parameters(this)
            use regime
            import ScatteringMode
            class(ScatteringMode) :: this
            real(knd) :: get_parameters
        end function
    end interface

    type, public, abstract, extends(ScatteringMode) :: AxisymmetricMode
    contains
        procedure, private :: calculate_matrices => calculate_axisymmetric_matrices
        procedure, public :: get_extinction => get_symmetric_extinction_factor
        procedure, public :: get_scattering => get_symmetric_scattering_factor
    end type AxisymmetricMode

    type, public, extends(AxisymmetricMode) :: AxisymmetricTE
    contains
        procedure, private :: calculate_initial => calculate_axisymmetric_initial_te
        final :: delete_axisymmetric_mode_te
    end type AxisymmetricTE

    type, public, extends(AxisymmetricMode) :: AxisymmetricTM
    contains
        procedure, private :: calculate_initial => calculate_axisymmetric_initial_tm
        final :: delete_axisymmetric_mode_tm
    end type AxisymmetricTM

    type, public, abstract, extends(ScatteringMode) :: NonAxisymmetricMode
    contains
        procedure, private :: calculate_matrices => calculate_nonaxisymmetric_matrices
        procedure, public :: get_extinction => get_non_symmetric_extinction_factor
        procedure, public :: get_scattering => get_non_symmetric_scattering_factor
    end type NonAxisymmetricMode

    type, public, extends(NonAxisymmetricMode) :: NonAxisymmetricTE
    contains
        procedure, private :: calculate_initial => calculate_nonaxisymmetric_initial_te
        final :: delete_nonaxisymmetric_mode_te
    end type NonAxisymmetricTE

    type, public, extends(NonAxisymmetricMode) :: NonAxisymmetricTM
    contains
        procedure, private :: calculate_initial => calculate_nonaxisymmetric_initial_tm
        final :: delete_nonaxisymmetric_mode_tm
    end type NonAxisymmetricTM

contains

    subroutine set_scattering_mode(this, scat, matrix_size, maxm)

        class(ScatteringMode) :: this
        type(Scatterer), target :: scat
        integer :: matrix_size, maxm

        this%scatterer => scat

        if (allocated(this%initial) .and. (matrix_size /= this%matrix_size .or. maxm /= this%maxm)) then
            deallocate(this%initial, this%tmatr, this%solution, this%A11, this%A31, this%A31inv)
        endif
        if (.not. allocated(this%initial)) then
            allocate(this%initial(maxm, matrix_size), this%tmatr(maxm, matrix_size, matrix_size), &
                    this%solution(maxm, matrix_size))
            allocate(this%A11(maxm, matrix_size, matrix_size), this%A31(maxm, matrix_size, matrix_size), &
                    this%A31inv(maxm, matrix_size, matrix_size))
        endif
        this%matrix_size = matrix_size
        this%maxm = maxm

        this%initial_calculated = .false.
        this%tmatrix_calculated = .false.
        this%solution_calculated = .false.
        this%matrices_calculated = .false.

    end subroutine set_scattering_mode

    function get_tmatrix(this) result(tmatr)
        class(ScatteringMode) :: this
        complex(knd) :: tmatr(this%maxm, this%matrix_size, this%matrix_size)
        integer :: m

        !        write(*,*) 'ms = ', this%matrix_size
        !        write(*,*) 'maxm = ', this%maxm

        if (.not. this%tmatrix_calculated) then
            call this%calculate_tmatrix
        endif

        tmatr = this%tmatr
        return
        !        do m = 1, this%maxm
        !            call multiply_by_diag_right(tmatr(m, :, :), this%matrix_size, this%scatterer%layer(m, 0)%r1)
        !            call multiply_by_diag_left(tmatr(m, :, :), this%matrix_size, 1q0 / this%scatterer%layer(m, 0)%r3)
        !        end do

    end function get_tmatrix

    subroutine calculate_solution(this)
        class(ScatteringMode) :: this
        integer :: i

        if (.not. this%initial_calculated) then
            call this%calculate_initial()
        endif

        if (.not. this%tmatrix_calculated) then
            call this%calculate_tmatrix()
        endif

        do i = 1, this%maxm
            this%solution(i, :) = matmul(this%tmatr(i, :, :), this%initial(i, :))
        enddo
        this%solution_calculated = .true.

    end subroutine calculate_solution

    function get_solution(this) result(solution)

        class(ScatteringMode) :: this
        complex(knd) :: solution(this%maxm, this%matrix_size)

        if (.not. this%solution_calculated) then
            call this%calculate_solution()
        endif

        solution = this%solution

    end function get_solution


    subroutine calculate_axisymmetric_matrices(this)

        class(AxisymmetricMode) :: this

        integer :: n
        complex(knd), allocatable, dimension(:, :) :: adder
        complex(knd), allocatable, dimension(:) :: R11, R31, R12, W1
        complex(knd) :: mu(0:1)

        !		write(*,*) 'in matrices thismatrsize = ', this%matrix_size

        allocate(adder(this%matrix_size, this%matrix_size), R11(this%matrix_size), &
                R31(this%matrix_size), R12(this%matrix_size), W1(this%matrix_size))

        if (.not.(this%scatterer%functions_calculated)) then
            call this%scatterer%calculate_functions()
        endif

        mu = this%scatterer%mu
        select type (this)
        class is (AxisymmetricTM)
            mu = this%scatterer%eps
        end select

        R11 = this%scatterer%layer(1, 0)%r1d / this%scatterer%layer(1, 0)%r1
        R31 = this%scatterer%layer(1, 0)%r3d / this%scatterer%layer(1, 0)%r3
        R12 = this%scatterer%layer(1, 1)%r1d / this%scatterer%layer(1, 1)%r1
        W1 = -1q0 / (R31 - R11)

        !		call calculate_delta(this%layer(0), this%layer(1), this%Delta, this%accuracy)

        !		write(*,*) 'delta'
        !		write(*,*) this%Delta(1:10, 1)

        !		write(*,*) 'R12'
        !		write(*,*) R12
        !		write(*,*) 'W1', size(W1)
        !		write(*,*) W1
        !		write(1,*) 'mode = ', mode, 'mu = ', mu
        this%A11(1, :, :) = (1q0 - mu(0) / mu(1)) * this%scatterer%ksi(1) / (this%scatterer%ksi(1) ** 2 - this%scatterer%f) * &
                this%scatterer%Delta(1, :, :)
        !		write(*,*) 'A11_1 = '
        !		write(*,*) this%A11(:, 1:1, mode)
        adder = -(mu(0) / mu(1)) * this%scatterer%Delta(1, :, :)
        call multiply_by_diag_right(adder, this%matrix_size, R12)
        this%A11(1, :, :) = this%A11(1, :, :) + adder
        !		write(*,*) 'A11_2 = '
        !		write(*,*) this%A11(:, 1:1, mode)
        !		write(1,*) 'A11_1 mode =', mode
        !		write(1,*) this%A11(1:10, 1:1, mode)
        this%A31(1, :, :) = this%A11(1, :, :)

        adder = this%scatterer%Delta(1, :, :)
        call multiply_by_diag_left(adder, this%matrix_size, R11)
        !		write(*,*) 'delta'
        !		write(*,*) this%Delta(1:10, 1)
        this%A11(1, :, :) = this%A11(1, :, :) + adder
        !		write(*,*) 'A11_3 = '
        !		write(*,*) this%A11(:, 1:1, mode)

        adder = this%scatterer%Delta(1, :, :)
        call multiply_by_diag_left(adder, this%matrix_size, R31)
        this%A31(1, :, :) = this%A31(1, :, :) + adder

        call multiply_by_diag_left(this%A11(1, :, :), this%matrix_size, W1)
        !		write(*,*) 'A11_4 = '
        !		write(*,*) this%A11(:, 1:1, mode)
        call multiply_by_diag_left(this%A31(1, :, :), this%matrix_size, W1)

        deallocate(adder, R11, R31, R12, W1)

        call multiply_by_diag_left(this%A31(1, :, :), this%matrix_size, 1q0 / this%scatterer%layer(1, 0)%r1)
        call multiply_by_diag_left(this%A11(1, :, :), this%matrix_size, 1q0 / this%scatterer%layer(1, 0)%r3)

        !		write(1,*) 'A11 mode =', mode
        !		write(1,*) this%A11(1:10, 1:1, mode)
        !		write(1,*) 'A31 mode = ',mode
        !		write(1,*) this%A31(1:10, 1:1, mode)
        this%matrices_calculated = .true.

    end subroutine calculate_axisymmetric_matrices

    function get_part_11(f, ksi, eps, mu, left_multiplier, right_multiplier, &
            Delta, Q01, Q11, Epsilon, matrix_size) result(res)
        complex(knd) :: eps(0:1), mu(0:1)
        real(knd) :: ksi
        integer :: matrix_size, f
        complex(knd) :: left_multiplier(matrix_size), right_multiplier(matrix_size), &
                adder(matrix_size, matrix_size), Delta(matrix_size, matrix_size), Q01(matrix_size, matrix_size), &
                Q11(matrix_size, matrix_size), Epsilon(matrix_size, matrix_size), identity(matrix_size, matrix_size), &
                res(matrix_size, matrix_size)

        call get_identity_matrix(identity, matrix_size)

        res = -(eps(1) / eps(0) - mu(0) / mu(1)) * f * ksi / (ksi ** 2 - f) * matmul(Q01, Epsilon) - &
                (eps(1) / eps(0) - 1q0) * ksi * matmul(Q01, identity - 2 * ksi**2 * Q11)

        !write(*, *) 'res = ', res
        !write(*, *) 'res = ', transpose(res)
        !write(*, *) 'ksi = ', ksi, 'ksi**2 = ', ksi**2, ksi**2q0, ksi*ksi
        !write(*, *)
        adder = (mu(0) / mu(1) - 1q0) * ksi**2 * Q01 - mu(0) / mu(1) * Delta
        call multiply_by_diag_right(adder, matrix_size, right_multiplier)
        res = res + adder

        adder = Delta + (eps(1) / eps(0) - 1q0) * ksi ** 2 * Q01
        call multiply_by_diag_left(adder, matrix_size, left_multiplier)
        res = res + adder

    end function get_part_11

    function get_part_12(f, ksi, eps, mu, left_multiplier, right_multiplier, &
            Delta, Q01, Q11, Kappa, Gamma11, matrix_size) result(result)

        complex(knd) :: eps(0:1), mu(0:1), identity(matrix_size, matrix_size)
        real(knd) :: ksi
        integer :: matrix_size, f
        complex(knd) :: left_multiplier(matrix_size), right_multiplier(matrix_size), &
                adder(matrix_size, matrix_size), Delta(matrix_size, matrix_size), Q01(matrix_size, matrix_size), &
                Q11(matrix_size, matrix_size), Kappa(matrix_size, matrix_size), Gamma11(matrix_size, matrix_size), &
                result(matrix_size, matrix_size)

        call get_identity_matrix(identity, matrix_size)

        result = -(eps(1) / eps(0) - mu(0) / mu(1)) * f / (ksi ** 2 - f) * &
                (matmul((ksi**2 * Q01 - Delta), Kappa) + matmul(Delta, Gamma11)) + &
                (eps(1) / eps(0) - 1q0) * f * ksi**2 * 2q0 * matmul(matmul(Q01, Q11), Gamma11)
        adder = (mu(0) / mu(1) - 1q0) * f * ksi * matmul(Q01, Gamma11)
        call multiply_by_diag_right(adder, matrix_size, right_multiplier)
        result = result + adder

        adder = (eps(1) / eps(0) - 1q0) * f * ksi * matmul(Q01, Gamma11)
        call multiply_by_diag_left(adder, matrix_size, left_multiplier)
        result = result + adder

    end function get_part_12

    function get_part_21(f, ksi, eps, mu, left_multiplier, right_multiplier, &
            Delta, Q01, Q11, Kappa, Gamma11, matrix_size) result(result)

        complex(knd) :: eps(0:1), mu(0:1)
        real(knd) :: ksi
        integer :: matrix_size, f
        complex(knd) :: left_multiplier(matrix_size), right_multiplier(matrix_size), &
                adder(matrix_size, matrix_size), Delta(matrix_size, matrix_size), Q01(matrix_size, matrix_size), &
                Q11(matrix_size, matrix_size), Kappa(matrix_size, matrix_size), Gamma11(matrix_size, matrix_size), &
                result(matrix_size, matrix_size)

        result = (eps(1) / eps(0) - mu(0) / mu(1)) * ksi**2 / (ksi ** 2 - f) * matmul(Q01, Kappa) - &
                (eps(1) / eps(0) - 1q0) * ksi**2 * 2q0 * matmul(matmul(Q01, Q11), Gamma11)
        adder = -(mu(0) / mu(1) - 1q0) * ksi * matmul(Q01, Gamma11)
        call multiply_by_diag_right(adder, matrix_size, right_multiplier)
        result = result + adder

        adder = -(eps(1) / eps(0) - 1q0) * ksi * matmul(Q01, Gamma11)
        call multiply_by_diag_left(adder, matrix_size, left_multiplier)
        result = result + adder

    end function get_part_21

    function get_part_22(f, ksi, eps, mu, left_multiplier, right_multiplier, &
            Delta, Q01, Q11, Epsilon, matrix_size) result(result)
        complex(knd) :: eps(0:1), mu(0:1)
        real(knd) :: ksi
        integer :: matrix_size, f
        complex(knd) :: left_multiplier(matrix_size), right_multiplier(matrix_size), &
                adder(matrix_size, matrix_size), Delta(matrix_size, matrix_size), Q01(matrix_size, matrix_size), &
                Q11(matrix_size, matrix_size), Epsilon(matrix_size, matrix_size), identity(matrix_size, matrix_size), &
                result(matrix_size, matrix_size)

        call get_identity_matrix(identity, matrix_size)

        result = (eps(1) / eps(0) - mu(0) / mu(1)) * ksi / (ksi ** 2 - f) * (f * matmul(Q01, Epsilon) + Delta) + &
                (eps(1) / eps(0) - 1q0) * ksi * matmul(Q01, identity - 2 * ksi**2 * Q11)

        adder = -(mu(0) / mu(1) - 1q0) * ksi**2 * Q01 - Delta
        call multiply_by_diag_right(adder, matrix_size, right_multiplier)
        result = result + adder

        adder = eps(1) / eps(0) * Delta - (eps(1) / eps(0) - 1q0) * ksi ** 2 * Q01
        call multiply_by_diag_left(adder, matrix_size, left_multiplier)
        result = result + adder

    end function get_part_22

    subroutine set_full_matrix(f, ksi, eps, mu, first_multiplier, left_multiplier, right_multiplier, &
            Delta, Q01, Q11, Kappa, Gamma11, Epsilon, matrix_size, result)

        complex(knd) :: eps(0:1), mu(0:1)
        real(knd) :: ksi
        integer :: matrix_size, f, i, j
        complex(knd) :: first_multiplier(matrix_size), left_multiplier(matrix_size), right_multiplier(matrix_size), &
                adder(matrix_size, matrix_size), Delta(matrix_size, matrix_size), Q01(matrix_size, matrix_size), &
                Q11(matrix_size, matrix_size), Kappa(matrix_size, matrix_size), Gamma11(matrix_size, matrix_size), &
                Epsilon(matrix_size, matrix_size), result(2 * matrix_size, 2 * matrix_size)

        !write(*,*) 'Delta = ', Delta
        !write(*,*) 'Gamma = ', Gamma11
        !write(*,*) 'Kappa = ', Kappa
        !write(*,*) 'Epsilon = ', Epsilon
        !write(*,*) 'Q01 = ', Q01
        !write(*,*) 'Q11 = ', Q11
        !write(*,*)
        ! write(*,*) 'ms = ', matrix_size
        result(1:matrix_size, 1:matrix_size) = &
                get_part_11(f, ksi, eps, mu, left_multiplier, right_multiplier, Delta, Q01, Q11, Epsilon, matrix_size)
        result(1:matrix_size, (matrix_size + 1):(2 * matrix_size)) = &
                get_part_12(f, ksi, eps, mu, left_multiplier, right_multiplier, Delta, Q01, Q11, Kappa, Gamma11, matrix_size)
        result((matrix_size + 1):(2 * matrix_size), 1:matrix_size) = &
                get_part_21(f, ksi, eps, mu, left_multiplier, right_multiplier, Delta, Q01, Q11, Kappa, Gamma11, matrix_size)
        result((matrix_size + 1):(2 * matrix_size), (matrix_size + 1):(2 * matrix_size)) = &
                get_part_22(f, ksi, eps, mu, left_multiplier, right_multiplier, Delta, Q01, Q11, Epsilon, matrix_size)
        ! write(*,*) 'befor matr = ', result(1:5,1:5)
        do i = 0, 1
            do j = 0, 1
                call multiply_by_diag_left(&
                        result((i * matrix_size + 1):((i + 1) * matrix_size), (j * matrix_size + 1):((j + 1) * matrix_size)), &
                        matrix_size, first_multiplier)
            end do
        end do

        ! write(*,*) 'after matr = ', result(1:5,1:5)
        ! write(*,*)
    end subroutine set_full_matrix

    subroutine calculate_nonaxisymmetric_matrices(this)

        class(NonAxisymmetricMode) :: this

        integer :: n, i, j, m
        complex(knd), allocatable, dimension(:, :) :: adder, identity
        complex(knd), allocatable, dimension(:) :: R11, R31, R12, W1
        complex(knd) :: mu(0:1), eps(0:1), initial_corrector(this%matrix_size), solution_corrector(this%matrix_size)
        complex(knd) :: k1, c1

        !		write(*,*) 'in matrices thismatrsize = ', this%matrix_size

        n = this%matrix_size / 2
        allocate(adder(n, n), identity(n, n), R11(n), R31(n), R12(n), W1(n))
        call get_identity_matrix(identity, n)

        if (.not.(this%scatterer%functions_calculated)) then
            call this%scatterer%calculate_functions()
        endif

        eps = this%scatterer%eps
        mu = this%scatterer%mu
        select type (this)
        class is (NonAxisymmetricTM)
            mu = this%scatterer%eps
            eps = this%scatterer%mu
        end select

        this%A11 = 0
        this%A31 = 0
        this%A31inv = 0
        !write(*,*) 'check_a11'
        !do j = 1, 3
        !    write(*,*) this%A31(1, j, 1:3)
        !end do
        !write(*,*)

        k1 = 2 * PI / this%scatterer%lambda
        do m = 1, this%maxm
            c1 = this%scatterer%layer(m, 0)%c
            write(*, *) 'k1 = ', k1, ' c1 = ', c1
            !write(*,*) 'check_a11'
            !do j = 1, 3
            !    write(*,*) this%A31(1, j, 1:3)
            !end do
            !write(*,*)

            R11 = this%scatterer%layer(m, 0)%r1d / this%scatterer%layer(m, 0)%r1
            R31 = this%scatterer%layer(m, 0)%r3d / this%scatterer%layer(m, 0)%r3
            R12 = this%scatterer%layer(m, 1)%r1d / this%scatterer%layer(m, 1)%r1
            W1 = -1q0 / (R31 - R11)
            !write(*,*) 'W1 =', qcmplx(0q0, 1q0) * c1 * (this%scatterer%ksi(1)**2 - this%scatterer%f) * &
            !        this%scatterer%layer(m, 0)%r1 * this%scatterer%layer(m, 0)%r3
            call set_full_matrix(this%scatterer%f, this%scatterer%ksi(1), eps, mu, &
                    W1, R11, R12, &
                    this%scatterer%Delta(m, :, :), this%scatterer%Q01(m, :, :), this%scatterer%Q11(m, :, :), &
                    this%scatterer%Kappa(m, :, :), this%scatterer%Gamma11(m, :, :), this%scatterer%Epsilon(m, :, :), &
                    n, this%A11(m, :, :))
            call set_full_matrix(this%scatterer%f, this%scatterer%ksi(1), eps, mu, &
                    W1, R31, R12, &
                    this%scatterer%Delta(m, :, :), this%scatterer%Q01(m, :, :), this%scatterer%Q11(m, :, :), &
                    this%scatterer%Kappa(m, :, :), this%scatterer%Gamma11(m, :, :), this%scatterer%Epsilon(m, :, :), &
                    n, this%A31(m, :, :))

            initial_corrector(1:n) = 1q0 / (this%scatterer%layer(m, 0)%r1 * k1)
            initial_corrector((n + 1):(2 * n)) = 1q0 / (this%scatterer%layer(m, 0)%r1 * c1)
            solution_corrector(1:n) = 1q0 / (this%scatterer%layer(m, 0)%r3 * k1)
            solution_corrector((n + 1):(2 * n)) = 1q0 / (this%scatterer%layer(m, 0)%r3 * c1)

            !write(*,*) 'sizes:', size(initial_corrector), size(this%scatterer%layer(m, 0)%r1), size(this%A11(m,1,:))
            call multiply_by_diag_left(this%A31(m, :, :), this%matrix_size, initial_corrector)
            call multiply_by_diag_left(this%A11(m, :, :), this%matrix_size, solution_corrector)
            !write(*,*) 'check_a11'
            !do j = 1, 3
            !    write(*,*) this%A31(1, j, 1:3)
            !end do
            !write(*,*)

        enddo

        deallocate(adder, R11, R31, R12, W1)

        this%matrices_calculated = .true.

    end subroutine calculate_nonaxisymmetric_matrices

    subroutine calculate_tmatrix(this)
        class(ScatteringMode) :: this
        integer :: i, j

        if (.not. this%matrices_calculated) then
            call this%calculate_matrices()
        endif

        do i = 1, this%maxm
            call inverse_matrix(this%A31(i, :, :), this%matrix_size, this%A31inv(i, :, :))
            this%tmatr(i, :, :) = -matmul(this%A11(i, :, :), this%A31inv(i, :, :))
        enddo

        this%tmatrix_calculated = .true.
    end subroutine calculate_tmatrix


    ! Initializers
    subroutine calculate_axisymmetric_initial_te(this)
        class(AxisymmetricTE) :: this
        integer :: i

        this%initial = 0

        do i = 1, this%matrix_size
            this%initial(1, i) = -2q0 * (qcmplx(0q0, 1q0) ** i) * &
                    this%scatterer%layer(1, 0)%s1(i, 1)
        enddo

        this%initial_calculated = .true.
    end subroutine calculate_axisymmetric_initial_te

    subroutine calculate_nonaxisymmetric_initial_te(this)
        class(NonAxisymmetricTE) :: this
        integer :: i, m
        real(knd) :: k1

        !        write(*,*) 'maxm = ', this%maxm
        !        write(*,*) 'ms = ', this%matrix_size
        !        write(*,*) 'sin = ', qsin(this%scatterer%alpha)
        !        write(*,*) 'this%initial(1, 1) = ', this%initial(1, 1)
        !        write(*,*) 'r1 = ', this%scatterer%layer(1, 0)%r1(1)
        !        write(*,*) 's11 = ', this%scatterer%layer(1, 0)%s1(1, 1)
        this%initial = 0
        k1 = 2 * PI / this%scatterer%lambda
        do m = 1, this%maxm
            do i = 1, this%matrix_size / 2
                this%initial(m, i) = 4q0 * (qcmplx(0q0, 1q0) ** (i - 1)) / k1 * &
                        this%scatterer%layer(m, 0)%s1(i, 1) / qsin(this%scatterer%alpha)
            enddo
        end do

        this%initial_calculated = .true.
    end subroutine calculate_nonaxisymmetric_initial_te

    subroutine calculate_axisymmetric_initial_tm(this)
        class(AxisymmetricTM) :: this
        integer :: i

        this%initial = 0

        do i = 1, this%matrix_size
            this%initial(1, i) = 2q0 * cqsqrt(this%scatterer%eps(0) / this%scatterer%mu(0)) * (qcmplx(0q0, 1q0) ** i) * &
                    this%scatterer%layer(1, 0)%s1(i, 1)
        enddo

        this%initial_calculated = .true.
    end subroutine calculate_axisymmetric_initial_tm

    subroutine calculate_nonaxisymmetric_initial_tm(this)
        class(NonAxisymmetricTM) :: this
        integer :: i, m

        this%initial = 0

        do m = 1, this%maxm
            do i = 1, this%matrix_size / 2
                this%initial(m, i) = -4q0 * (qcmplx(0q0, 1q0) ** (i - 1)) * &
                        cqsqrt(this%scatterer%eps(0) / this%scatterer%mu(0)) * &
                        this%scatterer%layer(m, 0)%s1(i, 1) / qsin(this%scatterer%alpha) * this%scatterer%layer(m, 0)%r1(i)
            enddo
        end do

        this%initial_calculated = .true.
    end subroutine calculate_nonaxisymmetric_initial_tm

    ! Destructors
    subroutine delete_scattering_mode(this)

        class(ScatteringMode), intent(inout) :: this

        if (allocated(this%tmatr)) then
            deallocate(this%tmatr)
        endif
        if (allocated(this%solution)) then
            deallocate(this%solution)
        endif
        if (allocated(this%initial)) then
            deallocate(this%initial)
        endif
        if (allocated(this%A11)) then
            deallocate(this%A11)
        endif
        if (allocated(this%A31)) then
            deallocate(this%A31)
        endif
        if (allocated(this%A31inv)) then
            deallocate(this%A31inv)
        endif

    end subroutine delete_scattering_mode

    subroutine delete_axisymmetric_mode_te(this)

        type(AxisymmetricTE), intent(inout) :: this

        call this%delete_mode()
    end subroutine delete_axisymmetric_mode_te

    subroutine delete_nonaxisymmetric_mode_te(this)

        type(NonAxisymmetricTE), intent(inout) :: this

        call this%delete_mode()
    end subroutine delete_nonaxisymmetric_mode_te

    subroutine delete_axisymmetric_mode_tm(this)

        type(AxisymmetricTM), intent(inout) :: this

        call this%delete_mode()
    end subroutine delete_axisymmetric_mode_tm

    subroutine delete_nonaxisymmetric_mode_tm(this)

        type(NonAxisymmetricTM), intent(inout) :: this

        call this%delete_mode()
    end subroutine delete_nonaxisymmetric_mode_tm

    function get_symmetric_extinction_factor(this) result(ext)
        class(AxisymmetricMode) :: this
        integer :: i
        complex(knd) :: ideg
        real(knd) :: ext

        if (.not. this%solution_calculated) then
            call this%calculate_solution()
        endif

        ext = 0q0

        ideg = qcmplx(0q0, -1q0)

        do i = 1, this%matrix_size
            ext = ext + real(this%solution(1, i) * ideg**i * this%scatterer%layer(1, 0)%s1(i, 1), knd)
        enddo

        ext = -ext * 4q0 * this%scatterer%common_factor

    end function get_symmetric_extinction_factor

    function get_non_symmetric_extinction_factor(this) result(ext)
        class(NonAxisymmetricMode) :: this
        integer :: i, m
        complex(knd) :: ideg
        real(knd) :: ext, k1

        if (.not. this%solution_calculated) then
            call this%calculate_solution()
        endif

        ext = 0q0
        k1 = 2 * PI / this%scatterer%lambda

        ideg = qcmplx(0q0, -1q0)

        write(*, *) 'maxm = ', this%maxm
        write(*, *) 'ms = ', this%matrix_size

        do m = 1, this%maxm
            !write(*,*) 'm = ', m
            !write(*,*) 'solution = ', this%solution(m, 1:5)
            !write(*,*) 'ext = ', ext
            !write(*,*) 's1 = ', this%scatterer%layer(m, 0)%s1(:, 1)
            !write(*,*) 's1d = ', this%scatterer%layer(m, 0)%s1d(:, 1)
            !write(*,*)
            do i = 1, this%matrix_size / 2
                !write(*,*) i, this%solution(m, i), this%scatterer%layer(m, 0)%s1(i, 1), &
                !        this%solution(m, i + this%matrix_size / 2), this%scatterer%layer(m, 0)%s1d(i, 1), &
                !        this%scatterer%layer(m, 0)%r3(i)
                !write(*,*) 'dext = ', ideg**(i - 1) * (k1 * this%solution(m, i) * this%scatterer%layer(m, 0)%s1(i, 1) - &
                !        this%solution(m, i + this%matrix_size / 2) * ideg * this%scatterer%layer(m, 0)%s1d(i, 1))
                ext = ext + real(ideg**(i - 1) * (k1 * this%solution(m, i) * this%scatterer%layer(m, 0)%s1(i, 1) - &
                        this%solution(m, i + this%matrix_size / 2) * ideg * this%scatterer%layer(m, 0)%s1d(i, 1)), knd)
            enddo
        end do

        ext = ext * 4q0 * this%scatterer%common_factor * qsin(this%scatterer%alpha)

    end function get_non_symmetric_extinction_factor

    function get_symmetric_scattering_factor(this) result(sca)
        class(AxisymmetricMode) :: this
        integer :: i
        complex(knd) :: ideg
        real(knd) :: sca

        if (.not. this%solution_calculated) then
            call this%calculate_solution()
        endif

        sca = sum(abs(this%solution(1, :))**2) * 2q0 * this%scatterer%common_factor

    end function get_symmetric_scattering_factor

    function get_non_symmetric_scattering_factor(this) result(sca)
        class(NonAxisymmetricMode) :: this
        integer :: i, j, m
        complex(knd) :: ideg, Omega(this%matrix_size / 2, this%matrix_size / 2), &
                Kappa(this%matrix_size / 2, this%matrix_size / 2), Tau(this%matrix_size / 2, this%matrix_size / 2)
        real(knd) :: sca, k1

        if (.not. this%solution_calculated) then
            call this%calculate_solution()
        endif

        k1 = 2 * PI / this%scatterer%lambda

        sca = 0
        ideg = qcmplx(0q0, 1q0)

        do m = 1, this%maxm
            call calculate_omega(this%scatterer%layer(m, 0), this%scatterer%layer(m, 0), Omega, &
                    this%matrix_size / 2, number_of_d_coeff(this%scatterer%c0))
            call calculate_kappa(this%scatterer%layer(m, 0), this%scatterer%layer(m, 0), Kappa, &
                    this%matrix_size / 2, number_of_d_coeff(this%scatterer%c0))
            call calculate_tau(this%scatterer%layer(m, 0), this%scatterer%layer(m, 0), Tau, &
                    this%matrix_size / 2, number_of_d_coeff(this%scatterer%c0))
            do i = 1, this%matrix_size / 2
                do j = 1, this%matrix_size / 2
                    sca = sca + real(ideg**(j - i) * (&
                    k1**2 * this%solution(m, i) * conjg(this%solution(m, j)) * Omega(i, j) + &
                    ideg * k1 * (&
                    this%solution(m, i + this%matrix_size / 2) * conjg(this%solution(m, j)) * Kappa(i, j) - &
                    conjg(this%solution(m, j + this%matrix_size / 2)) * this%solution(m, i) * Kappa(j, i)) + &
                    this%solution(m, i + this%matrix_size / 2) * conjg(this%solution(m, j + this%matrix_size / 2)) * &
                    Tau(i, j)), knd)
                end do
            end do
        end do
        sca = sca * this%scatterer%common_factor

    end function get_non_symmetric_scattering_factor
end module scattering_mode
