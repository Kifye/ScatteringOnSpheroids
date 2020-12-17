module scattering

    use regime
    use constants
    use matrix
    use integrals
    use spheroidal
    use elemfunc

    implicit none
    private

    type, public :: Scatterer

        type(SpheroidalCalculation), allocatable, dimension(:, :) :: layer

        real(knd) :: alpha, common_factor, lambda
        complex(knd) :: c0

        ! dependent on layers
        real(knd), allocatable, dimension(:) :: ksi   !  1..number_of_layers
        complex(knd), allocatable, dimension(:) :: eps, mu  !  0..number_of_layers

        complex(knd), allocatable, dimension(:, :, :) :: Delta, Kappa, Gamma01, Gamma11, Epsilon, Q01, Q11
        integer :: matrix_size, f, accuracy, maxm, number_of_layers

        logical functions_calculated

    contains

        procedure, public :: set

        procedure, public :: calculate_functions

        final :: delete_scatterer
    end type

contains
    subroutine set(this, f, rv, lambda, ab, alpha, eps, mu, matrix_size, number_of_layers, maxm, accuracy)

        class(Scatterer) :: this
        real(knd) :: rv, lambda, ab, alpha, arg(1)
        complex(knd) :: ni, mu(0:NUMBER_OF_LAYERS), eps(0:NUMBER_OF_LAYERS)
        integer :: matrix_size, f, i, number_of_layers, maxm, j
        integer, optional :: accuracy

        this%f = f

        ! set arrays dependant on number of layers
        ! layers, ksi, eps and mu are always allocated and deallocated together
        if (allocated(this%ksi) .and. number_of_layers /= this%number_of_layers) then
            deallocate(this%ksi, this%eps, this%mu)
        end if
        if (.not. allocated(this%ksi)) then
            allocate(this%ksi(1:number_of_layers), this%eps(0:number_of_layers), this%mu(0:number_of_layers))
        end if
        this%number_of_layers = number_of_layers
        this%mu = mu
        this%eps = eps

        this%alpha = alpha
        this%lambda = lambda
        if (f == 1) then
            this%ksi(1) = ab / qsqrt(ab * ab - 1)
            this%c0 = (1q0 / ab)**(1q0 / 3q0)
        else
            this%ksi(1) = 1q0 / qsqrt(ab * ab - 1)
            this%c0 = (1q0 / ab)**(2q0 / 3q0)
        endif
        this%c0 = 2q0 * PI * rv / lambda * qsqrt(ab**2q0 - 1q0) * this%c0

        this%accuracy = number_of_d_coeff(this%c0)
        if (present(accuracy)) then
            this%accuracy = accuracy
        endif

        arg(1) = qcos(alpha)

        ! set layers[1..maxm][0..number_of_layers]
        if (allocated(this%layer)) then
            deallocate(this%layer)
        end if
        this%maxm = maxm
        allocate(this%layer(1:this%maxm, 0:NUMBER_OF_LAYERS))
        do i = 1, this%maxm
            call this%layer(i, 0)%set(i, i + matrix_size - 1, this%c0, this%ksi(1), 1, arg, f)
            do j = 1, number_of_layers
                call this%layer(i, j)%set(i, i + matrix_size - 1, this%c0 * cqsqrt(eps(j) * mu(j)), this%ksi(j), 1, arg, f)
            enddo
        end do

        if (allocated(this%Delta)) then
            deallocate(this%Delta, this%Kappa, this%Gamma01, this%Gamma11, this%Epsilon, this%Q01, this%Q11)
        endif
        if (.not. allocated(this%Delta)) then
            allocate(this%Delta(this%maxm, matrix_size, matrix_size))
            allocate(this%Kappa(this%maxm, matrix_size, matrix_size))
            allocate(this%Gamma01(this%maxm, matrix_size, matrix_size))
            allocate(this%Gamma11(this%maxm, matrix_size, matrix_size))
            allocate(this%Epsilon(this%maxm, matrix_size, matrix_size))
            allocate(this%Q01(this%maxm, matrix_size, matrix_size))
            allocate(this%Q11(this%maxm, matrix_size, matrix_size))
        endif

        this%matrix_size = matrix_size

        this%common_factor = 1q0 / (cqabs(this%layer(1, 0)%c)**2q0 * &
                qsqrt((this%ksi(1)**2 - this%f) * (this%ksi(1)**2 - this%f * qcos(this%alpha)**2)))

        this%functions_calculated = .false.

    end subroutine set

    subroutine calculate_functions(this)
        class(Scatterer) :: this
        type(SpheroidalCalculation) :: calculation_for_q0, calculation_for_q1

        complex(knd), allocatable, dimension(:, :) :: tmp, identity, result, side
        integer :: i, j, size_q, accuracy_q

        do j = 1, this%maxm
            do i = 0, this%number_of_layers
                call this%layer(j, i)%calculate()
            enddo
            call calculate_delta(this%layer(j, 0), this%layer(j, 1), this%Delta(j, :, :), this%matrix_size, this%accuracy)
            call calculate_kappa(this%layer(j, 1), this%layer(j, 1), this%Kappa(j, :, :), this%matrix_size, this%accuracy)
            call calculate_gamma(this%layer(j, 0), this%layer(j, 1), this%Gamma01(j, :, :), this%matrix_size, this%accuracy)
            call calculate_gamma(this%layer(j, 1), this%layer(j, 1), this%Gamma11(j, :, :), this%matrix_size, this%accuracy)
            call calculate_epsilon(this%layer(j, 1), this%layer(j, 1), this%Epsilon(j, :, :), this%matrix_size, this%accuracy)

            size_q = 2 * this%matrix_size
            call calculation_for_q0%set(j, j + size_q - 1, this%layer(j, 0)%c, this%layer(j, 0)%ksi + 1q0, &
                    this%layer(j, 0)%narg, this%layer(j, 0)%arg, this%layer(j, 0)%spheroidal_type)
            call calculation_for_q0%calculate()
            call calculation_for_q1%set(j, j + size_q - 1, this%layer(j, 1)%c, this%layer(j, 1)%ksi + 1q0, &
                    this%layer(j, 1)%narg, this%layer(j, 1)%arg, this%layer(j, 1)%spheroidal_type)
            call calculation_for_q1%calculate()

            accuracy_q = 150
            allocate(tmp(size_q, size_q), identity(size_q, size_q), result(size_q, size_q), side(size_q, size_q))

            call calculate_gamma(calculation_for_q0, calculation_for_q1, tmp, size_q, accuracy_q)
            this%Gamma01(j,:,:) = tmp(1:this%matrix_size, 1:this%matrix_size)
            call calculate_gamma(calculation_for_q1, calculation_for_q1, tmp, size_q, accuracy_q)
            call calculate_delta(calculation_for_q0, calculation_for_q1, identity, size_q, accuracy_q)
            result = matmul(identity, tmp)
           ! write(*,*) 'must be equal 0'
           ! write(*,*) this%Gamma01(j, :, :)
           ! write(*,*) result(1:this%matrix_size, 1:this%matrix_size)
           ! write(*,*)

            call get_identity_matrix(identity, size_q)
            call calculate_gamma(calculation_for_q1, calculation_for_q1, tmp, size_q, accuracy_q)
            !do i = 1, size_q
            !    write(*,*) 'tmp = ', tmp(i,:)
            !end do
            tmp = this%ksi(1)**2 * identity - this%f * matmul(tmp, tmp)
            call inverse_matrix(tmp, size_q, result)
            this%Q01(j, :, :) = result(1:this%matrix_size, 1:this%matrix_size)

            call calculate_omega(calculation_for_q1, calculation_for_q1, tmp, size_q, accuracy_q)
            call get_identity_matrix(identity, size_q)
            tmp = (this%ksi(1)**2 - this%f) * identity + this%f * tmp
            call inverse_matrix(tmp, size_q, result)
            this%Q11(j, :, :) = result(1:this%matrix_size, 1:this%matrix_size)
            side = result
            result = matmul(identity, side)

            this%Q01(j, :, :) = result(1:this%matrix_size, 1:this%matrix_size)

            deallocate(tmp, identity, result, side)
        end do

        this%functions_calculated = .true.

    end subroutine calculate_functions


    subroutine delete_scatterer(this)

        type(Scatterer), intent(inout) :: this

        if (allocated(this%Delta)) then
            deallocate(this%Delta)
            deallocate(this%Kappa)
            deallocate(this%Gamma01)
            deallocate(this%Gamma11)
            deallocate(this%Epsilon)
            deallocate(this%Q01)
            deallocate(this%Q11)
        endif

        if (allocated(this%layer)) then
            deallocate(this%layer)
        endif
        if (allocated(this%ksi)) then
            deallocate(this%ksi, this%eps, this%mu)
        end if
    end subroutine delete_scatterer
end module scattering