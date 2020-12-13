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

        !		write(*,*) 'x = ', 2q0 * PI * rv / lambda
        !		write(*,*) 'c0 =', c0, 'c =', c

        !		write(*,*) 'thismatrsize = ', this%matrix_size
        !		write(*,*) 'matrsize = ', matrix_size

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
            !        if (allocated(this%Delta) .and. matrix_size /= this%matrix_size) then
            !			write(*,*) 'deallocating'
            deallocate(this%Delta, this%Kappa, this%Gamma01, this%Gamma11, this%Epsilon, this%Q01, this%Q11)
        endif
        !		write(*,*) 'ready to allocate scatterer'
        if (.not. allocated(this%Delta)) then
            allocate(this%Delta(this%maxm, matrix_size, matrix_size))
            allocate(this%Kappa(this%maxm, matrix_size, matrix_size))
            allocate(this%Gamma01(this%maxm, matrix_size, matrix_size))
            allocate(this%Gamma11(this%maxm, matrix_size, matrix_size))
            allocate(this%Epsilon(this%maxm, matrix_size, matrix_size))
            allocate(this%Q01(this%maxm, matrix_size, matrix_size))
            allocate(this%Q11(this%maxm, matrix_size, matrix_size))
        endif
        !		write(*,*) 'allocated scatterer'

        this%matrix_size = matrix_size

        this%common_factor = 1q0 / (cqabs(this%layer(1, 0)%c)**2q0 * &
                qsqrt((this%ksi(1)**2 - this%f) * (this%ksi(1)**2 - this%f * qcos(this%alpha)**2)))

        this%functions_calculated = .false.

    end subroutine set

    subroutine calculate_functions(this)
        class(Scatterer) :: this
        complex(knd), allocatable, dimension(:, :) :: tmp
        integer :: i, j
        complex(knd) :: identity(this%matrix_size, this%matrix_size)

        !		write(*,*) 'calculating functions'
        allocate(tmp(this%matrix_size, this%matrix_size))
        do j = 1, this%maxm
            do i = 0, this%number_of_layers
                call this%layer(j, i)%calculate()
            enddo
            call calculate_delta(this%layer(j, 0), this%layer(j, 1), this%Delta(j, :, :), this%matrix_size, this%accuracy)
            call calculate_kappa(this%layer(j, 1), this%layer(j, 1), this%Kappa(j, :, :), this%matrix_size, this%accuracy)
            call calculate_gamma(this%layer(j, 0), this%layer(j, 1), this%Gamma01(j, :, :), this%matrix_size, this%accuracy)
            call calculate_gamma(this%layer(j, 1), this%layer(j, 1), this%Gamma11(j, :, :), this%matrix_size, this%accuracy)
            call calculate_epsilon(this%layer(j, 1), this%layer(j, 1), this%Epsilon(j, :, :), this%matrix_size, this%accuracy)
            !call calculate_omega(this%layer(j, 1), this%layer(j, 1), tmp, this%matrix_size, this%accuracy)
            !this%Epsilon(j,:,:) = this%Epsilon(j,:,:) - tmp
            ! write(*,*) 'j = ', j
            ! write(*,*) 'Delta = ', this%Delta(j, 1:5, 1:5)
            ! write(*,*) 'Kappa = ', this%Kappa(j, 1:5, 1:5)
            ! write(*,*) 'Gamma11 = ', this%Gamma11(j, 1:5, 1:5)
            ! write(*,*) 'Epsilon = ', this%Epsilon(j, 1:5, 1:5)


            call get_identity_matrix(tmp, this%matrix_size)
            tmp = this%ksi(1)**2 * tmp - this%f * matmul(this%Gamma11(j, :, :), this%Gamma11(j, :, :))
            !write(*,*) 'tmp = ', tmp
            call inverse_matrix(tmp, this%matrix_size, this%Q01(j, :, :))

            call calculate_omega(this%layer(j, 1), this%layer(j, 1), tmp, this%matrix_size, this%accuracy)
            call get_identity_matrix(identity, this%matrix_size)
            tmp = (this%ksi(1)**2 - this%f) * identity + this%f * tmp
            call inverse_matrix(tmp, this%matrix_size, this%Q11(j, :, :))
            !write(*,*) 'must be equal'
            !write(*,*) this%Q11(j, :, :)
            !write(*,*) this%Q01(j, :, :)
            !write(*,*)



            !tmp = matmul(this%Delta(j, :, :), this%Gamma11(j,:,:))
            !write(*,*) 'must be equal'
            !write(*,*) tmp(:, :)
            !write(*,*) this%Gamma01(j, :, :)
            !write(*,*) 'must be identity', tmp(1:5,1:5)
            call calculate_omega(this%layer(j, 0), this%layer(j, 1), tmp, this%matrix_size, this%accuracy)
            tmp = (this%ksi(1)**2 - this%f) * this%Delta(j, :, :) + this%f * tmp
            call inverse_matrix(tmp, this%matrix_size, this%Q01(j, :, :))
            write(*,*) 'must be equal'
            write(*,*) this%Q01(j, :, :)
            write(*,*) matmul(this%Delta(j, :, :), this%Q11(j, :, :))
            this%Q01(j, :, :) = matmul(this%Delta(j, :, :), this%Q11(j, :, :))
        end do

        deallocate(tmp)
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