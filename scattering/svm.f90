! Created by drakosha on 20.12.2020.

module svm
    use scattering_mode
    use regime
    use scattering
    use matrix
    use integrals
    use constants


    implicit none

    type, public, extends(NonAxisymmetricTM) :: SvmTM
    contains
        procedure, private :: calculate_tmatrix => calculate_svm_tmatrix_tm
        procedure, private :: calculate_matrices => calculate_svm_matrices_tm
    end type SvmTM
contains

    subroutine calculate_svm_tmatrix_tm(this)
        class(SvmTM) :: this

        integer :: i, j, n
        complex(knd) :: mu(0:1), eps(0:1), initial_corrector(this%matrix_size), solution_corrector(this%matrix_size)
        complex(knd) :: k1, c1

        if (.not. this%matrices_calculated) then
            call this%calculate_matrices()
        endif

        n = this%matrix_size / 2
        write(*,*) 'n = ', n
        k1 = 2q0 * PI / this%scatterer%lambda

        do i = 1, this%maxm
            c1 = this%scatterer%layer(i, 0)%c
            !write(*,*) 'A31 = ', this%A31(i,1:10,1:10)
            !do j = 1, this%matrix_size
            !    write(*,*) this%A31(i,j,:)
            !end do

            !write(*,*) 'A11 = ', this%A11(i,1:10,1:10)
            !do j = 1, this%matrix_size
            !    write(*,*) this%A11(i,j,:)
            !end do

            call inverse_matrix(this%A31(i, :, :), this%matrix_size, this%A31inv(i, :, :))
            !write(*,*) 'A31inv = '
            !do j = 1, this%matrix_size
            !    write(*,*) this%A31inv(i,j,:)
            !end do
            this%tmatr(i, :, :) = -matmul(this%A31inv(i, :, :), this%A11(i, :, :))
            write(*,*) 'svm  = ', this%tmatr(i, 1:10, 1:10)
            !write(*,*) 'this%tmatr = '
            !do j = 1, this%matrix_size
            !    write(*,*) this%tmatr(i,j,:)
            !end do

            initial_corrector(1:n) = (this%scatterer%layer(i, 0)%r1 * k1)
            initial_corrector((n + 1):(2 * n)) = (this%scatterer%layer(i, 0)%r1 * c1)
            solution_corrector(1:n) = 1q0 / (this%scatterer%layer(i, 0)%r3 * k1)
            solution_corrector((n + 1):(2 * n)) = 1q0 / (this%scatterer%layer(i, 0)%r3 * c1)
            !write(*,*) 'n = ', n, 'ms = ', this%matrix_size, 'size(init) = ', size(initial_corrector), 'size(matr) = ',&
            !        size(this%tmatr(i,1,:)), 'size(calc) = ', size(this%scatterer%layer(i, 0)%r3)

            !write(*,*) 'sizes:', size(initial_corrector), size(this%scatterer%layer(m, 0)%r1), size(this%A11(m,1,:))
            call multiply_by_diag_right(this%tmatr(i, :, :), this%matrix_size, initial_corrector)
            call multiply_by_diag_left(this%tmatr(i, :, :), this%matrix_size, solution_corrector)

        enddo

        this%tmatrix_calculated = .true.

    end subroutine calculate_svm_tmatrix_tm

    function get_part_11(f, ksi, eps, mu, R2, R1, matrix_size) result(result)
        complex(knd) :: eps(0:1), mu(0:1)
        real(knd) :: ksi
        integer :: matrix_size, f, i
        complex(knd) :: R2(matrix_size, matrix_size), R1(matrix_size, matrix_size), tmp(matrix_size, matrix_size), &
                result(matrix_size, matrix_size)

        !write(*,*) 'ksi = ', ksi
        !write(*,*) 'R2 = '
        !do i = 1, matrix_size
        !    write(*,*) R2(i,:)
        !end do
        !write(*,*) 'R1 = '
        !do i = 1, matrix_size
        !    write(*,*) R1(i,:)
        !end do

        result = ksi * (R2 - R1)
        !write(*,*) 'result = '
        !do i = 1, matrix_size
        !    write(*,*) result(i,:)
        !end do


    end function get_part_11

    function get_part_12(f, ksi, eps, mu, R2, R1, Gamma, matrix_size) result(result)

        complex(knd) :: eps(0:1), mu(0:1), identity(matrix_size, matrix_size)
        real(knd) :: ksi
        integer :: matrix_size, f, i
        complex(knd) :: R2(matrix_size, matrix_size), R1(matrix_size, matrix_size), tmp(matrix_size, matrix_size), &
                Gamma(matrix_size, matrix_size), result(matrix_size, matrix_size)

        result = f * (matmul(Gamma, R2 - R1))
        !write(*,*) 'result = '
        !do i = 1, matrix_size
        !    write(*,*) result(i,:)
        !end do

    end function get_part_12

    function get_part_21(f, ksi, eps, mu, R2, R1, Gamma, Kappa, matrix_size) result(result)

        complex(knd) :: eps(0:1), mu(0:1)
        real(knd) :: ksi
        integer :: matrix_size, f, i
        complex(knd) :: R2(matrix_size, matrix_size), R1(matrix_size, matrix_size), tmp(matrix_size, matrix_size), &
                Gamma(matrix_size, matrix_size), Kappa(matrix_size, matrix_size), &
                result(matrix_size, matrix_size)

        !write(*,*) 'eps = ', eps(1)
        result = matmul(Gamma, 1q0 / eps(1) * R2 - R1) - (1q0 - 1q0 / eps(1)) * ksi / (ksi**2 - f) * Kappa
        !write(*,*) 'result = '
        !do i = 1, matrix_size
        !    write(*,*) result(i,:)
        !end do

    end function get_part_21

    function get_part_22(f, ksi, eps, mu, R2, R1, Gamma, Kappa, matrix_size) result(result)
        complex(knd) :: eps(0:1), mu(0:1)
        real(knd) :: ksi
        integer :: matrix_size, f, i
        complex(knd) :: R2(matrix_size, matrix_size), R1(matrix_size, matrix_size), tmp(matrix_size, matrix_size), &
                Gamma(matrix_size, matrix_size), Kappa(matrix_size, matrix_size), &
                result(matrix_size, matrix_size), identity(matrix_size, matrix_size)

        call get_identity_matrix(identity, matrix_size)

        result = 1q0 / eps(1) * (identity + ksi * R2) - (identity + ksi * R1) - &
                f * (1q0 - 1q0 / eps(1)) * ksi / (ksi**2 - f) * matmul(Kappa, Gamma)
        !write(*,*) 'result = '
        !do i = 1, matrix_size
        !    write(*,*) result(i,:)
        !end do

    end function get_part_22

    subroutine set_full_matrix(scat, R1diag, R2diag, m, result)

        class(Scatterer) :: scat
        real(knd) :: ksi
        integer :: i, j, m
        complex(knd) :: diag(scat%matrix_size), R1diag(scat%matrix_size), R2diag(scat%matrix_size), &
                R1(scat%matrix_size, scat%matrix_size), R2(scat%matrix_size, scat%matrix_size), &
                Kappa(scat%matrix_size, scat%matrix_size), Gamma(scat%matrix_size, scat%matrix_size), &
                result(2 * scat%matrix_size, 2 * scat%matrix_size), Delta(scat%matrix_size, scat%matrix_size), &
                tmp(scat%matrix_size, scat%matrix_size)

        call get_full_matrix_from_diag(R1diag, scat%matrix_size, R1)

        write(*,*) 'R1 = ', R1(1:10, 1:10)

        call get_full_matrix_from_diag(R2diag, scat%matrix_size, R2)
        write(*,*) 'R2 = ', R2(1:10, 1:10)
        call calculate_delta(scat%layer(m, 1), scat%layer(m, 0), Delta, scat%matrix_size, scat%accuracy)

        !write(*,*) 'Delta1 = ', Delta(1:10, 1:10)
        call inverse_matrix(Delta, scat%matrix_size, tmp)
        !R2 = matmul(tmp, matmul(R2, Delta))

        !write(*,*) 'invers = ', tmp(1:10, 1:10)
        !tmp = matmul(tmp, Delta)
        !write(*,*) 'check  = ', tmp(1:10, 1:10)
        R2 = matmul(R2, Delta)
        call calculate_delta(scat%layer(m, 0), scat%layer(m, 1), Delta, scat%matrix_size, scat%accuracy)
        !write(*,*) 'Delta2 = ', Delta(1:10, 1:10)
        !write(*,*)
        R2 = matmul(Delta, R2)

        write(*,*) 'acc = ', scat%accuracy

        call calculate_kappa(scat%layer(m, 0), scat%layer(m, 0), Kappa, scat%matrix_size, scat%accuracy)
        call calculate_gamma(scat%layer(m, 0), scat%layer(m, 0), Gamma, scat%matrix_size, scat%accuracy)


        !write(*,*) 'Delta = '
        !do i = 1, scat%matrix_size
        !    write(*,*) Delta(i,:)
        !end do
        !write(*,*) 'Gamma = '
        !do i = 1, scat%matrix_size
        !    write(*,*) Gamma(i,:)
        !end do
        !write(*,*) 'Kappa = '
        !do i = 1, scat%matrix_size
        !    write(*,*) Kappa(i,:)
        !end do

        !write(*,*) 'Epsilon = ', Epsilon
        !write(*,*) 'Q01 = ', Q01
        !write(*,*) 'Q11 = ', Q11
        !write(*,*)
        !write(*,*) 'ms = ', scat%matrix_size
        result(1:scat%matrix_size, 1:scat%matrix_size) = &
                get_part_11(scat%f, scat%ksi(1), scat%eps, scat%mu, R2, R1, scat%matrix_size)
        result(1:scat%matrix_size, (scat%matrix_size + 1):(2 * scat%matrix_size)) = &
                get_part_12(scat%f, scat%ksi(1), scat%eps, scat%mu, R2, R1, Gamma, scat%matrix_size)
        result((scat%matrix_size + 1):(2 * scat%matrix_size), 1:scat%matrix_size) = &
                get_part_21(scat%f, scat%ksi(1), scat%eps, scat%mu, R2, R1, Gamma, Kappa, scat%matrix_size)
        result((scat%matrix_size + 1):(2 * scat%matrix_size), (scat%matrix_size + 1):(2 * scat%matrix_size)) = &
                get_part_22(scat%f, scat%ksi(1), scat%eps, scat%mu, R2, R1, Gamma, Kappa, scat%matrix_size)

    end subroutine set_full_matrix


    subroutine calculate_svm_matrices_tm(this)
        class(SvmTM) :: this
        integer :: m
        complex(knd), allocatable, dimension(:) :: R0, R1, R2

        if (.not.(this%scatterer%functions_calculated)) then
            call this%scatterer%calculate_functions()
        endif

        allocate(R0(this%scatterer%matrix_size), R1(this%scatterer%matrix_size), R2(this%scatterer%matrix_size))

        this%A11 = 0
        this%A31 = 0
        this%A31inv = 0

        do m = 1, this%maxm

            R0 = this%scatterer%layer(m, 0)%r1d / this%scatterer%layer(m, 0)%r1
            R1 = this%scatterer%layer(m, 0)%r3d / this%scatterer%layer(m, 0)%r3
            R2 = this%scatterer%layer(m, 1)%r1d / this%scatterer%layer(m, 1)%r1
            !write(*, *) 'R0 = ', R0
            !write(*, *) 'R1 = ', R1

            call set_full_matrix(this%scatterer, R0, R2, m, this%A11(m, :, :))
            call set_full_matrix(this%scatterer, R1, R2, m, this%A31(m, :, :))

        enddo

        deallocate(R0, R1, R2)

        this%matrices_calculated = .true.

    end subroutine calculate_svm_matrices_tm
end module svm