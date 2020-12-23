! Created by drakosha on 21.12.2020.

module svm_calculation
    use scattering
    use scattering_mode
    use spheroidal
    use regime
    use constants
    use svm


    implicit none
    private

    type, public :: SvmCalculation
        type(Scatterer) :: scatterer
        type(SvmTM) :: svm_tm
        integer :: matrix_size
    contains
        procedure, public :: set, get_tmatrix, get_solution, get_extinction, get_scattering
        final :: delete_scattering
    end type SvmCalculation

contains

    subroutine set(this, rv, lambda, ab, alpha, matrix_size, eps, mu, f, number_of_layers, maxm, acc)

        class(SvmCalculation) :: this
        real(knd) :: rv, lambda, ab, alpha, arg(1), xv
        complex(knd) :: c0, c, mu(0:NUMBER_OF_LAYERS), eps(0:NUMBER_OF_LAYERS)
        integer :: matrix_size, f, i, number_of_layers, maxm
        integer, optional :: acc

        if (present(acc)) then
            call this%scatterer%set(f, rv, lambda, ab, alpha, eps, mu, matrix_size, number_of_layers, maxm, acc)
        else
            call this%scatterer%set(f, rv, lambda, ab, alpha, eps, mu, matrix_size, number_of_layers, maxm)
        end if

        this%matrix_size = matrix_size

        call this%svm_tm%set(this%scatterer, 2 * matrix_size, maxm)

    end subroutine set

    function get_tmatrix(this, mode, maxm, matrix_size) result(tmatr)
        class(SvmCalculation) :: this
        complex(knd) :: tmatr(maxm, matrix_size, matrix_size)
        integer :: i, mode, maxm, matrix_size

        tmatr = this%svm_tm%get_tmatrix()

    end function get_tmatrix

    function get_solution(this, mode, maxm, matrix_size) result(solution)
        class(SvmCalculation) :: this
        complex(knd) :: solution(maxm, matrix_size)
        integer :: i, mode, maxm, matrix_size

        solution = this%svm_tm%get_solution()

    end function get_solution

    function get_extinction(this, mode) result(ext)
        class(SvmCalculation) :: this
        real(knd) :: ext
        integer :: i, mode

        ext = this%svm_tm%get_extinction()

    end function get_extinction

    function get_scattering(this, mode) result(sca)
        class(SvmCalculation) :: this
        real(knd) :: sca
        integer :: i, mode

        sca = this%svm_tm%get_scattering()

    end function get_scattering

    subroutine delete_scattering(this)

        type(SvmCalculation), intent(inout) :: this

    end subroutine delete_scattering

end module svm_calculation