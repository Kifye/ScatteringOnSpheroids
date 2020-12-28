! Created by drakosha on 24.11.2020.

module scattering_calculation
    use scattering
    use scattering_mode
    use spheroidal
    use regime
    use constants

    implicit none
    private

    type, public :: ScatteringCalculation
        type(Scatterer) :: scatterer
        type(WavelengthPoint) :: calculation_point
        type(AxisymmetricTE) :: symte
        type(AxisymmetricTM) :: symtm
        type(NonAxisymmetricTE) :: nonsymte
        type(NonAxisymmetricTM) :: nonsymtm
        integer :: matrix_size, min_m, max_m
    contains
        procedure, public :: set, get_tmatrix, get_solution, get_extinction, get_scattering
        final :: delete_scattering
    end type ScatteringCalculation

contains

    subroutine set(this, rv, lambda, ab, alpha, matrix_size, eps, mu, f, number_of_layers, maxm, acc)

        class(ScatteringCalculation) :: this
        real(knd) :: rv, lambda, ab, alpha, arg(1), xv
        complex(knd) :: c0, c, mu(0:NUMBER_OF_LAYERS), eps(0:NUMBER_OF_LAYERS)
        integer :: matrix_size, f, i, number_of_layers, maxm
        integer, optional :: acc

        if (present(acc)) then
            call this%scatterer%set(f, rv, lambda, ab, alpha, eps, mu, matrix_size, number_of_layers, maxm, acc)
        else
            call this%scatterer%set(f, rv, lambda, ab, alpha, eps, mu, matrix_size, number_of_layers, maxm)
        end if

        call this%calculation_point%initialize(lambda, number_of_layers, eps, mu)

        this%matrix_size = matrix_size

        call this%symte%set(this%scatterer, matrix_size, 1)
        call this%symtm%set(this%scatterer, matrix_size, 1)
        call this%nonsymte%set(this%scatterer, 2 * matrix_size, maxm)
        call this%nonsymtm%set(this%scatterer, 2 * matrix_size, maxm)

    end subroutine set

    function get_tmatrix(this, mode, maxm, matrix_size) result(tmatr)
        class(ScatteringCalculation) :: this
        complex(knd) :: tmatr(maxm, matrix_size, matrix_size)
        integer :: i, mode, maxm, matrix_size

        tmatr = 0
        if (mode == 0) then
            tmatr = this%symte%get_tmatrix()
        elseif (mode == 1) then
            tmatr = this%symtm%get_tmatrix()
        elseif (mode == 2) then
            tmatr = this%nonsymte%get_tmatrix()
        elseif (mode == 3) then
            tmatr = this%nonsymtm%get_tmatrix()
        end if

    end function get_tmatrix

    function get_solution(this, mode, maxm, matrix_size) result(solution)
        class(ScatteringCalculation) :: this
        complex(knd) :: solution(maxm, matrix_size)
        integer :: i, mode, maxm, matrix_size

        solution = 0
        if (mode == 0) then
            solution = this%symte%get_solution()
        elseif (mode == 1) then
            solution = this%symtm%get_solution()
        elseif (mode == 2) then
            solution = this%nonsymte%get_solution()
        elseif (mode == 3) then
            solution = this%nonsymtm%get_solution()
        end if

    end function get_solution

    function get_extinction(this, mode) result(ext)
        class(ScatteringCalculation) :: this
        real(knd) :: ext
        integer :: i, mode

        ext = 0
        if (mode == 0) then
            ext = this%symte%get_extinction()
        elseif (mode == 1) then
            ext = this%symtm%get_extinction()
        elseif (mode == 2) then
            ext = this%nonsymte%get_extinction()
        elseif (mode == 3) then
            ext = this%nonsymtm%get_extinction()
        end if

    end function get_extinction

    function get_scattering(this, mode) result(sca)
        class(ScatteringCalculation) :: this
        real(knd) :: sca
        integer :: i, mode

        sca = 0
        if (mode == 0) then
            sca = this%symte%get_scattering()
        elseif (mode == 1) then
            sca = this%symtm%get_scattering()
        elseif (mode == 2) then
            sca = this%nonsymte%get_scattering()
        elseif (mode == 3) then
            sca = this%nonsymtm%get_scattering()
        end if

    end function get_scattering

    subroutine delete_scattering(this)

        type(ScatteringCalculation), intent(inout) :: this

    end subroutine delete_scattering

end module scattering_calculation