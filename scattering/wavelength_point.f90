! Created by drakosha on 26.12.2020.

module wavelength_point

    use regime
    use constants

    implicit none
    private

    type, public :: WavelengthPoint
        real(knd) :: k
        complex(knd), allocatable, dimension(:) :: eps, mu

    contains
        procedure :: initialize => initialize_wavelength_point
        final :: delete_wavelength_point
    end type WavelengthPoint

contains
    subroutine initialize_wavelength_point(this, lambda, number_of_layers, eps, mu)
        class(WavelengthPoint), intent(inout) :: this
        real(knd), intent(in) :: lambda
        integer, intent(in) :: number_of_layers
        complex(knd), intent(in) :: eps(0:number_of_layers), mu(0:number_of_layers)

        this%k = 2q0 * PI / lambda

        if (allocated(this%eps) .and. size(this%eps) /= number_of_layers + 1) then
            deallocate(this%eps, this%mu)
        end if

        if (.not.allocated(this%eps)) then
            allocate(this%eps(0:number_of_layers), this%mu(0:number_of_layers))
        end if

        this%eps = eps
        this%mu = mu
    end subroutine initialize_wavelength_point

    subroutine delete_wavelength_point(this)
        type(WavelengthPoint), intent(inout) :: this

        if (allocated(this%eps)) then
            deallocate(this%eps, this%mu)
        end if
    end subroutine delete_wavelength_point
end module wavelength_point