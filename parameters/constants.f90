module constants
    use regime
    use elemfunc

    implicit none
    real(knd), parameter :: PI = 4q0 * atan(1q0)
!    integer, parameter :: NUMBER_OF_LAYERS = 1
!    integer, parameter :: NUMBER_OF_MODES = 3
    integer, parameter :: ADDITIONAL_N = 10

contains

    integer function size_of_matrices(c)
        complex(knd) :: c
        integer :: ic
        ic = qint(cqabs(c))
        if (ic <= 10) then
            size_of_matrices = 2 * ic + 35
        elseif (ic <= 30) then
            size_of_matrices = 2 * ic + 25
        elseif (ic <= 70) then
            size_of_matrices = 2 * ic + 5
        else
            size_of_matrices = 3 * ic / 2 + 10
        endif
    end function size_of_matrices

!    integer function getmaxm(c)
!        complex(knd) :: c
!        integer :: ic
!        getmaxm = 1
!    end function getmaxm

    integer function number_of_d_coeff(c)
        complex(knd) :: c
        integer :: ic
        ic = qint(cqabs(c))
        number_of_d_coeff = 50 + 2 * ic
    end function number_of_d_coeff

end module constants