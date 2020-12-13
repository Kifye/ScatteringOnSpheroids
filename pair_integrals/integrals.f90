module integrals
    use regime
    use spheroidal

    implicit none

contains
    subroutine fill_common_multiplier(m, accuracy, common_multiplier)
        integer :: m, accuracy, r
        real(knd), allocatable, dimension(:) :: common_multiplier
        real(knd) :: factorial
        allocate(common_multiplier(0:accuracy))
        factorial = 1q0
        do r = 2, 2 * m
            factorial = factorial * r
        end do
        do r = 0, accuracy
            common_multiplier(r) = factorial * 2q0 / (2 * r + 2 * m + 1)
            factorial = factorial * (r + 2 * m + 1) / (r + 1)
        enddo
    end subroutine fill_common_multiplier

    ! sum from mod(n + 1, 2) because here n = n0 - m + 1 => mod(n0-m, 2) = mod(n0 - m + 1 + 1, 2) = mod(n + 1, 2)
    subroutine calculate_delta(first, second, Delta, matrix_size, accuracy)
        class(SpheroidalCalculation), intent(in) :: first, second
        real(knd), allocatable, dimension(:) :: common_multiplier
        complex(knd), intent(out) :: Delta(first%lnum, second%lnum)
        integer n, l, r, accuracy, matrix_size

        if (first%m /= second%m) then
            write(*, *) 'different m in delta'
        end if
        call fill_common_multiplier(first%m, accuracy, common_multiplier)
        Delta = 0q0

        !		write(*,*) 'n = ', first%lnum, 'l=', second%lnum
        !		write(*,*) 's1 = ', size(first%legendre(1,:)), 's2=', size(second%legendre(1,:))
        ! write(*,*) 'delta calculation'
        ! write(*,*) 'ms = ', matrix_size
        ! write(*,*) 'leg = ', size(first%legendre, 1), size(first%legendre, 2)
        ! write(*,*) first%legendre(2, 0:5)
        do n = 1, matrix_size
            do l = 1, matrix_size
                if (mod(abs(n - l), 2) == 1) then
                    cycle
                endif
                ! write(*,*) 'n = ', n, 'l = ', l
                ! write(*,*) 'leg1 = ', first%legendre(n, 0:5)
                ! write(*,*) 'leg2 = ', second%legendre(l, 0:5)
                ! write(*,*)
                do r = mod(n + 1, 2), accuracy, 2
                    Delta(n, l) = Delta(n, l) + first%legendre(n, r) * second%legendre(l, r) * &
                            common_multiplier(r);
                enddo
                ! write(*,*) 'Delta(n, l) =', Delta(n, l)
                ! write(*,*)
            enddo
        enddo

        !		write(*,*) 'delta = '
        !		write(*,*) Delta(1:4, 1:4)

        deallocate(common_multiplier)
    end subroutine calculate_delta

    real(knd) function gamma_c_lower(m, r)
        integer m, r
        gamma_c_lower = real(r, knd) / (2 * r + 2 * m - 1)
    end function gamma_c_lower
    real(knd) function gamma_c_upper(m, r)
        integer m, r
        gamma_c_upper = real(2 * m + 1 + r, knd) / (2 * r + 2 * m + 3)
    end function gamma_c_upper

    subroutine calculate_gamma(first, second, Gamma, matrix_size, accuracy)
        class(SpheroidalCalculation), intent(in) :: first, second
        real(knd), allocatable, dimension(:) :: common_multiplier
        complex(knd), intent(out) :: Gamma(matrix_size, matrix_size)
        integer n, l, r, accuracy, matrix_size, ix

        if (first%m /= second%m) then
            write(*, *) 'different m in gamma'
        end if
        call fill_common_multiplier(first%m, accuracy, common_multiplier)
        Gamma = 0q0

        !		write(*,*) 'n = ', first%lnum, 'l=', second%lnum
        !		write(*,*) 's1 = ', size(first%legendre(1,:)), 's2=', size(second%legendre(1,:))
        do n = 1, matrix_size
            do l = 1, matrix_size
                if (mod(abs(n - l), 2) == 0) then
                    cycle
                endif
                ix = mod(n + 1, 2)
                if (ix == 0) then
                    Gamma(n, l) = first%legendre(n, ix) * second%legendre(l, ix + 1) * gamma_c_upper(first%m, ix) * &
                            common_multiplier(ix)
                    ix = 2
                endif

                do r = ix, accuracy, 2
                    Gamma(n, l) = Gamma(n, l) + first%legendre(n, r) * common_multiplier(r) * &
                            (second%legendre(l, r - 1) * gamma_c_lower(first%m, r) + second%legendre(l, r + 1) * &
                                    gamma_c_upper(first%m, r))
                enddo
            enddo
        enddo

        !		write(*,*) 'delta = '
        !		write(*,*) Delta(1:4, 1:4)

        deallocate(common_multiplier)
    end subroutine calculate_gamma

    real(knd) function kappa_c_lower(m, r)
        integer m, r
        kappa_c_lower = -real(r * (r + m - 1), knd) / (2 * r + 2 * m - 1)
    end function kappa_c_lower
    real(knd) function kappa_c_upper(m, r)
        integer m, r
        kappa_c_upper = real((r + m + 2) * (2 * m + 1 + r), knd) / (2 * r + 2 * m + 3)
    end function kappa_c_upper

    subroutine calculate_kappa(first, second, Kappa, matrix_size, accuracy)
        class(SpheroidalCalculation), intent(in) :: first, second
        real(knd), allocatable, dimension(:) :: common_multiplier
        complex(knd), intent(out) :: Kappa(matrix_size, matrix_size)
        integer n, l, r, accuracy, matrix_size, ix

        if (first%m /= second%m) then
            write(*, *) 'different m in gamma'
        end if
        call fill_common_multiplier(first%m, accuracy, common_multiplier)
        Kappa = 0q0

        !		write(*,*) 'n = ', first%lnum, 'l=', second%lnum
        !		write(*,*) 's1 = ', size(first%legendre(1,:)), 's2=', size(second%legendre(1,:))
        do n = 1, matrix_size
            do l = 1, matrix_size
                if (mod(abs(n - l), 2) == 0) then
                    cycle
                endif
                ix = mod(n + 1, 2)
                if (ix == 0) then
                    Kappa(n, l) = first%legendre(n, ix) * second%legendre(l, ix + 1) * kappa_c_upper(first%m, ix) * &
                            common_multiplier(ix)
                    ix = 2
                endif

                do r = ix, accuracy, 2
                    Kappa(n, l) = Kappa(n, l) + first%legendre(n, r) * common_multiplier(r) * &
                            (second%legendre(l, r - 1) * kappa_c_lower(first%m, r) + second%legendre(l, r + 1) * &
                                    kappa_c_upper(first%m, r))
                enddo
            enddo
        enddo

        !		write(*,*) 'delta = '
        !		write(*,*) Delta(1:4, 1:4)

        deallocate(common_multiplier)
    end subroutine calculate_kappa
    ! Sigma
    real(knd) function sigma_c_lower(m, r)
        integer m, r
        sigma_c_lower = -real((r - 1) * r * (r + m - 1), knd) / &
                ((2 * r + 2 * m - 1) * (2 * r + 2 * m - 3))
    end function sigma_c_lower
    real(knd) function sigma_c_middle(m, r)
        integer m, r
        sigma_c_middle = real(3 * (r + m) * (r + m + 1) - m * m - 2, knd) / &
                ((2 * r + 2 * m - 1) * (2 * r + 2 * m + 3))
    end function sigma_c_middle
    real(knd) function sigma_c_upper(m, r)
        integer m, r
        sigma_c_upper = real((r + m + 2) * (r + 2 * m + 1) * (r + 2 * m + 2), knd) / &
                ((2 * r + 2 * m + 3) * (2 * r + 2 * m + 5))
    end function sigma_c_upper

    subroutine calculate_sigma(first, second, Sigma, matrix_size, accuracy)
        class(SpheroidalCalculation), intent(in) :: first, second
        real(knd), allocatable, dimension(:) :: common_multiplier
        complex(knd), intent(out) :: Sigma(matrix_size, matrix_size)
        integer n, l, r, accuracy, matrix_size, ix

        if (first%m /= second%m) then
            write(*, *) 'different m in epsilon'
        end if
        call fill_common_multiplier(first%m, accuracy, common_multiplier)
        Sigma = 0q0

        !		write(*,*) 'n = ', first%lnum, 'l=', second%lnum
        !		write(*,*) 's1 = ', size(first%legendre(1,:)), 's2=', size(second%legendre(1,:))
        do n = 1, matrix_size
            do l = 1, matrix_size
                if (mod(abs(n - l), 2) == 1) then
                    cycle
                endif
                ix = mod(n + 1, 2)

                Sigma(n, l) = first%legendre(n, ix) * common_multiplier(ix) * &
                        (second%legendre(l, ix) * sigma_c_middle(first%m, ix) + &
                                second%legendre(l, ix + 2) * sigma_c_upper(first%m, ix))

                do r = ix + 2, accuracy, 2
                    Sigma(n, l) = Sigma(n, l) + first%legendre(n, r) * common_multiplier(r) * &
                            (second%legendre(l, r - 2) * sigma_c_lower(first%m, r) + &
                                    second%legendre(l, r) * sigma_c_middle(first%m, r) + &
                                    second%legendre(l, r + 2) * sigma_c_upper(first%m, r))
                enddo
            enddo
        enddo

        !		write(*,*) 'delta = '
        !		write(*,*) Delta(1:4, 1:4)

        deallocate(common_multiplier)
    end subroutine calculate_sigma

    real(knd) function epsilon_c_lower(m, r)
        integer m, r
        epsilon_c_lower = -real((r - 1) * r * (r + m - 2), knd) / &
                ((2 * r + 2 * m - 1) * (2 * r + 2 * m - 3))
    end function epsilon_c_lower
    real(knd) function epsilon_c_middle(m, r)
        integer m, r
        epsilon_c_middle = real((r + m) * (r + m + 1) - 3 * m * m, knd) / &
                ((2 * r + 2 * m - 1) * (2 * r + 2 * m + 3))
    end function epsilon_c_middle
    real(knd) function epsilon_c_upper(m, r)
        integer m, r
        epsilon_c_upper = real((r + m + 3) * (r + 2 * m + 1) * (r + 2 * m + 2), knd) / &
                ((2 * r + 2 * m + 3) * (2 * r + 2 * m + 5))
    end function epsilon_c_upper

    subroutine calculate_epsilon(first, second, Epsilon, matrix_size, accuracy)
        class(SpheroidalCalculation), intent(in) :: first, second
        real(knd), allocatable, dimension(:) :: common_multiplier
        complex(knd), intent(out) :: Epsilon(matrix_size, matrix_size)
        integer n, l, r, accuracy, matrix_size, ix

        if (first%m /= second%m) then
            write(*, *) 'different m in epsilon'
        end if
        call fill_common_multiplier(first%m, accuracy, common_multiplier)
        Epsilon = 0q0

        !		write(*,*) 'n = ', first%lnum, 'l=', second%lnum
        !		write(*,*) 's1 = ', size(first%legendre(1,:)), 's2=', size(second%legendre(1,:))
        do n = 1, matrix_size
            do l = 1, matrix_size
                if (mod(abs(n - l), 2) == 1) then
                    cycle
                endif
                ix = mod(n + 1, 2)

                Epsilon(n, l) = first%legendre(n, ix) * common_multiplier(ix) * &
                        (second%legendre(l, ix) * epsilon_c_middle(first%m, ix) + &
                                second%legendre(l, ix + 2) * epsilon_c_upper(first%m, ix))

                do r = ix + 2, accuracy, 2
                    Epsilon(n, l) = Epsilon(n, l) + first%legendre(n, r) * common_multiplier(r) * &
                            (second%legendre(l, r - 2) * epsilon_c_lower(first%m, r) + &
                                    second%legendre(l, r) * epsilon_c_middle(first%m, r) + &
                                    second%legendre(l, r + 2) * epsilon_c_upper(first%m, r))
                enddo
            enddo
        enddo

        !		write(*,*) 'delta = '
        !		write(*,*) Delta(1:4, 1:4)

        deallocate(common_multiplier)
    end subroutine calculate_epsilon


    ! Omega
    real(knd) function omega_c_lower(m, r)
        integer m, r
        omega_c_lower = -real((r - 1) * r, knd) / &
                ((2 * r + 2 * m - 1) * (2 * r + 2 * m - 3))
    end function omega_c_lower
    real(knd) function omega_c_middle(m, r)
        integer m, r
        omega_c_middle = real(2 * ((r + m) * (r + m + 1) + m * m - 1), knd) / &
                ((2 * r + 2 * m - 1) * (2 * r + 2 * m + 3))
    end function omega_c_middle
    real(knd) function omega_c_upper(m, r)
        integer m, r
        omega_c_upper = -real((r + 2 * m + 1) * (r + 2 * m + 2), knd) / &
                ((2 * r + 2 * m + 3) * (2 * r + 2 * m + 5))
    end function omega_c_upper

    subroutine calculate_omega(first, second, Omega, matrix_size, accuracy)
        class(SpheroidalCalculation), intent(in) :: first, second
        real(knd), allocatable, dimension(:) :: common_multiplier
        complex(knd), intent(out) :: Omega(matrix_size, matrix_size)
        integer n, l, r, accuracy, matrix_size, ix

        if (first%m /= second%m) then
            write(*, *) 'different m in epsilon'
        end if
        call fill_common_multiplier(first%m, accuracy, common_multiplier)
        Omega = 0q0

        do n = 1, matrix_size
            do l = 1, matrix_size
                if (mod(abs(n - l), 2) == 1) then
                    cycle
                endif
                ix = mod(n + 1, 2)

                Omega(n, l) = first%legendre(n, ix) * common_multiplier(ix) * &
                        (second%legendre(l, ix) * omega_c_middle(first%m, ix) + &
                                second%legendre(l, ix + 2) * omega_c_upper(first%m, ix))

                do r = ix + 2, accuracy, 2
                    Omega(n, l) = Omega(n, l) + first%legendre(n, r) * common_multiplier(r) * &
                            (second%legendre(l, r - 2) * omega_c_lower(first%m, r) + &
                                    second%legendre(l, r) * omega_c_middle(first%m, r) + &
                                    second%legendre(l, r + 2) * omega_c_upper(first%m, r))
                enddo
            enddo
        enddo
        deallocate(common_multiplier)
    end subroutine calculate_omega


    real(knd) function tau_c_middle(m, r)
        integer m, r
        tau_c_middle = real((r + m) * (r + m + 1), knd)
    end function tau_c_middle

    subroutine calculate_tau(first, second, Tau, matrix_size, accuracy)
        class(SpheroidalCalculation), intent(in) :: first, second
        real(knd), allocatable, dimension(:) :: common_multiplier
        complex(knd), intent(out) :: Tau(first%lnum, second%lnum)
        integer n, l, r, accuracy, matrix_size

        if (first%m /= second%m) then
            write(*, *) 'different m in delta'
        end if
        call fill_common_multiplier(first%m, accuracy, common_multiplier)
        Tau = 0q0

        do n = 1, matrix_size
            do l = 1, matrix_size
                if (mod(abs(n - l), 2) == 1) then
                    cycle
                endif
                do r = mod(n + 1, 2), accuracy, 2
                    Tau(n, l) = Tau(n, l) + first%legendre(n, r) * second%legendre(l, r) * &
                            common_multiplier(r) * tau_c_middle(first%m, r);
                enddo
            enddo
        enddo

        deallocate(common_multiplier)
    end subroutine calculate_tau

end module integrals