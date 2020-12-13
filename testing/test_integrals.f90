! Created by drakosha on 03.12.2020.

program test_integrals
    use integrals
    use constants
    use spheroidal
    use matrix

    implicit none

    integer, parameter :: MAXN = 100
    integer, parameter :: DSIZE = 4
    integer, parameter :: ACC = 2000
    complex(16) c1, c2, res1, res2, d1(DSIZE, DSIZE), d2(DSIZE, DSIZE), d3(DSIZE, DSIZE), d4(DSIZE, DSIZE)
    real(16) cur, ksi, arg(1)
    integer n, m, l, i, j, numarg, cmod, carg, accuracy
    type(SpheroidalCalculation) :: first, second

    90	FORMAT(' ', 4I5, ' ', 5E40.25)
    100	FORMAT(' ', 1I5, ' ', 3E40.25)
    105	FORMAT('#', 1A5, ' ', 2A40)
    109	FORMAT('#', 4A5, ' ', 5A40)
    110	FORMAT(' ', 2I5, ' ', 4E40.25)
    open(1, file = 'output_kappa_gamma_10_10_20.txt')
    open(2, file = 'output_delta_modified.txt')

    c1 = qcmplx(5q0, 0.0q0)
    c2 = qcmplx(2.0q0, 0.0q0)
    ksi = 1.3d0
    numarg = 1
    arg(1) = qcos(0.5q0)
    accuracy = number_of_d_coeff(c1)
    do m = 1, 3
        n = m + DSIZE - 1
        write(*,*) 'm = ', m
        call first%set(m, n, c1, ksi, numarg, arg, 1)
        call first%calculate()
        call second%set(m, n, c2, ksi, numarg, arg, 1)
        call second%calculate()

        !do i = 1, DSIZE
        !    write(*,*) 'i = ', i, 'd1 = ', first%legendre(i, 0:5)
        !end do

        !call calculate_delta(first, first, d1, DSIZE, accuracy)
        !write(*,*) 'delta is identity'
        !write(*,*) d1

        !call calculate_delta(first, second, d1, DSIZE, number_of_d_coeff(c1))
        !call calculate_delta(second, first, d2, DSIZE, number_of_d_coeff(c1))
        !write(*,*) 'delta is symmetric'
        !write(*,*) d1
        !write(*,*) transpose(d2)

        !call calculate_gamma(first, second, d1, DSIZE, number_of_d_coeff(c1))
        !call calculate_gamma(second, first, d2, DSIZE, number_of_d_coeff(c1))
        !write(*,*) 'gamma is symmetric'
        !write(*,*) d1
        !write(*,*) transpose(d2)

        call calculate_gamma(first, second, d1, DSIZE, number_of_d_coeff(c1))
        call calculate_kappa(first, second, d2, DSIZE, number_of_d_coeff(c1))
        call calculate_kappa(second, first, d3, DSIZE, number_of_d_coeff(c1))
        write(*,*) '2gamma = kappa + kappa'
        write(*,*) 2 * d1
        write(*,*) d2 + transpose(d3)

        !call calculate_omega(first, second, d1, DSIZE, accuracy)
        !call calculate_omega(second, first, d2, DSIZE, accuracy)
        !write(*,*) 'omega is symmetric'
        !write(*,*) d1
        !write(*,*) transpose(d2)

        !call calculate_delta(first, second, d1, DSIZE, accuracy)
        !call calculate_omega(first, second, d4, DSIZE, accuracy)
        !call calculate_sigma(first, second, d2, DSIZE, accuracy)
        !call calculate_sigma(second, first, d3, DSIZE, accuracy)
        !write(*,*) 'sigma + sigma = 2delta - omega'
        !write(*,*) 2 * d1 - d4
        !write(*,*) d2 + transpose(d3)

        call calculate_delta(first, second, d1, DSIZE, accuracy)
        call calculate_omega(first, second, d4, DSIZE, accuracy)
        call calculate_epsilon(first, second, d2, DSIZE, accuracy)
        call calculate_epsilon(second, first, d3, DSIZE, accuracy)
        write(*,*) 'epsilon + epsilon = 2delta - 3omega'
        write(*,*) 2 * d1 - 3 * d4
        write(*,*) d2 + transpose(d3)

        call calculate_delta(first, second, d1, DSIZE, accuracy)
        call calculate_omega(first, first, d4, DSIZE, accuracy)
        call calculate_epsilon(first, first, d2, DSIZE, accuracy)
        call calculate_sigma(first, first, d3, DSIZE, accuracy)
        write(*,*) 'epsilon = sigma - omega'
        write(*,*) d3 - d4
        write(*,*) d2

        call calculate_delta(first, second, d1, DSIZE, accuracy)
        call calculate_omega(first, first, d4, DSIZE, accuracy)
        call calculate_gamma(first, first, d2, DSIZE, accuracy)
        call get_identity_matrix(d3, DSIZE)
        write(*,*) 'omega = identity - gamma^2'
        write(*,*) d4
        write(*,*) d3 - matmul(d2, d2)


        !call calculate_sigma(first, first, d1, DSIZE, accuracy)
        !call calculate_kappa(first, first, d2, DSIZE, accuracy)
        !call calculate_gamma(first, first, d3, DSIZE, accuracy)
        !write(*,*) 'sigma = kappa * gamma'
        !write(*,*) d1
        !write(*,*) matmul(d2, d3)

        call calculate_delta(first, second, d1, DSIZE, accuracy)
        call calculate_gamma(second, second, d2, DSIZE, accuracy)
        call calculate_gamma(first, second, d3, DSIZE, accuracy)
        write(*,*) 'gamma = delta * gamma'
        write(*,*) matmul(d1, d2)
        write(*,*) d3
    enddo
end program test_integrals