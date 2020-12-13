! Created by drakosha on 03.12.2020.

program test_symmetric
    use regime
    use spheroidal
    use matrix
    use integrals
    use scattering_calculation
    use constants
    use elemfunc

    implicit none

    integer, parameter :: MAXSIZE = 100
    !    type(Scatterer) :: scatter, scatter2
    type(SpheroidalCalculation) :: first, second
    complex(knd) :: ni, eps(0:1), mu(0:1), c0
    complex(knd), allocatable, dimension(:,:,:) :: tmatr
    complex(knd), allocatable, dimension(:,:) :: solution
    real(knd) :: lambda, rv, ab, alpha, ext_te, sca_te, dr, ext_tm, sca_tm, tmatr_abs_te, tmatr_abs_tm, &
            tmatr_rel_te, tmatr_rel_tm
    integer :: matrix_size, f, maxm, i
    type(ScatteringCalculation) :: scatter
    eps(0:1) = (/qcmplx(1q0, 0q0), qcmplx(2q0, 0q0)/)
    mu(0:1) = (/qcmplx(1q0, 0q0), qcmplx(1q0, 0q0)/)

    ! nm
    lambda = 1q0
    rv = 5q0
    ab = 2q0

    alpha = 0.5q0

    !	write(*,*) 'start'

    ni = qcmplx(qsqrt(2q0), 0.0q0)

    matrix_size = MAXSIZE
    f = 1

    if (f == 1) then
        write(*,*) 'prolate spheroids'
    else
        write(*,*) 'oblate spheroids'
    endif
    write(*,*) 'alpha = ', alpha, 'radian'
    write(*,*) 'lambda = ', lambda
    write(*,*) 'a/b = ', ab
    write(*,*) 'ri = ', ni
    write(*,*) 'eps = ', eps
    write(*,*) 'mu = ', mu
    write(*,*)

    if (f == 1) then
        dr = lambda / 2q0 / PI / qsqrt(ab**2 - 1) * ab**(1q0 / 3q0)
    else
        dr = lambda / 2q0 / PI / qsqrt(ab**2 - 1) * ab**(2q0 / 3q0)
    endif

    105	FORMAT('#', 3A8,' ', 2A6, ' ', 12A25)
    write(*,105) 'rv', 'x', 'c0', 'size', 'acc', &
            'te_tmatr_abs', 'te_tmatr_rel', 'te_ext', 'te_sca', 'te_abs', 'te_rel', &
            'tm_tmatr_abs', 'tm_tmatr_rel', 'tm_ext', 'tm_sca', 'tm_abs', 'tm_rel'

    maxm = 1
    do rv = dr, 100q0 * dr, 10q0 * dr
        if (f == 1) then
            c0 = (1q0 / ab)**(1q0 / 3q0)
        else
            c0 = (1q0 / ab)**(2q0 / 3q0)
        endif
        c0 = 2q0 * PI * rv / lambda * qsqrt(ab**2q0 - 1q0) * c0
        matrix_size = size_of_matrices(c0)
        if (allocated(tmatr)) then
            deallocate(tmatr)
        endif
        allocate(tmatr(maxm, matrix_size, matrix_size))
        allocate(solution(maxm, matrix_size))

        call scatter%set(rv, lambda, ab, alpha, matrix_size, eps, mu, f, 1, maxm)

        tmatr = scatter%get_tmatrix(0, maxm, matrix_size)
        !do i = 1, matrix_size
        !    write(*,*) tmatr(1, i, :)
        !end do
        !write(*,*) tmatr(1, 1:5, 1:5)
        call check_symmetric_result(tmatr(1,:,:), matrix_size, tmatr_abs_te, tmatr_rel_te)
        solution = scatter%get_solution(0, maxm, matrix_size)
        ! write(*,*) solution(1, :)
        ext_te = scatter%get_extinction(0)
        sca_te = scatter%get_scattering(0)

        tmatr = scatter%get_tmatrix(1, maxm, matrix_size)
        call check_symmetric_result(tmatr(1,:,:), matrix_size, tmatr_abs_tm, tmatr_rel_tm)

        ext_tm = scatter%get_extinction(1)
        sca_tm = scatter%get_scattering(1)

        100	FORMAT(' ', 3F8.4,' ', 2I6, ' ', 12E25.15)

        write(*,100) rv, 2q0 * PI * rv / lambda, cqabs(c0), matrix_size, number_of_d_coeff(c0), &
                tmatr_abs_te, tmatr_rel_te, &
                ext_te, sca_te, qabs(qabs(ext_te) - qabs(sca_te)), &
                qabs(qabs(ext_te) - qabs(sca_te)) / (qabs(ext_te) + qabs(sca_te)), &
                tmatr_abs_tm, tmatr_rel_tm, ext_tm, sca_tm, qabs(qabs(ext_tm) - qabs(sca_tm)), &
                qabs(qabs(ext_tm) - qabs(sca_tm)) / (qabs(ext_tm) + qabs(sca_tm))

        deallocate(tmatr)
        deallocate(solution)
        !exit
    enddo

end program test_symmetric