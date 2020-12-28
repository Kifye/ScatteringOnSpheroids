! Created by drakosha on 22.12.2020.

program test_ext_sca_full
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
    complex(knd) :: eps(0:1), mu(0:1), c0
    complex(knd), allocatable, dimension(:,:,:) :: tmatr
    complex(knd), allocatable, dimension(:,:) :: solution
    real(knd) :: lambda, rv, ab, alpha, ext_te, sca_te, dr, ext_tm, sca_tm, tmatr_abs_te, tmatr_abs_tm, &
            tmatr_rel_te, tmatr_rel_tm
    integer :: matrix_size, f, maxm, i
    type(ScatteringCalculation) :: scatter
    eps(0:1) = (/qcmplx(1q0, 0q0), qcmplx(1.7q0, 0.7q0)**2/)
    mu(0:1) = (/qcmplx(1q0, 0q0), qcmplx(1q0, 0q0)/)

    ! nm
    lambda = 1q0
    rv = 5q0
    ab = 1.1q0

    alpha = PI / 2q0

    !	write(*,*) 'start'

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
    write(*,*) 'ri = ', cqsqrt(eps(1) * mu(1))
    write(*,*) 'eps = ', eps
    write(*,*) 'mu = ', mu
    write(*,*)

    if (f == 1) then
        dr = lambda / 2q0 / PI / qsqrt(ab**2 - 1) * ab**(1q0 / 3q0)
    else
        dr = lambda / 2q0 / PI / qsqrt(ab**2 - 1) * ab**(2q0 / 3q0)
    endif
    maxm = 6
    write(*,*) 'maxm = ', maxm

    105	FORMAT('#', 3A8,' ', 2A6, ' ', 12A25)
    write(*,105) 'rv', 'x', 'c0', 'size', 'acc', &
            'te_tmatr_abs', 'te_tmatr_rel', 'te_ext', 'te_sca', 'te_abs', 'te_rel', &
            'tm_tmatr_abs', 'tm_tmatr_rel', 'tm_ext', 'tm_sca', 'tm_abs', 'tm_rel'

    dr = 3.0q0 / 2q0 / PI
    do rv = 1q0 * dr, 2q0 * dr, 10q0 * dr
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
        allocate(tmatr(maxm, 2 * matrix_size, 2 * matrix_size))
        allocate(solution(maxm, 2 * matrix_size))

        call scatter%set(rv, lambda, ab, alpha, matrix_size, eps, mu, f, 1, maxm)

        tmatr = scatter%get_tmatrix(2, maxm, 2 * matrix_size)
        !       write(*,*) 'TE:'
        !       do i = 1, 2 * matrix_size
        !           write(*,*) tmatr(1, i, :)
        !       end do
        !write(*,*) tmatr(1, 1:5, 1:5)
        !write(*,*) transpose(tmatr(1, 1:5, 1:5))
        call check_symmetric_result(tmatr(1,1:matrix_size,1:matrix_size), matrix_size, tmatr_abs_te, tmatr_rel_te)
        !write(*,*) 'up left symmetr', tmatr_abs_te, tmatr_rel_te
        call check_symmetric_result(tmatr(1,(matrix_size + 1):(2 * matrix_size),(matrix_size + 1):(2 * matrix_size)), &
                matrix_size, tmatr_abs_te, tmatr_rel_te)
        !write(*,*) 'low right symmetr', tmatr_abs_te, tmatr_rel_te
        !write(*,*) 'ur = ', tmatr(1,1:matrix_size,(matrix_size + 1):(2 * matrix_size))
        !write(*,*) 'ur = ', transpose(tmatr(1,1:matrix_size,(matrix_size + 1):(2 * matrix_size)))
        !write(*,*) 'll = ', tmatr(1,(matrix_size + 1):(2 * matrix_size),1:matrix_size)
        !write(*,*) 'll = ', transpose(tmatr(1,(matrix_size + 1):(2 * matrix_size),1:matrix_size))
        call check_symmetric_result(tmatr(1,:,:), 2 * matrix_size, tmatr_abs_te, tmatr_rel_te)
        solution = scatter%get_solution(2, maxm, 2 * matrix_size)
        ! write(*,*) solution(1, :)
        ext_te = scatter%get_extinction(0) + scatter%get_extinction(2)
        sca_te = scatter%get_scattering(0) + scatter%get_scattering(2)

        tmatr = scatter%get_tmatrix(3, maxm, 2 * matrix_size)
        !      write(*,*) 'TM:'
        !      do i = 1, 2 * matrix_size
        !          write(*,*) tmatr(1, i, :)
        !      end do
        call check_symmetric_result(tmatr(1,1:matrix_size,1:matrix_size), matrix_size, tmatr_abs_tm, tmatr_rel_tm)
        !write(*,*) 'up left symmetr', tmatr_abs_tm, tmatr_rel_tm
        call check_symmetric_result(tmatr(1,(matrix_size + 1):(2 * matrix_size),(matrix_size + 1):(2 * matrix_size)), &
                matrix_size, tmatr_abs_tm, tmatr_rel_tm)
        !write(*,*) 'low right symmetr', tmatr_abs_tm, tmatr_rel_tm
        !write(*,*) 'ur = ', tmatr(1,1:matrix_size,(matrix_size + 1):(2 * matrix_size))
        !write(*,*) 'ur = ', transpose(tmatr(1,1:matrix_size,(matrix_size + 1):(2 * matrix_size)))
        !write(*,*) 'll = ', tmatr(1,(matrix_size + 1):(2 * matrix_size),1:matrix_size)
        !write(*,*) 'll = ', transpose(tmatr(1,(matrix_size + 1):(2 * matrix_size),1:matrix_size))
        call check_symmetric_result(tmatr(1,:,:), 2 * matrix_size, tmatr_abs_tm, tmatr_rel_tm)

        ext_tm = scatter%get_extinction(1) + scatter%get_extinction(3)
        sca_tm = scatter%get_scattering(1) + scatter%get_scattering(3)

        100	FORMAT(' ', 3F8.4,' ', 2I6, ' ', 12E25.15)

        !write(*,*) 'rel = ', ext_te / sca_te

        ext_te = ext_te * qsqrt(ab**2 * qsin(alpha)**2 + qcos(alpha)**2) / ab ** (2q0 / 3q0)
        ext_tm = ext_tm * qsqrt(ab**2 * qsin(alpha)**2 + qcos(alpha)**2) / ab ** (2q0 / 3q0)
        sca_te = sca_te * qsqrt(ab**2 * qsin(alpha)**2 + qcos(alpha)**2) / ab ** (2q0 / 3q0)
        sca_tm = sca_tm * qsqrt(ab**2 * qsin(alpha)**2 + qcos(alpha)**2) / ab ** (2q0 / 3q0)

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
end program test_ext_sca_full