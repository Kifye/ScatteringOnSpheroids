! Created by drakosha on 21.12.2020.

program test_svm
    use regime
    use spheroidal
    use matrix
    use integrals
    use scattering_calculation
    use constants
    use elemfunc
    use svm
    use svm_calculation

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
    type(SvmCalculation) :: svmc
    eps(0:1) = (/qcmplx(1q0, 0q0), qcmplx(2q0, 0q0)/)
    mu(0:1) = (/qcmplx(1q0, 0q0), qcmplx(1q0, 0q0)/)

    ! nm
    lambda = 1q0
    rv = 5q0
    ab = 2q0

    alpha = 0.5q0

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

    105	FORMAT('#', 3A8,' ', 2A6, ' ', 12A25)
    write(*,105) 'rv', 'x', 'c0', 'size', 'acc', &
            'svm_tmatr_abs', 'svm_tmatr_rel', 'svm_ext', 'svm_sca', 'svm_abs', 'svm_rel', &
            'tm_tmatr_abs', 'tm_tmatr_rel', 'tm_ext', 'tm_sca', 'tm_abs', 'tm_rel'

    maxm = 1
    do rv = 1 * dr, 50q0 * dr, 10q0 * dr
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
        call svmc%set(rv, lambda, ab, alpha, matrix_size, eps, mu, f, 1, maxm)

        tmatr = svmc%get_tmatrix(3, maxm, 2 * matrix_size)
        !      write(*,*) 'svn TM:'
        !       do i = 1, 2 * matrix_size
        !           write(*,*) tmatr(1, i, :)
        !       end do
        write(*,*) 'tmatr svn  = ', transpose(tmatr(1, 1:10, 1:10))
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
        solution = svmc%get_solution(2, maxm, 2 * matrix_size)
        ! write(*,*) solution(1, :)
        ext_te = svmc%get_extinction(2)
        sca_te = svmc%get_scattering(2)

        tmatr = scatter%get_tmatrix(3, maxm, 2 * matrix_size)
        write(*,*) 'tmatr ebcm = ', tmatr(1, 1:10, 1:10)
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

        ext_tm = scatter%get_extinction(3)
        sca_tm = scatter%get_scattering(3)

        100	FORMAT(' ', 3F8.4,' ', 2I6, ' ', 12E25.15)

        !write(*,*) 'rel = ', ext_te / sca_te

        write(*,100) rv, 2q0 * PI * rv / lambda, cqabs(c0), matrix_size, number_of_d_coeff(c0), &
                tmatr_abs_te, tmatr_rel_te, &
                ext_te, sca_te, qabs(qabs(ext_te) - qabs(sca_te)), &
                qabs(qabs(ext_te) - qabs(sca_te)) / (qabs(ext_te) + qabs(sca_te)), &
                tmatr_abs_tm, tmatr_rel_tm, ext_tm, sca_tm, qabs(qabs(ext_tm) - qabs(sca_tm)), &
                qabs(qabs(ext_tm) - qabs(sca_tm)) / (qabs(ext_tm) + qabs(sca_tm))

        deallocate(tmatr)
        deallocate(solution)
        exit
    enddo

end program test_svm