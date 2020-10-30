program test_ext_sca_real
use regime
use spheroidal
use matrix
use integrals
use scattering

implicit none

	type(Scatterer) :: scatter, scatter2
	type(SpheroidalCalculation) :: first, second
	complex(knd) :: ni, tmatr(10, 10), res(10)
	real(knd) :: lambda, rv, ab, alpha, eps(0:1), mu(0:1), ext, sca, arg(1)
	integer :: matrix_size, f
	! nm
	lambda = 500q0
	rv = 250q0
	ab = 3q0
	
	alpha = 0.5q0
	
	write(*,*) 'start'
	
	ni = qcmplx(qsqrt(2q0), 0.0q0)
	data eps(0:1) /1q0, 2.0q0/ 
	data mu(0:1) /1q0, 1.0q0/ 
	
	matrix_size = 10
	f = 1
	
	arg(1) = 0.5q0
	write(*,*) 'check_delta'
	call first%set(1, matrix_size, qcmplx(2q0,0q0), 1.2q0, 1, arg, 1)
	call first%calculate
	call second%set(1, matrix_size, qcmplx(2q0,0q0), 1.2q0, 1, arg, 1)
	call second%calculate
	call calculate_delta(first, second, tmatr, 10)
	
	write(*,*) real(tmatr(1:5, 1:5),4)
	call scatter%set(ni, lambda, rv, ab, alpha, matrix_size, eps, mu, f)
	
	write(*,*) 'all set'
	tmatr = scatter%get_tmatrix(1)
	write(*,*) 'TE tmatrix'
	call check_symmetric(tmatr, matrix_size)

	tmatr = scatter%get_tmatrix(2)
	write(*,*) 'TM tmatrix'
	call check_symmetric(tmatr, matrix_size)

	
	ext = scatter%get_extinction_factor(1)
	sca = scatter%get_scattering_factor(1)
	
	write(*,*) 'TE: ext = ', ext, 'sca = ', sca


	ext = scatter%get_extinction_factor(2)
	sca = scatter%get_scattering_factor(2)
	
	write(*,*) 'TM: ext = ', ext, 'sca = ', sca
	
end program test_ext_sca_real
