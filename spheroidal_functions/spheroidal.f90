module spheroidal

use regime

use vb_prolate
use vb_oblate

implicit none
private

	type, public :: SpheroidalCalculation
		
		integer :: m, lnum, spheroidal_type
		integer :: maxd, narg
		complex(knd) :: c
		real(knd) :: ksi
		real(knd), allocatable, dimension(:) :: arg
		integer, allocatable, dimension(:) :: r1_exp, r1d_exp, r2_exp, r2d_exp
		integer, allocatable, dimension(:,:) :: s1_exp, s1d_exp
		complex(knd), allocatable, dimension(:) :: r1, r1d, r2, r2d, r3, r3d
		complex(knd), allocatable, dimension(:,:) :: s1, s1d
		complex(knd), allocatable, dimension(:,:) :: legendre
		
	contains
		procedure :: set, calculate
		final :: delete_calculation
	end type SpheroidalCalculation

	abstract interface
		subroutine  sphfunc(cc,m,lnum,ioprad,x1,r1c,ir1e,r1dc,ir1de,r2c, &
                            ir2e,r2dc,ir2de,iopang,iopnorm,narg,arg, &
                            s1,is1e,s1d,is1de, enr,maxd)
			integer :: m, lnum, ioprad, iopang, iopnorm, narg, maxd
			integer :: ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum)
			integer :: is1e(lnum, narg), is1de(lnum, narg)
			real(16) x1, arg(narg)
			complex(16) :: cc, r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum)
			complex(16) :: s1(lnum, narg), s1d(lnum, narg)
			complex(16), allocatable, dimension(:,:) :: enr
		end subroutine sphfunc
	end interface

contains
	subroutine set(this, m, n, c, ksi, narg, arg, f)
		class(SpheroidalCalculation) :: this
		integer, intent(in) :: m, n, narg
		real(knd), intent(in) :: ksi, arg(narg)
		complex(knd), intent(in) :: c
		integer, optional, intent(in) :: f
		
		if (present(f)) then
			this%spheroidal_type = f
		else
			this%spheroidal_type = 1
		endif
		
		this%m = m
		this%c = c
		if (this%spheroidal_type == 1) then
			this%ksi = ksi - 1q0
		else 
			this%ksi = ksi
		endif
		write(*,*) 'ksi =', this%ksi
		this%narg = narg

		if (.not.(allocated(this%r1) .and. (n - m + 1 /= this%lnum))) then
			if (allocated(this%r1)) then
				call delete_calculation(this)
			endif
			
			this%lnum = n - m + 1
			
			allocate(this%r1(this%lnum), this%r1d(this%lnum), &
			this%r2(this%lnum), this%r2d(this%lnum), this%r3(this%lnum), this%r3d(this%lnum))
			allocate(this%r1_exp(this%lnum), this%r1d_exp(this%lnum), &
			this%r2_exp(this%lnum), this%r2d_exp(this%lnum))
			allocate(this%s1(this%lnum, this%narg), this%s1d(this%lnum, this%narg))
			allocate(this%s1_exp(this%lnum, this%narg), this%s1d_exp(this%lnum, this%narg))
			allocate(this%arg(narg))
			
		endif
		write(*,*) 'allocated calculation'
		this%arg = arg
				
	end subroutine set
	
	subroutine calculate(this)
		class(SpheroidalCalculation), intent(inout) :: this
		complex(16), allocatable, dimension(:,:) :: enr
		procedure (sphfunc), pointer :: sph_ptr => null ()
		integer :: i, j

		if (.not. allocated(this%r1)) then
			write(*,*) 'Nothing to calculate!'
			return
		endif
		
		if (this%spheroidal_type == 1) then
			sph_ptr => cprofcn
		else 
			sph_ptr => coblfcn
		endif
		
!		write(*,*) 'calculateing for c =', this%c
				
		this%r1 = 0
		this%r2 = 0
		this%r3 = 0
		this%r1d = 0
		this%r2d = 0
		this%r3d = 0
		this%r1_exp = 0
		this%r2_exp = 0
		this%r1d_exp = 0
		this%r2d_exp = 0
		this%s1 = 0
		this%s1_exp = 0
		this%s1d = 0
		this%s1d_exp = 0
		call sph_ptr(this%c, this%m, this%lnum, 2, this%ksi, &
		this%r1, this%r1_exp, this%r1d, this%r1d_exp, this%r2, this%r2_exp, this%r2d, this%r2d_exp, &
		2, 1, this%narg, this%arg, this%s1, this%s1_exp, this%s1d, this%s1d_exp, enr, this%maxd)
		
		this%r1 = this%r1 * (10q0**this%r1_exp)
		this%r1d = this%r1d * (10q0**this%r1d_exp)
		this%r2 = this%r2 * (10q0**this%r2_exp)
		this%r2d = this%r2d * (10q0**this%r2d_exp)
		this%s1 = this%s1 * (10q0**this%s1_exp)
		this%s1d = this%s1d * (10q0**this%s1d_exp)
		
		this%r3 = this%r1 + qcmplx(0q0, 1q0) * this%r2
		this%r3d = this%r1d + qcmplx(0q0, 1q0) * this%r2d
		
		if (allocated(this%legendre)) then
			deallocate(this%legendre)
		endif
		this%maxd = 2 * this%maxd + 1
		allocate(this%legendre(this%lnum, 0:this%maxd))
		this%legendre = 0
		do i = 1, this%lnum
			this%legendre(i, mod(i - this%m, 2)) = qcmplx(1q0, 0q0)
			do j = mod(i - this%m, 2) + 2, this%maxd, 2
				this%legendre(i, j) = this%legendre(i, j - 2) * enr(i, j / 2)
			enddo
		enddo

		deallocate(enr)
		do i = 1, this%lnum
			call normalize(this%legendre(i,:), this%maxd)
		enddo
!		write(*,*) 'lnum=', this%lnum
!		write(*,*) 'narg=', this%narg
!		write(*,*) 'r1:', this%r1(1:5)
!		write(*,*) 'r1d:', this%r1d(1:5)
!		write(*,*) 's1:', this%s1(1:10,1)
!		write(*,*) 's1_exp:', this%s1_exp(1:10,1)
!		write(*,*) 'legendre:', this%legendre(1, 0:5)
		
	end subroutine calculate
	
	subroutine delete_calculation(this)
		type(SpheroidalCalculation), intent(inout) :: this
		if (allocated(this%r1)) then
			deallocate(this%r1_exp, this%r1d_exp, this%r2_exp, this%r2d_exp)
			deallocate(this%r1, this%r1d, this%r2, this%r2d, this%r3, this%r3d)
		endif
		if (allocated(this%s1)) then
			deallocate(this%s1, this%s1d)
			deallocate(this%s1_exp, this%s1d_exp)
		endif
		if (allocated(this%legendre)) then
			deallocate(this%legendre)
		endif
		if (allocated(this%arg)) then
			deallocate(this%arg)
		endif
	end subroutine delete_calculation
	
	subroutine normalize(d, maxd)
		integer :: n, r, maxd
		complex(knd) d(0:maxd), norm
		
		norm = 0q0
		do r = 0, maxd
			norm = norm + d(r) ** 2 * (r + 1) * (r + 2) * 2q0 / (2 * r + 3)
		enddo
		norm = cqsqrt(norm)
		d = d / norm
	end subroutine normalize
end module spheroidal