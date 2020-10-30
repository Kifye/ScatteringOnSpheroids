module integrals
use regime
use spheroidal

contains
	subroutine calculate_delta(first, second, Delta, accuracy)
		class(SpheroidalCalculation), intent(in) :: first, second
		real(knd), allocatable, dimension(:) :: common_multiplier
		complex(knd), intent(out) :: Delta(first%lnum, second%lnum)
		integer n, l, r, accuracy
		
		if (allocated(common_multiplier)) then 
			deallocate(common_multiplier)
		endif
		allocate(common_multiplier(0:accuracy))
		do r = 0, accuracy
			common_multiplier(r) = (r + 1) * (r + 2) * 2q0 / (2 * r + 2 * first%m + 1)
		enddo

		Delta = 0q0
				
		write(*,*) 'n = ', first%lnum, 'l=', second%lnum
		write(*,*) 's1 = ', size(first%legendre), 's2=', size(second%legendre)
		do n = first%m, first%lnum
			do l = first%m, second%lnum
				if (mod(abs(n-l), 2) == 1) then
					continue
				endif
				do r = mod(n - first%m, 2), accuracy, 2
					Delta(n, l) = Delta(n, l) + &
					first%legendre(n, r) * second%legendre(l, r) * & 
					common_multiplier(r);
				enddo
			enddo
		enddo
		
!		write(*,*) 'delta = '
!		write(*,*) Delta(1:4, 1:4)
		
		deallocate(common_multiplier)
	end subroutine calculate_delta

end module integrals