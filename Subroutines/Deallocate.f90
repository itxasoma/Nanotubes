	subroutine DEALLOCATOR
	use GLOBAL
	implicit none
	if (allocated(R)) deallocate(R)	
	if (allocated(Spin)) deallocate(Spin)
	if (allocated(Nneigh)) deallocate(Nneigh)
	if (allocated(Nn_num)) deallocate(Nn_num)
	if (allocated(Hexch)) deallocate(Hexch)
	if (allocated(Hdip)) deallocate(Hdip)
	if (allocated(Nani)) deallocate(Nani)
	if (allocated(tempe)) deallocate(tempe)
!	if (allocated(Labels)) deallocate(Labels)
	
	end subroutine DEALLOCATOR
