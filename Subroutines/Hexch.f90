		subroutine HEXCHANGE
!-----------------------------------------------------------------------------------------------------------------------------------------      
!		Computes the exchange fields using the neigbours computed  
!		@ beggining of the program in Nneigh.f90
!----------------------------------------------------------------------------------------------------------------------------------------- 
		 use GLOBAL
		 implicit none
		 integer :: i,j, nn

		 Hexch=0.0d0
		 do i=1,Nspin 
			  nn = Nn_num(i)
			  do j=1,nn ! loop over neighbours
					  Hexch(i,:) = Hexch(i,:) + Spin(Nneigh(i,j),:)            
			  enddo
		 enddo

		end subroutine HEXCHANGE
