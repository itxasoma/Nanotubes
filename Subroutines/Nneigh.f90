		subroutine NEIGHBOURS
!-----------------------------------------------------------------------------------------------------------------------------------------
!		Finds the nearest neighbours of all the spins
!		Matrix Nneigh(i,:) stores all the nn of spin i
!		Nn_num(i) stores the number of nn of spin i
!-----------------------------------------------------------------------------------------------------------------------------------------
		use GLOBAL
		implicit none

		integer :: i,j,nn
		real(8) :: r_ij


		Nneigh=0
		do i=1,Nspin
			 nn=0
			 do j=1,Nspin
				  r_ij=dsqrt(dot_product(R(i,:)-R(j,:),R(i,:)-R(j,:)))
				  if ((r_ij.gt.0d0).and.(r_ij.lt.1.01d0)) then
						nn=nn+1
						Nneigh(i,nn)=j
				  endif
			 enddo
			 Nn_num(i)=nn
		enddo

		if (Writeaux2(7).eq.'Yes') then
				open(1,file= trim(FileNneigh)//'.out')
				do i=1,Nspin
					 write(1,*) i, Nneigh(i,:) ! Writes j nn's of atom i
				enddo
				close(1)
		endif

		end subroutine NEIGHBOURS
