      subroutine UPDATEFIELDS(n, Dspin)
!-----------------------------------------------------------------------------------------------------------------------------------------
!		Updates the dipolar and exchange field when a MC trial is accepted
!-----------------------------------------------------------------------------------------------------------------------------------------
		use GLOBAL
		implicit none
		integer,intent(in) :: n
		integer :: m,Nei
		real(8),dimension(3) :: DSpin

		! Update dipolar fields
		do m=1,Nspin
			if (n.ne.m) then
				Hdip(m,1)= Hdip(m,1) -((W(n,m,1)*DSpin(1)+W(n,m,4)*DSpin(2)+W(n,m,5)*DSpin(3)))          
				Hdip(m,2)= Hdip(m,2) -((W(n,m,4)*DSpin(1)+W(n,m,2)*DSpin(2)+W(n,m,6)*DSpin(3)))         
				Hdip(m,3)= Hdip(m,3) -((W(n,m,5)*DSpin(1)+W(n,m,6)*DSpin(2)+W(n,m,3)*DSpin(3)))          
			endif
		enddo
		! Update exchange fields of the Nn_num(n) neigbours of n
		do m=1,Nn_num(n)
			Nei=Nneigh(n,m)									! Label of the mth neigh of n
			Hexch(Nei,:)=Hexch(Nei,:)+DSpin(:)	! Update exch field of Nei
		enddo     
      end subroutine UPDATEFIELDS 