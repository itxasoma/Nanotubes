		subroutine ENERGY
!--------------------------------------------------------------------------------------------------------
!		Computes the energies: dipolar (2), exchange (3), anisotropy (4), total (1)
!		Computes the magnetization
!--------------------------------------------------------------------------------------------------------
		use GLOBAL
		implicit none
		integer :: i

		Ene=0.d0
		Mag=0.d0
		do i=1,Nspin
				Ene(2)=Ene(2) - dot_product(Spin(i,:),Hdip(i,:))                     	! Dipolar energy
				Ene(3)=Ene(3) - dot_product(Spin(i,:),Hexch(i,:))                   	! Exchange energy
				Ene(4)=Ene(4) - (dot_product(Spin(i,:),Nani(i,:)))**2            	! Anisotropy energy
				Mag(:)=Mag(:)+Spin(i,:)
		enddo
		Ene(1)=gdip*Ene(2)+J_ex*Ene(3)+K_ani*Ene(4)

   	end subroutine ENERGY