		subroutine MONTE(Beta)
!-----------------------------------------------------------------------------------------------------------------------------------------      
!		Performs a MC step trying to change spins Nspin times per step
!----------------------------------------------------------------------------------------------------------------------------------------- 
		use GLOBAL
		implicit none

		integer :: i,j,k

		real(8) :: Beta
		real(8), dimension(DimeE) :: DeltaEne
		real(8), dimension(:) :: Snew(DimeS), Sold(DimeS), deltaS(DimeS)
		real(8) :: r1, r2


		acc=0.0d0	! Accepted trials
		Ene2=0.d0	
		do k=1,Nspin	! Loop over spins, selecting a spin i @random
			call RANDOM_NUMBER(r1)
			i=int(Nspin*r1)+1 
			Sold(:)=Spin(i,:)
			call NEWANGLE(Sold(:),Snew(:))  ! Selects a trial spin direction                                                               
			DeltaS (:) = Snew(:) - Sold(:)  ! Computes the spin change
			call DELTAE(i, Sold, Snew, DeltaEne)  ! Computes the energy change DeltaEne                        

			call RANDOM_NUMBER(r2)
			if (r2.lt.exp(-DeltaEne(1)*Beta).or.(DeltaEne(1).lt.0.0d0)) then ! Acceptance decision according to Boltzmann Distrib.
				acc=acc+1
				Spin(i,:)=Snew(:)
				call UPDATEFIELDS(i,deltaS)
				Ene(:)=Ene(:)+DeltaEne(:) 
				Mag(:)=Mag(:)+DeltaS(:)
			endif
		enddo	! Ends loop over spins
		Ene2 = Ene(1)**2
		Magmod=sqrt(dot_product(Mag,Mag))

		end subroutine MONTE