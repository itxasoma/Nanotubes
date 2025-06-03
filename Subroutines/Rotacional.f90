		subroutine ROTACIONAL
!-----------------------------------------------------------------------------------------------------------------------------------------
!		Finds the vorticity rho(:) and rot (to be completed and checked)
!-----------------------------------------------------------------------------------------------------------------------------------------
		use GLOBAL
		implicit none

		integer :: i,j,k,l
		real(8) :: dx,dy,dz
		real(8), dimension(Natoms,Nlayers+1) :: Sx,Sy,Sz
		real(8), dimension(Natoms,Nlayers+1) :: Rx,Ry,Rz
		integer,dimension(Natoms) :: hh
		integer,dimension(Nlayers+1) :: pp

		do i=1, Nspin
				k=mod(i-1,Natoms)+1     !num atom de la capa j
				l=(i-1)/Natoms+1    !Capa
				Sx(k,l)=Spin(i,1)
				Sy(k,l)=Spin(i,2)
				Sz(k,l)=Spin(i,3)
				Rx(k,l)=R(i,1)
				Ry(k,l)=R(i,2)
				Rz(k,l)=R(i,3)
!				print *,'i,k,l',i,k,l,Sx(k,l),Sy(k,l),Sz(k,l)
		enddo
		hh(:) = 0
		pp(:) = 0

		do k=1, Natoms-1
				hh(k) = k + 1
		enddo
		hh(Natoms) = 1

		do l=1, Nlayers
				pp(l)=l +1
		enddo

		Rho(:) = 0.0d0
		Delrho(:) = 0.0d0

		do i=1,Nspin
				k=mod(i-1,Natoms)+1     !num atom de la capa j
				l=(i-1)/Natoms+1  !Capa
!				print *,'************  ', k,l,hh(k),pp(l)
				dx=Rx(hh(k),l)-Rx(k,l)
				dy=Ry(hh(k),l)-Ry(k,l)
				dz=Rz(k,pp(l))-Rz(k,l)
				print *,'************  ', k,l,dx,dy,dz

				if (abs(dx) <= 1.0d-4) then ! 3,7 , even layers
						Delrho(1) =  (Sz(hh(k),l)-Sz(k,l))/dy - (Sy(k,pp(l))-(Sy(k,l)))/dz
						Delrho(2) =  - (Sx(k,pp(l))-(Sx(k,l)))/dz
						Delrho(3) =  - (Sx(hh(k),l)-(Sx(k,l)))/dy
						print *,'A', k, l, Delrho(:),dx,dy,dz
				elseif (abs(dy) <= 1.0e-4) then ! 1,5 , even layers
						print *,'BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBb'
						Delrho(1) =  - (Sy(k,pp(l))-(Sy(k,l)))/dz
						Delrho(2) =  (Sz(hh(k),l)-Sz(k,l))/dx - (Sx(k,hh(l))-(Sx(k,l)))/dz
						Delrho(3) =  (Sy(hh(k),l)-Sy(k,l))/dx
						print *,'B',k, l, Delrho(:),dx,dy,dz
				elseif (abs(dz) <= 1.0e-4) then ! None
						Delrho(1) =  (Sz(hh(k),l)-Sz(k,l))/dy
						Delrho(2) =  (Sz(hh(k),l)-Sz(k,l))/dx
						Delrho(3) =  (Sy(hh(k),l)-Sy(k,l))/dx - (Sx(hh(k),l)-(Sx(k,l)))/dy
						print *,'C' ,k, l, Delrho(:),dx,dy,dz
				else ! 2, 4, 6, 8 , even layers // all  , odd layers
						Delrho(1) = (Sz(hh(k),l)-Sz(k,l))/dy - (Sy(k,pp(l))-(Sy(k,l)))/dz
						Delrho(2) =  (Sz(hh(k),l)-Sz(k,l))/dx - (Sx(k,pp(l))-(Sx(k,l)))/dz
						Delrho(3) =  (Sy(hh(k),l)-Sy(k,l))/dx - (Sx(hh(k),l)-(Sx(k,l)))/dy
						print *,'D', k, l, Delrho(:),dx,dy,dz
				endif
				Rho(:) = Rho(:)+Delrho(:)
		enddo
		Rot = sqrt(rho(1)**2+rho(2)**2+rho(3)**2)/Nspin
		Rho(:) = Rho(:)/Nspin
		print *,'Rho= ',Rho(:),'   Rot= ', Rot
		end subroutine ROTACIONAL