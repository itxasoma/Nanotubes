		subroutine CONFIGURATION
		use GLOBAL
		implicit none
		integer :: i,j, SEED
		real(8) :: theta, phi

		PARAMETER (SEED=16469784)

		call srand(seed)

		do j=1,NSpin
			select case (SpinOrientation)
			case("Alignedz")
				Spin(j,1)=0.0d0
				Spin(j,2)=0.0d0
				Spin(j,3)=1.0d0
			case("Alignedx")
				Spin(j,1)=1.0d0
				Spin(j,2)=0.0d0
				Spin(j,3)=0.0d0
			case("Alignedy")
				Spin(j,1)=0.0d0
				Spin(j,2)=1.0d0
				Spin(j,3)=0.0d0
			case("RandomIs")
				Spin(j,1)=0
				Spin(j,2)=0
				if (rand().lt.0.5d0) then
					Spin(j,3)=1
				else
					Spin(j,3)=-1
				endif
			case("Random")
				theta=dacos(2.0d0*rand()-1.0d0)
				phi=2.0d0*pi*rand()       				!evenly distributed on a sphere surface, next i convert them to cartesian components
				Spin(j,1)=1.0*sin((theta))*cos((phi))
				Spin(j,2)=1.0*sin((theta))*sin((phi))
				Spin(j,3)=1.0*cos((theta))
			case("VortexZ")
				phi = atan2(R(j,2),R(j,1))
				Spin(j,1)=-dsin(phi)
				Spin(j,2)=dcos(phi)
				Spin(j,3)=0.0d0
				! print *,j, phi*180/pi,R(j,3)
			case("VortexX")
				phi = atan2(R(j,2),R(j,1))
				Spin(j,2)=-dsin(phi)
				Spin(j,3)=dcos(phi)
				Spin(j,1)=0.0d0
			case("Normal")
				phi = atan2(R(j,2),R(j,1))
				Spin(j,1)=dcos(phi)
				Spin(j,2)=dsin(phi)
				Spin(j,3)=0.0d0
			end select
		enddo

!		if (Writeaux2(6).eq.'Yes') then
!				open(1,file= trim(FileSpin)//'_ini.out')
!				do i=1,Nspin
			!!			write(1,*) i, Spin(i,:)
					phi = atan2(Spin(i,2),Spin(i,1))
					theta=acos(Spin(i,3))
		!			print *,'THETA= ',theta,' PHI= ',phi,Spin(i,:)
			!!		      write(1,*) R(i,:),Pos_t(i,:),0.0d0,Spin(i,:),dcos(theta)*dcos(phi),dcos(theta)*dsin(phi),dsin(theta)
					! First validated version
		!!			write(1,*) R(i,:),Pos_t(i,:),0.0d0,Spin(i,:),dsin(theta)*dsin(phi),dcos(theta),dsin(theta)*dcos(phi)
!					write(1,*) R(i,:),Spin(i,:), Theta, Phi
!				enddo
!				close(1)
!		endif


		end subroutine CONFIGURATION
