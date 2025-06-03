		include 'Parameters.f90'
		include 'Subroutines.f90'
!-----------------------------------------------------------------------------------------------------------------------------------------
!		Main program: Nanotubes
!		cd  Users/Oscar\ Iglesias/OneDrive\ -\ Universitat\ de\ Barcelona/TFG\ Oriol\ 2018/New_Code/
!		cd /mnt/c/Users/oscar/OneDrive\ -\ Universitat\ de\ Barcelona/TFG_MarcG_2020_Mine
!		cd /mnt/c/Users/oscar/OneDrive\ -\ Universitat\ de\ Barcelona/TFG_Itxaso_2024_Mine/Code
!    gfortran -O4 -o Nanotubes_main Nanotubes_main.f90
!															By Oscar Iglesias & Marc Gris , Started: Feb-June 2020
!-----------------------------------------------------------------------------------------------------------------------------------------
		program MAIN
		use GLOBAL
		use SUBROUTINES
		implicit none

		integer :: i,j,k,acc_total,preMC
		integer, allocatable, dimension(:) :: Labels
		real(8) :: phi,theta,maxz

		real(8) :: Norm,C_e,Susc
		IJ = 25347     !Llavors per als nombres aleatoris
		KL = 15122
		call sranmar(IJ,KL)

		call READDATA
		call FILES
		if (version.eq.'Zigzag') then
			call GEOMETRY
		else if (version.eq.'Angle') then
			call GEOMETRY_OSCAR
		endif
		print *,'Finished Geometry '


		allocate(Spin(Nspin,DimeS))
		allocate(Nneigh(Nspin,6))	! Depends on the lattice: 4 for a squared lattice
		allocate(Nn_num(Nspin))
		allocate(Hexch(Nspin,DimeS))
		allocate(Hdip(Nspin, DimeS))
		allocate(Nani(Nspin,Dime))

		call CONFIGURATION
		print *,'Finished Configuration '

!	This was added by Oscar 26th June 2020 to get spin positions ordered according to S_z component
! In this way, in Zigzag and Angle generated tubes spins are in the same order
!		if (version.eq.'Angle') then
!				allocate(Labels(Nspin))
!				call ORDER(Labels)
!				maxz = maxval(R(i:Nspin,3))
!				open(1,file= 'Config_ord.out')
!				do i=1,NSpin
!					write(1,*) R(i,:), Spin(i,:)
!					write(1,*) R(i,1), R(i,2), R(i,3) + 1.0d0* maxz, Spin(i,:)
!				enddo
!				close(1)
!		endif
!		print *,'Finished Order '
!		stop

		call NEIGHBOURS
		print *,'Finished Neighbours '
		if (Version .eq. 'Zigzag') then
				!call ROTACIONAL
				!print *,'Finished Rotacional '
		else
!				call ROTACIONAL2
		endif
		call Ws
		print *,'Finished Ws '
		call HDIPOLAR
		print *,'Finished HDipolar '
		call HEXCHANGE
		print *,'Finished Exchange '
		call ANISOTROPY
		print *,'Finished Anisotropy '
		call ENERGY
		print *,'Finished Energy'
		call TEMPERATURES
		print *,'Finished Temps '

		print *, 'Initial Energies= ',Ene/Nspin
		print *, 'Initial Magnetization= ',Mag/Nspin



		Norm=1.0d0/real(Nspin)/real(Mctot)

		! Main loop over temperatures stored in Temp(:)
		do i=1, Ntemp
					! Thermalization over preMC steps without averaging Ene nor Magn
					do j=1, preMC
							call MONTE(1/Tempe(i))
					enddo
					print *,'Themalization Finished'

					acc_total = 0
					Ene_av(:) = 0.0d0
					Mag_av(:) = 0.0d0
					Magmod_av = 0.0d0
					Ene2_av=0.0d0
					Mag2_av=0.0d0
					Rot_av(:) = 0.0d0

					! MC loop ver MCtot step: averages Ene, Ene2, Magn, Magn2
					do k=1,MCtot
							call MONTE(1/Tempe(i))
							acc_total=acc_total + acc
							Ene_av(:)=Ene_av(:) + Ene(:)
							Ene2_av=Ene2_av + Ene2
							Magmod_av=Magmod_av + Magmod
							Mag_av(:)=Mag_av(:) + Mag(:)
							Mag2_av=Mag2_av+Magmod**2
							if (Version.eq. 'Zigzag') then
!									call ROTACIONAL
!									Rot_av(:) = Rot_av(:) + Rho(:)
							endif
					enddo
					print *,'Acceptance rate= ',real(acc_total)/real(MCtot)/Nspin*100,' %'

					! Computes the specific heat and susceptibility at each temperature
					!	Cap_e=(Ene_av(1)*Norm)**2-Ene2_av*Norm2)/Tempe(i)**2
					C_e=(dabs(Ene2_av*Norm**2*real(Mctot)- (Ene_av(1)*Norm)**2) /Tempe(i)/Tempe(i))*Nspin
					Susc= (dabs(Mag2_av*Norm**2*real(Mctot)-(Magmod_av*Norm)**2))*Nspin/Tempe(i)
					print *,'T= ',Tempe(i), 'Ene= ', Ene_av*Norm,C_e
					print *, 'Mag= ',Mag_av(:)*Norm,sqrt(dot_product(Mag_av,Mag_av))*Norm
					if (Version.eq. 'Zigzag') then
!							print *, 'Rot= ', Rot_av(:)/real(Mctot)
					endif

					write(100,*) Tempe(i), Ene_av(:)*Norm, Mag_av(:)*Norm, &
					sqrt(dot_product(Mag_av,Mag_av))*Norm, C_e,Susc
		enddo


		do i=1,Nspin
				Phi = atan2(Spin(i,2),Spin(i,1))
				Theta=acos(Spin(i,3))
				write(200,*) R(i,:), Spin(i,:), Theta, Phi     ! Writes final configuration in files Conf_end_
		enddo
		! We do not really need this now (12 March 2025)
		 do i=1,Nspin		! Writes final unfolded configuration in files Conf_plane_
				Phi = atan2(Spin(i,2),Spin(i,1))
				Theta=acos(Spin(i,3))
				if (Version .eq. 'Zigzag') then
						write(201,*) R(i,:),Spin(i,:),Theta,Phi
				elseif (Version .eq. 'Angle') then
						write(201,*) R(i,:),Pos_t(i,:),0.0d0,Spin(i,:),Theta,Phi
				endif
		 enddo

		close(100)
		close(200)
		close(201)

		call DEALLOCATOR

		end program MAIN

! **********************************************
		integer function GCD(a,b)
		implicit none
		integer, intent(in) :: a,b
		integer :: temp
		if (a<b) then
				temp=a
		else
			temp=b
		endif
		do while ((mod(a,temp)/=0).or.(mod(b,temp)/=0))
			temp=temp-1
		enddo
		gcd=temp
		end function GCD

! **********************************************
! Funció ranmar:
! Rep els valors numèrics creats per la subrutina anterior.
! Torna un valor aleatori pertanyent a una distribució uniforme de 0 a 1
! **********************************************

		REAL*8 FUNCTION ranmar()
		IMPLICIT NONE
		REAL*8 RVAL, U(97), C, CD, CM
		INTEGER I97, J97
		COMMON /RANMIN1/ U, C, CD, CM, I97, J97

		RVAL = U(I97) - U(J97)
		IF ( RVAL .LT. 0.D0 ) RVAL = RVAL + 1.0
		U(I97) = RVAL
		I97 = I97 - 1
		IF (I97 .EQ. 0) I97 = 97
		J97 = J97 - 1
		IF (J97 .EQ. 0) J97 = 97
		C = C - CD
		IF ( C .LT. 0.D0 ) C = C + CM
		RVAL = RVAL - C
		IF ( RVAL .LT. 0.D0 ) RVAL = RVAL + 1.0
		IF ( RVAL .EQ. 0.D0 ) RVAL = 2.0D-38

		RANMAR = RVAL
		END function ranmar
