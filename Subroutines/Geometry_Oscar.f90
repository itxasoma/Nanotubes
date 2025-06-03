		subroutine GEOMETRY_OSCAR
!-----------------------------------------------------------------------------------------------------------------------------------------
!		This subroutine builts the nanotube by folding a planar lattice along a certain
!		direction given by the chiral
!-----------------------------------------------------------------------------------------------------------------------------------------
		use GLOBAL
		implicit none
		integer :: i,j,k,n
		integer:: count
		integer :: GCD

		integer :: ndiv
		integer :: Np   ! Number of points in the original lattice
!!		character(10) :: case

		real(8) :: a,b
		real(8) :: Radius,Phi
		real(8) :: Chmod,Theta,Tmod
		real(8),dimension(2) :: a1,a2,Ch,T
		real(8),allocatable,dimension(:,:) :: Pos ! Stores the point positions in the atlas

		! Input data:
		alpha=alpha/180.0d0*pi

		! Memory allocation for: original, rotated and trimmed lattices (atlas)
		allocate(Pos((2*nx+1)*(2*ny+1),2))
		allocate(Pos_r((2*nx+1)*(2*ny+1),2))
		allocate(Pos_t((2*nx+1)*(2*ny+1),2))

		! Basis vectors for the original lattice
		a1(1)=1.
		a1(2)=0.
		a2(1)=dcos(alpha)
		a2(2)=dsin(alpha)
		!	Work in progress (10th May 2020): forget for the moment
		!	if (trim(Lattice).eq.'Graphene') then
		!	   alpha=2.*pi/3.
		!		a1(1)=3.0*dcos(alpha)
		!		a1(2)=dsin(alpha)
		!		a2(1)=3.0*dcos(alpha)
		!		a2(2)=-dsin(alpha)
		!		Lattice= 'General'
		!	endif

!-----------------------------------------------------------------------------------------------------------------------------------------
!		Defining the parameters of the nanotube:
!		Chiral vector and length: Ch(:), Chmod
!		Radius: Radius
!		Translation vector: T(:), Tmod
!-----------------------------------------------------------------------------------------------------------------------------------------
! 	Chiral vector Ch and chiral length Chmod for a (n1,n2) nanotube:
! 	After folding, Ch will be directed along the circumf. of the tube
! 	and Chmod will be the perimeter of the tube (2*pi*R).
!-----------------------------------------------------------------------------------------------------------------------------------------
		Ch(:)=float(n1)*a1(:)+float(n2)*a2(:)
		Chmod=sqrt(dot_product(Ch,Ch))
!-----------------------------------------------------------------------------------------------------------------------------------------

		! Chiral angle for a (n1,n2) nanotube:
		if (n2.eq.0) then ! This gives a columnar (AA) tube  for alpha = 90
			Theta=0
			ndiv=n1
		else if (n1.eq.0) then
			Theta=pi/2.
			ndiv=n2
		else
			if (trim(Lattice).eq.'SC') then ! This is the case of Alpha= 90
				! Chiral angle
				Theta=atan(float(n2)/float(n1))
			else if (trim(Lattice).eq.'General') then ! This is the case of Alpha neq 90, For alpha= 60: for all n1 and  n2=0 gives (AB) tube
				Theta=acos(dot_product(Ch,a1)/Chmod)
			endif
			! Greatest common divisor
			ndiv=GCD(n1,n2)
		endif

		Radius=Chmod/pi
!-----------------------------------------------------------------------------------------------------------------------------------------
		! Translation vector
		if (trim(Lattice).eq.'SC') then
			T(1)= (-n2)/ndiv  ! For SC: Modified according to J. Phys.: Condens. Matter, 2009, 21 , 075301
			T(2)= (n1)/ndiv
		else if (trim(Lattice).eq.'General') then
			T(1)= (2*n2+n1)/ndiv  ! Original mine (this is for CNTs or triagular NT's)
			T(2)= -(2*n1+n2)/ndiv
!			T(1)= (-n2)/ndiv  ! For SC: Modified according to J. Phys.: Condens. Matter, 2009, 21 , 075301
!			T(2)= (n1)/ndiv
		endif
		! Tube length
		Tmod=sqrt(dot_product(T,T))
	!!	Tlength=Final_length/Tmod
		Tmod=Tmod*Tlength
!-----------------------------------------------------------------------------------------------------------------------------------------

		print '(a,2(1xf8.4),a,2(1xf8.4))',"Basis vectors, a1= ",a1,", a2= ",a2
		print '(a,2(1xf8.4),2(1xi2))',"Chiral vector= ", Ch,n1,n2
		print '(a,(1xf8.4),a,(1xf8.4))',"Chiral length=",Chmod,", Chiral angle(Theta)= ",Theta*180.0d0/pi
		print '(a,2(1xf8.4),a,1xf8.4,a,1xf8.4)',"Translation vector, T= ",T,", Tmod= ",Tmod,", Tlength= ",Tlength
		print '(a8,f8.4)',"Radius= ",Radius

		! 1) Original lattice
		k=0
		Pos=0.0d0
		if (trim(Lattice).eq.'General') then 	! This is for a general lattice (not really checked if it is ok)
			print *,'General'
			do i=-nx,nx,1
				do j=-ny,ny,1
					k=k+1
					Pos(k,:)=float(i)*a1(:)+float(j)*a2(:)
				enddo
			enddo
		else if(trim(Lattice).eq.'SC') then 	! This is for Squared original lattice
			print *,'SC'
			do i=-nx,nx,1
				do j=-ny,ny,1
					k=k+1
					Pos(k,1)=i
					Pos(k,2)=j
				enddo
			enddo
		endif
		Np=k
		if (Writeaux2(2).eq.'Yes') then
				open(1,file= trim("Original.out"))
					do i=1,Np
						write(1,*) i,Pos(i,:),0.0d0
					enddo
				close(1)
		endif


		! 2) Rotation of the lattice
		do k=1,Np
			Pos_r(k,1)=(Pos(k,1)*n1-Pos(k,2)*n2)/Chmod
			Pos_r(k,2)=(Pos(k,1)*n2+Pos(k,2)*n1)/Chmod
		enddo

		if (Writeaux2(3).eq.'Yes') then
				open(1,file= trim("Rotated.out"))
					do i=1,Np
						write(1,*) i,Pos_r(i,:),0.0d0
					enddo
				close(1)
		endif

		! 3) Trimming the lattice to fit it to the tube unit cell
		j=0
		do i=1,Np
			if (&
			((Pos_r(i,1).le.(Chmod*1.000001)).and.(Pos_r(i,1).ge.(-Chmod*0.9988)))&
			.and.&
			(abs(Pos_r(i,2)).le.(Tmod*1.000001111))&  !!!!!
			) then
				j=j+1
				Pos_t(j,:)=Pos_r(i,:)
			endif
		enddo
		NSpin=j

		if (Writeaux2(4).eq.'Yes') then
				open(1,file= trim("Trimmed.out"))
					do i=1,NSpin
						write(1,*) i,Pos_t(i,:),0.0d0
					enddo
				close(1)
		endif
		print *,"Final number of Spins= ",Nspin

		! 4) Finally, assign the 3D coordinates of the spins in the tube surface
		allocate (R(Nspin,DimeS))

		do i=1,NSpin
			Phi=Pos_t(i,1)/Radius
			R(i,1)=Radius*dsin(phi)
			R(i,2)=Radius*dcos(phi)
			R(i,3)=Pos_t(i,2)
		enddo

		if (Writeaux2(5).eq.'Yes') then
				open(1,file= trim(FilePos)//'_ini.out')
					do i=1, NSpin
						write(1,*) i,R(i,:)
					enddo
				close(1)
				open(1,file= trim("Geo.out")) ! This is only to plot lattice in Mathematica
					write(1,*) Radius
					write(1,*) Tmod
					write(1,*) Nspin
					write(1,*) maxval(R(:,3))
				close(1)
		endif
		print *,'Length= ',2.0d0*maxval(R(:,3))
	! 	deallocate(Pos_t,Pos_r)

		end subroutine GEOMETRY_OSCAR