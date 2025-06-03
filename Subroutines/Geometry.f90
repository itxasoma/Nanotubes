	subroutine GEOMETRY
!-----------------------------------------------------------------------------------------------------------------------------------------
!		This subroutine was the original one, that builts the zig-zag nanotube as a stack of
!		rings of radius Radius with Natoms stacked in Nlayers
!-----------------------------------------------------------------------------------------------------------------------------------------
		use GLOBAL
		implicit none
		integer :: i,j,k,n
		integer:: count
		real(8) :: a,b
		real(8) :: Radius

		a=1.
		b=sqrt(2.)
		Radius= Natoms*b*a/2.d0/pi
		count=0
		j=0

		Nspin=Natoms*Nlayers
		print *,"Final number of Spins= ",Nspin
		allocate(R(Nspin,DimeS))

		do k=1,Nlayers
			j=j+1
			i=0
			do n=1,Natoms
				i=i+1
				if (mod(j,2).eq.0) then
					count=count+1
					R(count,1)=Radius*cos((2*n+1)*pi/Natoms)
					R(count,2)=Radius*sin((2*n+1)*pi/Natoms)
					R(count,3)=b*a*k/2.
				else
					count=count+1
					R(count,1)=Radius*cos((2*n)*pi/Natoms)
					R(count,2)=Radius*sin((2*n)*pi/Natoms)
					R(count,3)=b*a*k/2.
				endif
			enddo
		enddo
		print *,'Length= ',maxval(R(:,3))-minval(R(:,3))
		if (Writeaux2(1).eq.'Yes') then
				open(1,file= trim(FilePos))
					do i=1,Nspin
						write(1,*) i, R(i,:)
					enddo
				close(1)
		endif
		end subroutine GEOMETRY
