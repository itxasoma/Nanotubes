		subroutine NEWANGLE (Sold, Snew)
!-----------------------------------------------------------------------------------------------------------------------------------------
!		Finds a new spin orientation according to different schemes: Random, Ising, Cone
!-----------------------------------------------------------------------------------------------------------------------------------------
		use GLOBAL
		implicit none

		real(8) :: r1,r2,r3,r4,xsq,x1,x2, xsqsq,ranmar
		real(8), dimension(3) :: ra
		real(8), intent(out),dimension(:) :: Snew(DimeS)
		real(8), intent(in), dimension(:) :: Sold(DimeS)
		real(8) :: theta, phi

		if (NewOrientation.eq."Random") then
			 call RANDOM_NUMBER(r1)
			 call RANDOM_NUMBER(r2)
			 theta=dacos(2.0d0*r1-1.0d0)
			 phi=2.0d0*pi*r2      
			 Snew(1)=dsin(theta)*dcos(phi)
			 Snew(2)=dsin(theta)*dsin(phi)
			 Snew(3)=dcos(theta)
		else if (NewOrientation.eq."Ising") then 
			Snew=Sold   	
			Snew(3)=-Sold(3)   	
		else if (NewOrientation.eq."Cone") then
			xsq = 1.0d0 + 1.d-6
			do while(xsq.ge.1.0d0)
				call RANDOM_NUMBER(r3)
				call RANDOM_NUMBER(r4)
	!           r1=ranmar()
	!           r2=ranmar()
				x1=(2.d0*r3-1.d0)
				x2=(2.d0*r4-1.d0)
				xsq=x1**2+x2**2
			enddo
			xsqsq=sqrt(1.d0-xsq)
			ra(1)=2.d0*x1*xsqsq
			ra(2)=2.d0*x2*xsqsq
			ra(3)=1.d0-2.d0*xsq
			Snew(:) = Sold(:)+rmax*ra(:)                        !New spin direction
			Snew(:) = Snew(:)/sqrt(dot_product(Snew,Snew))   !New normalized spin
	!        print *,'Cone... ',Snew
		 endif

		end subroutine NEWANGLE
