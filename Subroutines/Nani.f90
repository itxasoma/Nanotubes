		subroutine ANISOTROPY
!-----------------------------------------------------------------------------------------------------------------------------------------      
!		Assigns anisotropy directions Nani
!----------------------------------------------------------------------------------------------------------------------------------------- 
 		use GLOBAL
		implicit none

		integer :: j
		real(8) :: phi

		if (v_ani.eq."Uniaxial") then       ! Uniaxial along tube axis
			do j=1, Nspin
					Nani(j,1)=0.0d0
					Nani(j,2)=0.0d0	
					Nani(j,3)=1.0d0
			enddo
		else if (v_ani.eq."Radial") then  ! Radial, perp to tube surface
			do j=1, Nspin
            phi = atan(R(j,2)/R(j,1))
            Nani(j,1)=cos(phi)
            Nani(j,2)=sin(phi)
            Nani(j,3)=0.0d0
			enddo
		endif
    
		end subroutine ANISOTROPY