		subroutine TEMPERATURES
!-----------------------------------------------------------------------------------------------------------------------------------------      
!		Computes and stores the intermediate Ntemp temperatures in array Tempe 
!		Two possible protocols: 1) Linear, 2) Annealing
!----------------------------------------------------------------------------------------------------------------------------------------- 
		use global
		implicit none

		integer :: i
		real(8) :: deltaT, Temp

		if (protocol.eq."Linear") then
			deltaT=Dtemp
			Ntemp = nint((tempini-tempfin)/deltaT)
			print*, "Ntemp= ",Ntemp," Tini,Tfin= ",tempini,tempfin,Dtemp
			allocate(Tempe(Ntemp))
			do i=1,Ntemp
				Tempe(i) = tempini - i*deltaT
			enddo
		else if (protocol.eq."Annealing") then
			Ntemp = 1 + nint((log(tempfin/tempini))/log(Dtemp)) 
			print*,"Ntemp= ",Ntemp
			allocate(Tempe(Ntemp))
			Temp=tempini
			do i=1,Ntemp
				Tempe(i)= temp*dtemp
				Temp=Tempe(i)
			enddo
      endif

		end subroutine