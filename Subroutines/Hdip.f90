      subroutine HDIPOLAR
!-----------------------------------------------------------------------------------------------------------------------------------------      
!		Computes the dipolar fields using the dipolar interaction matrix computed 
!		@ beggining of the program in Ws.f90
!----------------------------------------------------------------------------------------------------------------------------------------- 
      use GLOBAL
      implicit none

      integer :: n, m

      Hdip=0.d0

      do n=1,Nspin
            do m=1,Nspin
                  Hdip(n,1)= Hdip(n,1) + (W(n,m,1)*Spin(m,1)+W(n,m,4)*Spin(m,2) &   
                  +W(n,m,5)*Spin(m,3))          !   COMPONENT X DEL CAMP DIPOLAR QUE PATEIX n  Wxx SX +Wxy Sy + Wxz Sz

                  Hdip(n,2)= Hdip(n,2) + (W(n,m,4)*Spin(m,1)+W(n,m,2)*Spin(m,2) &
                  +W(n,m,6)*Spin(m,3))          !   COMPONENT Y

                  Hdip(n,3)= Hdip(n,3) + (W(n,m,5)*Spin(m,1)+W(n,m,6)*Spin(m,2) &
                  +W(n,m,3)*Spin(m,3))          !   COMPONENT Z
            end do
      end do
		Hdip=-Hdip

      end subroutine HDIPOLAR