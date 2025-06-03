! -------------------------------------------------------------------------------------------
! 	Computes the energy difference of a spin n before and after the change
! -------------------------------------------------------------------------------------------  
      subroutine DELTAE(n, Sold, Snew, Dene)
      use GLOBAL
      implicit none

      integer, intent(in) :: n            ! Spin to be changed
      real(8), intent(in), dimension(DimeS) :: Snew, Sold
      real(8), intent(out), dimension(DimeE) :: Dene
      real(8), dimension(dimeS) :: Dspin

      Dspin= Snew - Sold
      Dene=0.0d0

      Dene(2)= -2.*dot_product(Dspin(:), Hdip(n,:))
      Dene(3)= -2*dot_product(Dspin(:), Hexch(n,:))
      Dene(4)= dot_product(Sold(:),Nani(n,:))**2 - dot_product(Snew(:),Nani(n,:))**2
      Dene(1)= gdip*Dene(2) + J_ex*Dene(3)+ K_ani*Dene(4) 
      
      end subroutine DELTAE