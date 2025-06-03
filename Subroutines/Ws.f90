!     ########################################################################################
!     SUBROUTINE WITHOUT PERIODIC BOUNDARY CONDITIONS
!     ########################################################################################
     subroutine Ws
      use GLOBAL
      implicit none

      integer :: n, m
      integer :: a, b
      real(8), dimension(Dime) :: rnm         ! Vectorial difference between 2 spins positions
      real(8) :: rmod,rmod2,rmod3          ! The actual distance, as dot product of rnm

!--------- This subroutine calculates and stores in a global W array all the dipolar matrices that will be used
!--------- to calculate the dipolar effective field  

      allocate(W(Nspin,Nspin,6))
      W=0d0


      do n=1, Nspin
          do m=1, Nspin
              if(n.eq.m)then
                  W(n,m,:) = 0.d0 !La interacci� dipolar no afecta a la pr�pia part�cula
              else
                  rnm(:) = R(m, :) - R(n, :) !La dist�ncia �s la difer�ncia entre els vectors posici�     
                  rmod = dsqrt(dot_product(rnm, rnm))
                  rmod3=1/rmod**3
                  rmod2=rmod**2
                        W(n,m,1)=rmod3*(1d0-(3d0*rnm(1)*rnm(1))/rmod2)     ! xx
                        W(n,m,2)=rmod3*(1d0-(3d0*rnm(2)*rnm(2))/rmod2)     ! yy
                        W(n,m,3)=rmod3*(1d0-(3d0*rnm(3)*rnm(3))/rmod2)     ! zz
                        W(n,m,4)=rmod3*(-3d0*rnm(1)*rnm(2))/rmod2     ! xy
                        W(n,m,5)=rmod3*(-3d0*rnm(1)*rnm(3))/rmod2     ! xz
                        W(n,m,6)=rmod3*(-3d0*rnm(2)*rnm(3))/rmod2     ! yz
!                        print *,n,m,W(n,m,3)
              end if
          end do
      end do
      end subroutine Ws