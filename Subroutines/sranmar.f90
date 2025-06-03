! **********************************************
! Subrutina sranmar:
! Rep dos nombres enters ij i kl.
! Torna una sèrie de valors numèrics per mitjà d'un common
! **********************************************
     
      SUBROUTINE sranmar(IJ,KL)
!0 <= IJ <= 31328; 0 <= KL <= 30081
      INTEGER, intent(in) :: IJ, KL
      integer :: I,J,K,L,II,JJ,Mst
      INTEGER :: I97, J97
      REAL*8 :: U(97), C, CD, CM, r, at
      COMMON /RANMIN1/ U, C, CD, CM, I97, J97

!Si ij o kl queda fora del domini el programa s'atura
      IF ( IJ.LT.0.OR.IJ.GT.31328.OR.KL.LT.0.OR.KL.GT.30081) THEN
         STOP
      ENDIF

      I = MOD(IJ/177, 177) + 2
      J = MOD(IJ    , 177) + 2
      K = MOD(KL/169, 178) + 1
      L = MOD(KL,     169)

      DO 2 II = 1, 97
         r = 0.D0
         at = 0.5D0
         DO 3 JJ = 1, 24
            Mst = MOD(MOD(I*J, 179)*K, 179)
            I = J
            J = K
            K = Mst
            L = MOD(53*L+1, 169)
            IF (MOD(L*Mst, 64) .GE. 32) THEN
               r = r + at
            ENDIF
            at = 0.5D0 * at
3        CONTINUE
         U(II) = r
2     CONTINUE

      C = 362436.D0 / 16777216.D0
      CD = 7654321.D0 / 16777216.D0
      CM = 16777213.D0 /16777216.D0

      I97 = 97
      J97 = 33

      END subroutine sranmar



