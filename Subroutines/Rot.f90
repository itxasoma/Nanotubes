SUBROUTINE ComputeCurl(spin, curl, Lx, Ly, Lz)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Lx, Ly, Lz  ! Lattice dimensions
  INTEGER :: N  ! Total number of spins
  REAL(8), INTENT(IN) :: spin(Lx*Ly*Lz, 3)  ! Spin array: (N, 3) where N = Lx*Ly*Lz
  REAL(8), INTENT(OUT) :: curl(Lx*Ly*Lz, 3)  ! Curl array

  INTEGER :: n, nxp, nxm, nyp, nym, nzp, nzm
  INTEGER :: i, j, k

  N = Lx * Ly * Lz

  DO n = 1, N
     ! Compute (i, j, k) from n
     i = MOD(n-1, Lx) + 1
     j = MOD((n-1)/Lx, Ly) + 1
     k = (n-1) / (Lx * Ly) + 1

     ! Identify neighbor indices (handling boundary conditions)
     nxp = n + 1;  IF (i == Lx) nxp = n - Lx + 1
     nxm = n - 1;  IF (i == 1)  nxm = n + Lx - 1

     nyp = n + Lx; IF (j == Ly) nyp = n
     nym = n - Lx; IF (j == 1)  nym = n

     nzp = n + Lx * Ly; IF (k == Lz) nzp = n
     nzm = n - Lx * Ly; IF (k == 1)  nzm = n

     ! Compute curl components
     curl(n,1) = (spin(nyp,3) - spin(nym,3)) / 2.0 - (spin(nzp,2) - spin(nzm,2)) / 2.0
     curl(n,2) = (spin(nzp,1) - spin(nzm,1)) / 2.0 - (spin(nxp,3) - spin(nxm,3)) / 2.0
     curl(n,3) = (spin(nxp,2) - spin(nxm,2)) / 2.0 - (spin(nyp,1) - spin(nym,1)) / 2.0
  END DO

END SUBROUTINE ComputeCurl

! *****************************************************************

MODULE SpinLattice
  IMPLICIT NONE
  CONTAINS

  ! Subroutine to compute the nearest neighbors of each spin
  SUBROUTINE ComputeNeighbors(neighbors, Lx, Ly, Lz)
    INTEGER, INTENT(IN) :: Lx, Ly, Lz  ! Lattice dimensions
    INTEGER, INTENT(OUT) :: neighbors(Lx*Ly*Lz, 6)  ! Neighbor indices

    INTEGER :: n, i, j, k, N

    N = Lx * Ly * Lz

    DO n = 1, N
       ! Compute (i, j, k) from n
       i = MOD(n-1, Lx) + 1
       j = MOD((n-1)/Lx, Ly) + 1
       k = (n-1) / (Lx * Ly) + 1

       ! Identify neighbor indices with periodic boundary conditions
       neighbors(n,1) = n + 1;  IF (i == Lx) neighbors(n,1) = n - Lx + 1
       neighbors(n,2) = n - 1;  IF (i == 1)  neighbors(n,2) = n + Lx - 1

       neighbors(n,3) = n + Lx; IF (j == Ly) neighbors(n,3) = n
       neighbors(n,4) = n - Lx; IF (j == 1)  neighbors(n,4) = n

       neighbors(n,5) = n + Lx * Ly; IF (k == Lz) neighbors(n,5) = n
       neighbors(n,6) = n - Lx * Ly; IF (k == 1)  neighbors(n,6) = n
    END DO
  END SUBROUTINE ComputeNeighbors

  ! Subroutine to compute the curl using the precomputed neighbors
  SUBROUTINE ComputeCurl(spin, curl, neighbors, Lx, Ly, Lz)
    INTEGER, INTENT(IN) :: Lx, Ly, Lz  ! Lattice dimensions
    INTEGER, INTENT(IN) :: neighbors(Lx*Ly*Lz, 6)  ! Neighbor indices
    REAL(8), INTENT(IN) :: spin(Lx*Ly*Lz, 3)  ! Spin array: (N, 3)
    REAL(8), INTENT(OUT) :: curl(Lx*Ly*Lz, 3)  ! Curl array

    INTEGER :: n, nxp, nxm, nyp, nym, nzp, nzm

    DO n = 1, Lx*Ly*Lz
       ! Get neighbor indices from the precomputed table
       nxp = neighbors(n,1)
       nxm = neighbors(n,2)
       nyp = neighbors(n,3)
       nym = neighbors(n,4)
       nzp = neighbors(n,5)
       nzm = neighbors(n,6)

       ! Compute curl components using central differences
       curl(n,1) = (spin(nyp,3) - spin(nym,3)) / 2.0 - (spin(nzp,2) - spin(nzm,2)) / 2.0
       curl(n,2) = (spin(nzp,1) - spin(nzm,1)) / 2.0 - (spin(nxp,3) - spin(nxm,3)) / 2.0
       curl(n,3) = (spin(nxp,2) - spin(nxm,2)) / 2.0 - (spin(nyp,1) - spin(nym,1)) / 2.0
    END DO
  END SUBROUTINE ComputeCurl

END MODULE SpinLattice
! *****************************************************************

PROGRAM Main
  USE SpinLattice
  IMPLICIT NONE

  INTEGER, PARAMETER :: Lx = 10, Ly = 10, Lz = 10  ! Lattice size
  INTEGER :: N
  INTEGER :: neighbors(Lx*Ly*Lz, 6)
  REAL(8) :: spin(Lx*Ly*Lz, 3), curl(Lx*Ly*Lz, 3)

  N = Lx * Ly * Lz

  ! Initialize spin field with some values (example)
  CALL RANDOM_NUMBER(spin)

  ! Compute neighbors once
  CALL ComputeNeighbors(neighbors, Lx, Ly, Lz)

  ! Compute curl using precomputed neighbors
  CALL ComputeCurl(spin, curl, neighbors, Lx, Ly, Lz)

END PROGRAM Main
