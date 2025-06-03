module Subroutines
		use GLOBAL
		implicit none

		contains
! Here comes the list of files to be included with the different Subroutines used throughout the code
		include "Subroutines/sranmar.f90"
!		include "Subroutines/Readdata.f90"		!llegeix data d'un arxiu i allotja memòria
        include "Subroutines/Readdata_2.f90"        !llegeix data d'un arxiu i allotja memòria
		include "Subroutines/Geometry.f90"		!posició dels àtoms
		include "Subroutines/Geometry_Oscar.f90"		!posició dels àtoms
		include "Subroutines/Config.f90"		!configuració spins
		include "Subroutines/Nneigh.f90"		!matriu primers veins
		include "Subroutines/Hexch.f90"			!camp exchange energia camps locals
		include "Subroutines/Ws.f90"			!matriu dipolar
		include "Subroutines/Hdip.f90"			!camp energia dipolar
		include "Subroutines/Energy.f90"
		include "Subroutines/Deallocate.f90"	!Desallotja memòria
		include "Subroutines/Monte.f90"
		include "Subroutines/Newangle.f90"
		include "Subroutines/Updatefields.f90"
		include "Subroutines/DeltaE.f90"
		include "Subroutines/Files.f90"
		include "Subroutines/Nani.f90"			!Vector anisotropia
		include "Subroutines/Tempes.f90"
		include "Subroutines/Rotacional.f90"
		include "Subroutines/Order.f90"

		end module Subroutines
