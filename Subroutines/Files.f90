		subroutine FILES
!------------------------------------------------------------------------------
!		Opens units to write files
!		Creates file names according to input geometric parameters fo the nanotubes
!------------------------------------------------------------------------------
		use GLOBAL
		implicit none
		character(35) ::  char_sysdim
		character(:), allocatable :: output_prefix
        character(10) :: char_Kani
!		char_sysdim = trim(adjustl(char_Nx))//'_'//trim(adjustl(char_Ny))//'_'//trim(adjustl(char_Nz))

		! Create a prefix string by concatenating different character variables and
		! tags (using //).
		! The function trim(adjustl()) is used to get rid of trailing whitespaces.
!		output_prefix = trim(adjustl(char_sysdim))//&
!				'_Hzee_'//trim(adjustl(char_Hzee))//&
!				'_Dij_'//trim(adjustl(char_Dij))//&
!				'_Jexc_'//trim(adjustl(char_Jexc))//&
!				'_Kani_'//trim(adjustl(char_Kani))//'_'

		write(chargdip,'(f6.3)') gdip
		write(charrad,'(i0.2)' ) Natoms
		write(charlay,'(i0.2)') Nlayers
!		write(chargdip,9) gdip
!		9   format("_g",f5.3)
!		write(charrad,10 )Natoms
!		10  format("R",i0.2)
!		write(charlay,11) Nlayers
!		11  format("_L",i0.2)

		write(charn1,'(i2.0)') n1
		write(charn2,'(i2.0)') n2
		write(charlength,'(f4.1)') Tlength
        write(char_Kani, '(f6.3)') K_ani

        output_prefix = trim(adjustl(char_sysdim)) // &
                    trim(adjustl(charrad)) // trim(adjustl(charlay)) // &
                    trim(adjustl(charlength)) // '_K_' // trim(adjustl(char_Kani))

		select case (Version)
		case ("Zigzag")
            char_sysdim = 'R_' // trim(adjustl(charrad)) // &
                            '_L_' // trim(adjustl(charlay)) // &
                            '_g_' // trim(adjustl(chargdip)) // &
                            '_K_' // trim(adjustl(char_Kani))

			! Files in different subroutines

				! File for the lattice positions
!				open(101,file= trim(FilePos)//'_'//char_sysdim//'.out')
				! File for the spin configuration
!				open(102,file= trim(FileSpin)//'_'//char_sysdim//'.out')
				! File for the neighbours (just for checking)
!				open(103,file= trim(FileNneigh)//'_'//char_sysdim//'.out')

			! Files in main program

				! File for the energies, magn's, C_e, Susc @ different temperatures
				open(100,file= "Ener_"//trim(adjustl(char_sysdim))//'.out')
				! File with the final configuration
				open(200,file= "Conf_end_"//trim(adjustl(char_sysdim))//'.out')
				! Files with the final configuration, but in unfolded coodinates
!!				open(201,file= "Conf_plane_"//trim(adjustl(char_sysdim))//'.out')

		case ("Angle")
            char_sysdim = 'N1_' // trim(adjustl(charn1)) // &
                            '_N2_' // trim(adjustl(charn2)) // &
                            '_Tl_' // trim(adjustl(charlength)) // &
                            '_g_' // trim(adjustl(chargdip)) // &
                            '_K_' // trim(adjustl(char_Kani))
				! Files in different subroutines
!				open(102,file= trim(FileSpin)//'_'//trim(adjustl(char_sysdim))//'.out')
!				open(103,file= trim(FileNneigh)//'_'//trim(adjustl(char_sysdim))//'.out')

				! Files in main program
				open(100,file= "Ener_"//trim(adjustl(char_sysdim))//'.out')
				open(200,file= "Conf_end_"//trim(adjustl(char_sysdim))//'.out')
!!				open(201,file= "Conf_plane_"//trim(adjustl(char_sysdim))//'.out')
		end select

		end subroutine FILES
