		subroutine READDATA
!-----------------------------------------------------------------------------------------------------------------------------------------
!		Reads external data from *.dat file
!-----------------------------------------------------------------------------------------------------------------------------------------
		use GLOBAL
		implicit none
		character(24) :: input
		integer :: iostat

		! Get input filename from the command line argument
		call get_command_argument(1, input)

		! Check if a filename was provided
		if (input == '') then
				print *, "Error: No input file provided."
				print *, "Usage: ./Nanotubes_main <input_file>"
				stop
		end if

		open (10,file=input, status='old', action='read', iostat=iostat)
			 if (iostat /= 0) then
				  print *, "Error: Unable to open file: ", trim(input)
				  stop
			 end if
			read(10,*) Version
			read(10,*) Lattice
			read(10,*) Natoms,Nlayers
			read(10,*) alpha
			read(10,*) n1,n2
			read(10,*) nx,ny
			read(10,*) Final_length
			read(10,*) Tlength
			read(10,*) J_ex
			read(10,*) gdip
			read(10,*) K_ani
			read(10,*) v_ani
			read(10,*) MCtot
			read(10,*) Protocol
			read(10,*) tempini
			read(10,*) tempfin
			read(10,*) Dtemp
			read(10,*) NewOrientation
			read(10,*) rmax
			read(10,*) SpinOrientation
			read(10,*) Writeaux2
			read(10,*) FilePos
			read(10,*) FileSpin
			read(10,*) FileNneigh
		close(10)
		rmax=dtan(rmax)
		end subroutine READDATA

