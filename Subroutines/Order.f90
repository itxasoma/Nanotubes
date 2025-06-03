			subroutine ORDER(Labels)
	! ***************************************************************************
	!	This orders the coordinates of lattice sites in order of increasing z, y, x  component (in this order).
	!		26th June 2020	
	!	Next goal is to keep track of the spin number label once the array is ordered.
	!		Completed 20th July 2020
	! 	Got inspiration from: https://stackoverflow.com/questions/35324153/sorting-arrays-by-rows-in-fortran
	!	
	!	It returns: 	- R(i,j) ordered
	!					- Labels(i=1,Nspins): original spin labels corresponding to the order array.
	!
	! ***************************************************************************
		use GLOBAL
		implicit none

		integer :: i,j,k,locmaxz
		real(8) :: maxz
		real(8) :: tmp(DimeS)
		real(8) :: R_ord(Nspin,DimeS)
		integer,intent(out) :: Labels(Nspin)

		do i=1, Nspin
			do j=1,3
				if ((abs(R(i,j)).gt.0).and.(abs(R(i,j)).lt.1.0d-10)) then
					R(i,j)=0.0d0
				endif
			enddo
		enddo
		goto 1

!	This part orders R only according to increasing 3rd column (NOT used anymore)
		do i=1,Nspin	
			maxz=minval(R(i:Nspin,3))
			locmaxz=minloc(R(i:Nspin,3),dim=1)+i-1
			print *,'MAXIMUM VALUE= ',locmaxz, maxz
			tmp(:)=R(i,:)
			R(i,:)=R(locmaxz,:)
			R(locmaxz,:)=tmp(:)
	!		if (R(i,3).lt.maxz) then
	!	 	else
	!		endif
		enddo
		goto 2

	1	num_cols=3
		num_rows=Nspin

		call SORTROW_MINE(R,Labels)
		open(33, file='labels.out')
			do k=1,Nspin
				write(33,*)k,labels(k)
			enddo
		close(33)
		!	   do k=1,Nspin
		!		   write(34,*)k,a(k,:)
		!	   enddo
		!	   write(34,*)'** ', i

	2	print*,'end'
		end subroutine ORDER

!	This subroutine orders R according to increasing 1st, 2nd, 3rd columns
		subroutine SORTROW_MINE(a,Labels)
		use GLOBAL
		implicit none
		real(8), intent(inout) :: a(num_rows, num_cols)
		integer :: i, j,k,itmp
		real(8) :: tmp(num_cols)
		integer,intent(out) :: Labels(Nspin)

		labels=(/(i,i=1,Nspin)/)

	!	open(34, file='Conford.out')

		do i = 1,num_rows
			do j = i+1, num_rows
			 if (islarger(a(i,:), a(j,:))) then
				itmp=labels(i)
				labels(i)=labels(j)
				labels(j)=itmp
				tmp(:) = a(i,:)
				a(i,:) = a(j,:)
				a(j,:) = tmp(:)
			 end if
			end do
		end do

	!	close(34)
		end subroutine SORTROW_MINE

		function ISLARGER(a, b)
		implicit none
		real(8), intent(in) :: a(num_cols), b(num_cols)
		logical :: islarger
		integer :: i
	!!	do i = 1, num_cols !! Uncomment if order is 1t,2nd,3rd columns
		do i = num_cols,1,-1 !! Uncomment if order is 3rd,2nd,1st columns
			if (a(i) > b(i)) then
			 islarger = .TRUE.
			 return
			end if
			if (b(i) > a(i)) then
			 islarger = .FALSE.
			 return
			end if
		end do
		islarger = .FALSE.
		return
		end function ISLARGER


	! 
	! From: https://stackoverflow.com/questions/49141698/how-can-you-sort-rows-in-an-array-based-on-ascending-column-values
	!
		subroutine sortrow(a)
		use GLOBAL
		implicit none
		real(8), intent(inout) :: a(num_cols, num_rows)
		integer :: i, j
		integer :: tmp(num_cols)

		do i = 1, num_rows
			do j = i+1, num_rows
			 if (islarger(a(:,i), a(:,j))) then
				  tmp(:) = a(:, i)
				  a(:, i) = a(:, j)
				  a(:, j) = tmp(:)
			 end if
			end do
		end do
		end subroutine sortrow

