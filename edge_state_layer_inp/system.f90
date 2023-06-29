module system
implicit none
	complex(kind=8), allocatable,private		::	hops(:,:,:),rspauli(:,:,:,:)	!real space Hamiltonian and spin
	integer,private					::	num_wann,rvecnum	! dimension of Hamiltonian
	integer,private,allocatable			::	irvec(:,:), nrpts(:)
	complex(kind=8),allocatable,private		::	momentum(:,:,:)

	public						::	getH, setH_hoppingrspauli, getspin


contains



subroutine setH_hoppingrspauli(B,return_num_wann,MPI_COMM_WORLD,my_id,master)
!set up the realspace hamiltonian and spin operator
	real(kind=8),intent(in)			::	B(3)	!exchangefield
	integer					::	num_lines,ix,iy,iz,band1, band2,dir
	real(kind=8)				::	rdum,idum
	integer,intent(in)			::	my_id,master
	
	integer, intent(in)			::	MPI_COMM_WORLD
	integer					::	ierr

	integer					::	i,j,ii
	integer					::	return_num_wann
	logical					::	fileexist


!get numwann
inquire(file='hopping.1',exist=fileexist)
if(.not.fileexist)then
	write(*,*)	'hopping.1 does not exists, calling MPI_ABORT'
	call MPI_Abort(MPI_COMM_WORLD, '1', ierr)
	return
endif
         open(200,file='hopping.1')
         num_lines=0
         num_wann=0
         do
            read(200,fmt=*,end=311)ix,iy,iz,band1,band2,rdum,idum
            num_lines=num_lines+1
            num_wann=max(num_wann,band1)
         enddo
 311     continue
	return_num_wann=num_wann

!read in hopping.1
         rvecnum=num_lines/(num_wann*num_wann)
	if(my_id.eq.master)	 write(*,*)"num_lines=",num_lines
        if(my_id.eq.master)	 write(*,*)"rvecnum=",rvecnum
        if(my_id.eq.master)	 write(*,*)"num_wann",num_wann
         allocate( hops(1:num_wann,1:num_wann,rvecnum) )
         allocate( irvec(3,rvecnum) )
         hops=0.0d0
         rewind(200)
         num_lines=0
         do
            read(200,fmt=*,end=300)ix,iy,iz,band1,band2,rdum,idum
            num_lines=num_lines+1
            rvecnum=(num_lines-1)/(num_wann*num_wann)+1
            irvec(1,rvecnum)=ix
            irvec(2,rvecnum)=iy
            irvec(3,rvecnum)=iz
            hops(band1,band2,rvecnum )=cmplx(rdum,idum)
         enddo
 300     continue
         close(200)

!read in rspauli.1
inquire(file='rspauli.1',exist=fileexist)
if(.not.fileexist)then
	write(*,*)	'rspauli.1 does not exists, calling MPI_ABORT'
	call MPI_Abort(MPI_COMM_WORLD, '1', ierr)
	return
endif
        allocate(rspauli(1:num_wann, 1:num_wann, 3, rvecnum))
        open(201,file='rspauli.1')
        num_lines=0
        Do
           read(201, fmt=*,end=500) ix,iy,iz,band1,band2,dir,rdum,idum
           num_lines=num_lines+1
           rvecnum=(num_lines-1)/(num_wann*num_wann*3)+1 !number of r vectors
           rspauli(band1, band2, dir, rvecnum)=cmplx(rdum,idum)
        End Do
 500    continue
        close(201)

	allocate(nrpts(rvecnum))
	open(14,file='nrpts_inp')
	do j=1,rvecnum/15
		read(14,'(15I5)')	(nrpts(15*(j-1)+i) ,i=1,15)
	enddo
	read(14,'(<mod(rvecnum,15)>I5)')	(nrpts(15*(rvecnum/15)+i) ,i=1,mod(rvecnum,15))
	close(14)


!add exchange field to real space hopping hamiltonian
	do ii=1,rvecnum
		do j=1,num_wann
			do i=1,num_wann
				hops(i,j,ii)=hops(i,j,ii)+0.5d0*(rspauli(i,j,1,ii)*B(1)+rspauli(i,j,2,ii)*B(2)+rspauli(i,j,3,ii)*B(3))
			enddo
		enddo
	enddo

	if(num_wann.eq.0)then
		write(*,*) 'trying to set Hamiltonian with zero dimension, calling MPI_ABORT'
		call MPI_Abort(MPI_COMM_WORLD, '1', ierr)
	endif



	write(*,*) 'setH_hoppingrspauli done'
end subroutine setH_hoppingrspauli

subroutine getH(Hamiltonian,kpt)
!gives back the fourier transformed Hamiltonian for one k vector
	real(kind=8), intent(in)			 	::	kpt(:)
	complex(kind=8), intent(out), allocatable		::	Hamiltonian(:,:)
	real(kind=8)						::	phas
	complex(kind=8)						::	fac


	integer							::	i,j,ii


	allocate(Hamiltonian(num_wann,num_wann))	
	Hamiltonian=cmplx(0.0d0,0.0d0)

         Do ii=1,rvecnum
		phas=irvec(1,ii)*kpt(1)+irvec(2,ii)*kpt(2)+irvec(3,ii)*kpt(3)

		fac=cmplx(cos(phas),sin(phas))
		do i=1,num_wann
		       do j=1,num_wann
				Hamiltonian(i,j)=Hamiltonian(i,j)+fac*hops(i,j,ii)/nrpts(ii)
			end do
		end do
	End Do



end subroutine getH

subroutine getspin(Spin,kpt)
!gives back the fourier transformed Spin operator for one k vector

	real(kind=8), intent(in)			 	::	kpt(:)
	complex(kind=8), intent(out), allocatable		::	Spin(:,:,:)
	real(kind=8)						::	phas
	complex(kind=8)						::	fac
	integer							::	ix,iy,iz

	integer							::	i,j,ii,jj


	allocate(Spin(num_wann,num_wann,3))	
	Spin=cmplx(0.0d0,0.0d0)

	Do ii=1,rvecnum
		ix=irvec(1,ii)
		iy=irvec(2,ii)
		iz=irvec(3,ii)

		phas=     ix*kpt(1)
		phas=phas+iy*kpt(2)
		phas=phas+iz*kpt(3)

		fac=cmplx(cos(phas),sin(phas))

		do i=1,num_wann
			do j=1,num_wann
				do jj=1,3
					Spin(i,j,jj)=Spin(i,j,jj)+fac*rspauli(i,j,jj,ii)/nrpts(ii)
				end do
			end do
		end do
	 End Do


end subroutine getspin

subroutine getmomentum(momentum,kpt,bravv)
	real(kind=8), intent(in)			 	::	kpt(:)
	real(kind=8),intent(in)					::	bravv(3,3)

	complex(kind=8), allocatable,intent(out)		::	momentum(:,:,:)

	real(kind=8)						::	phas
	complex(kind=8)						::	fac,fac2
	real(kind=8)						::	kder

	integer							::	ii,dir,n2,n1,ix,iy,iz

	allocate(momentum(num_wann,num_wann,3))
	momentum=cmplx(0d0,0d0)
	Do ii=1,rvecnum
		ix=irvec(1,ii)
		iy=irvec(2,ii)
		iz=irvec(3,ii)

		phas=     ix*kpt(1)
		phas=phas+iy*kpt(2)
		phas=phas+iz*kpt(3)

		fac=cmplx(-sin(phas),cos(phas))	!inculdes i of formula	3.11 master thesis

		do dir=1,3
			kder=bravv(dir,1)*ix+bravv(dir,2)*iy+bravv(dir,3)*iz 
			fac2=fac*kder
			do n2=1,num_wann
				do n1=1,num_wann
					momentum(n1,n2,dir)=momentum(n1,n2,dir)+fac2*hops(n1,n2,ii)/nrpts(ii)
				end do
			end do
		end do
	End Do
end subroutine getmomentum

end module system

