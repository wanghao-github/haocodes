module system_layer
implicit none
	complex(kind=8), allocatable,private		::	hops(:,:,:),rspauli(:,:,:,:)	!real space Hamiltonian and spin
	integer,private					::	num_wann,rvecnum	! dimension of Hamiltonian
	integer,private,allocatable			::	irvec(:,:),nrpts(:)		! real space vectors of hops and rspauli, multiples of it
	real(kind=8),allocatable,private		::	localisationpar(:)		!Parameter which defines where in the slab with respect to the layering direction the eigenvector is
	integer,private,allocatable			::	wannierfunctioninHam(:)		!assigns basisvektors of Hamiltonian to wannierfunctions
	integer, private				::	Hdim				!dimension of the Hamiltonian

	public						::	setH_layer,getlocalisationpar,getH_layer

	type 							:: atom
		integer						:: number	!position in array
		real(kind=8)					:: position(3)		!position in realspace
		integer,allocatable				:: wannierfunctions(:)	!numbers of wannierfunctions attributed to the atom
	end type atom

	type(atom), allocatable,private				::	atomarr(:)	!array of all atomtypes in the construction

	integer,private					::	fourdim,layerdir
	integer, private,allocatable			::	fourdir(:)
	integer,private					::	layerspread, layerspreadmin, layerspreadmax
	integer,allocatable,private			::	atomnumberarray(:)
	integer,allocatable,private			::	layerintarr(:)
contains

subroutine getlocalisationpar(return_localisationpar)
	real(kind=8),allocatable		::	return_localisationpar(:)
	
	allocate(return_localisationpar(size(localisationpar)))
	return_localisationpar=localisationpar
end subroutine getlocalisationpar

subroutine getH_layer(hamiltonian,sigma,k,MPI_COMM_WORLD)
	real(kind=8),intent(in)				::	k(:)
	integer,intent(in)				::	MPI_COMM_WORLD

	complex(kind=8),allocatable			::	fourHamilton(:,:,:)		!difference in layer, WF1, WF2
	complex(kind=8),allocatable			::	fourRspauli(:,:,:,:)		!difference in layer, WF1, WF2, spindir
	real(kind=8)					::	phase, phase2
	integer						::	fourdirection
	complex(kind=8)					::	fac

	integer						::	zvalue
	complex(kind=8), allocatable			::	hamiltonian(:,:),sigma(:,:,:)

	integer						::	ii,i,j,i1,ierr

	if(.not.allocated(hops))then
		write(*,*) 'hops not allocated, Hamiltonian for layer not set (setH_layer), CALLING MPI ABORT'
		call MPI_Abort(MPI_COMM_WORLD, '1', ierr)
		return
	endif

!setup of fouriertransformed part of Hamilontian
	allocate(fourHamilton(layerspreadmin:layerspreadmax,num_wann,num_wann))
	allocate(fourRspauli(layerspreadmin:layerspreadmax,num_wann,num_wann,3))
	fourHamilton=0d0
	fourRspauli=0d0
!	write(*,*) 'k1,2',k(1),k(2)
		do ii=1,rvecnum
			do i=1,num_wann
				do j=1,num_wann
					phase=0d0
					phase2=0d0
					do i1=1,fourdim
						fourdirection=fourdir(i1)
						!phase=phase+(-atomarr(atomnumberarray(i))%position(fourdirection)+atomarr(atomnumberarray(j))%position(fourdirection)+irvec(fourdirection,ii))*k(fourdirection)
						phase=phase+(irvec(fourdirection,ii))*k(fourdirection)
						
					enddo

					!phase=(irvec(1,ii))*k(1)+(irvec(2,ii))*k(2)
					fac=cmplx(cos(phase),sin(phase))
					fourHamilton(irvec(layerdir,ii),i,j)=fourHamilton(irvec(layerdir,ii),i,j)+fac*hops(i,j,ii)/nrpts(ii)
					!fourHamilton(irvec(ii,3),i,j)=fourHamilton(irvec(ii,3),i,j)+fac*hops(i,j,ii)
					do i1=1,3
						fourRspauli(irvec(layerdir,ii),i,j,i1)=fourRspauli(irvec(layerdir,ii),i,j,i1)+fac*rspauli(i,j,i1,ii)/nrpts(ii)
					enddo
				enddo!j over all second bands
			enddo!i over all first bands
		enddo! ii over all rvecnum


	if(.not.allocated(hamiltonian).and..not.allocated(sigma))then
		allocate(hamiltonian(Hdim,Hdim))
		allocate(sigma(Hdim,Hdim,3))
	endif

	hamiltonian=cmplx(0d0,0d0)
	sigma=cmplx(0d0,0d0)
	do i=1,Hdim
		do j=1,Hdim
			zvalue= layerintarr(i)-layerintarr(j)!richtigrum???
			!write(*,*) 'zvalue', i,j,zvalue
			if(zvalue.ge.layerspreadmin.and.zvalue.le.layerspreadmax)then
				hamiltonian(i,j)=hamiltonian(i,j)+fourHamilton(zvalue,wannierfunctioninham(i),wannierfunctioninham(j))
				do i1=1,3
					sigma(i,j,i1)=sigma(i,j,i1)+fourRspauli(zvalue,wannierfunctioninham(i),wannierfunctioninham(j),i1)
				enddo
			endif
		enddo !j 
	enddo !i

end subroutine getH_layer

!needs to give back final Hamiltonian dimension -> exclude schon hier
subroutine setH_layer(B,ndim,return_num_wann,MPI_COMM_WORLD,my_id,master)
	real(kind=8),intent(in)		::	B(3)
	integer,intent(inout)		::	ndim
	integer,intent(in)		::	MPI_COMM_WORLD,my_id,master
	integer, intent(out)		::	return_num_wann

	integer				::	ierr, i,j,ii

	integer					::	num_lines,ix,iy,iz,band1, band2,dir
	real(kind=8)				::	rdum,idum
	logical					::	fileexist


	integer					::	numberlayer,ndiffatom,nexcludeup,nexcludedown,nwannexup,nwannexdown
	integer,allocatable			::	wannperat(:)

	integer, allocatable			::	excludeup(:),excludedown(:),temp(:)



	!get numwann


	if(my_id.eq.master)	write(*,*) 'start setH_layer'

	inquire(file='layer_inp',exist=fileexist)
	if(.not.fileexist)then
		write(*,*)	'layer_inp does not exists, calling MPI_ABORT'
		call MPI_Abort(MPI_COMM_WORLD, '1', ierr)
		return
	endif

	open(13,file='layer_inp')
	read(13,*) numberlayer

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


!read in hopping.1
         rvecnum=num_lines/(num_wann*num_wann)
	 if(my_id.eq.master)	 write(*,*)"num_lines=",num_lines
         if(my_id.eq.master)	 write(*,*)"rvecnum=",rvecnum
         if(my_id.eq.master)	 write(*,*)"num_wann" , num_wann
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
            hops(band1,band2,rvecnum)=cmplx(rdum,idum)
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


!add exchange field
	do ii=1,rvecnum
		do j=1,num_wann
			do i=1,num_wann
				hops(i,j,ii)=hops(i,j,ii)+0.5d0*(rspauli(i,j,1,ii)*B(1)+rspauli(i,j,2,ii)*B(2)+rspauli(i,j,3,ii)*B(3))
			enddo
		enddo
	enddo

!setup of wannierfunctions and atomnumberarray used in final Hamiltonian
	allocate(atomnumberarray(num_wann))
	read(13,*)	ndiffatom
	allocate(atomarr(ndiffatom))
	atomnumberarray=0
	read(13,'(<num_wann>I3)') atomnumberarray(:)

	
	allocate(wannperat(ndiffatom))
	wannperat=0
	do i=1,num_wann
		wannperat(atomnumberarray(i))=wannperat(atomnumberarray(i))+1
	enddo
	do i=1,ndiffatom
		atomarr(i)%number=i
		allocate(atomarr(i)%wannierfunctions(wannperat(i)))
		atomarr(i)%wannierfunctions(wannperat(i))=0
	enddo
	allocate(temp(ndiffatom))
	temp=1
	do i=1,num_wann
		atomarr(atomnumberarray(i))%wannierfunctions(temp(atomnumberarray(i)))=i
		temp(atomnumberarray(i))=temp(atomnumberarray(i))+1
	enddo
	do i=1,ndiffatom
		if(ANY(atomarr(i)%wannierfunctions(:).eq.0).or.temp(i)-1.ne.wannperat(i))then
			write(*,*)	'mismatch assigning wannier functions to atoms, CALLING MPI ABORT'
			call MPI_Abort(MPI_COMM_WORLD, '1', ierr)
		endif
	enddo
	deallocate(temp)
	deallocate(wannperat)
	do i=1,ndiffatom
		read(13,*) atomarr(i)%position(:)
	enddo
	read(13,'(2I3)')	nexcludeup,nexcludedown
	if(nexcludeup.eq.0)then
		read(13,*)
		allocate(excludeup(0))
	else
		allocate(excludeup(nexcludeup))
		read(13,'(<nexcludeup>I3)')	excludeup(:)
	endif
	if(nexcludedown.eq.0)then
		read(13,*)
		allocate(excludedown(0))
	else
		allocate(excludedown(nexcludedown))
		read(13,'(<nexcludedown>I3)')	excludedown(:)
	endif
	read(13,'(I3)')	fourdim
	allocate(fourdir(fourdim))
	read(13,'(I3,<fourdim>I3)')	layerdir,fourdir(:) 
	close(13)
	
	return_num_wann=num_wann*numberlayer
	nwannexup=0
	do i=1,nexcludeup
		nwannexup=nwannexup+size(atomarr(excludeup(i))%wannierfunctions(:))
	enddo
	nwannexdown=0
	do i=1,nexcludedown
		nwannexdown=nwannexdown+size(atomarr(excludedown(i))%wannierfunctions(:))
	enddo
	return_num_wann=return_num_wann-nwannexup-nwannexdown
	Hdim=return_num_wann

	allocate(wannierfunctioninHam(Hdim))
	allocate(localisationpar(Hdim))
	Call connectWannFctWithHamil(nwannexdown,nwannexup,ndiffatom,num_wann,excludedown,excludeup,layerdir,numberlayer,wannierfunctioninHam)	!fills wannierfunctioninham with number of wannierfunction for each energy basisvector

!construct fourier transformed parts of Hamiltonian	
	layerspreadmin=0
	layerspreadmax=0
	do i=1,rvecnum
		if(irvec(layerdir,i).gt.layerspreadmax)then
			layerspreadmax=irvec(layerdir,i)
		endif
		if(irvec(layerdir,i).lt.layerspreadmin)then
			layerspreadmin=irvec(layerdir,i)
		endif
	enddo
	layerspread=layerspreadmax-layerspreadmin

	if(my_id.eq.master)	write(*,*)  'use ndim',ndim,'to get real ndim'

	if(my_id.eq.master)	write(*,*)	'setH_layer DONE'
end subroutine setH_layer


!fills wannierfunctioninHam with wannierfunction at dimension of the Hamiltonian and fills localisationpar with localisationsparameter from -1 to 1
subroutine connectWannFctWithHamil(nwannexdown,nwannexup,ndiffatom,num_wann,excludedown,excludeup,layerdir,numberlayer,wannierfunctioninHam)
	integer					::	i,j,i1,i2,numberlayer
	integer, intent(in)			::	nwannexup,nwannexdown,ndiffatom,num_wann,layerdir
	integer,intent(in),allocatable		::	excludedown(:),excludeup(:)
	integer,intent(inout),allocatable	::	wannierfunctioninHam(:)	
	real(kind=8)				::	locmax,locmin
	j=1
	i=1
	allocate(layerintarr(Hdim))
	do while(i.le.Hdim)
		if(i.le.num_wann-nwannexup)then
			i1=1
			i2=1
			do while (i.le.num_wann-nwannexup)
				if(.not.ANY(excludeup==i1))then
					wannierfunctioninHam(i)=atomarr(i1)%wannierfunctions(i2)
					localisationpar(i)=atomarr(i1)%position(layerdir)+numberlayer-1		!+numberlayer such that up is positive
					layerintarr(i)=1 

					if(i2.lt.size(atomarr(i1)%wannierfunctions(:)))then
						i2=i2+1
					else
						i2=1
						i1=i1+1
					endif
					i=i+1
				else
					i1=i1+1
				endif
			enddo
		elseif(i.le.Hdim+nwannexdown-num_wann)then
			do i1=1,ndiffatom
				do i2=1,size(atomarr(i1)%wannierfunctions(:))
					wannierfunctioninHam(i)=atomarr(i1)%wannierfunctions(i2)
					localisationpar(i)=atomarr(i1)%position(layerdir)+numberlayer-(i+nwannexup-1)/num_wann-1
					layerintarr(i)=(i+nwannexup-1)/num_wann+1
					i=i+1
				enddo
			enddo
		else
			i1=1
			i2=1
			do while (i.le.Hdim)
				if(.not.ANY(excludedown==i1))then
					wannierfunctioninHam(i)=atomarr(i1)%wannierfunctions(i2)
					localisationpar(i)=atomarr(i1)%position(layerdir)
					layerintarr(i)=numberlayer
					if(i2.lt.size(atomarr(i1)%wannierfunctions(:)))then
						i2=i2+1
					else
						i2=1
						i1=i1+1
					endif
					i=i+1
				else
					i1=i1+1
				endif
			enddo

		endif
	!	wannierfunctioninHam(i)
	enddo
	locmin=MINVAL(localisationpar,Hdim)
	locmax=MAXVAL(localisationpar,Hdim)
!	write(*,*)	locmin,locmax
	do i=1,Hdim
		localisationpar(i)=-1d0+2d0*(localisationpar(i)-locmin)/(locmax-locmin)
	enddo
end subroutine
end module system_layer

