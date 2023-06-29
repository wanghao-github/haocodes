module kmesh
implicit none
	public 		:: setkplane,getkmeshnumber,getkmesh,distributek
	private		:: neighpos_plane

!	integer		:: neighpos   !function for neighbours 
	real(kind=8),parameter				:: pi=3.141592653589793238462643383279502884197169399375105820974944d0


	type mesharr
		type(kpoint), pointer :: next
		type(kpoint), pointer :: last
	end type mesharr


	type 							:: kpoint
		integer						:: id		!position in array
		real(kind=8)					:: k(3)		!kpoints
		real(kind=8)					:: distance	!distance between line of kpoints to plot band structure 
		TYPE(mesharr), dimension(:), allocatable	:: meshpointer	!pointer to next kpoints in mesh in direction 1 and 2
		real(kind=8),allocatable			:: eigenvalues(:),Berrycurv(:,:), spinproj(:,:),localisation(:)
		complex(kind=8),allocatable			:: eigenvectors(:,:),projeigenvectors(:,:)

	end type kpoint
	type(kpoint),target, allocatable,private		::	karr(:)		!array containing kpointmesh
	integer, private					::	nkpt,meshdim			!number of entries in kmesh, dimension of kmesh(1 for highs and 2 for plane
	integer,private,allocatable				::	nkptdir(:)

contains

subroutine projectedbandstruct(projdir,Nproj)
	real(kind=8),intent(in)				::	projdir(:)
	integer,intent(in)				::	Nproj

	type(kpoint),target, allocatable		::	karrtemp(:)
	integer						::	i,j

	allocate(karrtemp(nkpt))
	do i=1,nkpt
		karrtemp(i)%id=karr(i)%id
		karrtemp(i)%k=karr(i)%k
		karrtemp(i)%distance=karr(i)%distance
	enddo	
	deallocate(karr)
	
	allocate(karr(nkpt*Nproj))
	do i=1,nkpt
		do j=1,Nproj
			karr((i-1)*Nproj+j)%id=(i-1)*Nproj+j
			karr((i-1)*Nproj+j)%distance=karrtemp(i)%distance
			karr((i-1)*Nproj+j)%k(1)=karrtemp(i)%k(1)+(j-1)*projdir(1)/Nproj*pi
			karr((i-1)*Nproj+j)%k(2)=karrtemp(i)%k(2)+(j-1)*projdir(2)/Nproj*pi
			karr((i-1)*Nproj+j)%k(3)=karrtemp(i)%k(3)+(j-1)*projdir(3)/Nproj*pi
		enddo
	enddo

	nkpt=nkpt*Nproj
end subroutine

!subroutine to set up 'karr' with kpointmesh in case of 'ktype'='plane' with input from 'plane_inp'
subroutine setkCuboid(MPI_COMM_WORLD)
	logical				::	fileexist
	integer, intent(in)		::	MPI_COMM_WORLD
	integer				::	ierr
	real(kind=8)			::	kbase(3), k1(3), k2(3), k3(3),kbasediv, k1div, k2div,k3div
	integer,allocatable		::	checkdistarr(:)

	integer				:: 	i1,i2,i3, karrrun,j
	
	write(*,*) 'starting setkCuboid'
	meshdim=3
	allocate(nkptdir(meshdim))
	inquire(file='cuboid_inp',exist=fileexist)
	if(.not.fileexist)then
		write(*,*)	'cuboid_inp does not exists, calling MPI_ABORT'
		call MPI_Abort(MPI_COMM_WORLD, '1', ierr)
		return
	endif
	Open(11,file='cuboid_inp')
	read(11,*)  nkptdir(1),nkptdir(2),nkptdir(3)
	read(11,*)  kbase(:), kbasediv
	read(11,*)  k1(:), k1div
	read(11,*)  k2(:), k2div
	read(11,*)  k3(:), k3div
	close(11)
	
	nkpt=nkptdir(1)*nkptdir(2)*nkptdir(3)
	kbase(:)=kbase(:)/kbasediv*pi
	k1(:)=k1(:)/k1div*pi
	k2(:)=k2(:)/k2div*pi
	k3(:)=k3(:)/k3div*pi
	write(*,*) 'k1',k1(:)
	write(*,*) 'k2',k2(:)
	write(*,*) 'k3',k3(:)
	
	write(*,*) 'number of kpoints in direction 1, 2 and 3: ', nkptdir(1), nkptdir(2), nkptdir(3)
	!write(*,'(A30)') 'basepoint of kplane rectangle :'
	!write(*,'(3F10.6)') kbase
	!write(*,'(A39)') 'k1 direction spanning the k-rectangle :'
	!write(*,'(3F10.6)') k1
	!write(*,'(A39)') 'k2 direction spanning the k-rectangle :'
	!write(*,'(3F10.6)') k2
	
	allocate(karr(nkpt))
	Call allocatenextlast(MPI_COMM_WORLD)

	karrrun=1
	do i3=1,nkptdir(3)
		do i2=1,nkptdir(2)
			do i1=1,nkptdir(1)
				karr(karrrun)%id=karrrun
				karr(karrrun)%k=(/(kbase(j)+(i1-1)*k1(j)/nkptdir(1)+(i2-1)*k2(j)/nkptdir(2)+(i3-1)*k3(j)/nkptdir(3),j=1,3)/)
				karr(karrrun)%meshpointer(1)%next=>karr(neighpos_plane(nkptdir,1,karrrun,MPI_COMM_WORLD))
				karr(karrrun)%meshpointer(1)%last=>karr(neighpos_plane(nkptdir,2,karrrun,MPI_COMM_WORLD))
				karr(karrrun)%meshpointer(2)%next=>karr(neighpos_plane(nkptdir,3,karrrun,MPI_COMM_WORLD))
				karr(karrrun)%meshpointer(2)%last=>karr(neighpos_plane(nkptdir,4,karrrun,MPI_COMM_WORLD))
				karr(karrrun)%meshpointer(3)%next=>karr(neighpos_plane(nkptdir,5,karrrun,MPI_COMM_WORLD))
				karr(karrrun)%meshpointer(3)%last=>karr(neighpos_plane(nkptdir,6,karrrun,MPI_COMM_WORLD))
				karrrun=karrrun+1
			enddo! i1
		enddo! i2
	enddo! i3
	write(*,*) 'outmost points spanning the cube'
	write(*,'(3f10.6)')	karr(1)%k(:)
	write(*,'(3f10.6)')	karr(nkptdir(1))%k(:)
	write(*,'(3f10.6)')	karr(1+(nkptdir(2)-1)*nkptdir(1))%k(:)
	write(*,'(3f10.6)')	karr(nkptdir(1)*nkptdir(2))%k(:)
	write(*,'(3f10.6)')	karr(1+nkptdir(1)*nkptdir(2)*(nkptdir(3)-1))%k(:)
	write(*,'(3f10.6)')	karr(nkptdir(1)+nkptdir(1)*nkptdir(2)*(nkptdir(3)-1))%k(:)
	write(*,'(3f10.6)')	karr(1+(nkptdir(2)-1)*nkptdir(1)+nkptdir(1)*nkptdir(2)*(nkptdir(3)-1))%k(:)
	write(*,'(3f10.6)')	karr(nkptdir(1)*nkptdir(2)*nkptdir(3))%k(:)
	write(*,*)


	allocate(checkdistarr(nkpt))
	checkdistarr=0
	do i1=1,nkpt
		checkdistarr(karr(i1)%meshpointer(1)%next%id)=checkdistarr(karr(i1)%meshpointer(1)%next%id)+1
		checkdistarr(karr(i1)%meshpointer(1)%last%id)=checkdistarr(karr(i1)%meshpointer(1)%last%id)+1
		checkdistarr(karr(i1)%meshpointer(2)%next%id)=checkdistarr(karr(i1)%meshpointer(2)%next%id)+1
		checkdistarr(karr(i1)%meshpointer(2)%last%id)=checkdistarr(karr(i1)%meshpointer(2)%last%id)+1
		checkdistarr(karr(i1)%meshpointer(3)%next%id)=checkdistarr(karr(i1)%meshpointer(3)%next%id)+1
		checkdistarr(karr(i1)%meshpointer(3)%last%id)=checkdistarr(karr(i1)%meshpointer(3)%last%id)+1
	enddo
	

	do i1=1,nkpt
		if(checkdistarr(i1).ne.2*meshdim)then
			write(*,*) 'ERROR checking the number of neighbouring kpoint'
			write(*,*) 'calling MPI_ABORT'
			call MPI_Abort(MPI_COMM_WORLD, '1', ierr)
			return
		endif
	enddo


endsubroutine setkCuboid

!subroutine to set up 'karr' with kpointmesh in case of 'ktype'='plane' with input from 'plane_inp'
subroutine setkplane(bravv,volume,MPI_COMM_WORLD)
	logical				::	fileexist
	integer, intent(in)		::	MPI_COMM_WORLD
	integer				::	ierr
	real(kind=8)			::	kbase(3), k1(3), k2(3), kbasediv, k1div, k2div
	real(kind=8)			::	k1ev(3),k2ev(3),temp(3),R1(3),R2(3)
	integer,allocatable		::	checkdistarr(:)
	real(kind=8),intent(out)	::	volume
	real(kind=8),intent(in)		::	bravv(3,3)

	integer				:: 	i1,i2, karrrun,j
	
	meshdim=2
	allocate(nkptdir(meshdim))
	inquire(file='plane_inp',exist=fileexist)
	if(.not.fileexist)then
		write(*,*)	'plane_inp does not exists, calling MPI_ABORT'
		call MPI_Abort(MPI_COMM_WORLD, '1', ierr)
		return
	endif
	Open(11,file='plane_inp')
	read(11,*)  nkptdir(1),nkptdir(2)
	read(11,*)  kbase(:), kbasediv
	read(11,*)  k1(:), k1div
	read(11,*)  k2(:), k2div
	close(11)	



	nkpt=nkptdir(1)*nkptdir(2)
	kbase(:)=kbase(:)/kbasediv*pi
	k1(:)=k1(:)/k1div*pi
	k2(:)=k2(:)/k2div*pi


	volume=0d0
	
	k1ev(:)=k1(:)/sqrt((k1(1))**2+(k1(2))**2+(k1(3))**2)
	k2ev(:)=k2(:)/sqrt((k2(1))**2+(k2(2))**2+(k2(3))**2)
	R1=0d0
	R2=0d0
	do i1=1,3
		R1(:)=R1(:)+k1ev(i1)*bravv(i1,:)
		R2(:)=R2(:)+k2ev(i1)*bravv(i1,:)
	enddo
	Call crossprod(temp,R1,R2)
	volume=sqrt(Dot_product(temp,temp))
	write(*,*) 'area of plane used for Chern number: ', volume

	

	write(*,*) 'number of kpoints in direction 1 and 2: ', nkptdir(1), nkptdir(2)
	!write(*,'(A30)') 'basepoint of kplane rectangle :'
	!write(*,'(3F10.6)') kbase
	!write(*,'(A39)') 'k1 direction spanning the k-rectangle :'
	!write(*,'(3F10.6)') k1
	!write(*,'(A39)') 'k2 direction spanning the k-rectangle :'
	!write(*,'(3F10.6)') k2
	
	allocate(karr(nkpt))
	Call allocatenextlast(MPI_COMM_WORLD)

	karrrun=1
	do i2=1,nkptdir(2)
		do i1=1,nkptdir(1)
			karr(karrrun)%id=karrrun
			karr(karrrun)%k=(/(kbase(j)+(i1-1)*k1(j)/nkptdir(1)+(i2-1)*k2(j)/nkptdir(2),j=1,3)/)
			karr(karrrun)%meshpointer(1)%next=>karr(neighpos_plane(nkptdir,1,karrrun,MPI_COMM_WORLD))
			karr(karrrun)%meshpointer(1)%last=>karr(neighpos_plane(nkptdir,2,karrrun,MPI_COMM_WORLD))
			karr(karrrun)%meshpointer(2)%next=>karr(neighpos_plane(nkptdir,3,karrrun,MPI_COMM_WORLD))
			karr(karrrun)%meshpointer(2)%last=>karr(neighpos_plane(nkptdir,4,karrrun,MPI_COMM_WORLD))
			karrrun=karrrun+1
		enddo! i1
	enddo! i2
	write(*,*) 'outmost points spanning the plane'
	write(*,'(3f10.6)')	karr(1)%k(:)
	write(*,'(3f10.6)')	karr(nkptdir(1))%k(:)
	write(*,'(3f10.6)')	karr((nkptdir(2)-1)*nkptdir(1)+1)%k(:)
	write(*,'(3f10.6)')	karr(nkptdir(1)*nkptdir(2))%k(:)
	write(*,*)

	allocate(checkdistarr(nkpt))
	checkdistarr=0
	do i1=1,nkpt
		checkdistarr(karr(i1)%meshpointer(1)%next%id)=checkdistarr(karr(i1)%meshpointer(1)%next%id)+1
		checkdistarr(karr(i1)%meshpointer(1)%last%id)=checkdistarr(karr(i1)%meshpointer(1)%last%id)+1
		checkdistarr(karr(i1)%meshpointer(2)%next%id)=checkdistarr(karr(i1)%meshpointer(2)%next%id)+1
		checkdistarr(karr(i1)%meshpointer(2)%last%id)=checkdistarr(karr(i1)%meshpointer(2)%last%id)+1
	enddo

	do i1=1,nkpt
		if(checkdistarr(i1).ne.2*meshdim)then
			write(*,*) 'ERROR checking the number of neighbouring kpoint'
			write(*,*) 'calling MPI_ABORT'
			call MPI_Abort(MPI_COMM_WORLD, '1', ierr)
			return
		endif
	enddo

endsubroutine setkplane

!sets kpoints along high symmetry lines defined in highsym_inp
subroutine setkhighsym(bravv,MPI_COMM_WORLD)
	integer, intent(in)		::	MPI_COMM_WORLD
	integer				::	ierr,nsym
	integer, allocatable		::	Npersym(:), divideby(:)
	real(kind=8),intent(in)		::	bravv(3,3)

	real(kind=8),allocatable	::	sympos(:,:)
	integer				::	i,j,ii,jj
	logical				::	fileexist
	real(kind=8)			::	recbravv(3,3),dummy(3),volume


	meshdim=1
	allocate(nkptdir(meshdim))
	!read in
	inquire(file='highsym_inp',exist=fileexist)
	if(.not.fileexist)then
		write(*,*)	'highsym_inp does not exists, calling MPI_ABORT'
		call MPI_Abort(MPI_COMM_WORLD, '1', ierr)
		return
	endif
	Open(12,file='highsym_inp')
	read(12,'(I3)')  nsym
	allocate(Npersym(nsym-1))
	read(12,'(<nsym-1>I6)') (Npersym(i),i=1,nsym-1)
	!if(my_id.eq.master)	write(*,*)	Npersym(:)

	allocate(sympos(nsym,3))
	allocate(divideby(nsym))
	do i=1,nsym
		read(12,'(3f16.8,I6)') sympos(i,:), divideby(i)
		sympos(i,:)=sympos(i,:)/divideby(i)*2*pi
		write(*,*)	'high symmetry position ', i, ' :'
		write(*,*)	sympos(i,:)
	enddo

	nkpt=1
	do i=1,nsym-1
		nkpt=nkpt+Npersym(i)
	enddo
	close(12)

	Call crossprod(dummy,bravv(2,:),bravv(3,:))
	volume=DOT_PRODUCT(bravv(1,:),dummy(:))
	recbravv(1,:)=2d0*pi*dummy(:)/volume
	Call crossprod(dummy,bravv(3,:),bravv(1,:))
	recbravv(2,:)=2d0*pi*dummy(:)/volume
	Call crossprod(dummy,bravv(1,:),bravv(2,:))
	recbravv(3,:)=2d0*pi*dummy(:)/volume
	write(*,*)	'recipriocal lattice vectors'
	write(*,*) recbravv(1,:)
	write(*,*) recbravv(2,:)
	write(*,*) recbravv(3,:)
	
!setting up of kpoints
	allocate(karr(nkpt))
	
	Call allocatenextlast(MPI_COMM_WORLD)
	karr(1)%id=1
	karr(1)%k=sympos(1,:)
	karr(1)%meshpointer(1)%next=>karr(2)
	karr(1)%meshpointer(1)%last=>karr(nkpt)
	karr(1)%distance=0.0d0
	ii=2
	do i=1,nsym-1
		do j=1,Npersym(i)
			karr(ii)%id=ii
			karr(ii)%k=(/(sympos(i,jj)+j*(sympos(i+1,jj)-sympos(i,jj))/Npersym(i),jj=1,3)/)
			dummy= MATMUL((karr(ii)%k-karr(ii-1)%k)/2/pi,recbravv)
			karr(ii)%distance=karr(ii-1)%distance+sqrt(dummy(1)**2+dummy(2)**2+dummy(3)**2)
			karr(ii)%meshpointer(1)%next=>karr(mod(ii,nkpt)+1)
			karr(ii)%meshpointer(1)%last=>karr(ii-1)
			ii=ii+1
		enddo
	enddo
	
	!writing of high symmetry points and their distances used as a x axis
	open(20,file='highsympoint.out')
	ii=1
	write(20,'(4F14.8)')	karr(ii)%distance,karr(ii)%k(:)
	do i=1,nsym-1
		ii=ii+Npersym(i)
		write(20,'(4F14.8)')	karr(ii)%distance,karr(ii)%k(:)
	enddo
	close(20)
	nkptdir(1)=nkpt
endsubroutine setkhighsym

subroutine crossprod(vecout,vec1,vec2)
	real(kind=8),intent(in)		::	vec1(:),vec2(:)
	real(kind=8),intent(out)	::	vecout(:)

	vecout(1)=vec1(2)*vec2(3)-vec1(3)*vec2(2)
	vecout(2)=vec1(3)*vec2(1)-vec1(1)*vec2(3)
	vecout(3)=vec1(1)*vec2(2)-vec1(2)*vec2(1)
end subroutine crossprod

subroutine getkmeshnumber(nkpt_out)
	integer,intent(out)	:: nkpt_out
	nkpt_out=nkpt
end subroutine getkmeshnumber


subroutine getkptdir(meshdim_out,nkptdir_out)
	integer,intent(out)			:: meshdim_out
	integer,intent(out),allocatable		:: nkptdir_out(:)
	meshdim_out=meshdim
	allocate(nkptdir_out(meshdim))
	nkptdir_out(:)=nkptdir(:)

endsubroutine getkptdir

!subroutine to get the kmesh from the main program
subroutine getkmesh(karr_return,MPI_COMM_WORLD)
	integer,intent(in)			:: MPI_COMM_WORLD
	integer					:: ierr
	type(kpoint),pointer, allocatable	::	karr_return(:)


	if(nkpt.gt.0)then
		allocate(karr_return(nkpt))
		karr_return=>karr
	else
		write(*,*) 'no kmesh set, CALL MPI_ABORT'
		call MPI_Abort(MPI_COMM_WORLD, '1', ierr)
		return
	endif
endsubroutine getkmesh

!function used to find next neighbours in the case of a plane
integer function neighpos_plane(nkptdir,dir,pos,MPI_COMM_WORLD)
	integer, intent(in)		:: MPI_COMM_WORLD
	integer				:: ierr
	integer				:: nkpt,i
	integer, intent(in)		:: dir,pos
	integer, intent(in),allocatable	::	nkptdir(:)
	
	nkpt=1
	do i=1,size(nkptdir)
		nkpt=nkpt*nkptdir(i)
	enddo
	select case(dir)
	case(1)
		if(mod(pos-1,nkptdir(1)).ge.nkptdir(1)-1)then
			neighpos_plane=pos-nkptdir(1)+1
		else
			neighpos_plane=pos+1
		endif
	case(2)
		if(mod(pos-1,nkptdir(1)).lt.1)then
			neighpos_plane=pos+nkptdir(1)-1
		else
			neighpos_plane=pos-1
		endif
	case(3)
		if(mod(pos-1,nkptdir(1)*nkptdir(2)).ge.nkptdir(1)*(nkptdir(2)-1))then
			neighpos_plane=pos-nkptdir(1)*(nkptdir(2)-1)
		else
			neighpos_plane=pos+nkptdir(1)
		endif
	case(4)
		if(mod(pos-1,nkptdir(1)*nkptdir(2)).lt.nkptdir(1))then
			neighpos_plane=pos+nkptdir(1)*(nkptdir(2)-1)
		else
			neighpos_plane=pos-nkptdir(1)
		endif
	case(5)
		if(pos.gt.nkptdir(1)*nkptdir(2)*(nkptdir(3)-1))then
			neighpos_plane=pos-nkptdir(1)*nkptdir(2)*(nkptdir(3)-1)
		else
			neighpos_plane=pos+nkptdir(1)*nkptdir(2)
		endif
	case(6)
		if(pos.le.nkptdir(1)*nkptdir(2))then
			neighpos_plane=pos+nkptdir(1)*nkptdir(2)*(nkptdir(3)-1)
		else
			neighpos_plane=pos-nkptdir(1)*nkptdir(2)
		endif
	case default
		neighpos_plane=0
		write(*,*) 'Failed to find neighbouring point in k-plane, wrong dir. CALLING MPI_ABORT'
		call MPI_Abort(MPI_COMM_WORLD, '1', ierr)
	endselect
	if(neighpos_plane.lt.1.or.neighpos_plane.gt.nkpt)then
		write(*,*) 'neighpos_plane',neighpos_plane,'nkpt',nkpt,'dir',dir,'pos',pos
		write(*,*) 'Failed to find neighbouring point in k-plane, outside of karr scope. CALLING MPI_ABORT'
		call MPI_Abort(MPI_COMM_WORLD, '1', ierr)
	endif
end function neighpos_plane


!subroutine used to get the k values of the kpoints and their positions in karr, which are supposed to be distributed to the different processors
subroutine distributek(num_proc, processid,distribution,kpointdistributed,MPI_COMM_WORLD)
	integer, intent(in)		::	num_proc,MPI_COMM_WORLD
	integer,allocatable		::	distribution(:)
	real(kind=8),allocatable	::	kpointdistributed(:,:)
	integer				::	perproc, remainder, processid
	integer				::	i,ierr

	if(processid.ge.num_proc)then
		write(*,*) 'trying to distribute kpoints to process', processid+1,'while only ' ,num_proc, ' processes are started. CALLING MPI_ABORT.'
		call MPI_Abort(MPI_COMM_WORLD, '1', ierr)
	endif
	
	if(allocated(distribution))then
	deallocate(distribution)
	endif
	
	if(allocated(kpointdistributed))then
	deallocate(kpointdistributed)
	endif

	perproc=nkpt/num_proc
	remainder=mod(nkpt,num_proc)
	if(processid.lt.remainder)then
		allocate(distribution(perproc+1))
		allocate(kpointdistributed(perproc+1,3))
		do i=1,perproc+1
			distribution(i)=i+(processid*(perproc+1))
			kpointdistributed(i,:)=karr(distribution(i))%k(:)
		enddo
	else
		allocate(distribution(perproc))
		allocate(kpointdistributed(perproc,3))
		do i=1,perproc
			distribution(i)=i+remainder*(perproc+1)+perproc*(processid-remainder)
			kpointdistributed(i,:)=karr(distribution(i))%k(:)
		enddo
	endif
end subroutine distributek

subroutine allocatenextlast(MPI_COMM_WORLD)
	integer, intent(in)		::	MPI_COMM_WORLD
	integer				::	ierr,i

	if(meshdim.gt.0.and.nkpt.ge.1)then
		do i=1,nkpt
			allocate(karr(i)%meshpointer(meshdim))
		enddo
	else
		write(*,*) 'meshdim .le.0 or nkpt .lt.1, failed to allocate next and last mesh, Calling MPI_ABORT'
		call MPI_Abort(MPI_COMM_WORLD, '1', ierr)
		return
	endif
end subroutine allocatenextlast

end module kmesh

