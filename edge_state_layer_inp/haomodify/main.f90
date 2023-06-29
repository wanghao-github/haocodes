program main
use system
use kmesh
use processH
use system_layer
implicit none
!mpiifort system_layer.f90 kmesh.f90 system.f90 processH.f90 main.f90 -L /tmp_mnt/local/intel/mkl/current/lib/em64t -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -o testHamiltonian.out
!mpiifort system_layer.f90 kmesh.f90 system.f90 processH.f90 main.f90 -L /tmp_mnt/local/intel/mkl/current/lib/em64t -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -g -check all -traceback -warn all -o testHamiltonian.out


!	11	plane_inp
!	12	highsym_inp
!	13	layer_inp
!	20	highsympoint.out
!	21	eigenval.out
!	22	kpts.out
!	23	Berrycurvchernno.out
!	24	time.out
!	25	4pointchernno.out
!	200	hopping.1
!	201	rspauli.1


!inputparameter
integer					::	ndim, Hdim	!ndim: number occupied states, Hdim: dimension of Hamiltonian
real(kind=8)				::	Bfield(3),bravv(3,3)
logical					::	calcferm
integer					::	cutofup,cutofdown	!highest and lowest included energy state in calculation of fermi energy

integer					::	Nproj
real(kind=8)				::	projdir(3)

!kpointmesh
character(len=5) 			::	ktype,Htype			!type of kpointset for input (plane,highs)
type(kpoint), pointer			::	kpoints(:)		!array of all kpoints
real(kind=8),allocatable		::	localkpts(:,:)
integer					::	nkpt, nlockpts		!nkpt: number of k-points, nlockpts: number of local kpts of processor
integer,allocatable			::	kdistribution(:)
real(kind=8)				::	volume

!time
real 					::	time(10)


!MPIvariables
integer					::	ierr, my_id, num_procs, master, passmpiworld
!integer					::	req1, req2, req3, req4, req5



complex(kind=8),allocatable		::	Hamiltonian(:,:), Spin(:,:,:)
real(kind=8),allocatable		::	localeigenvalues(:,:)
complex(kind=8),allocatable		::	localeigenvectors(:,:,:),localprojeigenvect(:,:,:)
real(kind=8),allocatable		::	localspinproj(:,:,:)
real(kind=8) ,allocatable		::	locallocalisations(:,:)
real(kind=8) ,allocatable		::	localBerrycurv(:,:,:)
real(kind=8)				::	fermienergy

real(kind=8)				::	spinchern

integer					::	percent,oldpercent

!logicals
logical					::	fileexist, success
!logical,allocatable			::	outputlogical(:)	!eigenval,eigenvect,spinproj,localisation
logical					::	percentoutput
logical					::	outputlogical(6)
logical					::	logcalc(5)	

!loop variables
integer					::	i,j

real(kind=8),allocatable		::	Chernno(:),chernno4p(:),projchernno4p(:)
!MPI start
Include 'mpif.h'
integer status(MPI_STATUS_SIZE)
Call MPI_INIT(ierr)
call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)
master=0
passmpiworld=MPI_COMM_WORLD


oldpercent=0
time=0d0

success=.true.

outputlogical=.false.
logcalc=.false.
percentoutput=.true.
outputlogical(1)=.true.
outputlogical(3)=.true.


!****************************************************************************
!input
!****************************************************************************
inquire(file='Chern_inp',exist=fileexist)
if(.not.fileexist.and.my_id.eq.0)then
	write(*,*)	'Chern_inp does not exists, calling MPI_ABORT'
	call MPI_Abort(MPI_COMM_WORLD, '1', ierr)
endif
Open(10,file='Chern_inp')
read(10,'(I3)')	ndim
read(10,*)  Bfield(1),  Bfield(2),  Bfield(3)
read(10,*)  bravv(1,1),bravv(1,2),bravv(1,3)
read(10,*)  bravv(2,1),bravv(2,2),bravv(2,3)
read(10,*)  bravv(3,1),bravv(3,2),bravv(3,3)
read(10,'(1X,A5)')  ktype
read(10,'(1X,A5)')  Htype
read(10,*)	calcferm
read(10,*)	cutofdown,cutofup
read(10,*)	logcalc(1)	!oldchernno
read(10,*)	logcalc(2)	!4pointchernno
read(10,*)	logcalc(3)	!4pointchernno
read(10,*)	logcalc(4),Nproj,projdir(1),projdir(2),projdir(3)	!surface projected bandstructure
close(10)

if(logcalc(2))	outputlogical(2)=.true.
if(logcalc(1))	outputlogical(5)=.true.
if(logcalc(3))	outputlogical(6)=.true.

if (my_id.eq.master)then
	write(*,*) 'number of occupied states: ', ndim
	write(*,'(A23,3(f10.6))') ' applied exchangefield: ', Bfield
	write(*,*) 'lattice vectors: '
	write(*,'(3f10.5)') bravv(1,1),bravv(1,2),bravv(1,3)
	write(*,'(3f10.5)') bravv(2,1),bravv(2,2),bravv(2,3)
	write(*,'(3f10.5)') bravv(3,1),bravv(3,2),bravv(3,3)
	if(calcferm)	write(*,*) 'Calculating Fermi energy considering energy states', cutofdown, 'to', cutofup
	write(*,*)
endif
if (my_id.eq.master)	CALL CPU_TIME(time(1))
!****************************************************************************************************
!reading in of realspace Hamiltonian and Spin and setting dimension of Hamiltonian in variable Hdim
!****************************************************************************************************
select case (Htype)
	case('hoprs')
		call setH_hoppingrspauli(Bfield,Hdim,passmpiworld,my_id,master)
		if (my_id.eq.master)	write(*,*)	'setH from hopping.1 and rspauli.1 done'
	case('layer')
		call setH_layer(Bfield,ndim,Hdim,passmpiworld,my_id,master)
		if (my_id.eq.master)	write(*,*)	'setH of layer constructed from hopping.1 and rspauli.1 DONE'
	case default
		write(*,*)	'no Hamiltonian selection found, calling MPI_ABORT'
		call MPI_Abort(MPI_COMM_WORLD, '1', ierr)
	endselect


if (my_id.eq.master)	write(*,*) 'dimension of Hamiltonian: ', Hdim
if (my_id.eq.master)	CALL CPU_TIME(time(2))

!****************************************************************************************************
!setup and distribution to different processors of the kpoints
!****************************************************************************************************
if(my_id.eq.master)then
	select case (ktype)
	case('plane')
		Call setkplane(bravv,volume,passmpiworld)
		Call getkmeshnumber(nkpt)
	case('cuboi')
		Call setkCuboid(passmpiworld)
		Call getkmeshnumber(nkpt)
	case('highs')
		Call setkhighsym(bravv,passmpiworld)
		write(*,*) 'setkhighsym DONE'
		Call getkmeshnumber(nkpt)
	case default
		write(*,*)	'no k-vector selection found, calling MPI_ABORT'
		call MPI_Abort(MPI_COMM_WORLD, '1', ierr)
	endselect
	write(*,*) 'number of kpoints:',nkpt
	if(logcalc(4))then
		 Call projectedbandstruct(projdir,Nproj)
		 Call getkmeshnumber(nkpt)
	endif
	
	Call getkmesh(kpoints,passmpiworld)
!distribute kpoints
	do i=1,num_procs-1
		Call distributek(num_procs,i,kdistribution,localkpts,passmpiworld)
		nlockpts=size(localkpts(:,1))
		call mpi_send(nlockpts,1, MPI_INTEGER , i, 1, MPI_COMM_WORLD, ierr)
		call mpi_send(localkpts,size(localkpts), MPI_real8, i, 2, MPI_COMM_WORLD, ierr)
	enddo
	Call distributek(num_procs,0,kdistribution,localkpts,passmpiworld)
	nlockpts=size(localkpts(:,1))
	write(*,*) 'distribution of kpoints to processors DONE'
	if (my_id.eq.master)	CALL CPU_TIME(time(3))
else  !my_id.eq.master
	call mpi_recv(nlockpts,1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, status,ierr)
	allocate(localkpts(nlockpts,3))
	call mpi_recv(localkpts,size(localkpts), MPI_real8, 0, 2, MPI_COMM_WORLD,status, ierr)
endif !my_id.eq.master

!****************************************************************************************************
!calculation block
!****************************************************************************************************
allocate(localeigenvalues(nlockpts,Hdim))
allocate(localeigenvectors(nlockpts,Hdim,Hdim))
allocate(localspinproj(nlockpts,Hdim,3))
if(outputlogical(6))	allocate(localprojeigenvect(nlockpts,Hdim,ndim))
if(outputlogical(5))	allocate(localBerrycurv(nlockpts,Hdim,3))
select case (Htype)
	case('hoprs')
		Call setprocessH(Hdim)

		do i=1,nlockpts
			Call getH(Hamiltonian,localkpts(i,:))
			Call getspin(Spin,localkpts(i,:))
			
			Call processHamiltonian(Hamiltonian,Hdim,localeigenvalues(i,:),localeigenvectors(i,:,:),Spin,localspinproj(i,:,:),MPI_COMM_WORLD)
			
			if(outputlogical(5))	Call getBerrycurv(Hdim,localeigenvalues(i,:),localeigenvectors(i,:,:),bravv,localkpts(i,:),localBerrycurv(i,:,:))
			if(outputlogical(6))	Call getprojeigenvector(ndim,Hdim,localprojeigenvect(i,:,:),localeigenvectors(i,:,:),localkpts(i,:),MPI_COMM_WORLD)
			if(percentoutput)then	!prozentanzeige
				percent=(i*100)/nlockpts
				if(percent.gt.oldpercent)then
					write(*,*) 'process', my_id, ' at ',percent, 'percent of Hamiltonian diagonalisation'
					oldpercent=percent
				endif
			endif
			
		enddo
		Call unsetprocessH()
	case('layer')
		Call setprocessH(Hdim)
		outputlogical(4)=.true.	!calculate localisation
		allocate(locallocalisations(nlockpts,Hdim))
		do i=1,nlockpts
			Call getH_layer(Hamiltonian,Spin,localkpts(i,:),MPI_COMM_WORLD)
			Call processHamiltonian_loc(Hamiltonian,Hdim,localeigenvalues(i,:),localeigenvectors(i,:,:),Spin,localspinproj(i,:,:),locallocalisations(i,:),MPI_COMM_WORLD)
			if(percentoutput)then	!prozentanzeige
				percent=(i*100)/nlockpts
				if(percent.gt.oldpercent)then
					write(*,*) 'process', my_id, ' at ',percent, 'percent of Hamiltonian diagonalisation'
					oldpercent=percent
				endif
			endif

		enddo
		Call unsetprocessH()
	case default
		write(*,*)	'no Hamiltonian selection found, calling MPI_ABORT'
		call MPI_Abort(MPI_COMM_WORLD, '1', ierr)
	endselect
if (my_id.eq.master)	CALL CPU_TIME(time(4))
call mpi_barrier(MPI_COMM_WORLD,ierr)
if (my_id.eq.master)	CALL CPU_TIME(time(9))
!****************************************************************************************************
!gather all scattered data from different processors
!****************************************************************************************************

call mpi_barrier(MPI_COMM_WORLD,ierr)
write(*,*) my_id, 'start gathering'
if(my_id.eq.master)then 
	nlockpts=size(localkpts(:,1))
	!allocate everything seperately -.- some better way?
	if(outputlogical(1))then
		do i=1,nkpt	
			allocate(kpoints(i)%eigenvalues(Hdim))
		enddo
		do j=1,nlockpts
			kpoints(kdistribution(j))%eigenvalues(:)=localeigenvalues(j,:)
		enddo
	endif
	if(outputlogical(2))then
		do i=1,nkpt		
			allocate(kpoints(i)%eigenvectors(Hdim,Hdim))
		enddo
		do j=1,nlockpts
			kpoints(kdistribution(j))%eigenvectors(:,:)=localeigenvectors(j,:,:)
		enddo
	endif
	if(outputlogical(3))then
		do i=1,nkpt		
			allocate(kpoints(i)%spinproj(Hdim,3))
		enddo
		do j=1,nlockpts
			kpoints(kdistribution(j))%spinproj(:,:)=localspinproj(j,:,:)
		enddo
	endif
	if(outputlogical(4))then
		do i=1,nkpt		
			allocate(kpoints(i)%localisation(Hdim))
		enddo
		do j=1,nlockpts
			kpoints(kdistribution(j))%localisation(:)=locallocalisations(j,:)
		enddo
	endif
	if(outputlogical(5))then
		do i=1,nkpt		
			allocate(kpoints(i)%Berrycurv(Hdim,3))
		enddo
		do j=1,nlockpts
			kpoints(kdistribution(j))%Berrycurv(:,:)=localBerrycurv(j,:,:)
		enddo
	endif
	if(outputlogical(6))then
		do i=1,nkpt		
			allocate(kpoints(i)%projeigenvectors(Hdim,ndim))
		enddo
		do j=1,nlockpts
			kpoints(kdistribution(j))%projeigenvectors(:,:)=localprojeigenvect(j,:,:)
		enddo
	endif
		
	write(*,*) 'done allocate all'
	do i=1,num_procs-1
		Call distributek(num_procs,i,kdistribution,localkpts,passmpiworld)
		nlockpts=size(kdistribution)

		if(outputlogical(1))then		
			deallocate(localeigenvalues)
			allocate(localeigenvalues(nlockpts,Hdim))
			call mpi_recv(localeigenvalues,size(localeigenvalues), MPI_real8, i, 3, MPI_COMM_WORLD,status, ierr)
			do j=1,nlockpts
				kpoints(kdistribution(j))%eigenvalues(:)=localeigenvalues(j,:)
			enddo
		endif
		if(outputlogical(2))then		
			deallocate(localeigenvectors)
			allocate(localeigenvectors(nlockpts,Hdim,Hdim))
			call mpi_recv(localeigenvectors,size(localeigenvectors), MPI_DOUBLE_COMPLEX, i, 4, MPI_COMM_WORLD,status,ierr)
			do j=1,nlockpts
				kpoints(kdistribution(j))%eigenvectors(:,:)=localeigenvectors(j,:,:)
			enddo
		endif
		if(outputlogical(3))then
			deallocate(localspinproj)
			allocate(localspinproj(nlockpts,Hdim,3))
			call mpi_recv(localspinproj,size(localspinproj), MPI_real8, i, 5, MPI_COMM_WORLD,status,ierr)
			do j=1,nlockpts

				kpoints(kdistribution(j))%spinproj(:,:)=localspinproj(j,:,:)
			enddo
		endif
		if(outputlogical(4))then		
			deallocate(locallocalisations)
			allocate(locallocalisations(nlockpts,Hdim))
			call mpi_recv(locallocalisations,size(locallocalisations), MPI_real8, i, 6, MPI_COMM_WORLD,status,ierr)
			do j=1,nlockpts
				kpoints(kdistribution(j))%localisation(:)=locallocalisations(j,:)
			enddo
		endif
		if(outputlogical(5))then		
			deallocate(localBerrycurv)
			allocate(localBerrycurv(nlockpts,Hdim,3))
			call mpi_recv(localBerrycurv,size(localBerrycurv), MPI_real8, i, 7, MPI_COMM_WORLD,status,ierr)
			do j=1,nlockpts
				kpoints(kdistribution(j))%Berrycurv(:,:)=localBerrycurv(j,:,:)
			enddo
		endif
		if(outputlogical(6))then		
			deallocate(localprojeigenvect)
			allocate(localprojeigenvect(nlockpts,Hdim,ndim))
			call mpi_recv(localprojeigenvect,size(localprojeigenvect), MPI_DOUBLE_COMPLEX, i, 8, MPI_COMM_WORLD,status,ierr)
			do j=1,nlockpts
				kpoints(kdistribution(j))%projeigenvectors(:,:)=localprojeigenvect(j,:,:)
			enddo
		endif
		write(*,*) i,'th receive done'
	enddo
else
	if(outputlogical(1))	call mpi_send(localeigenvalues,nlockpts*Hdim, MPI_real8, 0, 3, MPI_COMM_WORLD,ierr)
	if(outputlogical(2))	call mpi_send(localeigenvectors,nlockpts*Hdim*Hdim, MPI_DOUBLE_COMPLEX, 0, 4, MPI_COMM_WORLD, ierr)
	if(outputlogical(3))	call mpi_send(localspinproj,nlockpts*Hdim*3, MPI_real8, 0, 5, MPI_COMM_WORLD, ierr)
	if(outputlogical(4))	call mpi_send(locallocalisations,nlockpts*Hdim, MPI_real8, 0, 6, MPI_COMM_WORLD, ierr)
	if(outputlogical(5))	call mpi_send(localBerrycurv,nlockpts*Hdim*3, MPI_real8, 0, 7, MPI_COMM_WORLD, ierr)
	if(outputlogical(6))	call mpi_send(localprojeigenvect,nlockpts*Hdim*ndim, MPI_DOUBLE_COMPLEX, 0, 8, MPI_COMM_WORLD, ierr)
endif


!if(num_procs.gt.1.and.outputlogical(1))		Call MPI_WAIT(req1, status, ierr)
!if(num_procs.gt.1.and.outputlogical(2))		Call MPI_WAIT(req2, status, ierr)
!if(num_procs.gt.1.and.outputlogical(3))		Call MPI_WAIT(req3, status, ierr)
!if(num_procs.gt.1.and.outputlogical(4))		Call MPI_WAIT(req4, status, ierr)
!if(num_procs.gt.1.and.outputlogical(5))		Call MPI_WAIT(req5, status, ierr)
call mpi_barrier(MPI_COMM_WORLD,ierr)

if (my_id.eq.master)	CALL CPU_TIME(time(5))

!if(my_id.ne.master)then
if(allocated(localeigenvalues))		deallocate(localeigenvalues)
if(allocated(localeigenvectors))	deallocate(localeigenvectors)
if(allocated(localspinproj))		deallocate(localspinproj)
if(allocated(locallocalisations))	deallocate(locallocalisations)
if(allocated(localBerrycurv))		deallocate(localBerrycurv)
!endif


if(my_id.eq.master)	write(*,*)	'gather all results DONE'



if(my_id.eq.master.and.calcferm)then
	if(Htype.eq.'layer')	write(*,*) 'WARNING,ndim not adjusted for layer calculations to calculated fermi energy'
		!implement with counting how many occupied states per atom (deduct if excluded)
	Call getfermi(Hdim,ndim,cutofup,cutofdown,nkpt,kpoints,fermienergy,passmpiworld)
	if (my_id.eq.master)	CALL CPU_TIME(time(6))
else
	time(6)=time(5)
endif
!****************************************************************************************************
!output
!****************************************************************************************************


if(my_id.eq.master)then
	write(*,*) 'open eigenvals'
	open(21,file='eigenval.out')
	do i=1,nkpt
		do j=1,Hdim
			if(.not.outputlogical(4)) write(21,'(4f9.6,f12.6,3f10.5)') kpoints(i)%k(:), kpoints(i)%distance, kpoints(i)%eigenvalues(j), kpoints(i)%spinproj(j,1), kpoints(i)%spinproj(j,2), kpoints(i)%spinproj(j,3)
			if(outputlogical(4)) write(21,'(f9.6,f12.6,4f10.5)') kpoints(i)%distance,kpoints(i)%eigenvalues(j),kpoints(i)%spinproj(j,1),kpoints(i)%spinproj(j,2),kpoints(i)%spinproj(j,3), kpoints(i)%localisation(j)
		enddo
	enddo
	close(21)
	write(*,*)	'write eigenval.out DONE'
	CALL CPU_TIME(time(7))
	open(22,file='kpts.out')
	do i=1,nkpt
		write(22,'(3f16.10)')	kpoints(i)%k(:)
	enddo
	close(22)
	write(*,*) 'write kpts.out DONE'
endif

if(my_id.eq.master.and.logcalc(1))then
	allocate(Chernno(Hdim))
	Chernno=0d0

	do i=1,nkpt
		do j=1,Hdim
			Chernno(j)=Chernno(j)+kpoints(i)%Berrycurv(j,3)
		enddo
	enddo
	open(23,file='Berrycurvchernno.out')
	do i=1,Hdim
		write(23,*)	i, Chernno(i)/volume/nkpt*2.d0*pi
	enddo
	close(23)
	write(*,*)	'write Berrycurvchernno.out DONE'
endif

if(my_id.eq.master.and.logcalc(2))then
	allocate(chernno4p(Hdim))
	Call	get4chernno(kpoints,Hdim,nkpt,chernno4p,MPI_COMM_WORLD)
	open(25,file='4pointchernno.out')
	do i=1,Hdim
		write(25,'(I6,f20.14)') i,chernno4p(i)
	enddo
	close(25)
	write(*,*)	'write 4pointchernno.out DONE'
endif

if(my_id.eq.master.and.logcalc(3))then
	allocate(projchernno4p(Hdim))
	Call	getproj4chernno(kpoints,ndim,nkpt,projchernno4p,spinchern,MPI_COMM_WORLD)
	open(26,file='proj4pointchernno.out')
	do i=1,ndim
		write(26,'(I6,f20.14)') i,projchernno4p(i)
	enddo
	close(26)
	write(*,*) 'spinchernnumber is ', spinchern
	write(*,*)	'write proj4pointchernno.out DONE'
endif
if(my_id.eq.master)then

	CALL CPU_TIME(time(8))

	Open(24,file='time.out')
	write(24,*)	'read first input', time(1)
	write(24,*)	'Hamiltonian set', time(2)
	write(24,*)	'kpoints distributed', time(3)
	write(24,*)	'local diagonalisation done', time(4)
	write(24,*)	'all diagonalisations done', time(9)
	write(24,*)	'gathering done', time(5)
	write(24,*)	'getfermienergy done', time(6)
	write(24,*)	'write eigenvalues done', time(7)
	write(24,*)	'finished program', time(8)
endif
	write(*,*) my_id,'done'
call mpi_barrier(MPI_COMM_WORLD,ierr)
call mpi_finalize(ierr)

end program

