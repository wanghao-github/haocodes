
!module process Hamiltonian comes after sorting module


!module found in the internet
module qsort_c_module
! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing

! Made F conformant by Walt Brainerd

implicit none
public :: QsortC
private :: Partition

contains

recursive subroutine QsortC(A)
  real(kind=8), intent(in out), dimension(:) :: A
  integer :: iq

  if(size(A) > 1) then
     call Partition(A, iq)
     call QsortC(A(:iq-1))
     call QsortC(A(iq:))
  endif
end subroutine QsortC

subroutine Partition(A, marker)
  real(kind=8), intent(in out), dimension(:) :: A
  integer, intent(out) :: marker
  integer :: i, j
  real(kind=8) :: temp
  real(kind=8) :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition

end module qsort_c_module













module processH
use system_layer, only : getlocalisationpar
use system, only : getmomentum,getSpin
use kmesh
use qsort_c_module
implicit none
public		::	processHamiltonian,getfermi

!zheevx stuff
  real(kind=8),private                  	:: vl, vu
  real(kind=8),parameter,private 		:: abstol=1d-12		!einstellen
  Integer,private 				:: lwork, info 
  Integer,PARAMETER,private 			:: lwmax=180000
  Complex(kind=8),private 			:: work(lwmax)
  Real(kind=8),ALLOCATABLE,private 		:: rwork(:)
  integer,allocatable,private 			:: ifail(:),iwork(:)

real(kind=8) ,allocatable			::	processeigval(:)
complex(kind=8), allocatable			::	processeigvec(:,:)

external	::	zheevx
contains



subroutine processHamiltonian(H,Hdim,eigenvals,eigenvects,sigma,spinproj,MPI_COMM_WORLD)
complex(kind=8),allocatable,intent(in)	::	H(:,:)
integer,intent(in)			::	Hdim

real(kind=8) ,intent(out)		::	eigenvals(:)
real(kind=8) ,intent(out)		::	spinproj(:,:)
complex(kind=8), intent(in)		::	sigma(Hdim,Hdim,3)
complex(kind=8), intent(out)		::	eigenvects(:,:)

integer,intent(in)			::	MPI_COMM_WORLD
integer					::	ne,ierr

integer					::	i,j,ii,jj


	call zheevx('V','A','U',Hdim,H,Hdim,vl,vu,1,Hdim,abstol,ne,processeigval,processeigvec,Hdim,work,lwork,rwork,iwork, ifail,info)
	if(info.ne.0)then
		call zheevx('V','A','U',Hdim,H,Hdim,vl,vu,1,Hdim,abstol,ne,processeigval,processeigvec,Hdim,work,lwork,rwork,iwork, ifail,info)
		write(*,*) 'second try'
		if(info.eq.0)	GOTO 100
		write(*,*) 'zheevx of Hamiltonian diagonalisation info not 0, calling MPI abort'
		call MPI_Abort(MPI_COMM_WORLD, '1', ierr)
		100	write(*,*) 'second try succeded'	
	endif

	do i=1,Hdim
		do j=1,3
			do ii=1,Hdim
				do jj=1,Hdim
					spinproj(i,j)=spinproj(i,j)+conjg(processeigvec(jj,i))*sigma(jj,ii,j)*processeigvec(ii,i)
				enddo
			enddo
		enddo
	enddo

	eigenvals(:)=processeigval(:)
	eigenvects(:,:)=processeigvec(:,:)
end subroutine processHamiltonian

subroutine setprocessH(Hdim)
integer,intent(in)		::	Hdim
lwork=12*Hdim
allocate(rwork(7*Hdim))
allocate(iwork(5*Hdim))
allocate(ifail(Hdim))
allocate(processeigval(Hdim))
allocate(processeigvec(Hdim,Hdim))
end subroutine setprocessH




subroutine unsetprocessH
	deallocate(rwork)
	deallocate(iwork)
	deallocate(ifail)
	deallocate(processeigval)
	deallocate(processeigvec)
end subroutine unsetprocessH





subroutine processHamiltonian_loc(H,Hdim,eigenvals,eigenvects,sigma,spinproj,localisations,MPI_COMM_WORLD)
complex(kind=8),allocatable		::	H(:,:)
integer					::	MPI_COMM_WORLD
integer					::	Hdim,ne,ierr
real(kind=8) ,intent(out)		::	eigenvals(:)
real(kind=8) ,intent(out)		::	localisations(:)
real(kind=8) ,intent(out)		::	spinproj(:,:)
complex(kind=8), intent(in)		::	sigma(Hdim,Hdim,3)
complex(kind=8), intent(out)		::	eigenvects(:,:)
real(kind=8),allocatable		::	localisationpar(:)


integer					::	i,j,ii,jj



	Call getlocalisationpar(localisationpar)
	call zheevx('V','A','U',Hdim,H,Hdim,vl,vu,1,Hdim,abstol,ne,processeigval,processeigvec,Hdim,work,lwork,rwork,iwork, ifail,info)
	if(info.ne.0)then
		write(*,*)	'info: ', info
		write(*,*) 'zheevx of Hamiltonian diagonalisation info not 0, calling MPI abort'
		call MPI_Abort(MPI_COMM_WORLD, '1', ierr)	
	endif
	spinproj=0d0
	localisations=0d0
	do i=1,Hdim
		do ii=1,Hdim
			do jj=1,Hdim	
				do j=1,3
					spinproj(i,j)=spinproj(i,j)+conjg(processeigvec(jj,i))*sigma(jj,ii,j)*processeigvec(ii,i)
				enddo
			enddo
			localisations(i)=localisations(i)+(conjg(processeigvec(ii,i))*processeigvec(ii,i))*localisationpar(ii)	
		enddo
	enddo
	eigenvals(:)=processeigval(:)
	eigenvects(:,:)=processeigvec(:,:)

end subroutine processHamiltonian_loc









subroutine getfermi(Hdim,ndim,cutofup,cutofdown,nkpt,kpoints,fermienergy,MPI_COMM_WORLD)
	integer,intent(in)				::	Hdim,ndim,nkpt,cutofup,cutofdown,MPI_COMM_WORLD
	type(kpoint), pointer, intent(in)		::	kpoints(:)
	real(kind=8),intent(out)			::	fermienergy

	real(kind=8),allocatable			::	tempeigenval(:)
	integer						:: 	ierr
	integer						::	i,j
	integer						::	maxocc, Nstates

	if(cutofdown.lt.1.or.cutofup.gt.Hdim.or.cutofdown.gt.cutofup.or.cutofup.eq.ndim)then
		write(*,*)	'Hdim',Hdim,'cutofdown',cutofdown,'cutofup',cutofup
		write(*,*)	'Problem with cutoffs for Calculating the Fermi energy, calling MPI_ABORT'
		call MPI_Abort(MPI_COMM_WORLD, '1', ierr)
		return
	endif

	maxocc=nkpt*(ndim-cutofdown+1)
	Nstates=cutofup-cutofdown+1
	allocate(tempeigenval(Nstates*nkpt))
	do i=1,Nstates
		do j=1,nkpt
			tempeigenval((i-1)*nkpt+j)=kpoints(j)%eigenvalues(i+cutofdown-1)
		enddo
	enddo
	write(*,*)	'start sorting of the ', size(tempeigenval), 'entries'
	call QsortC(tempeigenval)
	fermienergy=(tempeigenval(maxocc)+tempeigenval(maxocc+1))/2d0
!check cutoffs

	do i=1,nkpt
		if(kpoints(i)%eigenvalues(cutofdown).gt.fermienergy.or.kpoints(i)%eigenvalues(cutofup).lt.fermienergy) then
			write(*,*)	'Hdim',Hdim,'cutofdown',cutofdown,'cutofup',cutofup
			write(*,*)	'cutoffs for calculating Fermi energy chosen to tight'
				if(kpoints(i)%eigenvalues(cutofdown).gt.fermienergy)	write(*,*)	i,'th kpoint component of lower cutoff higher than calculated Fermi energy'
				if(kpoints(i)%eigenvalues(cutofup).lt.fermienergy)	write(*,*)	i,'th kpoint component of upper cutoff lesser than calculated Fermi energy'
			write(*,*) 'calling MPI_ABORT'
			call MPI_Abort(MPI_COMM_WORLD, '1', ierr)	
		endif
	enddo

	write(*,*) 'fermienergy', fermienergy

end subroutine getfermi











subroutine getBerrycurv(Hdim,eigenvals,eigenvects,bravv,kpoint,Berrycurvout)
	integer,intent(in)			::	Hdim
	real(kind=8),intent(in)			::	eigenvals(:),kpoint(:),bravv(:,:)
	complex(kind=8),intent(in)		::	eigenvects(:,:)
	real(kind=8),intent(out)		::	Berrycurvout(:,:)

	complex(kind=8),allocatable		::	momentum(:,:,:)
	complex(kind=8)				::	momentum2(Hdim,Hdim,3)
	complex(kind=8)				::	eigendk(Hdim,Hdim,3)
	
	integer				::	dir,m,n
	

	Call getmomentum(momentum,kpoint,bravv)		!get Hderivatives


	Do dir=1,3
		Do m=1,Hdim
			Do n=1,Hdim
				IF ( n .ne. m) momentum2(m,n,dir)=Dot_product(eigenvects(:,m),MATMUL(momentum(:,:,dir),eigenvects(:,n)))/(eigenvals(m)-eigenvals(n))
			End Do
		End Do
	End Do

	deallocate(momentum)
	eigendk=cmplx(0.d0,0.d0)
	Do dir=1,3
		Do m=1,Hdim
			Do n=1,Hdim
				eigendk(:,m,dir)=eigendk(:,m,dir)+momentum2(n,m,dir)*eigenvects(:,n)	
			end Do
		End Do
	End Do


	Do m=1,Hdim
		Berrycurvout(m,1)=-dimag(Dot_product(eigendk(:,m,2),eigendk(:,m,3)))+dimag(Dot_product(eigendk(:,m,3),eigendk(:,m,2)))
		Berrycurvout(m,2)=-dimag(Dot_product(eigendk(:,m,3),eigendk(:,m,1)))+dimag(Dot_product(eigendk(:,m,1),eigendk(:,m,3)))
		Berrycurvout(m,3)=-dimag(Dot_product(eigendk(:,m,1),eigendk(:,m,2)))+dimag(Dot_product(eigendk(:,m,2),eigendk(:,m,1)))
         End Do
	
	

end subroutine getBerrycurv






subroutine getprojeigenvector(Nocc,Hdim,projeigvec,eigenvects,kpt,MPI_COMM_WORLD)
integer,intent(in)		::	Nocc,Hdim, MPI_COMM_WORLD
complex(kind=8),intent(out)	::	projeigvec(:,:)
complex(kind=8),intent(in)	::	eigenvects(:,:)
real(kind=8),intent(in)		::	kpt(:)


complex(kind=8)			::	spineigenvec(Nocc,Nocc)
real(kind=8)			::	spineigenval(Nocc)
complex(kind=8)			::	diagmatrix(Nocc,Nocc)
complex(kind=8),allocatable	::	Spin(:,:,:)
integer				::	ne,ierr

integer				::	i,j,ii


	allocate(Spin(Hdim,Hdim,3))
	Call getspin(Spin,kpt)
	
	diagmatrix=cmplx(0.d0,0.d0)
         Do i=1,Nocc
            Do j=1,Nocc
               diagmatrix(i,j)=Dot_product(eigenvects(:,i),MATMUL(Spin(:,:,3),eigenvects(:,j)))
            End Do
         End Do

         call zheevx('V','A','U',Nocc,diagmatrix,Nocc,vl,vu,1,Nocc,abstol,ne,spineigenval,spineigenvec,Nocc,work,lwork,rwork,iwork, ifail,info)

	if(-spineigenval(Nocc/2)+spineigenval(Nocc/2+1).lt.1d0)then
		write(*,*) 'difference of inner projected eigenvalues lesser than 1, calling MPI_ABORT'
		call MPI_Abort(MPI_COMM_WORLD, '1', ierr)
	endif
	projeigvec=cmplx(0.d0,0.d0)
       	Do i=1,Nocc
       		Do j=1,Nocc
			do ii=1,Hdim
        			projeigvec(ii,i)=projeigvec(ii,i)+spineigenvec(j,i)*eigenvects(ii,j)
			end do
		End Do
       	End Do

end subroutine getprojeigenvector




subroutine getproj4chernno(kpoints,ndim,nkpt,chernnoout,spinchernout,MPI_COMM_WORLD)
	real(kind=8),intent(out)		::	chernnoout(:)
	type(kpoint), pointer			::	kpoints(:)
	integer,intent(in)			::	ndim,nkpt,MPI_COMM_WORLD
	real(kind=8),intent(out)		::	spinchernout

	integer					:: 	ierr,meshdim
	integer,allocatable			::	nkptdir(:)
	complex(kind=8)				::	pointU(2,nkpt)
	
	
	integer					::	m,i


	Call getkptdir(meshdim,nkptdir)	!not necessary?
	if(meshdim.ne.2)then
		write(*,*) 'dimension of kmesh not 2 => 4 point Chern number calculation not implemented'
		write(*,*) 'calling MPI_ABORT'
		call MPI_Abort(MPI_COMM_WORLD, '1', ierr)
	endif
		
	chernnoout=0.0d0
	do m=1,ndim
		pointU=0.0d0
		do i=1,nkpt
			pointU(1,i)=Dot_product(kpoints(i)%projeigenvectors(:,m),kpoints(i)%meshpointer(1)%next%projeigenvectors(:,m))
			pointU(1,i)=pointU(1,i)/abs(pointU(1,i))

			pointU(2,i)=Dot_product(kpoints(i)%projeigenvectors(:,m),kpoints(i)%meshpointer(2)%next%projeigenvectors(:,m))
			pointU(2,i)=pointU(2,i)/abs(pointU(2,i))

		enddo
		do i=1,nkpt
			chernnoout(m)=chernnoout(m)+dble(log(pointU(1,i)*pointU(2,kpoints(i)%meshpointer(1)%next%id)*(pointU(1,kpoints(i)%meshpointer(2)%next%id)**-1)*(pointU(2,i)**-1))/2/pi/cmplx(0.0,1.0))
		enddo
	enddo
	spinchernout=0d0
	do i=1,ndim/2
		spinchernout=spinchernout-chernnoout(i)/2d0	!spin down anteile
	enddo
	do i=ndim/2+1,ndim
		spinchernout=spinchernout+chernnoout(i)/2d0	!spin up alteile
	enddo

end subroutine getproj4chernno


subroutine get4chernno(kpoints,Hdim,nkpt,chernnoout,MPI_COMM_WORLD)
	real(kind=8),intent(out)		::	chernnoout(:)
	type(kpoint), pointer			::	kpoints(:)
	integer,intent(in)			::	Hdim,nkpt,MPI_COMM_WORLD

	integer					:: 	ierr,meshdim
	integer,allocatable			::	nkptdir(:)
	complex(kind=8)				::	pointU(2,nkpt)
	
	
	integer					::	m,i


	Call getkptdir(meshdim,nkptdir)	!not necessary?
	if(meshdim.ne.2)then
		write(*,*) 'dimension of kmesh not 2 => 4 point Chern number calculation not implemented'
		write(*,*) 'calling MPI_ABORT'
		call MPI_Abort(MPI_COMM_WORLD, '1', ierr)
	endif
		
	chernnoout=0.0d0
	do m=1,Hdim
		pointU=0.0d0
		do i=1,nkpt
			pointU(1,i)=Dot_product(kpoints(i)%eigenvectors(:,m),kpoints(i)%meshpointer(1)%next%eigenvectors(:,m))
			pointU(1,i)=pointU(1,i)/abs(pointU(1,i))

			pointU(2,i)=Dot_product(kpoints(i)%eigenvectors(:,m),kpoints(i)%meshpointer(2)%next%eigenvectors(:,m))
			pointU(2,i)=pointU(2,i)/abs(pointU(2,i))

		enddo
		do i=1,nkpt
			chernnoout(m)=chernnoout(m)+dble(log(pointU(1,i)*pointU(2,kpoints(i)%meshpointer(1)%next%id)*(pointU(1,kpoints(i)%meshpointer(2)%next%id)**-1)*(pointU(2,i)**-1))/2/pi/cmplx(0.0,1.0))
		enddo
	enddo

end subroutine get4chernno

end module processH




