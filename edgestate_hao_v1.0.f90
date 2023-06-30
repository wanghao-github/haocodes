program hao_edgestates
    !Hao Wang rewrite the legacy code for the calculation of edgestate in Juelich Forchungsentrum, 2022/10/27
    !  
    !mpiifort -CB -r8 -qmkl edgestate_hao.f90 -o edge_states_hao.x
    implicit none

    integer              :: nwannexup,nwannexdown,ix,iy,iz,band1,band2,rvecnum,nslab1,nslab2,ii,zvalue,ib
    real                 :: fermi,rdum,idum,pi,phase,lattice_vec(3,3)
    integer,allocatable  :: nrpts(:),irvec(:,:)
    complex,allocatable  :: hops(:,:,:),ham(:,:),fourHamilton(:,:,:),hamiltonian(:,:)
    real,allocatable     :: k(:),localisationpar(:)
    integer              :: num_lines,num_wann,i,j,ierr,irank,isize,ik,numkpts,fourdim,fourdirection
    integer              :: Hdim,numberlayer,layerspreadmax,layerspreadmin,layerspread,layerdir
    complex              :: fac
	integer,allocatable	 :: wannperat(:),atomnumberarray(:),fourdir(:),layerintarr(:)
	integer,allocatable  :: excludeup(:),excludedown(:),temp(:),wannierfunctioninHam(:)
    integer              :: ndiffatom,nexcludeup,nexcludedown,return_num_wann,locmin,locmax,i1,i2
    real                 :: vl,vu,abstol,ik1
    integer              :: length,length1,length2,length3,length4
    real,allocatable     :: eigvals(:),eigvals_per_k(:,:)
    complex,allocatable  :: eigvecs(:,:)
    integer              :: ne,info,lwork,ik_cpu
    complex,allocatable  :: work(:)
    integer,allocatable  :: iwork(:)
    real,allocatable     :: rwork(:)
    integer,allocatable  :: ifail(:)

    type                 :: atom
    integer              :: number
    real(kind=8)         :: position(3)		    !position in realspace
    integer,allocatable  :: wannierfunctions(:)	!numbers of wannierfunctions attributed to the atom
    end type atom

    type(atom), allocatable ::   atomarr(:)

    INCLUDE 'mpif.h'
    integer :: stt(MPI_STATUS_SIZE)
    CALL MPI_INIT(ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,isize,ierr)

    abstol=2.0*tiny(abstol)
    pi = 3.14159265
    numkpts = 400
    if(irank.eq.0)then
        write(*,*) "isize=", isize
        open(200,file='hopping.1')
        num_lines = 0
        num_wann  = 0
        do
            read(200,fmt=*,end=201)ix,iy,iz,band1,band2,rdum,idum
            num_lines=num_lines+1
            num_wann=max(num_wann,band1)
        enddo
201     continue

        rvecnum=num_lines/(num_wann*num_wann)    
        allocate(hops(num_wann,num_wann,rvecnum))
        allocate(irvec(3,rvecnum))
        hops=0.0
        rewind(200) 
        num_lines=0
        do
            read(200,fmt=*,end=211)ix,iy,iz,band1,band2,rdum,idum
            num_lines=num_lines+1
            rvecnum=(num_lines-1)/(num_wann*num_wann)+1
            irvec(1,rvecnum)=ix
            irvec(2,rvecnum)=iy
            irvec(3,rvecnum)=iz
            hops(band1,band2,rvecnum)=cmplx(rdum,idum)
        enddo
211     continue
        close(200)
    endif

    if(irank.eq.0)then
        write(*,*) "irank=", irank
        write(*,*)"rvecnum=",rvecnum
        write(*,*)"num_wann=",num_wann
        write(*,*) size(irvec)
    endif

    call mpi_bcast(num_wann,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    call mpi_bcast(rvecnum,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    if(.not.allocated(hops)) then 
        allocate(hops(num_wann,num_wann,rvecnum))
    endif 
    length1=num_wann*num_wann*rvecnum
    call mpi_bcast(hops,length1,MPI_DOUBLE_COMPLEX,0,mpi_comm_world,ierr)

    if(irank == 0)then
        write(*,*) "irank=", irank
        open(100,file='layer_inp')
        read(100,*) numberlayer
        allocate(atomnumberarray(num_wann))
        read(100,*)	ndiffatom
        allocate(atomarr(ndiffatom))
        atomnumberarray=0
        read(100,'(<num_wann>I3)') atomnumberarray(:)
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
            if(ANY(atomarr(i)%wannierfunctions(:).eq.0).or.temp(i)-1.ne.wannperat(i)) then
                write(*,*)	'mismatch assigning wannier functions to atoms, CALLING MPI ABORT'
            endif
        enddo
        deallocate(temp)
        deallocate(wannperat)
        do i=1,ndiffatom
            read(100,*) atomarr(i)%position(:)
        enddo
        read(100,'(2I3)')	nexcludeup,nexcludedown
        if(nexcludeup.eq.0)then
            read(100,*)
            allocate(excludeup(0))
        else
            allocate(excludeup(nexcludeup))
            read(100,'(<nexcludeup>I3)')	    excludeup(:)
        endif
        if(nexcludedown.eq.0)then
            read(100,*)
            allocate(excludedown(0))
        else
            allocate(excludedown(nexcludedown))
            read(100,'(<nexcludedown>I3)')      excludedown(:)
        endif
        read(100,'(I3)')                        fourdim
        allocate(fourdir(fourdim))
        read(100,'(I3,<fourdim>I3)')            layerdir,fourdir(:) 
        close(100)
    endif

    allocate(nrpts(rvecnum))

    if(irank.eq.0)then 
        open(400,file='nrpts_inp')
        do j=1,rvecnum/15
            read(400,'(15I5)') (nrpts(15*(j-1)+i),i=1,15)
        enddo
        read(400,'(<mod(rvecnum,15)>I5)') (nrpts(15*(rvecnum/15)+i),i=1,mod(rvecnum,15))
        close(400)
    endif

    call mpi_bcast(nrpts,rvecnum,MPI_INTEGER,0,mpi_comm_world,ierr)

    if(irank.eq.0)then
        write(*,*) "irank=", irank
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

        write(*,*) "Hdim is", Hdim
        allocate(wannierfunctioninHam(Hdim))
        allocate(localisationpar(Hdim))

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
                        localisationpar(i)=atomarr(i1)%position(layerdir)+numberlayer-1
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
        enddo
        locmin=MINVAL(localisationpar,Hdim)
        locmax=MAXVAL(localisationpar,Hdim)
    	write(*,*)	"locmin,locmax", locmin,locmax
        do i=1,Hdim
            localisationpar(i)=-1d0+2d0*(localisationpar(i)-locmin)/(locmax-locmin)
        enddo

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

    endif

    allocate(eigvals_per_k(numkpts,Hdim))
    eigvals_per_k=0.0
    call mpi_bcast(numberlayer,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    call mpi_bcast(return_num_wann,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    call mpi_bcast(Hdim,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    call mpi_bcast(nexcludeup,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    call mpi_bcast(nexcludedown,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    call mpi_bcast(nwannexup,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    call mpi_bcast(nwannexdown,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    call mpi_bcast(layerspreadmin,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    call mpi_bcast(layerspreadmax,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    call mpi_bcast(layerspread,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    call mpi_bcast(locmax,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    call mpi_bcast(locmin,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    call mpi_bcast(ndiffatom,1,MPI_INTEGER,0,mpi_comm_world,ierr)

    if(.not.allocated(atomarr))then 
        allocate(atomarr(ndiffatom))
    endif
    call mpi_bcast(atomarr,ndiffatom,MPI_INTEGER,0,mpi_comm_world,ierr)

    if(.not.allocated(excludeup))then 
        allocate(excludeup(nexcludeup))
    endif
    call mpi_bcast(excludeup,nexcludeup,MPI_INTEGER,0,mpi_comm_world,ierr)

    if(.not.allocated(excludedown))then 
        allocate(excludedown(nexcludedown))
    endif
    call mpi_bcast(excludedown,nexcludedown,MPI_INTEGER,0,mpi_comm_world,ierr)
    if(.not.allocated(layerintarr))then 
        allocate(layerintarr(Hdim))
    endif
    call mpi_bcast(layerintarr,Hdim,MPI_INTEGER,0,mpi_comm_world,ierr)

    
        allocate(eigvals(Hdim))
        allocate(eigvecs(Hdim,Hdim))
        lwork=12.0*Hdim
        allocate(work(lwork) )
        allocate(rwork(7*Hdim) )
        allocate(iwork(5*Hdim) )
        allocate(ifail(Hdim) )
        allocate(k(numkpts))
        k=0.0
        allocate(fourHamilton(layerspreadmin:layerspreadmax,num_wann,num_wann))
        allocate(hamiltonian(Hdim,Hdim))
        call mpi_bcast(fourdirection,1,MPI_INTEGER,0,mpi_comm_world,ierr)
        call mpi_bcast(fourdim,1,MPI_INTEGER,0,mpi_comm_world,ierr)
        call mpi_bcast(numkpts,1,MPI_INTEGER,0,mpi_comm_world,ierr)

        call mpi_bcast(k,numkpts,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)
        
    if(irank.eq.0)then 
        do ik=1,numkpts
            k(ik) = ik*3*pi/numkpts
        enddo
    endif
    ik_cpu = 0
    do ik=1,numkpts
        ik_cpu=ik_cpu+1
        if(mod(ik_cpu-1,isize).ne.irank) cycle

            write(*,*) "ik", ik, "k(ik)=", k(ik),"irank=",irank
            eigvals_per_k(ik,:)=0
            fourHamilton=0d0
            do ii=1,rvecnum   
                
                do i=1,num_wann
                    do j=1,num_wann
                        phase=0d0
                        do i1=1,fourdim
                            fourdirection=fourdir(i1)
                            phase=phase+(irvec(fourdirection,ii))*k(ik)
                        enddo
                        fac=cmplx(cos(phase),sin(phase))
                        fourHamilton(irvec(layerdir,ii),i,j)=fourHamilton(irvec(layerdir,ii),i,j)+fac*hops(i,j,ii)/nrpts(ii)
                    enddo       
                enddo 
            enddo

        hamiltonian=cmplx(0d0,0d0)
        do i=1,Hdim
            do j=1,Hdim
                zvalue = layerintarr(i)-layerintarr(j)
                if (zvalue.ge.layerspreadmin.and.zvalue.le.layerspreadmax) then
                    hamiltonian(i,j)=hamiltonian(i,j)+fourHamilton(zvalue,wannierfunctioninham(i),wannierfunctioninham(j))
                endif
            enddo
        enddo
        call zheevx('V','A','U',Hdim,hamiltonian,Hdim,vl,vu,1,Hdim,abstol,ne,eigvals,eigvecs,Hdim,work,lwork,rwork,iwork, ifail,info)
        do ib=1,Hdim
            eigvals_per_k(ik,ib) = eigvals(ib)
        enddo
    enddo
    if(irank.eq.0)then   
        open(222,file='kpts.out',recl=10000)
        do ik=1,numkpts
            write(222,*),k(ik)
        enddo
        
        open(234,file='layeri',recl=10000)
        do i=1,Hdim
            write(234,*),layerintarr(i)
        enddo

        open(123,file='output_bands',recl=10000)
        do ib=1,Hdim     
            do ik=1,numkpts        
                write(123,*), k(ik),eigvals_per_k(ik,ib)
            enddo
            write(123,*)
        enddo
    endif
    write(*,*) "size of hamiltonian", size(hamiltonian)
    call mpi_barrier(mpi_comm_world,ierr)
    call MPI_Finalize(ierr)
end program hao_edgestates