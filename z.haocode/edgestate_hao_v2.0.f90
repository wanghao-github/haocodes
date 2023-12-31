program hao_edgestates
    !Hao Wang rewrite the legacy code for the calculation of edgestate in Juelich Forchungsentrum, 2022/10/27
    !mpiifort -CB -r8 -qmkl edgestate_hao.f90 -o edge_states_hao.x
    !  mpiifort -CB -r8 edgestate_hao.f90   -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -o edge_states_hao.x
    implicit none

    integer              :: nwannexup,nwannexdown,ix,iy,iz,band1,band2,rvecnum,nslab1,nslab2,ii,zvalue,ib,sendcount,Np,ijmax,Ndim,io
    real                 :: fermi,rdum,idum,pi,phase,lattice_vec(3,3),w,emin,emax
    integer,allocatable  :: nrpts(:),irvec(:,:)
    complex,allocatable  :: hops(:,:,:),ham(:,:),fourHamilton(:,:,:),hamiltonian(:,:)
    real,allocatable     :: k(:),localisationpar(:)
    integer              :: num_lines,num_wann,i,j,ierr,irank,status,request,isize,ik,numkpts,fourdim,fourdirection
    integer              :: Hdim,numberlayer,layerspreadmax,layerspreadmin,layerspread,layerdir
    complex              :: fac
    integer,allocatable  :: wannperat(:),atomnumberarray(:),fourdir(:),layerintarr(:)
    integer,allocatable  :: excludeup(:),excludedown(:),temp(:),wannierfunctioninHam(:)
    integer              :: ndiffatom,nexcludeup,nexcludedown,return_num_wann,locmin,locmax,i1,i2
    real                 :: vl,vu,abstol,ik1,time_start,time_end,eta
    integer              :: length,length1,length2,length3,length4,length5
    real,allocatable     :: eigvals(:),eigvals_per_k(:,:),temp_array(:),eigvals_per_k_mpi(:,:)
    complex,allocatable  :: eigvecs(:,:)
    integer              :: ne,info,lwork,ik_cpu,omeganum, omegamax,omegamin,nw_half
    complex,allocatable  :: work(:)
    integer,allocatable  :: iwork(:)
    real,allocatable     :: rwork(:)
    integer,allocatable  :: ifail(:)

    real(kind(1.0d0)), allocatable :: omega(:)
    real(kind(1.0d0)), allocatable :: dos_l(:,:)
    real(kind(1.0d0)), allocatable :: dos_r(:,:)
    real(kind(1.0d0)), allocatable :: dos_l_only(:,:)
    real(kind(1.0d0)), allocatable :: dos_r_only(:,:)
    real(kind(1.0d0)), allocatable :: dos_l_mpi(:,:)
    real(kind(1.0d0)), allocatable :: dos_r_mpi(:,:)
    real(kind(1.0d0)), allocatable :: dos_bulk(:,:)
    real(kind(1.0d0)), allocatable :: dos_bulk_mpi(:,:)
    ! real(kind(1.0d0)), allocatable :: dos_r_only(:,:)
    ! real(kind(1.0d0)), allocatable :: dos_l_only(:,:)


    complex(kind(1.0d0)), allocatable :: GLL(:,:)
    complex(kind(1.0d0)), allocatable :: GRR(:,:)
    complex(kind(1.0d0)), allocatable :: GB (:,:)
    complex(kind(1.0d0)), allocatable :: H00(:,:)
    complex(kind(1.0d0)), allocatable :: H01(:,:)
    complex(kind(1.0d0)), allocatable :: ones(:,:)


    REAL(kind(1.0d0)),    ALLOCATABLE  :: sx_l(:, :), sy_l(:, :), sz_l(:, :)
    REAL(kind(1.0d0)),    ALLOCATABLE  :: sx_r(:, :), sy_r(:, :), sz_r(:, :)
    REAL(kind(1.0d0)),    ALLOCATABLE  :: sx_l_mpi(:, :), sy_l_mpi(:, :), sz_l_mpi(:, :)
    REAL(kind(1.0d0)),    ALLOCATABLE  :: sx_r_mpi(:, :), sy_r_mpi(:, :), sz_r_mpi(:, :)
    COMPLEX(kind(1.0d0)), ALLOCATABLE  :: sigma_x(:,:), sigma_y(:,:), sigma_z(:,:)
    COMPLEX(kind(1.0d0)), ALLOCATABLE  :: ctemp(:,:)

    type                 :: atom
    integer              :: number
    real(kind=8)         :: position(3)         !position in realspace
    integer,allocatable  :: wannierfunctions(:) !numbers of wannierfunctions attributed to the atom
    end type atom

    type(atom), allocatable ::   atomarr(:)

    INCLUDE 'mpif.h'
    integer :: stt(MPI_STATUS_SIZE)
    CALL MPI_INIT(ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,isize,ierr)

    abstol=2.0*tiny(abstol)
    pi = 3.14159265
    Np = 2
    ijmax = 10
    ! numkpts = 1280
    omeganum = 500
    omegamax = -3
    omegamin = -4
    
    eta=(omegamax- omegamin)/dble(omeganum)*1.5d0
    
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


    call MPI_Barrier(mpi_comm_world, ierr)
    Ndim = Np*num_wann


    if(.not.allocated(irvec)) then 
        allocate(irvec(3,rvecnum))
    endif 
    call mpi_bcast(irvec,size(irvec),MPI_INTEGER,0,mpi_comm_world,ierr)

    if(.not.allocated(hops)) then 
        allocate(hops(num_wann,num_wann,rvecnum))
    endif 
    length1=num_wann*num_wann*rvecnum
    call mpi_bcast(hops,length1,MPI_DOUBLE_COMPLEX,0,mpi_comm_world,ierr)
    write(*,*) irank
    if(irank == 0)then
        write(*,*) "irank=", irank
        open(100,file='layer_inp')
        read(100,*) numkpts
        read(100,*) numberlayer
        allocate(atomnumberarray(num_wann))
        read(100,*)	ndiffatom
        allocate(atomarr(ndiffatom))
        atomnumberarray=0
        read(100,'(<num_wann>I3)') atomnumberarray(:)
        allocate(wannperat(ndiffatom))
        wannperat=0
        do i=1,num_wann                                       !wannperatom存了每个原子上的轨道个数
            wannperat(atomnumberarray(i))=wannperat(atomnumberarray(i))+1
        enddo
        do i=1,ndiffatom
            atomarr(i)%number=i
            allocate(atomarr(i)%wannierfunctions(wannperat(i)))     !每个原子的wannierfunctions的维度是每个原子上轨道的个数
            atomarr(i)%wannierfunctions(wannperat(i))=0
        enddo
        allocate(temp(ndiffatom))
        temp=1          !temp 维度是原子个数 初始化为全1
        do i=1,num_wann
            atomarr(atomnumberarray(i))%wannierfunctions(temp(atomnumberarray(i)))=i  ! 前面是确定了原子的编号  对应的wannie轨道是i
            temp(atomnumberarray(i))=temp(atomnumberarray(i))+1
        enddo
      !!! 比如 我的atomnumberarray是  1 1 1 1 2 2 2 2 3 3 4 4 4 5 5 5的话 
      !!! temp一开始 是 1 1 1 1 1
      !!  num_wann是16 循环第一句前半截获得了原子编号 temp循环变化 21111 ,31111,4111第一个原子的wannierfunctions 1 2 3 4
      !!  到i=5的时候就是第二个原子了   这样temp(2)=1 然后temp变化为4 2 1 1 1,4 3 1 1 1, 4 4 1 1 1, 第二个原子的wannierfunctions 为i就是5 6 7 8
      !!! 到i=9的时候  temp(3) 从一开始  第三个原子上的wannierfunctions按i变化是9 10
      !!!! 总的来说每个原子的wannierfunctions获得在num_wann里面的编号
        write(*,*) "temp ok"
        do i=1,ndiffatom
            if(ANY(atomarr(i)%wannierfunctions(:).eq.0).or.temp(i)-1.ne.wannperat(i)) then
                write(*,*) 'mismatch assigning wannier functions to atoms, CALLING MPI ABORT'
            endif
        enddo
        deallocate(temp)
        deallocate(wannperat)
        do i=1,ndiffatom
            read(100,*) atomarr(i)%position(:)
        enddo
        read(100,'(2I3)') nexcludeup,nexcludedown
        if(nexcludeup.eq.0)then
            read(100,*)
            allocate(excludeup(0))
        else
            allocate(excludeup(nexcludeup))
            read(100,'(<nexcludeup>I3)') excludeup(:)
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
        write(*,*) "nrpts_ok"
        write(*,*) "numkpts" ,numkpts
        write(*,*) "fourdim" ,fourdim
        write(*,*) "fourdir" ,fourdir
        write(*,*) "layerdir" ,layerdir
    endif

    call mpi_bcast(numkpts,1,MPI_INTEGER,0,mpi_comm_world,ierr)

    if(.not.allocated(nrpts))then 
        allocate(nrpts(rvecnum))
    endif
    call mpi_bcast(nrpts,rvecnum,MPI_INTEGER,0,mpi_comm_world,ierr)
    
    call mpi_bcast(numberlayer,1,MPI_INTEGER,0,mpi_comm_world,ierr)

    call mpi_bcast(ndiffatom,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    write(*,*) "ndiffatom = ",ndiffatom,irank 

    call mpi_bcast(fourdim,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    write(*,*) "fourdim=",fourdim,irank


    call MPI_Barrier(mpi_comm_world, ierr)

    if(.not.allocated(fourdir))then 
        allocate(fourdir(fourdim))
    endif
    call mpi_bcast(fourdir,fourdim,MPI_INTEGER,0,mpi_comm_world,ierr)
    write(*,*) "fourdir=",fourdir,irank 

    call mpi_bcast(layerdir,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    write(*,*) "layerdir =",layerdir,irank 

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
            if(i.le.num_wann-nwannexup)then ! 轨道编号在第一个原胞之内
                i1=1  !计数原子在超胞中的编号
                i2=1  !记录每原子里面的wannier编号
                do while (i.le.num_wann-nwannexup)
                    if(.not.ANY(excludeup==i1))then   !如果i1不在任何排除的原子之内
                        wannierfunctioninHam(i)=atomarr(i1)%wannierfunctions(i2) !把这个原子的每个wannier轨道编号存到超胞Ham里
                        localisationpar(i)=atomarr(i1)%position(layerdir)+numberlayer-1
                        layerintarr(i)=1 
    
                        if(i2.lt.size(atomarr(i1)%wannierfunctions(:)))then !i2是每个原子上的wannier编号
                            i2=i2+1                                            
                        else
                            i2=1         !超过了每个原子上的wannier轨道个数后重新计数
                            i1=i1+1      !切换到下一个原子
                        endif
                        i=i+1     !切换判断下个超胞wannier轨道 判断他在单胞中的编号
                    else
                        i1=i1+1   !如果这个原子被排除了
                    endif
                enddo
            elseif(i.le.Hdim+nwannexdown-num_wann)then   !除了第一个原胞外的其他轨道
                do i1=1,ndiffatom                         ! 遍历所有原子
                    do i2=1,size(atomarr(i1)%wannierfunctions(:))        !i2仍然是每个原子上的wannier轨道个数
                        wannierfunctioninHam(i)=atomarr(i1)%wannierfunctions(i2)!把超胞wannierfunctioninHam里面每个轨道对应到单胞中的轨道编号上
                        localisationpar(i)=atomarr(i1)%position(layerdir)+numberlayer-(i+nwannexup-1)/num_wann-1
                        !每个原子在超胞中的位置  相比较于上一个在最边缘的多减了(i+nwannexup-1)/num_wann这一项 假如i是48 num_wann是44的话这样就多减了1 因为
                        !Fortran 整型相除舍弃了小数部分 相当于把layerdir那个方向的加了一个整数 直到-2 -3 -4 ....
                        layerintarr(i)=(i+nwannexup-1)/num_wann+1   !这个地方和前面一样相当于计算出了i的层数
                        i=i+1                                      !计算下一个i
                    enddo
                enddo
            else                   !剩下的最下面不在excludedown里面的那半截 画图的话大概就是到最后Hdim那一点
                i1=1
                i2=1
                do while (i.le.Hdim)                  !最后到HDIM那一点的轨道
                    if(.not.ANY(excludedown==i1))then    !如果这剩下一点的i1没在排除轨道里面
                        wannierfunctioninHam(i)=atomarr(i1)%wannierfunctions(i2)  !判断超胞ham里每个轨道对应于单胞的哪个轨道
                        localisationpar(i)=atomarr(i1)%position(layerdir)!每个原子的位置 可以看出是从down开始计数的
                        layerintarr(i)=numberlayer           !这是最后一部分
                        if(i2.lt.size(atomarr(i1)%wannierfunctions(:)))then  !跟上面一样 如果每个原子上的wannier轨道用完了
                            i2=i2+1                                          !就开始循环下一个原子
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
        locmin=MINVAL(localisationpar,Hdim)                              ! 获取在这个超胞中最小的层数编号
        locmax=MAXVAL(localisationpar,Hdim)                              ! 最大
        write(*,*) "locmin,locmax", locmin,locmax
        do i=1,Hdim
            localisationpar(i)=-1d0+2d0*(localisationpar(i)-locmin)/(locmax-locmin) !把层数编号替换为分数坐标
        enddo

        layerspreadmin=0
        layerspreadmax=0
        do i=1,rvecnum
            if(irvec(layerdir,i).gt.layerspreadmax)then             !获取recrum在layerdir方向上的最大值和最小值
                layerspreadmax=irvec(layerdir,i)             
            endif
            if(irvec(layerdir,i).lt.layerspreadmin)then
                layerspreadmin=irvec(layerdir,i)
            endif
        enddo
        layerspread=layerspreadmax-layerspreadmin                   !获得原Hr在这个方向上的差值

    endif

    if (irank.eq.0) then    
        write(*,*)"here is no problem2"
      endif

    call mpi_bcast(layerspreadmin,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    call mpi_bcast(layerspreadmax,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    call mpi_bcast(layerspread,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    call mpi_bcast(locmax,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    call mpi_bcast(locmin,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    
    if (irank.eq.0) then
        write(*,*)"layerspreadmin",layerspreadmin
        write(*,*)"layerspreadmax",layerspreadmax
        write(*,*)"here is no problem"
    endif

    call mpi_bcast(Hdim,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    call mpi_bcast(return_num_wann,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    write(*,*)"here is no problem2132131"
call MPI_Barrier(mpi_comm_world, ierr)

        allocate(eigvals_per_k(numkpts,Hdim))
        allocate(eigvals_per_k_mpi(numkpts,Hdim))
    ! endif
    ! call mpi_bcast(eigvals_per_k,size(eigvals_per_k),MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)   
    eigvals_per_k=0.0
    eigvals_per_k_mpi=0.0
    if (irank.eq.0) then
        write(*,*)"here is no problem1"
    endif

    if(.not.allocated(atomarr))then 
        allocate(atomarr(ndiffatom))
    endif
    call mpi_bcast(atomarr,ndiffatom,MPI_INTEGER,0,mpi_comm_world,ierr)

    if (irank.eq.0) then    
        write(*,*)"here is no problem3"
      endif

    if(.not.allocated(excludeup))then 
        allocate(excludeup(nexcludeup))
    endif
    call mpi_bcast(excludeup,nexcludeup,MPI_INTEGER,0,mpi_comm_world,ierr)

    if (irank.eq.0) then    
        write(*,*)"here is no problem4"
      endif

    if(.not.allocated(excludedown))then 
        allocate(excludedown(nexcludedown))
    endif
    call mpi_bcast(excludedown,nexcludedown,MPI_INTEGER,0,mpi_comm_world,ierr)

    if (irank.eq.0) then    
        write(*,*)"here is no problem5"
      endif

    if(.not.allocated(layerintarr))then 
        allocate(layerintarr(Hdim))
    endif
    call mpi_bcast(layerintarr,Hdim,MPI_INTEGER,0,mpi_comm_world,ierr)

    if (irank.eq.0) then    
        write(*,*)"here is no problem6"
      endif

        allocate(eigvals(Hdim))
        allocate(eigvecs(Hdim,Hdim))
        allocate(temp_array(Hdim))
        lwork=12.0*Hdim
        allocate(work(lwork) )
        allocate(rwork(7*Hdim) )
        allocate(iwork(5*Hdim) )
        allocate(ifail(Hdim) )
        allocate(k(numkpts))
        allocate(omega(omeganum))
        allocate(ones(Ndim,Ndim))
        allocate(GLL(Ndim, Ndim))
        allocate(GRR(Ndim, Ndim))
        allocate(GB (Ndim, Ndim))
        allocate(H00(Ndim, Ndim))
        allocate(H01(Ndim, Ndim))

        GLL= 0d0
        GRR= 0d0
        GB = 0d0
        H00= 0d0
        H01= 0d0
        ones= 0d0
        k=0.0

        allocate(dos_l(numkpts, omeganum))
        allocate(dos_r(numkpts, omeganum))
        allocate(dos_l_only(numkpts, omeganum))
        allocate(dos_r_only(numkpts, omeganum))
        allocate(dos_l_mpi(numkpts, omeganum))
        allocate(dos_r_mpi(numkpts, omeganum))
        allocate(dos_bulk(numkpts, omeganum))
        allocate(dos_bulk_mpi(numkpts, omeganum))

        omega=0d0
        dos_l=0d0
        dos_r=0d0
        dos_l_only=0d0
        dos_r_only=0d0
        dos_l_mpi=0d0
        dos_r_mpi=0d0
        dos_bulk=0d0
        dos_bulk_mpi=0d0

        
        ALLOCATE( sx_l(numkpts, omeganum), sy_l(numkpts, omeganum), sz_l(numkpts, omeganum))
        ALLOCATE( sx_l_mpi(numkpts, omeganum), sy_l_mpi(numkpts, omeganum), sz_l_mpi(numkpts, omeganum))
        ALLOCATE( sx_r(numkpts, omeganum), sy_r(numkpts, omeganum), sz_r(numkpts, omeganum))
        ALLOCATE( sx_r_mpi(numkpts, omeganum), sy_r_mpi(numkpts, omeganum), sz_r_mpi(numkpts, omeganum))
        ALLOCATE( sigma_x(ndim,ndim), sigma_y(ndim,ndim), sigma_z(ndim,ndim), ctemp(ndim,ndim))
        sigma_x      = 0d0;      sigma_y      = 0d0;      sigma_z      = 0d0
        sx_l         = 0d0;      sy_l         = 0d0;      sz_l         = 0d0
        sx_r         = 0d0;      sy_r         = 0d0;      sz_r         = 0d0
        sx_l_mpi     = 0d0;      sy_l_mpi     = 0d0;      sz_l_mpi     = 0d0
        sx_r_mpi     = 0d0;      sy_r_mpi     = 0d0;      sz_r_mpi     = 0d0


        if (irank.eq.0) then    
            write(*,*)"here is no problem7"
        endif

        allocate(fourHamilton(layerspreadmin:layerspreadmax,num_wann,num_wann))
        allocate(hamiltonian(Hdim,Hdim))

        if(.not.allocated(wannierfunctioninham))then 
            allocate(wannierfunctioninham(Hdim))
        endif
        call mpi_bcast(wannierfunctioninham,Hdim,MPI_INTEGER,0,mpi_comm_world,ierr)

        do ik=1,numkpts
            ! k(ik) = ik*3/numkpts
            k(ik) = ik*3*pi/numkpts
        enddo

        do i= 1, omeganum
            omega(i)=omegamin+(i-1)*(omegamax-(omegamin))/dble(omeganum)
        enddo
        
        do i=1,Ndim
            ones(i,i)=1.0d0
        enddo
        
        nw_half = Num_wann/2
        do i=1, Np
           do j=1, nw_half
              sigma_x( Num_wann*(i-1)+j        , Num_wann*(i-1)+j+nw_half ) =  1.0d0
              sigma_x( Num_wann*(i-1)+j+nw_half, Num_wann*(i-1)+j         ) =  1.0d0
              sigma_y( Num_wann*(i-1)+j        , Num_wann*(i-1)+j+nw_half ) = -cmplx(0,1)
              sigma_y( Num_wann*(i-1)+j+nw_half, Num_wann*(i-1)+j         ) =  cmplx(0,1)
              sigma_z( Num_wann*(i-1)+j        , Num_wann*(i-1)+j         ) =  1.0d0
              sigma_z( Num_wann*(i-1)+j+nw_half, Num_wann*(i-1)+j+nw_half ) = -1.0d0
           enddo 
        enddo

        call mpi_barrier(mpi_comm_world,ierr)

    ik_cpu = 0
    do ik=1,numkpts
        ik_cpu=ik_cpu+1
        if(mod(ik_cpu-1,isize).ne.irank) cycle
            fourHamilton=0d0
            do ii=1,rvecnum   
                do i=1,num_wann
                    do j=1,num_wann
                        phase=0d0
                        do i1=1,fourdim
                            fourdirection=fourdir(i1)
                            phase=phase+(irvec(fourdirection,ii))*k(ik)
                        enddo
                        fac=cmplx(cos(phase),sin(phase))                    ! (原始的单胞Hr只在fourdir上做傅里叶变换)
                        fourHamilton(irvec(layerdir,ii),i,j)=fourHamilton(irvec(layerdir,ii),i,j)+fac*hops(i,j,ii)/nrpts(ii)
                    enddo       
                enddo 
            enddo
        write(*,*) "fourHam no problem irank=",irank
        hamiltonian=cmplx(0d0,0d0)
        do i=1,Hdim
            do j=1,Hdim
                zvalue = layerintarr(i)-layerintarr(j)      !计算超胞里面所有轨道之间的层号差别
                if (zvalue.ge.layerspreadmin.and.zvalue.le.layerspreadmax) then   !如果这个层号差别的数在单胞的里面
                    hamiltonian(i,j)=hamiltonian(i,j)+fourHamilton(zvalue,wannierfunctioninham(i),wannierfunctioninham(j))
                    !那么就把这些进行过傅里叶变换的矩阵元加在一起  超胞中的轨道编号与单方向傅里叶变换后的单胞wannier编号的对应
                    !注意这个是差值相同的都会加起来  比如3和1  4和2 ....
                endif
            enddo
        enddo

        write(*,*) "Ham no problem, irank =" ,irank

     call zheevx('V','A','U',Hdim,hamiltonian,Hdim,vl,vu,1,Hdim,abstol,ne,eigvals,eigvecs,Hdim,work,lwork,rwork,iwork, ifail,info)
       !对角化这个超胞哈密顿量
     write(*,*) "zheevx no problem, irank =" ,irank

       eigvals_per_k(ik, :) = eigvals(:)
       

       write(*,*) "write eigvals OK, irank =" ,irank

       !!!!接下来抄的WT的
        do i=1,Np
            do j=1,Np
                if (abs(i-j).le.(ijmax)) then
                    H00(num_wann*(i-1)+1:num_wann*i,num_wann*(j-1)+1:num_wann*j)=fourHamilton(j-i,:,:)
                    !!!这个是([[0 1],[-1,0]])的大块矩阵
                endif
            enddo
        enddo
        write(*,*) "H00 ok" 

        !!! H01new
        do i=1,Np
            do j=Np+1,Np*2
                if (j-i.le.ijmax) then
                    H01(num_wann*(i-1)+1:num_wann*i,num_wann*(j-1-Np)+1:num_wann*(j-Np))=fourHamilton(j-i,:,:)
                endif
            enddo
        enddo
        write(*,*) "H01 ok"
        ! !!!H01new 是([[2 3],[1,2]])的大块矩阵
        
        do j = 1, omeganum
            w=omega(j)
            call surfgreen_1985(Ndim,eta,w,GLL,GRR,GB,H00,H01,ones)
            ! write(*,*) "GRR GLL ok"
            do i= 1,num_wann
                dos_l(ik, j)=dos_l(ik,j)- aimag(GLL(i,i))
            enddo ! i
            ! write(*,*) "dos_l ok"
            do i= 1, num_wann
                io= Ndim- num_wann+i
                dos_r(ik, j)=dos_r(ik,j)- aimag(GRR(io,io))
            enddo ! i
            ! write(*,*) "dos_r ok"
            do i= 1, Ndim
                dos_bulk(ik, j)=dos_bulk(ik,j)- aimag(GB(i,i))
            enddo ! i
            ! write(*,*) "dos_bulk ok"



            !!! spinresolved
            call mat_mul(ndim,gll,sigma_x,ctemp)
            do i = 1, num_wann
                sx_l_mpi(ik, j) = sx_l_mpi(ik, j)- AIMAG(ctemp(i,i))
            enddo ! i
            call mat_mul(ndim,gll,sigma_y,ctemp)
            do i = 1, num_wann
                sy_l_mpi(ik, j) = sy_l_mpi(ik, j)- AIMAG(ctemp(i,i))
            enddo !
            call mat_mul(ndim,gll,sigma_z,ctemp)
            do i = 1, num_wann
                sz_l_mpi(ik, j) = sz_l_mpi(ik, j)- AIMAG(ctemp(i,i))
            enddo ! i
            call mat_mul(ndim,grr,sigma_x,ctemp)
            do i = 1, num_wann
                sx_r_mpi(ik, j) = sx_r_mpi(ik, j)- AIMAG(ctemp(i,i))
            enddo ! i
            call mat_mul(ndim,grr,sigma_y,ctemp)
            do i = 1, num_wann
                sy_r_mpi(ik, j) = sy_r_mpi(ik, j)- AIMAG(ctemp(i,i))
            enddo !
            call mat_mul(ndim,grr,sigma_z,ctemp)
            do i = 1, num_wann
                sz_r_mpi(ik, j) = sz_r_mpi(ik, j)- AIMAG(ctemp(i,i))
            enddo ! i
        enddo ! j
!
    

!    call mpi_allreduce(dos_l, dos_l_mpi, size(dos_l),MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
!    call mpi_allreduce(dos_r, dos_r_mpi, size(dos_r),MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
!    call mpi_allreduce(dos_bulk, dos_bulk_mpi, size(dos_bulk),MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    
    write(*,*) "AllREDUCE DONE"

enddo

!!! ALLREDUCE spin
call mpi_allreduce(sx_l_mpi, sx_l, SIZE(sx_l), mpi_double_precision,mpi_sum, mpi_comm_world, ierr)
call mpi_allreduce(sy_l_mpi, sy_l, SIZE(sy_l), mpi_double_precision,mpi_sum, mpi_comm_world, ierr)
call mpi_allreduce(sz_l_mpi, sz_l, SIZE(sz_l), mpi_double_precision,mpi_sum, mpi_comm_world, ierr)
call mpi_allreduce(sx_r_mpi, sx_r, SIZE(sx_r), mpi_double_precision,mpi_sum, mpi_comm_world, ierr)
call mpi_allreduce(sy_r_mpi, sy_r, SIZE(sy_r), mpi_double_precision,mpi_sum, mpi_comm_world, ierr)
call mpi_allreduce(sz_r_mpi, sz_r, SIZE(sz_r), mpi_double_precision,mpi_sum, mpi_comm_world, ierr)



call MPI_ALLREDUCE(eigvals_per_k,eigvals_per_k_mpi,size(eigvals_per_k),MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
call mpi_reduce(dos_l, dos_l_mpi, size(dos_l),MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
call mpi_reduce(dos_r, dos_r_mpi, size(dos_r),MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
call mpi_reduce(dos_bulk, dos_bulk_mpi, size(dos_bulk),MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)

! dos_l_mpi= dos_l
! dos_r_mpi= dos_r
! dos_bulk_mpi= dos_bulk

dos_l=log(abs(dos_l_mpi))
dos_r=log(abs(dos_r_mpi))
dos_bulk=log(abs(dos_bulk_mpi)+0.000000001)
write(*,*) "dos_r ok"

do ik=1, numkpts
    do j=1, omeganum
        dos_l_only(ik, j)= dos_l_mpi(ik, j)- dos_bulk_mpi(ik, j)
        if (dos_l_only(ik, j)<0) dos_l_only(ik, j)=0.0000000001
        dos_r_only(ik, j)= dos_r_mpi(ik, j)- dos_bulk_mpi(ik, j)
        if (dos_r_only(ik, j)<0) dos_r_only(ik, j)=0.0000000001
    enddo
enddo
write(*,*) "dos_r_only ok"

if (irank.eq.0)then
    open(369, file='dos.dat_l')
    open(370, file='dos.dat_r')
    open(371, file='dos.dat_bulk')
    open(372, file='spindos.dat_l')
    open(373, file='spindos.dat_r')
    do ik=1, numkpts
        do j=1, omeganum
            WRITE(369, '(30f16.8)')k(ik), omega(j), dos_l(ik, j), log(dos_l_only(ik, j))
            WRITE(370, '(30f16.8)')k(ik), omega(j), dos_r(ik, j), log(dos_r_only(ik, j))
            WRITE(371, '(30f16.8)')k(ik), omega(j), dos_bulk(ik, j)

            ! s0(1)=sx_l(ik, j); s0(2)=sy_l(ik, j); s0(3)=sz_l(ik, j); 
            ! call rotate(s0, s1)
            write(372, '(30f16.8)')k(ik),omega(j),sx_l(ik, j),sy_l(ik, j),sz_l(ik, j); 
            ! s0(1)=sx_r(ikp, j); s0(2)=sy_r(ikp, j); s0(3)=sz_r(ikp, j); 
            ! call rotate(s0, s1)
            write(373, '(30f16.8)')k(ik),omega(j),sx_r(ik, j),sy_r(ik, j),sz_r(ik, j); 

        enddo
        WRITE(369, *) ' '
        WRITE(370, *) ' '
        WRITE(371, *) ' '
        WRITE(372, *) ' '
        WRITE(373, *) ' '
    enddo
    CLOSE(369)
    CLOSE(370)
    CLOSE(371)
    CLOSE(372)
    CLOSE(373)

    WRITE(*,*)'ndim',ndim
    WRITE(*,*) 'numkpts,omeganum,eta',numkpts, omeganum,eta
    WRITE(*,*)'calculate density of state successfully'
endif

emin= minval(omega)
emax= maxval(omega)

write(*,*) "emin=",emin,"emax=",emax

!> write script for gnuplot
if (irank==0) then
    open(374, file='surfdos_l.gnu')
    write(374, '(a)')"set encoding iso_8859_1"
    write(374, '(a)')'#set terminal  postscript enhanced color'
    write(374, '(a)')"#set output 'surfdos_l.eps'"
    write(374, '(3a)')'#set terminal  pngcairo truecolor enhanced', &
       '  font ", 60" size 1920, 1680'
    write(374, '(3a)')'set terminal  png truecolor enhanced', &
       ' font ", 60" size 1920, 1680'
    write(374, '(a)')"set output 'surfdos_l.png'"
    write(374,'(2a)') 'set palette defined (-10 "#194eff", ', &
       '0 "white", 10 "red" )'
    write(374, '(a)')'#set palette rgbformulae 33,13,10'
    write(374, '(a)')'set style data linespoints'
    write(374, '(a)')'set size 0.8, 1'
    write(374, '(a)')'set origin 0.1, 0'
    write(374, '(a)')'unset ztics'
    write(374, '(a)')'unset key'
    write(374, '(a)')'set pointsize 0.8'
    write(374, '(a)')'set pm3d'
    write(374, '(a)')'#set view equal xyz'
    write(374, '(a)')'set view map'
    write(374, '(a)')'set border lw 3'
    write(374, '(a)')'#set cbtics font ",48"'
    write(374, '(a)')'#set xtics font ",48"'
    write(374, '(a)')'#set ytics font ",48"'
    write(374, '(a)')'#set ylabel font ",48"'
    write(374, '(a)')'set ylabel "Energy (eV)"'
    write(374, '(a)')'#set xtics offset 0, -1'
    write(374, '(a)')'#set ylabel offset -6, 0 '
    write(374, '(a, f18.5, a, f18.5, a)')'set yrange [', emin, ':', emax, ']'
    write(374, '(a)')'set pm3d interpolate 2,2'
    write(374, '(2a)')"splot 'dos.dat_l' u 1:2:3 w pm3d"

    write(374, '(3a)')'set terminal png truecolor enhanced',' font ", 30" size 1920, 1680'
    write(374, '(a)')"set output 'spindos_l.png'"
    write(374, '(a)')"set multiplot layout 3, 1"
    write(374, '(a)')"set title 'sx'"
    write(374, '(2a)')"splot 'spindos.dat_l' u 1:2:3 w pm3d "
    write(374, '(a)')"set title 'sy'"
    write(374, '(2a)')"splot 'spindos.dat_l' u 1:2:4 w pm3d"
    write(374, '(a)')"set title 'sz'"
    write(374, '(2a)')"splot 'spindos.dat_l' u 1:2:5 w pm3d"

    CLOSE(374)

endif

    if(irank.eq.0)then
        open(222,file='kpts.out',recl=10000)
        do ik=1,numkpts
            write(222,*),k(ik)
        enddo
        close(222) 
       
        open(234,file='layeri',recl=10000)
        do i=1,Hdim
            write(234,*),layerintarr(i)
        enddo
        close(234) 
  !      
    !    open(123,file='dos_r',recl=10000)
    !         write(123,*), dos_r_mpi
    !   close(123)
      
      open(777,file='output_bands',recl=10000)
      do ib=1,Hdim
         do ik=1,numkpts
              write(777,*), k(ik),eigvals_per_k_mpi(ik,ib)
          enddo
      enddo
      close(777)
    !  enddo
    endif

    call mpi_barrier(mpi_comm_world,ierr)
    call MPI_Finalize(ierr)
end program hao_edgestates
subroutine now(time_now)

    implicit none
    integer   :: time_new(8)
    real      :: time_now
    call Date_and_time(values=time_new)
    time_now= time_new(3)*24*3600+time_new(5)*3600+&
              time_new(6)*60+time_new(7)+time_new(8)/1000d0  
    return
 end subroutine now

 subroutine surfgreen_1985(Ndim,eta,omega,GLL,GRR,GB,H00,H01,ones)
    implicit none

    ! inout variables     
    ! the factor 2 is induced by spin
    ! energy hbar omega
    real(kind(1.0d0)),intent(in) :: omega
    ! real(kind(1.0d0)),intent(in) :: omegamax,omegamin,omeganum  
    integer, intent(in) :: Ndim
    ! H00 Hamiltonian between nearest neighbour-quintuple-layers
    complex(kind(1.0d0)),intent(in) :: H00(Ndim,Ndim)

    ! H01 Hamiltonian between next-nearest neighbour-quintuple-layers
    complex(kind(1.0d0)),intent(in) :: H01(Ndim,Ndim)

    ! temp hamiltonian

    complex(kind(1.0d0)),intent(in)   :: ones(Ndim,Ndim)

    ! surface green function
    complex(kind(1.0d0)),intent(inout)  :: GLL(Ndim,Ndim)
    complex(kind(1.0d0)),intent(inout)  :: GRR(Ndim,Ndim)

    !> bulk green's function
    complex(kind(1.0d0)),intent(inout)  :: GB(Ndim,Ndim)

    ! >> local variables
    ! iteration number
    integer :: iter

    ! maximun iteration 
    integer ,parameter:: itermax=100

    ! accuracy control
    real(kind(1.0d0)) :: accuracy=1e-16

    ! a real type temp variable
    real(kind(1.0d0)) :: real_temp
    real(kind(1.0d0)) :: eta
    ! omegac=omega(i)+I * eta
    complex(kind(1.0d0)) :: omegac 


    ! some variables in Eq.(11)
    complex(kind(1.0d0)), allocatable :: alphai(:, :) 
    complex(kind(1.0d0)), allocatable :: betai(:, :) 
    complex(kind(1.0d0)), allocatable :: epsiloni(:, :) 
    complex(kind(1.0d0)), allocatable :: epsilons(:, :) 
    complex(kind(1.0d0)), allocatable :: epsilons_t(:, :) 

    complex(kind(1.0d0)), allocatable :: mat1 (:, :) 
    complex(kind(1.0d0)), allocatable :: mat2 (:, :) 

    ! g0= inv(w-e_i)
    complex(kind(1.0d0)), allocatable :: g0 (:, :) 

    ! eta = (omegamax- omegamin)/omeganum*2d0
    ! allocate some variables
    allocate(alphai(Ndim, Ndim)) 
    allocate(betai (Ndim, Ndim)) 
    allocate(epsiloni (Ndim, Ndim)) 
    allocate(epsilons (Ndim, Ndim)) 
    allocate(epsilons_t(Ndim, Ndim)) 
    allocate(mat1(Ndim, Ndim)) 
    allocate(mat2(Ndim, Ndim)) 
    allocate(g0(Ndim, Ndim)) 

    epsiloni= H00
    epsilons= H00
    epsilons_t= H00
    alphai  = H01
    betai   = conjg(transpose(H01))
   !print *, sqrt(sum(abs(H00)**2)), 'H00'

    ! w+i*0^+
    omegac= cmplx(omega, eta)
   !print *, omegac

    ! begin iteration
    do iter=1, itermax

       g0= omegac*ones- epsiloni
       call inv(Ndim, g0)

       ! a_i-1*(w-e_i-1)^-1
       call mat_mul(Ndim, alphai, g0, mat1 )
       
       ! b_i-1*(w-e_i-1)^-1
       call mat_mul(Ndim, betai, g0, mat2 )

       ! a_i-1*(w-e_i-1)^-1*b_i-1
       call mat_mul(Ndim, mat1, betai, g0)
       epsiloni= epsiloni+ g0
      !print *, sqrt(sum(abs(epsiloni)**2)), 'ei'
       ! es_i= es_i-1 + a_i-1*(w-e_i-1)^-1*b_i-1
       epsilons= epsilons+ g0
      !print *, sqrt(sum(abs(epsilons)**2)), 'es'
      !pause

       ! b_i-1*(w-e_i-1)^-1*a_i-1
       call mat_mul(Ndim, mat2, alphai, g0)
       epsiloni= epsiloni+ g0
       ! es_i= es_i-1 + a_i-1*(w-e_i-1)^-1*b_i-1
       epsilons_t= epsilons_t+ g0

       ! a_i= a_i-1*(w-e_i-1)^-1*a_i-1 
       call mat_mul(Ndim, mat1, alphai, g0)
       alphai= g0
       ! b_i= b_i-1*(w-e_i-1)^-1*b_i-1 
       call mat_mul(Ndim, mat2, betai, g0)
       betai= g0

      !real_temp=maxval(abs(alphai))   
       real_temp=sum(abs(alphai))   
      !if (cpuid.eq.0) print *, iter, real_temp
       if (real_temp.le.accuracy) exit

    enddo ! end of iteration

    ! calculate surface green's function
    GLL= omegac*ones- epsilons
    call inv(Ndim, GLL)

    GRR= omegac*ones- epsilons_t
    call inv(Ndim, GRR)

    GB = omegac*ones- epsiloni
    call inv(Ndim, GB)

    return
 end subroutine surfgreen_1985



subroutine inv(ndim,Amat)

    implicit none

    integer,parameter :: dp=8
    integer           :: i
    integer           :: info
    integer,intent(in):: ndim
    integer,allocatable   :: ipiv(:)
    complex(kind(1.0d0)),parameter :: zone=(1.0d0,0.0d0)
    complex(kind(1.0d0)),intent(inout):: Amat(ndim,ndim)
    complex(kind(1.0d0)),allocatable :: Bmat(:,:)

    allocate(ipiv(ndim))
    allocate(Bmat(ndim,ndim))

    ipiv=0

    Bmat= (0d0, 0d0)
    do i=1,ndim
       Bmat(i,i)= zone
    enddo

    call zgesv(ndim,ndim,Amat,ndim,ipiv,Bmat,ndim,info)

    if(info.ne.0)print *,'something wrong with zgesv'

    Amat=Bmat
    
    return

end subroutine inv 



subroutine mat_mul(nmatdim,A,B,C)
      
    implicit none
    
    integer,intent(in) :: nmatdim    
    complex(kind(1.0d0)) :: ALPHA
    complex(kind(1.0d0)) :: BETA 
    complex(kind(1.0d0)), intent(in)  :: A(nmatdim ,nmatdim)
    complex(kind(1.0d0)), intent(in)  :: B(nmatdim ,nmatdim)
    !complex(kind(1.0d0)) :: mat_mul(nmatdim,nmatdim)
    complex(kind(1.0d0)), intent(out) :: C(nmatdim,nmatdim)
 
    ALPHA=1.0d0 
    BETA=0.0D0

    C(:,:)=(0.0d0,0.0d0)
 
    call ZGEMM('N','N',nmatdim,nmatdim,nmatdim,ALPHA,A,nmatdim,B,nmatdim,BETA,C,nmatdim)
    return

end subroutine mat_mul


! subroutine rotate(R1, R2)
!     use para, only : dp, Urot
!     implicit none
!     real(dp), intent(in) :: R1(3)
!     real(dp), intent(inout) :: R2(3)
 
!     R2(1)= Urot(1, 1)*R1(1)+ Urot(1, 2)*R1(2)+ Urot(1, 3)*R1(3)
!     R2(2)= Urot(2, 1)*R1(1)+ Urot(2, 2)*R1(2)+ Urot(2, 3)*R1(3)
!     R2(3)= Urot(3, 1)*R1(1)+ Urot(3, 2)*R1(2)+ Urot(3, 3)*R1(3)
 
!     return
!  end subroutine rotate