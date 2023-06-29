program magvec_rspauli
! mpiifort -CB -r8 magvec_rspauli.f90 -qmkl -o magvec_rspauli.x
    implicit none

    complex,allocatable:: hops(:,:,:),hops_mag_temp(:,:,:,:,:)
    complex,allocatable:: rspauli(:,:,:,:)
    complex,allocatable:: pauli(:,:,:)
    real               :: rdum,idum
    integer            :: ix,iy,iz,band1,band2,num_wann,num_lines
    integer            :: ik1,ik2,ik3,kp1,kp2,kp3
    real               :: twopi
    real               :: phas
    complex            :: fac,fac2
    complex,allocatable:: ham(:,:) 
    real,allocatable   :: eigvals(:)
    complex,allocatable:: eigvecs(:,:)
    real               :: kder,amat(3,3),cross(3)
    real               :: fermi_min,fermi_max
    real               :: occupation_number
    integer,allocatable:: irvec(:,:),MAG_WANN_ORBS_INDEX(:)
    integer,allocatable:: nrpts(:)
    real               :: bohrincm,condq,abstol,pi,volume
    integer            :: i,j,ii,num_steps,length,length2,ne,length4
    integer            :: rvecnum
    integer            :: grid
    integer            :: maxdim
    integer            :: ierr,isize,irank
    real,allocatable   :: conductivity12(:),conductivity_temp(:)
    real,allocatable   :: conductivity13(:),conductivity23(:)
    real,allocatable   :: conductivity12_ahe(:)
    real,allocatable   :: conductivity13_ahe(:),conductivity23_ahe(:)
    real,allocatable   :: fermienergy(:)
    real,allocatable   :: occupation(:),occupation_temp(:)
    complex,allocatable:: momentum(:,:,:)
    complex,allocatable:: momentum2(:,:,:)
    complex,allocatable:: spinmomentum(:,:,:)
    complex,allocatable:: spinmomentum2(:,:,:)
    real               :: kpoints(3)
    real               :: vl,vu
    integer            :: info
    complex,allocatable:: work(:)
    integer            :: lwork
    integer,allocatable:: iwork(:)
    real,allocatable   :: rwork(:)
    integer,allocatable:: ifail(:)
    integer            :: n1,n2,n3,n4,dir,step,orb,num_occ
    character(len=8)   :: tempt
    integer            :: theta_num,theta_i,phi_j,cycle_count,mag_wann_num
    real               :: mag_strength,mag_theta,mag_phi,mag_field_x,mag_field_y,mag_field_z
    real,allocatable   :: theta(:)
    real               :: phi(8)
    real, allocatable  :: mag_field(:,:,:)
    logical            :: l_mag_vec
    ! real, allocatable :: hops_in_field(:,:,:,:,:)

    INCLUDE 'mpif.h'
    integer stt(MPI_STATUS_SIZE) 
    CALL MPI_INIT(ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,isize,ierr)

    abstol    = 2.0*tiny(abstol)
    twopi     = 2*3.141592654
    bohrincm  = 0.529177*1.e-8
    condq     = 38.74*1.e-6
    pi        = 3.141592654

    ! theta_num = 2
    ! mag_strength = 1.0
    ! allocate(theta(theta_num))
    ! allocate(mag_field(3,8,theta_num))
    ! if(irank.eq.0)then
    !     open(555,file='angle',recl=10000)
    !     open(666,file='magfield',recl=10000)
    !     theta = 0.0
    !     do i=1,theta_num
    !         theta(i) =  (i-1) * (pi/2)/theta_num
    !         write(555,*),"theta=",theta(i)
    !     end do
    !     do j=1,8
    !         phi(j) = (j-1) * (pi/4)
    !         write(555,*),"phi=",phi(j)
    !     end do
    !     do theta_i = 1,theta_num
    !         do phi_j = 1,8
    !             mag_field(1,phi_j,theta_i) = mag_strength * cos(theta(theta_i))
    !             mag_field(2,phi_j,theta_i) = mag_strength * sin(theta(theta_i)) * cos(phi(phi_j))
    !             mag_field(3,phi_j,theta_i) = mag_strength * sin(theta(theta_i)) * sin(phi(phi_j))
    !             write(666,*),"mag=",mag_field(1,phi_j,theta_i),mag_field(2,phi_j,theta_i),mag_field(3,phi_j,theta_i)
    !         enddo
    !     enddo
    !     close(555)
    !     close(666)
    ! endif
!!

if(irank == 0)then
    open(100,file='ahe_inp')
    read(100,*) amat(1,:)
    read(100,*) amat(2,:)
    read(100,*) amat(3,:)
    read(100,*) fermi_min,fermi_max,num_steps
    read(100,*) grid
    read(100,*) maxdim
    read(100,*) occupation_number
    read(100,*) l_mag_vec                        
    read(100,*) mag_strength                     
    read(100,*) mag_theta                        
    read(100,*) mag_phi                          
    read(100,*) mag_wann_num                     
    allocate(mag_wann_orbs_index(mag_wann_num))  
    read(100,*) mag_wann_orbs_index(mag_wann_num)
    close(100)
endif

    call mpi_bcast(amat,9,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)
    call mpi_bcast(fermi_min,1,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)
    call mpi_bcast(occupation_number,1,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)
    call mpi_bcast(fermi_max,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
    call mpi_bcast(num_steps,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    call mpi_bcast(maxdim,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    call mpi_bcast(grid,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    call mpi_bcast(l_mag_vec,1,MPI_LOGICAL,0,mpi_comm_world,ierr)
    call mpi_bcast(mag_strength,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
    call mpi_bcast(mag_theta,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
    call mpi_bcast(mag_phi,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
    call mpi_bcast(mag_wann_num,1,MPI_INTEGER, 0,mpi_comm_world,ierr)
    call mpi_bcast(mag_wann_orbs_index,mag_wann_num,MPI_INTEGER,0,mpi_comm_world,ierr)
    
    if(irank == 0)then  
        mag_field_x = mag_strength * cos(mag_theta)
        mag_field_y = mag_strength * sin(mag_theta) * cos(mag_phi)
        mag_field_z = mag_strength * sin(mag_theta) * sin(mag_phi)

        print *, mag_field_x, mag_field_y, mag_field_z
        print *, mag_wann_orbs_index
    endif

    call mpi_bcast(mag_field_x,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
    call mpi_bcast(mag_field_y,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
    call mpi_bcast(mag_field_z,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)

    open(200,file='hopping.1')
    num_lines=0
    num_wann=0
    do
        read(200,fmt=*,end=311)ix,iy,iz,band1,band2,rdum,idum
        num_lines=num_lines+1
        num_wann=max(num_wann,band1)
    enddo
311 continue
    rvecnum=num_lines/(num_wann*num_wann)
    allocate(hops(1:num_wann,1:num_wann,rvecnum))
    allocate(irvec(3,rvecnum))
    hops=0.0
    rewind(200)
    num_lines=0
    do
        read(200,fmt=*,end=300)ix,iy,iz,band1,band2,rdum,idum
        num_lines=num_lines+1
        rvecnum=(num_lines-1)/(num_wann*num_wann)+1
        irvec(1,rvecnum)=ix
        irvec(2,rvecnum)=iy
        irvec(3,rvecnum)=iz
        hops( band1,band2,rvecnum )=cmplx(rdum,idum)
    enddo
300 continue
    close(200)

    allocate(rspauli(1:num_wann, 1:num_wann, 3, rvecnum))
    open(400,file='./rspauli.1')
    num_lines=0
    Do
       read(400, fmt=*,end=500) ix,iy,iz,band1,band2,dir,rdum,idum
       num_lines=num_lines+1
       rvecnum=(num_lines-1)/(num_wann*num_wann*3)+1
       rspauli(band1, band2, dir, rvecnum)=cmplx(rdum,idum)
    End Do
500 continue

    close(400)
    allocate(nrpts(rvecnum))
    open(14,file='nrpts_inp')
    do j=1,rvecnum/15
        read(14,'(15I5)') (nrpts(15*(j-1)+i) ,i=1,15)
    enddo
    read(14,'(<mod(rvecnum,15)>I5)') (nrpts(15*(rvecnum/15)+i),i=1,mod(rvecnum,15))
    close(14)

    cross(1)=amat(1,2)*amat(2,3)-amat(1,3)*amat(2,2)
    cross(2)=amat(1,3)*amat(2,1)-amat(1,1)*amat(2,3)
    cross(3)=amat(1,1)*amat(2,2)-amat(1,2)*amat(2,1)

    volume=cross(1)*amat(3,1)+cross(2)*amat(3,2)+cross(3)*amat(3,3)
    if (maxdim.eq.2) then
        volume=volume/amat(3,3)/(0.529177e-8)
    endif   

    allocate(occupation(num_steps))
    allocate(occupation_temp(num_steps))
    allocate(conductivity12(num_steps))
    allocate(conductivity_temp(num_steps))
    allocate(conductivity13(num_steps))
    allocate(conductivity23(num_steps))
    allocate(conductivity12_ahe(num_steps))
    allocate(conductivity13_ahe(num_steps))
    allocate(conductivity23_ahe(num_steps))
    allocate(fermienergy(num_steps))

    occupation=0.0
    conductivity12=0.0
    conductivity13=0.0
    conductivity23=0.0
    conductivity12_ahe=0.0
    conductivity13_ahe=0.0
    conductivity23_ahe=0.0

    allocate(ham(num_wann,num_wann))
    allocate(pauli(num_wann,num_wann,3))
    allocate(momentum(num_wann,num_wann,3))
    allocate(momentum2(num_wann,num_wann,3))
    allocate(spinmomentum (num_wann,num_wann,3))
    allocate(spinmomentum2(num_wann,num_wann,3))

    allocate(eigvals(num_wann))
    allocate(eigvecs(num_wann,num_wann))
    lwork=12.0*num_wann
    allocate(work(lwork))
    allocate(rwork(17*num_wann))
    allocate(iwork(15*num_wann))
    allocate(ifail(15*num_wann))

    call mpi_bcast(num_wann,1,MPI_INTEGER,0,mpi_comm_world,ierr)
    call mpi_bcast(rvecnum,1,MPI_INTEGER,0,mpi_comm_world,ierr)

    if(.not.allocated(nrpts))allocate(nrpts(rvecnum))
    call mpi_bcast(nrpts,rvecnum,MPI_INTEGER,0,mpi_comm_world,ierr)

    if(.not.allocated(irvec))allocate(irvec(3,rvecnum))
    length=3*rvecnum
    call mpi_bcast(irvec,length,MPI_INTEGER,0,mpi_comm_world,ierr)

    if(.not.allocated(hops))allocate(hops(num_wann,num_wann,rvecnum))
    length2=num_wann*num_wann*rvecnum
    call mpi_bcast(hops,length2,mpi_complex,0,mpi_comm_world,ierr)

    ! if(.not.allocated(hops_mag_temp))allocate(hops_mag_temp(theta_num,8,num_wann,num_wann,rvecnum))
    ! length4 = theta_num*8*num_wann*num_wann*rvecnum
    ! call mpi_bcast(hops,length4,mpi_complex,0,mpi_comm_world,ierr)

    ! hops_mag_temp=cmplx(0.d0,0.d0)
    ! if(irank.eq.0)then
    !     do theta_i = 1,theta_num
    !         do phi_j = 1,8
    !             do ii=1,rvecnum
    !                 do j=1,num_wann
    !                     do i=1,num_wann
    !                         hops_mag_temp(theta_i,phi_j,j,i,ii)=hops(j,i,ii)
    !                     enddo
    !                 enddo
    !             enddo
    !         enddo
    !     enddo
    ! endif
    ! if(irank.eq.0)then
    !     open(886,file='hops_mag_temp1',recl=10000)
    !     write(886,*) hops_mag_temp
    !     close(886)
    ! endif
    ! if(irank.eq.0)then
    !     ! do theta_i = 1,theta_num
    !     !     do phi_j = 1,8
    !             do ii=1,rvecnum
    !                 do j=1,num_wann
    !                     do i=1,num_wann
    !                     ! if ((j .gt. 2) .and. (j .lt. 4) .and. (i .gt. 2) .and. (i .lt. 4)) then
    !                         hops_mag_temp(theta_i,phi_j,j,i,ii)=hops_mag_temp(theta_i,phi_j,j,i,ii)+0.5d0*(rspauli(j,i,1,ii)*mag_field(1,phi_j,theta_i)&
    !                         +rspauli(j,i,2,ii)*mag_field(2,phi_j,theta_i)+rspauli(j,i,3,ii)*mag_field(3,phi_j,theta_i))
    !                     ! endif
    !     	            enddo
    ! 	            enddo
    !         !     enddo
    !         ! enddo
    !     enddo 
    ! endif
! if(irank.eq.0)then
    if(l_mag_vec)then
        do ii=1,rvecnum
            do i=1,num_wann
                do j=1,num_wann
                    ! if(ANY(i==mag_wann_orbs_index) .and. ANY(j==mag_wann_orbs_index))then
                        hops(j,i,ii) = hops(j,i,ii) + mag_field_x * rspauli(j,i,1,ii)                  &
                        + mag_field_y * rspauli(j,i,2,ii) + mag_field_z * rspauli(j,i,3,ii)
                    ! endif
                enddo !j
            enddo !i
        enddo
    endif
! endif
    ! if(irank.eq.0)then
    !     open(888,file='hops_mag_temp2',recl=10000)
    !     write(888,*) hops_mag_temp
    !     close(888)
    ! endif

! cycle_count = 1000
! do theta_i = 1,theta_num
!     print*, "theta_i=", theta_i
!     do phi_j = 1,8
!         print*, "phj=", theta_i
!         cycle_count = cycle_count + 1
        
    ik1=0
    do kp1=0,grid-1
        kpoints(1)=-0.5+real(kp1)/real(grid)
            do kp2=0,grid-1
                kpoints(2)=-0.5+real(kp2)/real(grid)
                do kp3=0,(grid-1)*(maxdim-2)
                    ik1=ik1+1
                    if(mod(ik1-1,isize).ne.irank)cycle
                        ! kpoints(1)=-0.5+real(kp1)/real(grid)
                        ! kpoints(2)=-0.5+real(kp2)/real(grid)
                    kpoints(3)=-0.5+real(kp3)/real(grid)
                    ham=0.0  
                    pauli=cmplx(0.d0,0.d0)
                    do ii=1,rvecnum
                        ix=irvec(1,ii)
                        iy=irvec(2,ii)
                        iz=irvec(3,ii)
                        phas=     iy*kpoints(2)
                        phas=phas+iz*kpoints(3)
                        phas=phas+ix*kpoints(1)
                        phas=phas*twopi
                        fac=cmplx(cos(phas),sin(phas))
                        do i=1,num_wann
                            do j=1,num_wann
                                ! ham(j,i)=ham(j,i)+fac*hops(j,i,ii)/nrpts(ii)
                                ham(j,i)=ham(j,i)+fac*hops(j,i,ii)/nrpts(ii)
                                do dir=1,3
                                    pauli(j,i,dir)=pauli(j,i,dir)+fac*rspauli(j,i,dir,ii)
                                enddo
                            enddo
                        enddo 
                    enddo
                    call zheevx('V','A','U',num_wann,ham,num_wann,&
                                vl,vu,1,num_wann,abstol,ne,eigvals,eigvecs,num_wann,&
                                work,lwork,rwork,iwork,ifail,info)
                    if(info.ne.0)stop 'zheevx'

                    momentum=0.0
                    do ii=1,rvecnum
                        ix=irvec(1,ii)
                        iy=irvec(2,ii)
                        iz=irvec(3,ii)
                        phas=     iy*kpoints(2)
                        phas=phas+iz*kpoints(3)
                        phas=phas+ix*kpoints(1)
                        phas=phas*twopi
                        fac=cmplx(-sin(phas),cos(phas))
                        do dir=1,3
                            kder=amat(dir,1)*ix+amat(dir,2)*iy+amat(dir,3)*iz
                            fac2=fac*kder                  
                            do n2=1,num_wann
                                do n1=1,num_wann
                                    ! momentum(n1,n2,dir)= momentum(n1,n2,dir)+fac2*hops(n1,n2,ii)
                                    momentum(n1,n2,dir)= momentum(n1,n2,dir)+fac2*hops(n1,n2,ii)
                                enddo
                            enddo
                        enddo
                    enddo
                    spinmomentum=0.0
                    if (.true.) then
                        do dir=1,3
                            spinmomentum(:,:,dir)= (MATMUL(pauli(:,:,3),momentum(:,:,dir))+MATMUL(momentum(:,:,dir),pauli(:,:,3)))/2.d0
                        enddo
                    else
                        do ii=1,rvecnum
                            ix=irvec(1,ii)
                            iy=irvec(2,ii)
                            iz=irvec(3,ii)
                            phas=     iy*kpoints(2)
                            phas=phas+iz*kpoints(3)
                            phas=phas+ix*kpoints(1)
                            phas=phas*twopi
                            fac=cmplx(-sin(phas),cos(phas))
                            do dir=1,3
                                kder=amat(dir,1)*ix+amat(dir,2)*iy+amat(dir,3)*iz
                                fac2=fac*kder                  
                                do n2=1,num_wann
                                    do n1=1,num_wann                   
                                        ! spinmomentum(n1,n2,dir)= spinmomentum(n1,n2,dir)+fac2*hops(n1,n2,ii)
                                        spinmomentum(n1,n2,dir)= spinmomentum(n1,n2,dir)+fac2*hops(n1,n2,ii)
                                    enddo
                                enddo
                            enddo
                        enddo
                    endif
                    momentum2=0.0
                    do dir=1,3
                        do n2=1,num_wann
                            do n4=1,num_wann
                                do n1=1,num_wann
                                    do n3=1,num_wann
                                        momentum2(n4,n2,dir)=momentum2(n4,n2,dir)+ momentum(n3,n1,dir)*eigvecs(n1,n2)*conjg(eigvecs(n3,n4))
                                    enddo !n4
                                enddo !n3
                            enddo !n2
                        enddo !n1
                    enddo  !dir
                    spinmomentum2=0.0
                    do dir=1,3
                        do n2=1,num_wann
                            do n4=1,num_wann
                                do n1=1,num_wann
                                    do n3=1,num_wann
                                        spinmomentum2(n4,n2,dir)=spinmomentum2(n4,n2,dir)+spinmomentum(n3,n1,dir)*eigvecs(n1,n2)*conjg(eigvecs(n3,n4))
                                    enddo
                                enddo
                            enddo
                        enddo
                    enddo

                    do step=1,num_steps
                        fermienergy(step) = fermi_min+(fermi_max-fermi_min)*real(step)/real(num_steps)
                        iy=0
                        do ix=1,num_wann
                            if(eigvals(ix).le.fermienergy(step))then
                                iy=iy+1
                            endif
                        enddo
                        num_occ=iy
                        occupation(step)=occupation(step)+num_occ
  
                        do ik2=num_occ+1,num_wann
                            do orb=1,num_occ
                                if (abs(eigvals(orb)-eigvals(ik2)) .gt. 0.0000001d0) then
                                    conductivity12(step)=conductivity12(step)-aimag(momentum2(orb,ik2,1)*conjg(spinmomentum2(orb,ik2,2)))/(eigvals(orb)-eigvals(ik2))**2
                                    conductivity13(step)=conductivity13(step)-aimag(momentum2(orb,ik2,1)*conjg(spinmomentum2(orb,ik2,3)))/(eigvals(orb)-eigvals(ik2))**2
                                    conductivity23(step)=conductivity23(step)-aimag(momentum2(orb,ik2,2)*conjg(spinmomentum2(orb,ik2,3)))/(eigvals(orb)-eigvals(ik2))**2
                                end if
                            enddo
                        enddo 
                        do ik2=num_occ+1,num_wann
                            do orb=1,num_occ
                                conductivity12_ahe(step)=conductivity12_ahe(step)-aimag(momentum2(orb,ik2,1)*conjg(momentum2(orb,ik2,2)))/(eigvals(orb)-eigvals(ik2))**2
                                conductivity13_ahe(step)=conductivity13_ahe(step)-aimag(momentum2(orb,ik2,1)*conjg(momentum2(orb,ik2,3)))/(eigvals(orb)-eigvals(ik2))**2
                                conductivity23_ahe(step)=conductivity23_ahe(step)-aimag(momentum2(orb,ik2,2)*conjg(momentum2(orb,ik2,3)))/(eigvals(orb)-eigvals(ik2))**2
                            enddo
                        enddo
                    enddo
            enddo
        enddo
    enddo

    conductivity12_ahe=conductivity12_ahe/grid**maxdim/volume/bohrincm*twopi*condq*2.0
    conductivity13_ahe=conductivity13_ahe/grid**maxdim/volume/bohrincm*twopi*condq*2.0
    conductivity23_ahe=conductivity23_ahe/grid**maxdim/volume/bohrincm*twopi*condq*2.0
        
    conductivity12=conductivity12/grid**maxdim/volume/bohrincm*twopi*condq*2.0
    conductivity13=conductivity13/grid**maxdim/volume/bohrincm*twopi*condq*2.0
    conductivity23=conductivity23/grid**maxdim/volume/bohrincm*twopi*condq*2.0
        
    occupation=occupation/grid**maxdim
        
    call MPI_REDUCE(occupation,occupation_temp,num_steps,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)      
    occupation=occupation_temp
        
    call MPI_REDUCE(conductivity12,conductivity_temp,num_steps,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)      
    conductivity12=conductivity_temp
        
    call MPI_REDUCE(conductivity13,conductivity_temp,num_steps,MPI_DOUBLE_PRECISION,MPI_SUM,0, mpi_comm_world,ierr)      
    conductivity13=conductivity_temp
        
    call MPI_REDUCE(conductivity23,conductivity_temp,num_steps,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)      
    conductivity23=conductivity_temp
        
    call MPI_REDUCE(conductivity12_ahe,conductivity_temp,num_steps,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)      
    conductivity12_ahe=conductivity_temp
        
    call MPI_REDUCE(conductivity13_ahe,conductivity_temp,num_steps,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)      
    conductivity13_ahe=conductivity_temp
        
    call MPI_REDUCE(conductivity23_ahe,conductivity_temp,num_steps,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)      
    conductivity23_ahe=conductivity_temp
        
    if(irank.eq.0)then
        ! print*, 'cycle_count=', cycle_count
        ! write(tempt,'(i8)') cycle_count
        open(123,file='output_ahe_condquant',recl=10000)
        do step=1,num_steps
            write(123,*)"fermienergy=",fermienergy(step)
            write(123,*)"occupation=",occupation(step)
            write(123,*)"conductivity12=",conductivity12_ahe(step)/(0.5*77.48e-6)
            write(123,*)"conductivity13=",conductivity13_ahe(step)/(0.5*77.48e-6)
            write(123,*)"conductivity23=",conductivity23_ahe(step)/(0.5*77.48e-6)
            write(123,*)"************************************"
        enddo !step
        close(123)

        ! write(tempt,'(i8)') cycle_count
        open(123,file='output_she_condquant',recl=10000)
        do step=1,num_steps
            write(123,*)"fermienergy=",fermienergy(step)
            write(123,*)"occupation=",occupation(step)
            write(123,*)"conductivity12=",conductivity12(step)/(0.5*77.48e-6)
            write(123,*)"conductivity13=",conductivity13(step)/(0.5*77.48e-6)
            write(123,*)"conductivity23=",conductivity23(step)/(0.5*77.48e-6)
            write(123,*)"************************************"
        enddo !step
        close(123)
    endif

call mpi_barrier(mpi_comm_world,ierr)
call MPI_Finalize(ierr)

end program magvec_rspauli