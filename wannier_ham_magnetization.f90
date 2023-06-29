program wannier_ham_magnetization

    implicit none
                                                    
    complex,allocatable :: hops(:,:,:) 
    complex,allocatable :: rsnabla(:,:,:,:) 
    complex,allocatable :: rspauli(:,:,:,:) 
    complex,allocatable :: pauli(:,:,:) 
    complex,allocatable :: paulifft(:,:,:) 
    complex,allocatable :: paulifft2(:,:,:) 
    real                :: rdum,idum 
    integer             :: ix,iy,iz,band1,band2,h,num_wann,m,n
    integer             :: ik1,ik2,ik3 
    real                :: twopi 
    real                :: phas 
    complex             :: fac,fac2 
    complex,allocatable :: ham(:,:) 
    real                :: vl,vu 
    integer             :: ne,j 
    real                :: abstol 
    real,allocatable    :: eigvals(:) 
    complex,allocatable :: eigvecs(:,:) 
    integer             :: info 
    complex,allocatable :: work(:) 
    integer             :: lwork 
    integer,allocatable :: iwork(:) 
    real,allocatable    :: rwork(:) 
    integer,allocatable :: ifail(:) 
    real                :: kpoints(3)
    real                :: scale ,Beta_fake,mu
    integer             :: maxhopx2,maxhopy2,maxhopz2,dire 
    real,allocatable    :: fermienergy(:) 
    real,allocatable    :: deviation(:) 
    integer             :: grid,i1,i2,i3,i4,orb,Nk1,Nk2,Nk3,knv3,ikx,iky,ikz,ik
    real,allocatable    :: conductivity(:),conductivity2(:),sigma_tensor_ahc_x(:)
    real,allocatable    :: conductivity13(:),conductivity23(:) 
    real,allocatable    :: conductivity_ahe(:) ,sigma_tensor_ahc_y(:),sigma_tensor_ahc_z(:)
    real,allocatable    :: conductivity13_ahe(:),conductivity23_ahe(:) 
    real,allocatable    :: conductivity_fsur(:) 
    real,allocatable    :: conductivity13_fsur(:) 
    real,allocatable    :: conductivity23_fsur(:),sigma_tensor_ahc_mpi(:,:),sigma_tensor_ahc_mpi2(:,:)
    integer             :: ierr,isize,irank,kp1,kp2,kp3,num_steps_tot
    integer             :: ix1,ix2,ix3,num_occ 
    complex,allocatable :: momentum(:,:,:) 
    complex,allocatable :: momentum2(:,:,:)                                                                   
    complex,allocatable :: spinmomentum(:,:,:) 
    complex,allocatable :: spinmomentum2(:,:,:) 
 
    complex             :: berry 
 
    integer,allocatable :: sortarray(:) 
    integer             :: n1,n2,n3,n4,dir 
    real,allocatable    :: occupation(:),occupation2(:) 
    integer,allocatable :: nrpts(:) 
    real                :: occupation_number 
    integer             :: step,i,ii,num_steps,jj
    real                :: fermi_min,fermi_max 
    logical             :: l_tb,l_nabla,l_bfield 
    real                :: bfield(3) 
    integer             :: rvecnum,num_lines,length 
    integer,allocatable :: irvec(:,:) 
    real,allocatable    :: kpts(:,:) 
    real                :: kder,amat(3,3) 
    real                :: volume,cross(3) 
    real                :: bohrincm,condq 
    integer             :: maxdim,num_kpts 
    real                :: minenerg,maxenerg,gamma,kbT
    logical             :: l_bandstruc,l_fermisurf 
    real,allocatable    :: magnetic(:,:) 
    complex(kind(1.0d0)), allocatable :: Omega_x(:), Omega_y(:), Omega_z(:)
    complex(kind(1.0d0)), allocatable :: Omega_x_t(:), Omega_y_t(:), Omega_z_t(:)
    real(kind(1.0d0)) :: fermi

    real(kind(1.0d0)) :: K3D_start_cube(3)
    real(kind(1.0d0)) :: K3D_vec1_cube(3)
    real(kind(1.0d0)) :: K3D_vec2_cube(3)
    real(kind(1.0d0)) :: K3D_vec3_cube(3)

    INCLUDE 'mpif.h' 
    integer stt(MPI_STATUS_SIZE) 

    CALL MPI_INIT(ierr)                 
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,isize,ierr)

    abstol=2.0*tiny(abstol) 
    twopi=2*3.141592654          
    bohrincm=0.5291772*1.e-8 
    condq=38.7405*1.e-6 
    kbT = 0.0861733*1e-3*300    
    if(irank.eq.0)then 
        open(300,file='ahe_inp') 
        read(300,*) amat(1,:) 
        read(300,*) amat(2,:) 
        read(300,*) amat(3,:) 
        read(300,*) fermi_min,fermi_max,num_steps
        read(300,*) Nk1,Nk2,Nk3
        read(100,*) occupation_number
        read(100,*) l_mag_vec     
        read(100,*) mag_strength  
        read(100,*) mag_theta
        read(100,*) mag_phi
        read(100,*) mag_wann_num
        allocate(mag_wann_orbs_index(mag_wann_num))
        read(100,*) mag_wann_orbs_index(mag_wann_num)
        close(300) 
    endif 

    grid = Nk1

    call mpi_bcast(amat,9,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)                             
    call mpi_bcast(fermi_min,1,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)                             
    call mpi_bcast(occupation_number,1,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)                             
    call mpi_bcast(fermi_max,1,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)                             
    call mpi_bcast(num_steps,1,MPI_INTEGER,0,mpi_comm_world,ierr)                             
    call mpi_bcast(maxdim,1,MPI_INTEGER,0,mpi_comm_world,ierr)                             
    call mpi_bcast(grid,1,MPI_INTEGER,0,mpi_comm_world,ierr)    
    call mpi_bcast(Nk1,1,MPI_INTEGER,0,mpi_comm_world,ierr)   
    call mpi_bcast(Nk2,1,MPI_INTEGER,0,mpi_comm_world,ierr)   
    call mpi_bcast(Nk3,1,MPI_INTEGER,0,mpi_comm_world,ierr)                            
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
        
        ! print *, mag_field_x, mag_field_y, mag_field_z
        ! print *, mag_wann_orbs_index
    endif
     
    call mpi_bcast(mag_field_x,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
    call mpi_bcast(mag_field_y,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
    call mpi_bcast(mag_field_z,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)      
                                                                            
    if(irank.eq.0)then 
        open(200,file='hopping.1') 
        num_lines=0 
        num_wann=0 
        do 
            read(200,fmt=*,end=311)ix,iy,iz,band1,band2,rdum,idum 
            num_lines=num_lines+1 
            num_wann=max(num_wann,band1) 
        enddo 
311     continue 
        rvecnum=num_lines/(num_wann*num_wann)
        write(*,*)"num_lines=",num_lines
        write(*,*)"rvecnum=",rvecnum
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
300     continue 
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
500     continue 
        close(400) 
        write(*,*) 'rvecnum',rvecnum 
        allocate(nrpts(rvecnum)) 
        open(14,file='nrpts_inp')

        do j=1,rvecnum/15 
            read(14,'(15I5)') (nrpts(15*(j-1)+i) ,i=1,15) 
        enddo

        read(14,'(<mod(rvecnum,15)>I5)') (nrpts(15*(rvecnum/15)+i),i=1,mod(rvecnum,15))           
        close(14) 
        write(*,*) nrpts 
        write(*,*) 'size',size(nrpts)                                                                       
    endif 
                                                                            
    if(irank.eq.0)then 
        write(*,*)"rvecnum=",rvecnum 
        write(*,*)"num_wann=",num_wann 
    endif 
                                                                            
    call mpi_bcast(num_wann,1,MPI_INTEGER,0,mpi_comm_world,ierr)                             
    call mpi_bcast(rvecnum,1,MPI_INTEGER,0,mpi_comm_world,ierr)                             
                                                                                                                                           
    if(.not.allocated(nrpts))allocate(nrpts(rvecnum)) 
    call mpi_bcast(nrpts,rvecnum,MPI_INTEGER,0,mpi_comm_world,ierr)                             
                                                                            
    if(.not.allocated(irvec))allocate(irvec(3,rvecnum)) 
    length=3*rvecnum 
    call mpi_bcast(irvec,length,MPI_INTEGER,0,mpi_comm_world,ierr)                             
                                                                            
    if(.not.allocated(hops))allocate(hops(num_wann,num_wann,rvecnum)) 
    length=num_wann*num_wann*rvecnum 
    call mpi_bcast(hops,length,MPI_DOUBLE_COMPLEX,0,mpi_comm_world,ierr)                             
                                                                            
                                                                                                                                             
        cross(1)=amat(1,2)*amat(2,3)-amat(1,3)*amat(2,2) 
        cross(2)=amat(1,3)*amat(2,1)-amat(1,1)*amat(2,3) 
        cross(3)=amat(1,1)*amat(2,2)-amat(1,2)*amat(2,1) 
                                                                            
        volume=cross(1)*amat(3,1)+cross(2)*amat(3,2)+cross(3)*amat(3,3) 
        ! if(maxdim.eq.2)then 
        !    volume=volume/amat(3,3)/(0.529e-8) !!!   unit is cm
        ! endif 

        K3D_start_cube= (/ 0.0,  0.0,   0.0/)
        K3D_vec1_cube = (/ 1.0,  0.0,   0.0/)
        K3D_vec2_cube = (/ 0.0,  1.0,   0.0/)
        K3D_vec3_cube = (/ 0.0,  0.0,   1.0/)    
        
        allocate( sortarray(num_steps) ) 
        allocate( deviation(num_steps) ) 
        allocate( occupation(num_steps) ) 
        allocate( occupation2(num_steps) ) 
        allocate( conductivity(num_steps) ) 
        allocate( conductivity2(num_steps) ) 
        allocate( conductivity13(num_steps) ) 
        allocate( conductivity23(num_steps) ) 
        allocate( conductivity_ahe(num_steps) ) 
        allocate( conductivity13_ahe(num_steps) ) 
        allocate( conductivity23_ahe(num_steps) ) 
        allocate( fermienergy(num_steps) ) 
        
        ! allocate(sigma_tensor_ahc_x(num_steps))
        ! allocate(sigma_tensor_ahc_y(num_steps))
        ! allocate(sigma_tensor_ahc_z(num_steps))
        allocate(sigma_tensor_ahc_mpi(3,num_steps))
        allocate(sigma_tensor_ahc_mpi2(3,num_steps))
        allocate(Omega_x(num_wann), Omega_y(num_wann), Omega_z(num_wann))
        allocate(Omega_x_t(num_wann), Omega_y_t(num_wann), Omega_z_t(num_wann))
    
    
        magnetic=0.0 
        occupation=0.0 
        conductivity=0.0 
        conductivity13=0.0 
        conductivity23=0.0 
        conductivity_ahe=0.0 
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
        allocate( work(lwork) ) 
        allocate( rwork(17*num_wann) ) 
        allocate( iwork(15*num_wann) ) 
        allocate( ifail(15*num_wann) ) 

        if(l_mag_vec)then        
            if(.not.allocated(rspauli)) allocate(rspauli(num_wann,num_wann,3,rvecnum))              
                length=num_wann*num_wann*rvecnum*3 
                call mpi_bcast(rspauli,length,MPI_DOUBLE_COMPLEX,0,mpi_comm_world,ierr) 
            do ii=1,rvecnum
                do i=1,num_wann
                    do j=1,num_wann
                        do jj =1,mag_wann_num
                            if((i==mag_wann_orbs_index(jj)) .or. (j==mag_wann_orbs_index(jj)))then
                                hops(j,i,ii) = hops(j,i,ii) + mag_field_x * rspauli(j,i,1,ii)     &
                                + mag_field_y * rspauli(j,i,2,ii) + mag_field_z * rspauli(j,i,3,ii)
                            endif
                        enddo
                    enddo !j
                enddo !i
            enddo
        endif        

        knv3= Nk1*Nk2*Nk3
        do ik=1+ irank,knv3,isize
            ! if(mod(ik,isize).ne.irank)cycle
            ikx= (ik-1)/(nk2*nk3)+1
            iky= ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
            ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
            kpoints= K3D_start_cube + K3D_vec1_cube*(ikx-1)/dble(nk1)  &
                                    + K3D_vec2_cube*(iky-1)/dble(nk2)  &
                                    + K3D_vec3_cube*(ikz-1)/dble(nk3)
          ! ik1=0
          ! do kp1=0,grid-1 
          !  kpoints(1)=-0.5+real(kp1)/real(grid) 
          !  do kp2=0,grid-1 
          !   kpoints(2)=-0.5+real(kp2)/real(grid) 
          !   do kp3=0,(grid-1)*(maxdim-2) 
          !    ik1=ik1+1 
          !    if(mod(ik1-1,isize).ne.irank)cycle                     
          !    kpoints(3)=-0.5+real(kp3)/real(grid)        
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
                        ham(j,i)=ham(j,i)+fac*hops(j,i,ii)/nrpts(ii) 
                        do dir=1,3 
                            pauli(j,i,dir)=pauli(j,i,dir)+fac*rspauli(j,i,dir,ii) 
                        enddo                                               
                    enddo                                                 
                enddo
            enddo
        call zheevx('V','A','U',num_wann,ham,num_wann,                      &
         &        vl,vu,1,num_wann,abstol,ne,eigvals,eigvecs,num_wann,      &
         &        work,lwork,rwork,iwork,ifail,info) 
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
                        momentum(n1,n2,dir)=momentum(n1,n2,dir)+fac2*hops(n1,n2,ii)                                 
                    enddo
                enddo
            enddo
        enddo

        momentum2=0.0 
        do dir=1,3
            do n2=1,num_wann
                do n4=1,num_wann
                    do n1=1,num_wann
                        do n3=1,num_wann
                            momentum2(n4,n2,dir)=momentum2(n4,n2,dir)+momentum(n3,n1,dir)*eigvecs(n1,n2)*conjg(eigvecs(n3,n4))                                             
                        enddo
                    enddo
                enddo
            enddo
        enddo
        !  print*,"momentum2",momentum2                       
        Beta_fake = 1.0/kbT
    
        Omega_x=cmplx(0.d0,0.d0);Omega_y=cmplx(0.d0,0.d0); Omega_z=cmplx(0.d0,0.d0)
        do m=1,num_wann
            do n=1,num_wann
                if (abs(eigvals(m)-eigvals(n))<1e-6) cycle 
                    Omega_z(m)= Omega_z(m)-aimag(momentum2(m,n,1)*conjg(momentum2(m,n,2)))/(eigvals(m)-eigvals(n))**2
                    Omega_y(m)= Omega_y(m)-aimag(momentum2(m,n,1)*conjg(momentum2(m,n,3)))/(eigvals(m)-eigvals(n))**2 
                    Omega_x(m)= Omega_x(m)-aimag(momentum2(m,n,2)*conjg(momentum2(m,n,3)))/(eigvals(m)-eigvals(n))**2   
            enddo !n
        enddo  !m
                ! print*,"omega_x ok"
    
        Omega_x_t=cmplx(0.d0,0.d0);Omega_y_t=cmplx(0.d0,0.d0); Omega_z_t=cmplx(0.d0,0.d0)
        do step=1,num_steps                                                             
            fermienergy(step) = fermi_min + (fermi_max-fermi_min)*real(step)/real(num_steps)        
            mu = fermienergy(step)
                ! print*,"mu=",mu
                ! print*,"step=",step
            iy=0 
            do ix=1,num_wann 
                if(eigvals(ix).le.fermienergy(step))then 
                      iy=iy+1 
                endif 
            enddo 
                                                                               
            num_occ=iy                                                          
            occupation(step)=occupation(step)+num_occ 
            do m= 1, num_wann
                Omega_x_t(m)= Omega_x(m)*(-1.0/(1*300))*((eigvals(m)-mu)*fermi(eigvals(m)-mu, Beta_fake)+(kbT*log(1+exp(-(eigvals(m)-mu)/kbT))))
                Omega_y_t(m)= Omega_y(m)*(-1.0/(1*300))*((eigvals(m)-mu)*fermi(eigvals(m)-mu, Beta_fake)+(kbT*log(1+exp(-(eigvals(m)-mu)/kbT))))
                Omega_z_t(m)= Omega_z(m)*(-1.0/(1*300))*((eigvals(m)-mu)*fermi(eigvals(m)-mu, Beta_fake)+(kbT*log(1+exp(-(eigvals(m)-mu)/kbT))))
                ! Omega_x_t(m)= Omega_x(m)*fermi(eigvals(m)-mu, Beta_fake)
                ! Omega_y_t(m)= Omega_y(m)*fermi(eigvals(m)-mu, Beta_fake)
                ! Omega_z_t(m)= Omega_z(m)*fermi(eigvals(m)-mu, Beta_fake)
            enddo
            ! print*,"Omega_x_t OK"
            sigma_tensor_ahc_mpi(1, step)= sigma_tensor_ahc_mpi(1, step)+ real(sum(Omega_x_t))
            sigma_tensor_ahc_mpi(2, step)= sigma_tensor_ahc_mpi(2, step)+ real(sum(Omega_y_t))
            sigma_tensor_ahc_mpi(3, step)= sigma_tensor_ahc_mpi(3, step)+ real(sum(Omega_z_t))
            ! print*,"sigma_tensor_ahc_mpi(1, step)=",sigma_tensor_ahc_mpi(1, step)
        enddo
          ! if(irank.eq.0)then 
          ! print*,"kpoint=",kpoints
          ! print*,"irank=",irank
          ! print*,"isize=",isize
          ! endif
       ! enddo
       ! enddo
    enddo!end k
    
    occupation=occupation/knv3
    call MPI_REDUCE(occupation,occupation2,num_steps,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)                                      
    occupation=occupation2    
    ! num_steps_tot =num_steps*3                                                        
    call MPI_REDUCE(sigma_tensor_ahc_mpi,sigma_tensor_ahc_mpi2,size(sigma_tensor_ahc_mpi2),MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
    sigma_tensor_ahc_mpi=sigma_tensor_ahc_mpi2
    ! if (maxdim == 2)then
       ! sigma_tensor_ahc_mpi=sigma_tensor_ahc_mpi/knv3/(volume/amat(3,3)/(0.529117e-8))/bohrincm*condq*2.0*twopi
    ! endif

    ! if (maxdim == 3)then
       sigma_tensor_ahc_mpi=sigma_tensor_ahc_mpi/knv3/volume*condq*2.0*twopi
    ! endif
    conductivity_ahe=sigma_tensor_ahc_mpi(3,:) 
    conductivity13_ahe=sigma_tensor_ahc_mpi(2,:)   
    conductivity23_ahe=sigma_tensor_ahc_mpi(1,:)

    if(irank.eq.0)then 
        open(123,file='output_ahe_condquant',recl=10000) 
        do step=1,num_steps 
            write(123,*)"fermienergy=",fermienergy(step)                                                          
            write(123,*)"occupation=",occupation(step)                                
            ! write(123,*)"conductivity=",conductivity_ahe(step)     /(0.5*77.48e-6)*amat(3,3)                   
            ! write(123,*)"conductivity13=",conductivity13_ahe(step) /(0.5*77.48e-6)*amat(3,3)                 
            ! write(123,*)"conductivity23=",conductivity23_ahe(step) /(0.5*77.48e-6)*amat(3,3)   
            write(123,*)"conductivity=",conductivity_ahe(step)     /bohrincm *100               
            write(123,*)"conductivity13=",conductivity13_ahe(step) /bohrincm *100          
            write(123,*)"conductivity23=",conductivity23_ahe(step) /bohrincm *100                       
            write(123,*)"************************************"                                                  
        enddo 
        close(123)                                                                                                                        
    endif

    call mpi_barrier(mpi_comm_world,ierr) 
    call MPI_Finalize(ierr) 
                                                                                                                                 
    end program wannier_ham_magnetization
         
function fermi(omega, Beta_fake) result(value)
    implicit none
    real(kind(1.0d0)), intent(in) :: omega
    real(kind(1.0d0)), intent(in) :: Beta_fake
    real(kind(1.0d0)) :: value 
    if (Beta_fake*omega .ge. 20d0) then
       value = 0.0
    elseif (Beta_fake*omega .le. -20d0)then
       value = 1.0
    else
       value= 1.0/(1.0+exp(Beta_fake*omega))
    endif
    return
end function fermi    