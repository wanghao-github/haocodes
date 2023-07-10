program anomalous_nernst_effect
!*************************************************                      
!rewrite by HaoWang March 7 2023 replace momentum2 by wt
!*************************************************                      
    implicit none                                                      
    complex,allocatable:: hops   (:,:,:),spin_dir(:,:),spin_dir_mpi(:,:),spin_dir_mpi2(:,:)
    complex,allocatable:: rsnabla(:,:,:,:),rspauli_final(:,:,:,:),rspauli_ori(:,:,:,:)
    complex,allocatable:: rspauli(:,:,:,:),rspauli1(:,:,:,:),rspauli2(:,:,:,:),rspauli3(:,:,:,:),rspauli4(:,:,:,:)
    complex,allocatable:: pauli(:,:,:),rspauli5(:,:,:,:),rspauli6(:,:,:,:),rspauli7(:,:,:,:),rspauli8(:,:,:,:)
    complex,allocatable:: paulifft(:,:,:) 
    complex,allocatable:: paulifft2(:,:,:) 
    
    complex,allocatable:: spin_texture(:,:,:), spin_texture_mpi(:,:,:)
    real               :: rdum,idum 
    integer            :: ix,iy,iz,band1,band2,h,num_wann,num_wann_2,m,n
    integer            :: ik1,ik2,ik3 
    real               :: twopi,temperature
    real               :: phas 
    complex            :: fac,fac2,zi
    complex,allocatable:: ham(:,:),spin_sigma_x(:,:),spin_sigma_y(:,:),spin_sigma_z(:,:),spin_sigma_temp(:,:)
    complex,allocatable:: spin_sigma_x_comp(:,:),spin_sigma_y_comp(:,:),spin_sigma_z_comp(:,:),spin_sigma_t_mpi(:,:),spin_sigma_t_mpi2(:,:)
    real               :: vl,vu 
    integer            :: ne,j 
    real               :: abstol,time_start,time_end ,time_start1,time_end1
    real,allocatable   :: eigvals(:),eigvals_x(:),eigvals_y(:),eigvals_z(:)
    complex,allocatable:: eigvecs(:,:),eigvecs_dag(:,:),mat_temp(:,:),eigvecs_x(:,:)
    integer            :: info 
    complex,allocatable:: work(:),eigvecs_y(:,:),eigvecs_z(:,:)
    integer            :: lwork 
    integer,allocatable:: iwork(:) 
    real,allocatable   :: rwork(:) 
    integer,allocatable:: ifail(:),select_atom1(:,:),select_atom2(:,:),select_atom3(:,:),select_atom4(:,:)
    integer,allocatable:: select_atom5(:,:),select_atom6(:,:),select_atom7(:,:),select_atom8(:,:)
    real               :: kpoints(3)
    real               :: scale ,Beta_fake,mu
    integer            :: maxhopx2,maxhopy2,maxhopz2,dire,ik
    real,allocatable   :: fermienergy(:),spindirx(:),spindiry(:),spindirz(:)
    real,allocatable   :: deviation(:) 
    integer            :: grid,i1,i2,i3,i4,orb,Nk1,Nk2,Nk3,knv3,ikx,iky,ikz
    real,allocatable   :: conductivity(:),conductivity2(:),sigma_tensor_ahc_x(:)
    real,allocatable   :: conductivity13(:),conductivity23(:) 
    real,allocatable   :: conductivity_ahe(:) ,sigma_tensor_ahc_y(:),sigma_tensor_ahc_z(:)
    real,allocatable   :: conductivity13_ahe(:),conductivity23_ahe(:) 
    real,allocatable   :: conductivity_fsur(:) 
    real,allocatable   :: conductivity13_fsur(:) 
    real,allocatable   :: conductivity23_fsur(:),sigma_tensor_ahc_mpi(:,:),sigma_tensor_ahc_mpi2(:,:)
    integer            :: ierr,isize,irank,kp1,kp2,kp3,num_steps_tot
    integer            :: ix1,ix2,ix3,num_occ 
    
    complex,allocatable:: momentum(:,:,:) 
    complex,allocatable:: momentum2(:,:,:)                                                                 
    complex,allocatable:: spinmomentum(:,:,:) 
    complex,allocatable:: spinmomentum2(:,:,:) 
    
    integer,parameter  :: Dp=kind(1.0d0)      
    complex            :: berry 

    integer,allocatable:: sortarray(:),mag_wann_orbs_index1(:),mag_wann_orbs_index2(:),mag_wann_orbs_index3(:),mag_wann_orbs_index8(:)
    integer            :: n1,n2,n3,n4,dir,mag_wann_num1,mag_wann_num2,mag_wann_num5,mag_wann_num6,mag_wann_num7,mag_wann_num8
    real,allocatable   :: occupation(:),occupation2(:)
    integer,allocatable:: nrpts(:),mag_wann_orbs_index4(:),mag_wann_orbs_index5(:),mag_wann_orbs_index6(:),mag_wann_orbs_index7(:)
    real               :: occupation_number
    integer            :: step,i,ii,num_steps,mag_wann_num3,mag_wann_num4
    real               :: fermi_min,fermi_max
    logical            :: l_tb,l_nabla,l_bfield,l_mag_vec
    real               :: bfield(3)
    integer            :: rvecnum,num_lines,length 
    integer,allocatable:: irvec(:,:)
    real,allocatable   :: kpts(:,:) ,crvec(:,:)
    real               :: kder,amat(3,3) 
    real               :: volume,cross(3),mag_field_x1,mag_field_y1,mag_field_z1,mag_field_x2,mag_field_y2,mag_field_z2
    real               :: mag_field_x3,mag_field_y3,mag_field_z3,mag_field_x4,mag_field_y4,mag_field_z4
    real               :: mag_field_x5,mag_field_y5,mag_field_z5,mag_field_x6,mag_field_y6,mag_field_z6
    real               :: mag_field_x7,mag_field_y7,mag_field_z7,mag_field_x8,mag_field_y8,mag_field_z8
    real               :: bohrincm,condq,mag_strength,mag_theta1,mag_phi1,mag_theta2,mag_phi2
    real               :: mag_theta3,mag_phi3,mag_theta4,mag_phi4,mag_theta5,mag_phi5,mag_theta6,mag_phi6
    real               :: mag_theta7,mag_phi7,mag_theta8,mag_phi8
    integer            :: maxdim,num_kpts,num_spin
    real               :: minenerg,maxenerg,gamma,kbT
    logical            :: l_bandstruc,l_fermisurf 
    real,allocatable   :: magnetic(:,:) 
    complex(kind(1.0d0)), allocatable :: Omega_x(:), Omega_y(:), Omega_z(:)
    complex(kind(1.0d0)), allocatable :: Omega_x_t(:), Omega_y_t(:), Omega_z_t(:),spin_sigma_x_t(:),spin_sigma_y_t(:),spin_sigma_z_t(:)
    real(kind(1.0d0)) :: fermi
    
    real(kind(1.0d0)) :: K3D_start_cube(3)
    real(kind(1.0d0)) :: K3D_vec1_cube(3)
    real(kind(1.0d0)) :: K3D_vec2_cube(3)
    real(kind(1.0d0)) :: K3D_vec3_cube(3)

!!!以上全部都是变量类型声明
    INCLUDE 'mpif.h' 
    integer stt(MPI_STATUS_SIZE) 
                                                                      
    CALL MPI_INIT(ierr)                                                            
    CALL MPI_COMM_RANK (MPI_COMM_WORLD,irank,ierr) 
    CALL MPI_COMM_SIZE (MPI_COMM_WORLD,isize,ierr) 
                                                                                                                                
    abstol=2.0*tiny(abstol) 
    twopi=2*3.141592654 
    zi = cmplx(0.,1.)                                                            
    bohrincm=0.5291772*1.e-8 
    condq=38.7405*1.e-6 
!!!在所有CPU上定义的常数

    if(irank.eq.0)then 
        open(300,file='ahe_inp') 
        read(300,*)amat(1,:) 
        read(300,*)amat(2,:) 
        read(300,*)amat(3,:) 
        read(300,*)fermi_min,fermi_max,num_steps 
        read(300,*)Nk1,Nk2,Nk3 
        read(300,*)maxdim 
        read(300,*)occupation_number 
        read(300,*)temperature
        read(300,*)l_mag_vec         
        read(300,*)mag_strength                
        read(300,*)mag_theta1                  
        read(300,*)mag_phi1
        read(300,*)mag_theta2
        read(300,*)mag_phi2
        read(300,*)mag_theta3
        read(300,*)mag_phi3
        read(300,*)mag_theta4                 
        read(300,*)mag_phi4
        read(300,*)mag_theta5
        read(300,*)mag_phi5
        read(300,*)mag_theta6
        read(300,*)mag_phi6
        read(300,*)mag_theta7
        read(300,*)mag_phi7
        read(300,*)mag_theta8
        read(300,*)mag_phi8       
        read(300,*)mag_wann_num1
        read(300,*)mag_wann_num2
        read(300,*)mag_wann_num3
        read(300,*)mag_wann_num4
        read(300,*)mag_wann_num5
        read(300,*)mag_wann_num6
        read(300,*)mag_wann_num7
        read(300,*)mag_wann_num8
        allocate(mag_wann_orbs_index1(mag_wann_num1))
        allocate(mag_wann_orbs_index2(mag_wann_num2))
        allocate(mag_wann_orbs_index3(mag_wann_num3))
        allocate(mag_wann_orbs_index4(mag_wann_num4))
        allocate(mag_wann_orbs_index5(mag_wann_num5))
        allocate(mag_wann_orbs_index6(mag_wann_num6))
        allocate(mag_wann_orbs_index7(mag_wann_num7))
        allocate(mag_wann_orbs_index8(mag_wann_num8))
        read(300,*)mag_wann_orbs_index1
        read(300,*)mag_wann_orbs_index2
        read(300,*)mag_wann_orbs_index3
        read(300,*)mag_wann_orbs_index4
        read(300,*)mag_wann_orbs_index5
        read(300,*)mag_wann_orbs_index6
        read(300,*)mag_wann_orbs_index7
        read(300,*)mag_wann_orbs_index8
        close(300) 
        write(*,*) mag_wann_orbs_index1
        write(*,*) mag_wann_orbs_index2
        write(*,*) mag_wann_orbs_index3
        write(*,*) mag_wann_orbs_index4
        write(*,*) mag_wann_orbs_index5
        write(*,*) mag_wann_orbs_index6
        write(*,*) mag_wann_orbs_index7
        write(*,*) mag_wann_orbs_index8
    endif 

!!! 读入所有输入的轨道

    grid = Nk1

    call mpi_bcast(amat,9,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)                             
    call mpi_bcast(fermi_min,1,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)                             
    call mpi_bcast(occupation_number,1,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)                             
    call mpi_bcast(fermi_max,1,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)                             
    call mpi_bcast(num_steps,1,MPI_INTEGER, 0,mpi_comm_world,ierr)                             
    call mpi_bcast(maxdim,1,MPI_INTEGER, 0,mpi_comm_world,ierr)                             
    call mpi_bcast(grid,1,MPI_INTEGER,0,mpi_comm_world,ierr)    
    call mpi_bcast(Nk1,1,MPI_INTEGER, 0,mpi_comm_world,ierr)   
    call mpi_bcast(Nk2,1,MPI_INTEGER, 0,mpi_comm_world,ierr)   
    call mpi_bcast(Nk3,1,MPI_INTEGER, 0,mpi_comm_world,ierr)                            
    !call mpi_bcast(l_tb,1,MPI_LOGICAL,0,mpi_comm_world,ierr)                             
    !call mpi_bcast(l_nabla,1,MPI_LOGICAL,0,mpi_comm_world,ierr)                             
    !call mpi_bcast(l_bfield,1,MPI_LOGICAL, 0,mpi_comm_world,ierr)                             
    !call mpi_bcast(l_fermisurf,1,MPI_LOGICAL,0,mpi_comm_world,ierr)                             
    !call mpi_bcast(bfield,3,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)                             
    !call mpi_bcast(gamma,1,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)                             
    call mpi_bcast(bfield,3,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)        
    call mpi_bcast(l_mag_vec,1,MPI_LOGICAL,0,mpi_comm_world,ierr)
    call mpi_bcast(mag_strength,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
    call mpi_bcast(mag_theta1,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
    call mpi_bcast(mag_phi1,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
    call mpi_bcast(mag_theta2,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
    call mpi_bcast(mag_phi2,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
    call mpi_bcast(mag_theta3,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
    call mpi_bcast(mag_phi3,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
    call mpi_bcast(mag_theta4,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
    call mpi_bcast(mag_phi4,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
    call mpi_bcast(mag_theta5,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
    call mpi_bcast(mag_phi5,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
    call mpi_bcast(mag_theta6,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
    call mpi_bcast(mag_phi6,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
    call mpi_bcast(mag_theta7,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
    call mpi_bcast(mag_phi7,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
    call mpi_bcast(mag_theta8,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
    call mpi_bcast(mag_phi8,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)

    ! call mpi_bcast(mag_wann_num1,1,MPI_INTEGER, 0,mpi_comm_world,ierr)
    ! call mpi_bcast(mag_wann_orbs_index1,mag_wann_num1,MPI_INTEGER,0,mpi_comm_world,ierr)
    ! call mpi_bcast(mag_wann_num2,1,MPI_INTEGER, 0,mpi_comm_world,ierr)
    ! call mpi_bcast(mag_wann_orbs_index2,mag_wann_num2,MPI_INTEGER,0,mpi_comm_world,ierr)
    call mpi_bcast(temperature,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
    call mpi_bcast(kbT,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)

    kbT = 0.0861733*1e-3*temperature
    amat= amat*0.5291772

    mag_field_x1 = mag_strength * cos(mag_theta1)
    mag_field_y1 = mag_strength * sin(mag_theta1) * cos(mag_phi1)
    mag_field_z1 = mag_strength * sin(mag_theta1) * sin(mag_phi1)

    mag_field_x2 = mag_strength * cos(mag_theta2)
    mag_field_y2 = mag_strength * sin(mag_theta2) * cos(mag_phi2)
    mag_field_z2 = mag_strength * sin(mag_theta2) * sin(mag_phi2)
        
    mag_field_x3 = mag_strength * cos(mag_theta3)
    mag_field_y3 = mag_strength * sin(mag_theta3) * cos(mag_phi3)
    mag_field_z3 = mag_strength * sin(mag_theta3) * sin(mag_phi3)

    mag_field_x4 = mag_strength * cos(mag_theta4)
    mag_field_y4 = mag_strength * sin(mag_theta4) * cos(mag_phi4)
    mag_field_z4 = mag_strength * sin(mag_theta4) * sin(mag_phi4)

    mag_field_x5 = mag_strength * cos(mag_theta5)
    mag_field_y5 = mag_strength * sin(mag_theta5) * cos(mag_phi5)
    mag_field_z5 = mag_strength * sin(mag_theta5) * sin(mag_phi5)

    mag_field_x6 = mag_strength * cos(mag_theta6)
    mag_field_y6 = mag_strength * sin(mag_theta6) * cos(mag_phi6)
    mag_field_z6 = mag_strength * sin(mag_theta6) * sin(mag_phi6)

    mag_field_x7 = mag_strength * cos(mag_theta7)
    mag_field_y7 = mag_strength * sin(mag_theta7) * cos(mag_phi7)
    mag_field_z7 = mag_strength * sin(mag_theta7) * sin(mag_phi7)

    mag_field_x8 = mag_strength * cos(mag_theta8)
    mag_field_y8 = mag_strength * sin(mag_theta8) * cos(mag_phi8)
    mag_field_z8 = mag_strength * sin(mag_theta8) * sin(mag_phi8)

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
        allocate( hops(1:num_wann,1:num_wann,rvecnum) ) 
        allocate( irvec(3,rvecnum) ) 
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

        allocate(crvec(3,rvecnum))
        
        do ii=1,rvecnum
            crvec(:,ii)= amat(1,:)*irvec(1,ii)+amat(2,:)* irvec(2,ii)+amat(3,:)* irvec(3,ii)
        enddo
        
        allocate(rspauli(1:num_wann, 1:num_wann, 3, rvecnum)) 
        allocate(rspauli1(1:num_wann, 1:num_wann, 3, rvecnum))
        allocate(rspauli2(1:num_wann, 1:num_wann, 3, rvecnum)) 
        allocate(rspauli3(1:num_wann, 1:num_wann, 3, rvecnum))
        allocate(rspauli4(1:num_wann, 1:num_wann, 3, rvecnum))
        allocate(rspauli5(1:num_wann, 1:num_wann, 3, rvecnum))
        allocate(rspauli6(1:num_wann, 1:num_wann, 3, rvecnum))
        allocate(rspauli7(1:num_wann, 1:num_wann, 3, rvecnum))
        allocate(rspauli8(1:num_wann, 1:num_wann, 3, rvecnum))
        allocate(rspauli_ori(1:num_wann, 1:num_wann, 3, rvecnum))
        allocate(rspauli_final(1:num_wann, 1:num_wann, 3, rvecnum))
        allocate(spin_sigma_x(num_wann,num_wann))
        allocate(spin_sigma_y(num_wann,num_wann))
        allocate(spin_sigma_z(num_wann,num_wann))
        
        spin_sigma_x=0.0d0;spin_sigma_y=0.0d0;spin_sigma_z=0.0d0

        num_wann_2 = num_wann/2
        write(*,*) "num_wann_2=", num_wann_2
        do j=1, num_wann_2
            spin_sigma_x(j, num_wann_2+j)=1.0d0
            spin_sigma_x(j+num_wann_2, j)=1.0d0
            spin_sigma_y(j, num_wann_2+j)=-zi
            spin_sigma_y(j+num_wann_2, j)=zi
            spin_sigma_z(j, j)= 1.0d0
            spin_sigma_z(j+num_wann_2, j+num_wann_2)=-1.0d0
        enddo

        write(*,*) "spin_sigma_x="
        write(*,*) spin_sigma_x
        write(*,*) "spin_sigma_y="
        write(*,*) spin_sigma_y
        write(*,*) "spin_sigma_z="
        write(*,*) spin_sigma_z
      
        open(400,file='./rspauli.1') 
        num_lines=0

        do 
            read(400, fmt=*,end=500) ix,iy,iz,band1,band2,dir,rdum,idum 
            num_lines=num_lines+1 
            rvecnum=(num_lines-1)/(num_wann*num_wann*3)+1 
            rspauli(band1, band2, dir, rvecnum)=cmplx(rdum,idum) 
        enddo

500     continue 
        
        close(400)  
        
        allocate(nrpts(rvecnum)) 
        
        open(14,file='nrpts_inp') 
        do j=1,rvecnum/15 
                read(14,'(15I5)') (nrpts(15*(j-1)+i) ,i=1,15) 
        enddo 
        read(14,'(<mod(rvecnum,15)>I5)') (nrpts(15*(rvecnum/15)+i),i=1,mod(rvecnum,15))           
        close(14)                                   
        
!         open(444,file='./rspauli.2')
!         num_lines=0
!         Do
!             read(444, fmt=*,end=555) ix,iy,iz,band1,band2,dir,rdum,idum
!             num_lines=num_lines+1
!             rvecnum=(num_lines-1)/(num_wann*num_wann*3)+1
!             rspauli_ori(band1, band2, dir, rvecnum)=cmplx(rdum,idum)
!         End Do
! 555     continue
!         close(444)
      
    endif 
      
        call mpi_bcast(num_wann,1,MPI_INTEGER,0,mpi_comm_world,ierr)                             
        call mpi_bcast(rvecnum,1,MPI_INTEGER,0,mpi_comm_world,ierr) 

    if(irank == 0)then

        allocate(select_atom1(num_wann,num_wann))
        allocate(select_atom2(num_wann,num_wann))
        allocate(select_atom3(num_wann,num_wann))
        allocate(select_atom4(num_wann,num_wann))
        allocate(select_atom5(num_wann,num_wann))
        allocate(select_atom6(num_wann,num_wann))
        allocate(select_atom7(num_wann,num_wann))
        allocate(select_atom8(num_wann,num_wann))

        select_atom1=0
        select_atom2=0
        select_atom3=0
        select_atom4=0
        select_atom5=0
        select_atom6=0
        select_atom7=0
        select_atom8=0

        write(*,*) "mag_wann_orbs_index1", mag_wann_orbs_index1
        write(*,*) "mag_wann_orbs_index2", mag_wann_orbs_index2
        write(*,*) "mag_wann_orbs_index3", mag_wann_orbs_index3
        write(*,*) "mag_wann_orbs_index4", mag_wann_orbs_index4
        write(*,*) "mag_wann_orbs_index5", mag_wann_orbs_index5
        write(*,*) "mag_wann_orbs_index6", mag_wann_orbs_index6
        write(*,*) "mag_wann_orbs_index7", mag_wann_orbs_index7
        write(*,*) "mag_wann_orbs_index8", mag_wann_orbs_index8

        do m = 1,num_wann/2
        
        if (ANY(mag_wann_orbs_index1 == m))then
            write(*,*) "satisified 1", m
            select_atom1(m,m) = 1
            select_atom1(m+num_wann/2,m+num_wann/2) = 1 
        endif

        if (ANY(mag_wann_orbs_index2 == m))then
            write(*,*) "satisified 2", m
            select_atom2(m,m) = 1
            select_atom2(m+num_wann/2,m+num_wann/2) = 1
        endif

        if (ANY(mag_wann_orbs_index3 == m))then
            write(*,*) "satisified 3", m
            select_atom3(m,m) = 1
            select_atom3(m+num_wann/2,m+num_wann/2) = 1
        endif

        if (ANY(mag_wann_orbs_index4 == m))then
            write(*,*) "satisified 4", m
            select_atom4(m,m) = 1
            select_atom4(m+num_wann/2,m+num_wann/2) = 1
        endif

        if (ANY(mag_wann_orbs_index5 == m))then
            write(*,*) "satisified 5", m
            select_atom5(m,m) = 1
            select_atom5(m+num_wann/2,m+num_wann/2) = 1
        endif

        if (ANY(mag_wann_orbs_index6 == m))then
            write(*,*) "satisified 6", m
            select_atom6(m,m) = 1
            select_atom6(m+num_wann/2,m+num_wann/2) = 1
        endif

        if (ANY(mag_wann_orbs_index7 == m))then
            write(*,*) "satisified 7", m
            select_atom7(m,m) = 1
            select_atom7(m+num_wann/2,m+num_wann/2) = 1
        endif

        if (ANY(mag_wann_orbs_index8 == m))then
            write(*,*) "satisified 8", m
            select_atom8(m,m) = 1
            select_atom8(m+num_wann/2,m+num_wann/2) = 1
        endif

        enddo

        ! select_atom2 = reshape(select_atom2,[num_wann,num_wann])
        ! do n = 1, num_wann
        !     write(*, '(108I2)') (select_atom2(n, m),m=1,n)
        ! end do

        write(*,*) mag_field_x1,mag_field_y1,mag_field_z1,mag_field_x2,mag_field_y2,mag_field_z2
        write(*,*) mag_field_x3,mag_field_y3,mag_field_z3,mag_field_x4,mag_field_y4,mag_field_z4
        write(*,*) mag_field_x5,mag_field_y5,mag_field_z5,mag_field_x6,mag_field_y6,mag_field_z6
        write(*,*) mag_field_x7,mag_field_y7,mag_field_z7,mag_field_x8,mag_field_y8,mag_field_z8

        rspauli_final=0.0d0

        if (l_mag_vec .eq. .true.)then
            do ii = 1,rvecnum
            ! for m in 1,num_wann
                ! for n in 1,num_wann
                do dir = 1,3

                    rspauli1(:,:,dir,ii) = MATMUL(MATMUL(select_atom1,rspauli(:,:,dir,ii)),select_atom1)
                    rspauli2(:,:,dir,ii) = MATMUL(MATMUL(select_atom2,rspauli(:,:,dir,ii)),select_atom2)
                    rspauli3(:,:,dir,ii) = MATMUL(MATMUL(select_atom3,rspauli(:,:,dir,ii)),select_atom3)
                    rspauli4(:,:,dir,ii) = MATMUL(MATMUL(select_atom4,rspauli(:,:,dir,ii)),select_atom4)
                    rspauli5(:,:,dir,ii) = MATMUL(MATMUL(select_atom5,rspauli(:,:,dir,ii)),select_atom5)
                    rspauli6(:,:,dir,ii) = MATMUL(MATMUL(select_atom6,rspauli(:,:,dir,ii)),select_atom6)
                    rspauli7(:,:,dir,ii) = MATMUL(MATMUL(select_atom7,rspauli(:,:,dir,ii)),select_atom7)
                    rspauli8(:,:,dir,ii) = MATMUL(MATMUL(select_atom8,rspauli(:,:,dir,ii)),select_atom8)

                enddo
                !!! 下面是传统的直接把rspauli加到hop上的方法
                hops(:,:,ii) = hops(:,:,ii) + mag_field_x1 * rspauli1(:,:,1,ii)/2.0d0 + &
                mag_field_y1 * rspauli1(:,:,2,ii)/2.0d0 + mag_field_z1 * rspauli1(:,:,3,ii)/2.0d0 + &
                mag_field_x2 * rspauli2(:,:,1,ii)/2.0d0 + mag_field_y2 * rspauli2(:,:,2,ii)/2.0d0 + & 
                mag_field_z2 * rspauli2(:,:,3,ii)/2.0d0 + mag_field_x3 * rspauli3(:,:,1,ii)/2.0d0 + &
                mag_field_y3 * rspauli3(:,:,2,ii)/2.0d0 + mag_field_z3 * rspauli3(:,:,3,ii)/2.0d0 + &
                mag_field_x4 * rspauli4(:,:,1,ii)/2.0d0 + mag_field_y4 * rspauli4(:,:,2,ii)/2.0d0 + &
                mag_field_z4 * rspauli4(:,:,3,ii)/2.0d0 + mag_field_x5 * rspauli5(:,:,1,ii)/2.0d0 + &
                mag_field_y5 * rspauli5(:,:,2,ii)/2.0d0 + mag_field_z5 * rspauli5(:,:,3,ii)/2.0d0 + &
                mag_field_x6 * rspauli6(:,:,1,ii)/2.0d0 + mag_field_y6 * rspauli6(:,:,2,ii)/2.0d0 + &
                mag_field_z6 * rspauli6(:,:,3,ii)/2.0d0 + mag_field_x7 * rspauli7(:,:,1,ii)/2.0d0 + &
                mag_field_y7 * rspauli7(:,:,2,ii)/2.0d0 + mag_field_z7 * rspauli7(:,:,3,ii)/2.0d0 + &
                mag_field_x8 * rspauli8(:,:,1,ii)/2.0d0 + mag_field_y8 * rspauli8(:,:,2,ii)/2.0d0 + &
                mag_field_z8 * rspauli8(:,:,3,ii)/2.0d0 
            
                !!! 这个是我又想出来的制造完美spiral然后再加 原始rspauli的方法2

                ! rspauli_final(:,:,1,ii) = rspauli_final(:,:,1,ii) + mag_field_x1 * rspauli1(:,:,1,ii)/2.0d0 + &
                ! mag_field_x2 * rspauli2(:,:,1,ii)/2.0d0 + mag_field_x3 * rspauli3(:,:,1,ii)/2.0d0 + &
                ! mag_field_x4 * rspauli4(:,:,1,ii)/2.0d0 + mag_field_x5 * rspauli5(:,:,1,ii)/2.0d0 + &
                ! mag_field_x6 * rspauli6(:,:,1,ii)/2.0d0 + mag_field_x7 * rspauli7(:,:,1,ii)/2.0d0 + &
                ! mag_field_x8 * rspauli8(:,:,1,ii)/2.0d0

                ! rspauli_final(:,:,2,ii) = rspauli_final(:,:,2,ii) + mag_field_y1 * rspauli1(:,:,2,ii)/2.0d0 + &
                ! mag_field_y2 * rspauli2(:,:,2,ii)/2.0d0 + mag_field_y3 * rspauli3(:,:,2,ii)/2.0d0 + &
                ! mag_field_y4 * rspauli4(:,:,2,ii)/2.0d0 + mag_field_y5 * rspauli5(:,:,2,ii)/2.0d0 + &
                ! mag_field_y6 * rspauli6(:,:,2,ii)/2.0d0 + mag_field_y7 * rspauli7(:,:,2,ii)/2.0d0 + &
                ! mag_field_y8 * rspauli8(:,:,2,ii)/2.0d0

                ! rspauli_final(:,:,3,ii) = rspauli_final(:,:,3,ii) + mag_field_z1 * rspauli1(:,:,3,ii)/2.0d0 + &
                ! mag_field_z2 * rspauli2(:,:,3,ii)/2.0d0 + mag_field_z3 * rspauli3(:,:,3,ii)/2.0d0 + &
                ! mag_field_z4 * rspauli4(:,:,3,ii)/2.0d0 + mag_field_z5 * rspauli5(:,:,3,ii)/2.0d0 + &
                ! mag_field_z6 * rspauli6(:,:,3,ii)/2.0d0 + mag_field_z7 * rspauli7(:,:,3,ii)/2.0d0 + &
                ! mag_field_z8 * rspauli8(:,:,3,ii)/2.0d0

            enddo
            !!! 和方法2是一起的
            ! do ii = 1,rvecnum
            !     hops(:,:,ii)=hops(:,:,ii)+rspauli_final(:,:,1,ii)+rspauli_final(:,:,2,ii)+rspauli_final(:,:,3,ii) + &
            !     rspauli_ori(:,:,1,ii)+rspauli_ori(:,:,2,ii)+rspauli_ori(:,:,3,ii)
            ! enddo
        endif
    endif
     ! rspauli_final = 

    if(irank.eq.0)then
    open(555,file='new_rspauli1')
    ! open(987,file='rspauli_ori')
    ! open(567,file='rspauli_final')
    do ii = 1, rvecnum
        do i = 1, num_wann
            do j= 1, num_wann
                do dir =1,3
                    write(444,'(6I5,3F16.8)') irvec(1,ii) ,irvec(2,ii), irvec(3,ii), j,i, dir, rspauli1(j,i,dir,ii)
                    ! write(987,'(6I5,3F16.8)') irvec(1,ii) ,irvec(2,ii),irvec(3,ii), j,i, dir, rspauli_ori(j,i,dir,ii)
                    ! write(567,'(6I5,3F16.8)') irvec(1,ii) ,irvec(2,ii), irvec(3,ii), j,i, dir, rspauli_final(j,i,dir,ii)
                enddo
            enddo
        enddo
    enddo
    
    ! close(567)
    ! close(987)
    close(555)
    endif  

    if(irank.eq.0)then
    open(666,file='hopping2')
        do ii = 1, rvecnum
            do i = 1, num_wann
                do j= 1, num_wann
                    write(666,'(5I5,3F16.8)') irvec(1,ii) ,irvec(2,ii), irvec(3,ii), j,i,hops(j,i,ii)
                enddo
            enddo
        enddo
    close(666)
    endif

    !!! 在所有CPU上 allocate数组
    if(.not.allocated(nrpts))allocate(nrpts(rvecnum)) 
    call mpi_bcast(nrpts,rvecnum,MPI_INTEGER,0,mpi_comm_world,ierr)                             

    if(.not.allocated(irvec))allocate(irvec(3,rvecnum)) 
    length=3*rvecnum 
    call mpi_bcast(irvec,length,MPI_INTEGER,0,mpi_comm_world,ierr)                             

    if(.not.allocated(crvec))allocate(crvec(3,rvecnum))
    call mpi_bcast(crvec,length,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr) 
    
    if(.not.allocated(hops))allocate(hops(num_wann,num_wann,rvecnum)) 
    length=num_wann*num_wann*rvecnum 
    call mpi_bcast(hops,length,MPI_DOUBLE_COMPLEX,0,mpi_comm_world,ierr)                             

    if(.not.allocated(rspauli))allocate(rspauli(num_wann,num_wann,3,rvecnum))
    length=num_wann*num_wann*rvecnum*3
    call mpi_bcast(rspauli,length,MPI_DOUBLE_COMPLEX,0,mpi_comm_world,ierr) 
    
    if(.not.allocated(rspauli_ori))allocate(rspauli_ori(num_wann,num_wann,3,rvecnum))
    length=num_wann*num_wann*rvecnum*3
    call mpi_bcast(rspauli_ori,length,MPI_DOUBLE_COMPLEX,0,mpi_comm_world,ierr) 
    
    if(.not.allocated(rspauli_final))allocate(rspauli_final(num_wann,num_wann,3,rvecnum))
    length=num_wann*num_wann*rvecnum*3
    call mpi_bcast(rspauli_final,length,MPI_DOUBLE_COMPLEX,0,mpi_comm_world,ierr)

    if(.not.allocated(spin_sigma_x))allocate(spin_sigma_x(num_wann,num_wann))
    call mpi_bcast(spin_sigma_x,size(spin_sigma_x),MPI_DOUBLE_COMPLEX,0,mpi_comm_world,ierr)

    if(.not.allocated(spin_sigma_y))allocate(spin_sigma_y(num_wann,num_wann))
    call mpi_bcast(spin_sigma_y,size(spin_sigma_y),MPI_DOUBLE_COMPLEX,0,mpi_comm_world,ierr)

    if(.not.allocated(spin_sigma_z))allocate(spin_sigma_z(num_wann,num_wann))
    call mpi_bcast(spin_sigma_z,size(spin_sigma_z),MPI_DOUBLE_COMPLEX,0,mpi_comm_world,ierr)
                                                                                                                        
    cross(1)=amat(1,2)*amat(2,3)-amat(1,3)*amat(2,2) 
    cross(2)=amat(1,3)*amat(2,1)-amat(1,1)*amat(2,3) 
    cross(3)=amat(1,1)*amat(2,2)-amat(1,2)*amat(2,1) 
                                                                        
    volume=cross(1)*amat(3,1)+cross(2)*amat(3,2)+cross(3)*amat(3,3) 

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
    
    ! allocate(pauli(num_wann,num_wann,3))
    allocate(sigma_tensor_ahc_x(num_steps))
    allocate(sigma_tensor_ahc_y(num_steps))
    allocate(sigma_tensor_ahc_z(num_steps))
    allocate(sigma_tensor_ahc_mpi(3,num_steps))
    allocate(sigma_tensor_ahc_mpi2(3,num_steps))
    allocate(Omega_x(num_wann), Omega_y(num_wann), Omega_z(num_wann))
    allocate(Omega_x_t(num_wann), Omega_y_t(num_wann), Omega_z_t(num_wann))

    allocate(spin_texture(Nk1*Nk2*Nk3,num_wann,3))
    allocate(spin_texture_mpi(Nk1*Nk2*Nk3,num_wann,3))   
    allocate(spin_dir(3,num_wann))
    allocate(spin_dir_mpi(3,num_wann))
    allocate(spin_dir_mpi2(3,num_wann))
    ! allocate(spin_dir_y(num_wann))
    ! allocate(spin_dir_z(num_wann))    
    allocate(spin_sigma_x_t(num_wann))
    allocate(spin_sigma_y_t(num_wann))
    allocate(spin_sigma_z_t(num_wann))

    allocate(spindirx(num_wann))
    allocate(spindiry(num_wann))
    allocate(spindirz(num_wann))    
    allocate(spin_sigma_t_mpi(3,num_steps))
    allocate(spin_sigma_t_mpi2(3,num_steps))
    magnetic=0.0 
    occupation=0.0 
    conductivity=0.0 
    conductivity13=0.0 
    conductivity23=0.0 
    conductivity_ahe=0.0 
    conductivity13_ahe=0.0 
    conductivity23_ahe=0.0 
    sigma_tensor_ahc_mpi = 0.0               
    sigma_tensor_ahc_mpi2 = 0.0                               

    spin_dir = 0.0
    spin_dir_mpi = 0.0
    spin_dir_mpi2 = 0.0
    allocate(ham(num_wann,num_wann)) 
    allocate(pauli(num_wann,num_wann,3)) 
    allocate(momentum (num_wann,num_wann,3)) 
    allocate(momentum2(num_wann,num_wann,3)) 

    allocate(eigvals(num_wann)) 
    allocate(eigvals_x(num_wann))
    allocate(eigvals_y(num_wann))
    allocate(eigvals_z(num_wann))
    allocate(eigvecs(num_wann,num_wann)) 
    allocate(eigvecs_dag(num_wann,num_wann))
    allocate(eigvecs_x(num_wann,num_wann))
    allocate(eigvecs_y(num_wann,num_wann))
    allocate(eigvecs_z(num_wann,num_wann))  
    allocate(mat_temp(num_wann,num_wann))
    allocate(spin_sigma_temp(num_wann,num_wann))    
    allocate(spin_sigma_x_comp(num_wann,num_wann))
    allocate(spin_sigma_y_comp(num_wann,num_wann))
    allocate(spin_sigma_z_comp(num_wann,num_wann))  

    lwork=12.0*num_wann 

    allocate( work(lwork) ) 
    allocate( rwork(17*num_wann) ) 
    allocate( iwork(15*num_wann) ) 
    allocate( ifail(15*num_wann) ) 

    spin_texture=0.0d0
    spin_texture_mpi=0.0d0
    time_start = 0.0
    knv3= Nk1*Nk2*Nk3
    call now(time_start)
    do ik= 1+ irank, knv3, isize
        if (irank .eq. 0 .and. mod(ik/isize, 1) .eq. 0) then
            call now(time_end) 
            write(*, '(a, i18, "/", i18, a, f10.2, "s")') 'ik/knv3', &
            ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/isize
            time_start= time_end
        endif

        ikx= (ik-1)/(nk2*nk3)+1
        iky= ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
        ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
        kpoints = K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1)  &
                + K3D_vec2_cube*(iky-1)/dble(nk2) + K3D_vec3_cube*(ikz-1)/dble(nk3)

        ham=0.0
        pauli=0.0d0

        do ii=1,rvecnum 
            ix=irvec(1,ii) 
            iy=irvec(2,ii) 
            iz=irvec(3,ii) 
            phas=     iy*kpoints(2) 
            phas=phas+iz*kpoints(3) 
            phas=phas+ix*kpoints(1) 
            phas=phas*twopi 
            fac=cmplx(cos(phas),sin(phas)) 
            ham(:,:)=ham(:,:)+fac*hops(:,:,ii)/nrpts(ii)                                            
            
            ! pauli(:,:,1)=pauli(:,:,1)+fac*rspauli_final(:,:,1,ii)
            ! pauli(:,:,2)=pauli(:,:,2)+fac*rspauli_final(:,:,2,ii)
            ! pauli(:,:,3)=pauli(:,:,3)+fac*rspauli_final(:,:,3,ii)
        enddo 
 
        call zheevx('V','A','U',num_wann,ham,num_wann,vl,vu,1,num_wann,abstol,ne,&
                    eigvals,eigvecs,num_wann,work,lwork,rwork,iwork,ifail,info)                        
        if(info.ne.0)stop 'zheevx'                              

        eigvecs_dag=conjg(transpose(eigvecs))

        ! call zheevx('V','A','U',num_wann,pauli(:,:,1),num_wann,vl,vu,1,num_wann,abstol,ne,&
            ! eigvals_x,eigvecs_x,num_wann,work,lwork,rwork,iwork,ifail,info)
        ! call zheevx('V','A','U',num_wann,pauli(:,:,2),num_wann,vl,vu,1,num_wann,abstol,ne,&
            ! eigvals_y,eigvecs_y,num_wann,work,lwork,rwork,iwork,ifail,info)
        ! call zheevx('V','A','U',num_wann,pauli(:,:,3),num_wann,vl,vu,1,num_wann,abstol,ne,&
            ! eigvals_z,eigvecs_z,num_wann,work,lwork,rwork,iwork,ifail,info)

        spin_sigma_x_comp = 0.0d0
        spin_sigma_y_comp = 0.0d0
        spin_sigma_z_comp = 0.0d0

        spin_sigma_x_comp = MATMUL(eigvecs_dag(:,:),MATMUL(spin_sigma_x,eigvecs(:,:)))
        spin_sigma_y_comp = MATMUL(eigvecs_dag(:,:),MATMUL(spin_sigma_y,eigvecs(:,:)))
        spin_sigma_z_comp = MATMUL(eigvecs_dag(:,:),MATMUL(spin_sigma_z,eigvecs(:,:)))
  
        ! spin_sigma_x_comp=(MATMUL(MATMUL(eigvecs,spin_sigma_x),eigvecs_dag)+MATMUL(MATMUL(eigvecs_dag,spin_sigma_x),eigvecs))/2
        ! spin_sigma_y_comp=(MATMUL(MATMUL(eigvecs,spin_sigma_y),eigvecs_dag)+MATMUL(MATMUL(eigvecs_dag,spin_sigma_y),eigvecs))/2
        ! spin_sigma_z_comp=(MATMUL(MATMUL(eigvecs,spin_sigma_z),eigvecs_dag)+MATMUL(MATMUL(eigvecs_dag,spin_sigma_z),eigvecs))/2

        ! spin_sigma_x_comp = MATMUL(MATMUL(eigvecs_x,spin_sigma_x),conjg(transpose(eigvecs_x))
        ! spin_sigma_y_comp = MATMUL(MATMUL(eigvecs_y,spin_sigma_y),conjg(transpose(eigvecs_y))
        ! spin_sigma_z_comp = MATMUL(MATMUL(eigvecs_z,spin_sigma_z),conjg(transpose(eigvecs_z))

        ! spin_sigma_x_comp = 0.0d0
        ! spin_sigma_y_comp = 0.0d0
        ! spin_sigma_z_comp = 0.0d0
        ! spin_sigma_temp = 0.0d0        
        ! spin_sigma_x_comp = MATMUL(MATMUL(eigvecs,pauli(:,:,1)),eigvecs_dag)
        ! spin_sigma_y_comp = MATMUL(MATMUL(eigvecs,pauli(:,:,2)),eigvecs_dag)
        ! spin_sigma_z_comp = MATMUL(MATMUL(eigvecs,pauli(:,:,3)),eigvecs_dag)
        
        ! call mat_mul(num_wann,spin_sigma_x,eigvecs,spin_sigma_temp)
        ! call mat_mul(num_wann,eigvecs_dag,spin_sigma_temp,spin_sigma_x_comp)
        ! call mat_mul(num_wann,spin_sigma_y(:,:),eigvecs,spin_sigma_temp)
        ! call mat_mul(num_wann,eigvecs_dag,spin_sigma_temp,spin_sigma_y_comp(:,:))
        ! call mat_mul(num_wann,spin_sigma_z(:,:),eigvecs,spin_sigma_temp)
        ! call mat_mul(num_wann,eigvecs_dag,spin_sigma_temp,spin_sigma_z_comp(:,:))

        call zheevx('V','A','U',num_wann,spin_sigma_x_comp,num_wann,vl,vu,1,num_wann,abstol,ne,&
            eigvals_x,eigvecs_x,num_wann,work,lwork,rwork,iwork,ifail,info)
        call zheevx('V','A','U',num_wann,spin_sigma_y_comp,num_wann,vl,vu,1,num_wann,abstol,ne,&
            eigvals_y,eigvecs_y,num_wann,work,lwork,rwork,iwork,ifail,info)
        call zheevx('V','A','U',num_wann,spin_sigma_z_comp,num_wann,vl,vu,1,num_wann,abstol,ne,&
            eigvals_z,eigvecs_z,num_wann,work,lwork,rwork,iwork,ifail,info)


        spin_texture(ik,:,1)=eigvals_x
        spin_texture(ik,:,2)=eigvals_y
        spin_texture(ik,:,3)=eigvals_z

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
            momentum(:,:,1)= momentum(:,:,1)+fac*crvec(1,ii)*hops(:,:,ii)/nrpts(ii) 
            momentum(:,:,2)= momentum(:,:,2)+fac*crvec(2,ii)*hops(:,:,ii)/nrpts(ii)
            momentum(:,:,3)= momentum(:,:,3)+fac*crvec(3,ii)*hops(:,:,ii)/nrpts(ii)
        enddo      
        momentum2= 0.0
        call mat_mul(num_wann,momentum(:,:,1),eigvecs,mat_temp)
        call mat_mul(num_wann,eigvecs_dag,mat_temp,momentum2(:,:,1))
        call mat_mul(num_wann,momentum(:,:,2),eigvecs,mat_temp)
        call mat_mul(num_wann,eigvecs_dag,mat_temp,momentum2(:,:,2))
        call mat_mul(num_wann,momentum(:,:,3),eigvecs,mat_temp)
        call mat_mul(num_wann,eigvecs_dag,mat_temp,momentum2(:,:,3))       

        ! if (irank == 0) then
            ! write(*,*) "momentum2 = ", momentum2
        ! endif
        ! if (kbT == 0) then
            ! Beta_fake = 10000
        !else
        Beta_fake = 1.0/kbT
        ! endif      
        Omega_x=0d0;Omega_y=0d0; Omega_z=0d0
        ! spin_dir = 0.0d0
        ! spin_dir_mpi = 0.0d0
        ! spin_dir_y = 0.0d0
        ! spin_dir_z = 0.0d0
        do m=1,num_wann
            do n=1,num_wann
                if (abs(eigvals(m)-eigvals(n)) .gt. 0.00000001d0) then
                    Omega_z(m)= Omega_z(m)-aimag(momentum2(m,n,1)*conjg(momentum2(m,n,2)))/(eigvals(m)-eigvals(n))**2
                    Omega_y(m)= Omega_y(m)-aimag(momentum2(m,n,1)*conjg(momentum2(m,n,3)))/(eigvals(m)-eigvals(n))**2
                    Omega_x(m)= Omega_x(m)-aimag(momentum2(m,n,2)*conjg(momentum2(m,n,3)))/(eigvals(m)-eigvals(n))**2
                endif
                ! if (m.eq.n) then
                    ! spin_dir(1,m) = spin_sigma_x_comp(m,m)
                    ! spin_dir(2,m) = spin_sigma_y_comp(m,m)
                    ! spin_dir(3,m) = spin_sigma_z_comp(m,m)
                ! endif
            enddo !n
        enddo  !m

        ! do m=1,num_wann
            ! spin_dir(1,m) = eigvals_x(m)
            ! spin_dir(2,m) = eigvals_y(m)
            ! spin_dir(3,m) = eigvals_z(m)
        ! enddo
        
         
        ! spin_dir_mpi = spin_dir
         
        ! if (irank == 0)then
        ! write(*,*) "spindir_x",spin_dir(1,:)
        ! write(*,*) "spindir_y",spin_dir(2,:)
        ! write(*,*) "spindir_z",spin_dir(3,:)
        ! endif
        Omega_x = Omega_x*2.0
        Omega_y = Omega_y*2.0
        Omega_z = Omega_z*2.0   
         
        do step=1,num_steps                                                             
            fermienergy(step) = fermi_min + (fermi_max-fermi_min)*real(step)/real(num_steps)        
            mu = fermienergy(step)
            do m= 1, num_wann
                ! Omega_x_t(m)= Omega_x(m)*(-1.0/(1*300))*((eigvals(m)-mu)*fermi(eigvals(m)-mu, Beta_fake)+(kbT*log(1+exp(-(eigvals(m)-mu)/kbT))))
                ! Omega_y_t(m)= Omega_y(m)*(-1.0/(1*300))*((eigvals(m)-mu)*fermi(eigvals(m)-mu, Beta_fake)+(kbT*log(1+exp(-(eigvals(m)-mu)/kbT))))
                ! Omega_z_t(m)= Omega_z(m)*(-1.0/(1*300))*((eigvals(m)-mu)*fermi(eigvals(m)-mu, Beta_fake)+(kbT*log(1+exp(-(eigvals(m)-mu)/kbT))))
                Omega_x_t(m) = Omega_x(m)*fermi(eigvals(m)-mu, Beta_fake)
                Omega_y_t(m) = Omega_y(m)*fermi(eigvals(m)-mu, Beta_fake)
                Omega_z_t(m) = Omega_z(m)*fermi(eigvals(m)-mu, Beta_fake)

                ! spin_sigma_x_t(m)=spin_sigma_x_comp(m,m)*fermi(eigvals(m)-mu, Beta_fake)
                ! spin_sigma_y_t(m)=spin_sigma_y_comp(m,m)*fermi(eigvals(m)-mu, Beta_fake)
                ! spin_sigma_z_t(m)=spin_sigma_z_comp(m,m)*fermi(eigvals(m)-mu, Beta_fake)

            enddo
                sigma_tensor_ahc_mpi(1, step)= sigma_tensor_ahc_mpi(1, step) + real(sum(Omega_x_t))
                sigma_tensor_ahc_mpi(2, step)= sigma_tensor_ahc_mpi(2, step) + real(sum(Omega_y_t))
                sigma_tensor_ahc_mpi(3, step)= sigma_tensor_ahc_mpi(3, step) + real(sum(Omega_z_t))

                ! spin_sigma_t_mpi(1,step)= spin_sigma_t_mpi(1, step) + real(sum(spin_sigma_x_t))
                ! spin_sigma_t_mpi(2,step)= spin_sigma_t_mpi(2, step) + real(sum(spin_sigma_y_t))
                ! spin_sigma_t_mpi(3,step)= spin_sigma_t_mpi(3, step) + real(sum(spin_sigma_z_t))
        enddo ! end step
    enddo ! end kpt

    num_steps_tot =num_steps*3
    num_spin = num_wann*3                                             
    call MPI_REDUCE(sigma_tensor_ahc_mpi,sigma_tensor_ahc_mpi2,num_steps_tot,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
    ! call MPI_REDUCE(spin_dir_mpi,spin_dir_mpi2,num_spin,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
    ! call MPI_REDUCE(spin_sigma_t_mpi,spin_sigma_tmpi2,num_wann*3,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
    call MPI_REDUCE(spin_texture,spin_texture_mpi,size(spin_texture),MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
    
    sigma_tensor_ahc_mpi=sigma_tensor_ahc_mpi2
    sigma_tensor_ahc_mpi=sigma_tensor_ahc_mpi/knv3/volume*condq*1.0e8*twopi

    conductivity_ahe=sigma_tensor_ahc_mpi(3,:) 
    conductivity13_ahe=sigma_tensor_ahc_mpi(2,:)   
    conductivity23_ahe=sigma_tensor_ahc_mpi(1,:)
   
    ! spindirx = spin_dir_mpi2(1,:)
    ! spindiry = spin_dir_mpi2(2,:)
    ! spindirz = spin_dir_mpi2(3,:)

    if(irank.eq.0)then 
        open(123,file='output_ahe_condquant',recl=10000) 
        do step=1,num_steps 
            write(123,*)"fermienergy=",fermienergy(step)                                                          
            write(123,*)"occupation=",occupation(step)                                
            ! write(123,*)"conductivity=",conductivity_ahe(step)/condq*amat(3,3)*1.0e-8  !2D direct convert to integer
            ! write(123,*)"conductivity13=",conductivity13_ahe(step)/condq*amat(3,3)*1.0e-8                 
            ! write(123,*)"conductivity23=",conductivity23_ahe(step)/condq*amat(3,3)*1.0e-8
            write(123,*)"conductivity=",conductivity_ahe(step)   !unit is same with wannier_tool  S/cm
            write(123,*)"conductivity13=",conductivity13_ahe(step)
            write(123,*)"conductivity23=",conductivity23_ahe(step)                      
            write(123,*)"************************************"                                                  
        enddo 
        close(123)      
         
        open(456,file='spin_texture',recl=10000) 
        ! write(456,*) spin_dir_mpi
        do ik = 1,knv3
            ikx=(ik-1)/(nk2*nk3)+1
            iky=((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
            ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)  
            do m=1,num_wann
                write(456,*)"kx= ",(ikx-1)/dble(nk1),"ky= ",(iky-1)/dble(nk2),"kz= ",(ikz-1)/dble(nk3)
                write(456,*)"m= ", m
                write(456,*)"spindir_x=",spin_texture_mpi(ik,m,1)
                write(456,*)"spindir_y=",spin_texture_mpi(ik,m,2)
                write(456,*)"spindir_z=",spin_texture_mpi(ik,m,3)
                write(456,*)"************************************" 
            enddo
!            write(456,*)"band_index",m
!            write(456,*)"spindir_x",spin_dir_mpi2(1,m)
!           write(456,*)"spindir_y",spin_dir_mpi2(2,m)
!           write(456,*)"spindir_z",spin_dir_mpi2(3,m)
!           write(456,*)"************************************"
        enddo
        close(456)
    endif 
    call mpi_barrier(mpi_comm_world,ierr) 
    call MPI_Finalize(ierr) 

end program anomalous_nernst_effect


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

subroutine now(time_now)

    implicit none
    integer   :: time_new(8)
    real      :: time_now
    call Date_and_time(values=time_new)
    time_now= time_new(3)*24*3600+time_new(5)*3600+&
            time_new(6)*60+time_new(7)+time_new(8)/1000d0  
    return
end subroutine now

subroutine mat_mul(nmatdim,A,B,C)

    implicit none
    integer,intent(in) :: nmatdim

    complex(kind(1.0d0)) :: ALPHA
    complex(kind(1.0d0)) :: BETA

    complex(kind(1.0d0)), intent(in)  :: A(nmatdim ,nmatdim)
    complex(kind(1.0d0)), intent(in)  :: B(nmatdim ,nmatdim)
    complex(kind(1.0d0)), intent(out) :: C(nmatdim,nmatdim)

    ALPHA=1.0d0
    BETA=0.0D0

    C(:,:)=(0.0d0,0.0d0)

    call ZGEMM('N','N',nmatdim,nmatdim,nmatdim,ALPHA,A,nmatdim,B,nmatdim,BETA,C,nmatdim)

    return
end subroutine mat_mul


