program bidebug2

         implicit none

         complex,allocatable:: hops   (:,:,:)
         complex,allocatable:: rsnabla(:,:,:,:)
         complex,allocatable:: rspauli(:,:,:,:)
         complex,allocatable:: pauli(:,:,:)
         complex,allocatable:: paulifft(:,:,:)
         complex,allocatable:: paulifft2(:,:,:)
         real               :: rdum,idum
         integer            :: ix,iy,iz,band1,band2,h,num_wann
         integer            :: ik1,ik2,ik3
         real               :: twopi
         real               :: phas
         complex            :: fac,fac2
         complex,allocatable:: ham(:,:)
         real               :: vl,vu
         integer            :: ne,j
         real               :: abstol
         real,allocatable   :: eigvals(:)
         complex,allocatable:: eigvecs(:,:)
         integer            :: info
         complex,allocatable:: work(:)
         integer            :: lwork
         integer,allocatable:: iwork(:)
         real,allocatable   :: rwork(:)
         integer,allocatable:: ifail(:)
         real               :: kpoints(3)
         real               :: scale
         integer            :: maxhopx2,maxhopy2,maxhopz2,dire
         real,allocatable   :: fermienergy(:)
         real,allocatable   :: deviation(:)
         integer            :: grid,i1,i2,i3,i4,orb
         real,allocatable   :: conductivity(:),conductivity2(:)
         real,allocatable   :: conductivity13(:),conductivity23(:)
         real,allocatable   :: conductivity_ahe(:)
         real,allocatable   :: conductivity13_ahe(:),conductivity23_ahe(:)
         real,allocatable   :: conductivity_fsur(:)
         real,allocatable   :: conductivity13_fsur(:)
         real,allocatable   :: conductivity23_fsur(:)
         integer            :: ierr,isize,irank,kp1,kp2,kp3
         integer            :: ix1,ix2,ix3,num_occ
         complex,allocatable:: momentum(:,:,:)
         complex,allocatable:: momentum2(:,:,:)
         complex,allocatable:: spinmomentum(:,:,:)
         complex,allocatable:: spinmomentum2(:,:,:)
         complex            :: berry

         integer,allocatable:: sortarray(:)
         integer            :: n1,n2,n3,n4,dir
         real,allocatable   :: occupation(:),occupation2(:)
         integer,allocatable   :: nrpts(:)
         real               :: occupation_number
         integer            :: step,i,ii,num_steps
         real               :: fermi_min,fermi_max
         logical            :: l_tb,l_nabla,l_bfield
         real               :: bfield(3)
         integer            :: rvecnum,num_lines,length
         integer,allocatable:: irvec(:,:)
         real,allocatable   :: kpts(:,:)
         real               :: kder,amat(3,3)
         real               :: volume,cross(3)
         real               :: bohrincm,condq
         integer            :: maxdim,num_kpts
         real               :: minenerg,maxenerg,gamma
         logical            :: l_bandstruc,l_fermisurf
         real,allocatable   :: magnetic(:,:)

         abstol=2.0*tiny(abstol)
         twopi=2*3.141592654
         bohrincm=0.529*1.e-8
         condq=38.74*1.e-6
         pi = 3.141592654

            open(300,file='ahe_inp')
            read(300,*)amat(1,:)
            read(300,*)amat(2,:)
            read(300,*)amat(3,:)
            read(300,*)fermi_min,fermi_max,num_steps
            read(300,*)grid
            read(300,*)maxdim
            read(300,*)l_tb,l_nabla
            read(300,*)l_bfield,bfield
            read(300,*)occupation_number
            read(300,*)l_fermisurf,gamma
            close(300)

            open(200,file='hopping.1')
            num_lines=0
            num_wann=0
            do
               read(200,fmt=*,end=311)ix,iy,iz,band1,band2,rdum,idum
               num_lines=num_lines+1
               num_wann=max(num_wann,band1)
            enddo
     311    continue
            rvecnum=num_lines/(num_wann*num_wann)
            write(*,*)"num_lines=",num_lines
            write(*,*)"rvecnum=",rvecnum
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
     300    continue
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
     500   continue
           close(400)
           write(*,*) 'rvecnum',rvecnum
           allocate(nrpts(rvecnum))
           open(14,file='nrpts_inp')
           do j=1,rvecnum/15
                   read(14,'(15I5)') (nrpts(15*(j-1)+i) ,i=1,15)
           enddo
           read(14,'(<mod(rvecnum,15)>I5)') (nrpts(15*(rvecnum/15)+i),     &
        &                                   i=1,mod(rvecnum,15))
           close(14)
         write(*,*) nrpts
         write(*,*) 'size',size(nrpts)



         cross(1)=amat(1,2)*amat(2,3)-amat(1,3)*amat(2,2)
         cross(2)=amat(1,3)*amat(2,1)-amat(1,1)*amat(2,3)
         cross(3)=amat(1,1)*amat(2,2)-amat(1,2)*amat(2,1)

         volume=cross(1)*amat(3,1)+cross(2)*amat(3,2)+cross(3)*amat(3,3)
         if(maxdim.eq.2)then
            volume=volume/amat(3,3)/(0.529e-8)
         endif

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

         magnetic=0.0
         occupation=0.0
         conductivity=0.0
         conductivity13=0.0
         conductivity23=0.0
         conductivity_ahe=0.0
         conductivity13_ahe=0.0
         conductivity23_ahe=0.0

         allocate(        ham(num_wann,num_wann))
         allocate(    pauli(num_wann,num_wann,3))
         allocate(momentum (num_wann,num_wann,3))
         allocate(momentum2(num_wann,num_wann,3))

         allocate(spinmomentum (num_wann,num_wann,3))
         allocate(spinmomentum2(num_wann,num_wann,3))

         allocate(eigvals(num_wann))
         allocate(eigvecs(num_wann,num_wann))
         print*,"num_wann=",num_wann

         allocate( work(lwork) )
         allocate( rwork(17*num_wann) )
         allocate( iwork(15*num_wann) )
         allocate( ifail(15*num_wann) )

         integer           :: theta_num
         real              :: mag_strength
         theta_num = 10
         real,allocatable  :: theta(:)
         real              :: phi(8)
         allocate(theta(theta_num))
         theta = 0.0

         do i=0,theta_num
            theta(i) =  i * (pi/2)/theta_num

         end do

         do j=1,8
            phi(j) = (j-1) * (pi/4)
         end do

         real, allocatable :: mag_field(:,:,:)
         allocate(mag_field(3,8,theta_num))

         do phi_j = 1,8
            do theta_i = 0,theta_num
               mag_field(1,phi_j,theta_i) = mag_strength * cos(theta(theta_i))
               mag_field(2,phi_j,theta_i) = mag_strength * sin(theta(theta_i)) * cos(phi(phi_j))
               mag_field(3,phi_j,theta_i) = mag_strength * sin(theta(theta_i)) * sin(phi(phi_j))
            end do
         end do

         real, allocatable :: hops_in_field(:,:,:,:,:)

         allocate(hops_in_field(theta_num, 8,num_wann,num_wann,rvecnum))

         do num_mag_atom=1,8
            do ii=1,rvecnum
               do i=1,num_wann
               do j=1,num_wann
                  hops_in_field(num_mag_atom,j,i,ii)=hops(num_mag_atom,j,i,ii)+
        &                   bfield(1)*rspauli(j,i,1,ii)*0.5   &
        &                      +bfield(2)*rspauli(j,i,2,ii)*0.5               &
        &                      +bfield(3)*rspauli(j,i,3,ii)*0.5

         end program bidebug2
