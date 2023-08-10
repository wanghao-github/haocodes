      program hoti_version_hao
!*************************************************
! mpiifort -CB -r8 ho_hamtricell2.F -qmkl -o 2tricell2.x
!*************************************************
      implicit none

      complex,allocatable   :: hops(:,:,:),ham(:,:),hij(:,:,:)
      integer,allocatable   :: irvec(:,:),numbers(:)
      integer,allocatable   :: nrpts(:),excludeup(:),excludedown(:)
      real,allocatable      :: eigvals(:),btompositions(:,:)
      real,allocatable      :: wannier_center(:,:),atompos(:,:)
      real,allocatable      :: atompositions(:,:),testat(:,:)
      real,allocatable      :: wannier_center_real(:,:)
      real                  :: rdum, idum,fermi,amat(3,3)
      integer               :: ierr, isize, irank, nwannexup,nwannexdown
      integer               :: nslab1,nslab2,num_lines,num_wann,ix,iy,iz
      integer               :: band1,band2,rvecnum,i,j,num_wann_y,preadmax
      integer               :: preadmin,hdim,up,down
      integer               :: ijmax,numa,numb,h_dim,htotal
      integer,allocatable   :: atomnumbertot(:),atomorbitaltot(:)
      integer,allocatable   :: atomnumbertot2(:),atomorbitaltot2(:)     
      integer,allocatable   :: atomnumber(:),atomorbital(:)
      integer,allocatable   :: layera(:),layerb(:)
      real,allocatable      :: posreal(:,:)
      real,allocatable      :: positions(:,:)
      real                  :: t1,t2,t3,t4,t5,t6
      integer               :: kmin,kmax
      logical               :: l_vec

      INCLUDE 'mpif.h'
      integer :: stt(MPI_STATUS_SIZE)
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_RANK (MPI_COMM_WORLD,irank,ierr)
      CALL MPI_COMM_SIZE (MPI_COMM_WORLD,isize,ierr)
     
      if(irank.eq.0)then
         open(300,file='inp')
         read(300,*) fermi
         read(300,*)nslab1,nslab2

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
         write(*,*)"num_wann=",num_wann
         write(*,*)"rvecnum=",rvecnum
         allocate( hops(num_wann,num_wann,rvecnum) )
         allocate( irvec(3,rvecnum) )
         allocate( numbers(num_wann))

         read(300,*) numbers(:)  ! atomnumberarray(:)
         read(300,*) amat(1,:)
         read(300,*) amat(2,:)
         read(300,*) amat(3,:)
         read(300,*) t1,t2,t3,t4,t5,t6
         read(300,*) kmin,kmax,l_vec
         close(300)

         write(*,*) 'da,db,d1,d2,d3:',t1,t2,t3,t4,t5,t6
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
            hops(band1,band2,rvecnum)=cmplx(rdum,idum)
         enddo
 300     continue
         close(200)
   
         write(*,*) 'rvecnum',rvecnum
         allocate(nrpts(rvecnum))
         open(14,file='nrpts_inp')
         do j=1,rvecnum/15
            read(14,'(15I5)') (nrpts(15*(j-1)+i) ,i=1,15)
         enddo
         read(14,'(<mod(rvecnum,15)>I5)') (nrpts(15*(rvecnum/15)+i),i=1,mod(rvecnum,15))
         close(14)

         write(*,*) 'size',size(nrpts)
         allocate(wannier_center(3,num_wann))
         wannier_center  = 0.d0
         open(400,file='position')
         do j=1,num_wann/2
            read(400, *) wannier_center(:,j)
         enddo
         close(400)
         do j=1,num_wann/2
            wannier_center(:,j+num_wann/2)=wannier_center(:,j)
         enddo
         numa = nslab1
         numb = nslab2
         ijmax = 12
         allocate(wannier_center_real(3,num_wann))
         wannier_center_real=0.0
         do i =1, num_wann
            wannier_center_real(1, num_wann)=                &
               amat(1,1)*wannier_center(1,num_wann)+         &
               amat(2,1)*wannier_center(2,num_wann)+         &
               amat(3,1)*wannier_center(3,num_wann)
            wannier_center_real(2, num_wann)=                 &
               amat(1,2)*wannier_center(1,num_wann)+          &
               amat(2,2)*wannier_center(2,num_wann)+          &
               amat(3,2)*wannier_center(3,num_wann)           
            wannier_center_real(3, num_wann)=                 &
               amat(1,3)*wannier_center(1,num_wann)+          &
               amat(2,3)*wannier_center(2,num_wann)+          &
               amat(3,3)*wannier_center(3,num_wann)           
         enddo
          
         write(*,*) 'num_wann= ', num_wann, 'nslab1 = ',nslab1
         write(*,*) 'nwannexup=',nwannexup,'nwannexdown=',nwannexdown
         write(*,*) 'excludeup=',excludeup,'excludedown',excludedown
         write(*,*)'numa,numb=',numa,numb

         htotal = num_wann*numa*numb
         hdim= num_wann*numa*numb
          
         allocate(atomnumber(hdim))
         ! allocate(hij( -4:4,hdim,hdim))
         write(*,*) 'its time for get dim'
         write(*,*) 'its time for a tri_cut program'
         allocate(atompositions(3,h_dim))
         allocate(ham(h_dim,h_dim))
         allocate(atomorbital(h_dim))

         call tri_cut(nwannexup,nwannexdown,num_wann,nslab1,rvecnum,   &
      excludeup,excludedown,numbers,irvec,hops,nrpts,wannier_center,   &
      atompositions,atomnumber,preadmax,preadmin,htotal,up,down,       &
      amat,numb,h_dim,ham,atomorbital,t1,t2,t3,t4,t5,t6)

         if(allocated(excludeup))       deallocate(excludeup)
         if(allocated(excludedown))     deallocate(excludedown)

         open(99, file='orbital.dat',recl = 10000)
         do i = 1, h_dim
            write(99,*) i, atompositions(1,i),atompositions(2,i),atompositions(3,i)
         enddo
         close(99)

         open(119,file='atomorbitals.dat')
         do i=1,h_dim
            write(119,*) i, atomorbital(i)
         enddo
         close(119)
         write(*,*) "atom dat yes!"

         allocate(posreal(3,h_dim))
         open(199, file='orbitalreal.dat',recl = 10000)
         do i = 1, h_dim
            posreal(1, i)=                         &
               amat(1,1)*atompositions(1,i)+       &
               amat(2,1)*atompositions(2,i)+       &
               amat(3,1)*atompositions(3,i)
            posreal(2, i)=                         &
               amat(1,2)*atompositions(1,i)+       &
               amat(2,2)*atompositions(2,i)+       &
               amat(3,2)*atompositions(3,i)
            posreal(3, i)=                         &
               amat(1,3)*atompositions(1,i)+       &
               amat(2,3)*atompositions(2,i)+       &
               amat(3,3)*atompositions(3,i)  
            write(199,*) i, posreal(1,i),posreal(2,i),posreal(3,i)
         enddo
         close(199)
         write(*,*) "orbitalreal dat yes!"

         if(allocated(atomnumber))      deallocate(atomnumber)
         if(allocated(atompositions))   deallocate(atompositions)
         if(allocated(hij))             deallocate(hij)
         if(allocated(hops))            deallocate(hops)
         if(allocated(irvec))           deallocate(irvec)
         if(allocated(nrpts))           deallocate(nrpts)
         if(allocated(wannier_center))  deallocate(wannier_center)
         if(allocated(btompositions))   deallocate(btompositions)
         if(allocated(posreal))         deallocate(posreal)
 
         if(allocated(atomnumbertot))   deallocate(atomnumbertot)
         if(allocated(atomnumbertot2))  deallocate(atomnumbertot2)
         if(allocated(atomorbitaltot))  deallocate(atomorbitaltot)
         if(allocated(atomorbitaltot2)) deallocate(atomorbitaltot2)
         if(allocated(layera))          deallocate(layera)
         if(allocated(layerb))          deallocate(layerb)
         if(allocated(positions))       deallocate(positions)
         
        allocate(eigvals(h_dim))
        write(*,*) 'fermi level is set as ', fermi
        if(l_vec) then
          call solve_ham_vec(h_dim,ham,eigvals,fermi,kmin,kmax)
        else
          call solve_ham(h_dim,ham,eigvals,fermi)
        endif

        if(allocated(ham))             deallocate(ham)
        if(allocated(eigvals))         deallocate(eigvals)
        if(allocated(excludeup))       deallocate(excludeup)
        if(allocated(excludedown))     deallocate(excludedown)
        if(allocated(numbers))         deallocate( numbers)

        endif

      call mpi_barrier(mpi_comm_world,ierr)
      call MPI_Finalize(ierr)

      end program hoti_version_hao

      subroutine tri_cut(nwannexup,nwannexdown,num_wann,numa,rvecnum,  &
       excludeup,excludedown,numbers,irvec,hops,nrpts,wanniercenters,  &
       atompos,atomnumber,preadmax,preadmin,htotal,up,down,            &
       amat,numb,h_dim,ham,atomorbital,t1,t2,t3,t4,t5,t6)
        implicit none
        integer,intent(in)    :: nwannexup,nwannexdown,num_wann,numa,numb
        integer,intent(in)    :: rvecnum,up,down
        integer,intent(in)    :: excludeup(up),excludedown(down)
        integer,intent(in)    :: numbers(num_wann),irvec(3,rvecnum) 
        integer,intent(in)    :: nrpts(rvecnum)
        integer               :: htotal
        integer               :: h_dim,i,i1,ii,ix,iy,iz,zvalue,hi,hj,j,i2
        integer               :: spreadmin,spreadmax,preadmax,preadmin
        integer               :: ijmax,m,n
        integer               :: atomnumbertot(htotal),atomorbitaltot(htotal)
        integer               :: atomnumbertot2(htotal),atomorbitaltot2(htotal)
        integer               :: layerdir(htotal),cut_num,layeratot(htotal)
        integer               :: layerbtot(htotal),atomnumber(htotal)
        integer               :: atomorbital(h_dim)
        integer,allocatable   :: layera(:),layerb(:)
        real,intent(in)       :: amat(3,3)
        integer               :: layerx,layery
        real                  :: ra,tol
        real,intent(in)       :: wanniercenters(3,num_wann)
        complex,intent(in)    :: hops(num_wann,num_wann,rvecnum)
        complex,allocatable   :: tmp(:,:,:),tmpp(:,:,:,:)
        real,allocatable      :: testat(:,:)
        real                  :: atompos(3,h_dim)
        real,allocatable      :: positions(:,:)
        real                  :: positionmed(3,htotal)
        real                  :: t1,t2,t3,t4,t5,t6
        complex               :: ham(h_dim,h_dim)

        tol=1.0e-8
        cut_num=0
        atomnumbertot2=0
        atomorbitaltot2=0
        layeratot=0
        layerbtot=0
        ijmax=12
        write(*,*) ' in tri_cut htotal = ' , htotal
        allocate(positions(3,htotal))
        write(*,*) 'allocate array positions yes! ' 
        positions=0.0
        ra=sqrt(amat(1,1)**2+amat(1,2)**2+amat(1,3)**2)
        write(*,*) 'Ra = ' , amat(1,1),amat(1,2),amat(1,3)
        write(*,*) 'Ra = ' , amat(2,1),amat(2,2),amat(2,3)
        write(*,*) 'Ra = ' , amat(3,1),amat(3,2),amat(3,3)
        write(*,*) 'abs(Ra) = ' , ra
        positionmed=0.0
         i=1
         i2=1
        do m=1,numa
           do n=1,numb
              i1=1
           do i=1,num_wann
              atomnumbertot(i) = numbers(i1)
              atomorbitaltot(i) = i1
              positionmed(1,i) = wanniercenters(1,i1)+m-1
              positionmed(2,i) = wanniercenters(2,i1)+n-1
              positionmed(3,i) = wanniercenters(3,i1)
              layerdir(i) = (i-1)/num_wann+1
              if(positionmed(2,i) > -1*positionmed(1,i)+t1    &
      .and. positionmed(2,i) < 0.5 * positionmed(1,i)+t2      &
      .and. positionmed(2,i) >   2 * positionmed(1,i)+t3      &
      ) then  

            positions(1,i2) = positionmed(1,i)
            positions(2,i2) = positionmed(2,i)
            positions(3,i2) = positionmed(3,i)

            atomnumbertot2(i2)  = atomnumbertot(i)
            atomorbitaltot2(i2) = atomorbitaltot(i)
            layeratot(i2) = m
            layerbtot(i2) = n

            i2=i2+1           
            cut_num = cut_num+1 
            endif 
            i1=i1+1
          enddo
         enddo
        enddo
        write(*,*)  "total_dim:",i2-1

        h_dim = i2-1
        allocate(testat(3,h_dim))
        write(*,*) 'allocate array testat yes!'
        atompos=0.0
        write(*,*) 'allocate array atompos yes!', h_dim
        ham=0.0
        write(*,*) 'allocate array ham yes! ' , h_dim
        do i =1,h_dim
          atompos(:,i)=positions(:,i)
        enddo

        atomorbital(:)=atomorbitaltot2(1:h_dim)
        allocate(layera(h_dim))
        layera(:)=layeratot(1:h_dim)
        allocate(layerb(h_dim))
        layerb(:)=layerbtot(1:h_dim)
        deallocate( positions)
        write(*,*) 'atomorbital yes! ' , h_dim

         do i=1, h_dim
            atomorbital(i)=mod(atomorbital(i),num_wann)
         enddo
         do i=1,h_dim
            if ( atomorbital(i).eq.0) then
            atomorbital(i)=num_wann
         endif
         enddo
         write(*,*) 'atomorbital mod yes! ' , h_dim
         write(*,*) 'h_dim = ' , h_dim
         write(*,*) 'excludeup=', excludeup
         write(*,*) 'excludedown=', excludedown
         
        spreadmin=0
        spreadmax=0
        do i=1,rvecnum
         if(irvec(1,i).gt.spreadmax)then
            spreadmax=irvec(1,i)
         endif
         if(irvec(1,i).lt.spreadmin)then
            spreadmin=irvec(1,i)
         endif
        enddo

        preadmin=0
        preadmax=0
        do i=1,rvecnum
         if(irvec(2,i).gt.preadmax)then
            preadmax=irvec(2,i)
         endif
         if(irvec(2,i).lt.preadmin)then
            preadmin=irvec(2,i)
         endif
        enddo
        write(*,*) 'spreadmin=',spreadmin,'spreadmax',spreadmax
        write(*,*) 'preadmin=',preadmin,'preadmax',preadmax

        allocate(tmp(spreadmin:spreadmax,num_wann,num_wann))
        allocate(tmpp(spreadmin:spreadmax,preadmin:preadmax,num_wann,num_wann))

        write(*,*) 'read hopping'  
        tmpp=0.0d0
        tmp=0.0d0
        do ii=1,rvecnum
            ix = irvec(1,ii)
            iy = irvec(2,ii)
            iz = irvec(3,ii)
            if (abs(ix).le.ijmax)then
               if (abs(iy).le.ijmax)then
               do i=1,num_wann
                  do j=1,num_wann             
                     tmp(ix,i,j)=tmp(ix,i,j)+hops(i,j,ii)/nrpts(ii)
                     tmpp(ix,iy,i,j)=tmpp(ix,iy,i,j)+ hops(i,j,ii)/nrpts(ii)
                  enddo
               enddo
               endif
            endif
        enddo
         
        write(*,*) 'begin to creat matrix'
        i2=1
        do i =1,h_dim
            do j =1,h_dim
            layerx = layera(j) - layera(i)
            layery = layerb(j) - layerb(i)

          if(layerx.ge.spreadmin.and.layerx.le.spreadmax.and.&
            layery.ge.preadmin.and.layery.le.preadmax)then
           
            i2 =i2+1
            if(atompos(2,i)> -1*atompos(1,i)+t1    &
         .and. atompos(2,i)<0.5*atompos(1,i)+t2    &
         .and. atompos(2,i)>  2*atompos(1,i)+t3    &
            )then
            if(atompos(2,j)> -1*atompos(1,i)+t1    &
         .and. atompos(2,j)<0.5*atompos(1,i)+t2    &
         .and. atompos(2,j)>  2*atompos(1,i)+t3    &
            )then

      !       if(atompos(2,i)<t1                   &
      !  .and. atompos(2,i)<atompos(1,i)+t2        &
      !  .and. atompos(2,i)>t3                     &
      !  .and. atompos(2,i)>atompos(1,i)+t4        &
      !  .and. atompos(1,i)>t5                     &
      !  .and. atompos(1,i)<t6                     &
      !  ) then
   
      !       if(atompos(2,j)<t1                   &
      !  .and. atompos(2,j)<atompos(1,j)+t2        &
      !  .and. atompos(2,j)>t3                     &
      !  .and. atompos(2,j)>atompos(1,j)+t4        &
      !  .and. atompos(1,j)>t5                     &
      !  .and. atompos(1,j)<t6                     &
      !  ) then
              
            hi = atomorbital(i)
            hj = atomorbital(j)

            ham(i,j)=ham(i,j)+tmpp(layerx,layery,hi,hj)
           else
            write(*,*) 'j out bound'
            endif
           else 
            write(*,*) 'i out bound'
            endif
          endif
         enddo
        enddo

        write(*,*) 'all matrix element is ', i2
        write(*,*) 'the Triangle cut is gone'
        
        if(allocated(tmp))         deallocate(tmp) 
        if(allocated(tmpp))        deallocate(tmpp)
        deallocate(layera)
        deallocate(layerb)
       end subroutine

      subroutine solve_ham(NN,A,eigvals,fermi )
        implicit none
        integer,intent(in)  :: NN
        complex,intent(in)  :: A(NN,NN)
        real,intent(in)     :: fermi
        real,intent(inout)  :: eigvals(NN)
        integer             :: info,kp1
        integer             :: lwork
        complex,allocatable :: work(:)
        real,allocatable    :: rwork(:)
        write(*,*) 'slove hdim begin! '

        lwork=16*NN
        allocate(  work(lwork) )
        allocate( rwork(lwork) )
        info = 0
        call zheev('N','U',NN,A,NN,eigvals,work,lwork,rwork,info)
        write(*,*) 'slove matrix yes! '
        if(info.ne.0)then
          write(*,*) 'wrong about zheev'
          stop 'zheevx'
        endif

         do kp1 = 1, NN
            eigvals(kp1) = eigvals(kp1) - fermi      
         enddo
         write(*,*) 'begin to write state! '
         open(111, file='state.dat',recl = 10000)
         do kp1 = 1, NN
            write(111,*) kp1, eigvals(kp1) 
         enddo
         close(111)
         write(*,*) ' write state yes! '
        if(allocated(work))    deallocate(work)
        if(allocated(rwork))   deallocate(rwork)
        write(*,*) 'hamiltonian without eigenvector has been given'
      end subroutine

      subroutine solve_ham_vec(NN,A,eigvals,i,j,m,n,fermi,kmin,kmax )
        implicit none
        integer,intent(in)  :: NN,i,j,m,n
        complex,intent(in)  :: A(NN,NN)
        real,intent(inout)  :: eigvals(NN)
        real,intent(in)     :: fermi
        integer             :: info,kp1,k
        integer             :: lwork
        complex,allocatable :: work(:)
        real,allocatable    :: rwork(:)
        integer,intent(in)  :: kmin,kmax
        character( len = 8 )         :: tempt

         lwork=16*NN
         allocate(  work(lwork) )
         allocate( rwork(lwork) )
         info = 0
         call zheev('V','U',NN,A,NN,eigvals,work,lwork,rwork,info)
     
         if(info.ne.0)then
           write(*,*) 'wrong about zheev'
           stop 'zheevx'
         endif
         do kp1 = 1, NN
            eigvals(kp1) = eigvals(kp1) - fermi
         enddo
         open(111, file='state.dat',recl = 10000)
         do kp1 = 1, NN
            write(111,*) kp1, eigvals(kp1)
         enddo
         close(111)
         
        do k = kmin,kmax 
         write( tempt,'(i8)' ) k
         open(k,file= trim(adjustl(tempt)) //'.dat',recl = 10000)
         do kp1 = 1,NN
             write(k,*) kp1, abs(A(kp1,k))**2
         enddo
         close(k)
        enddo
  
        if(allocated(work))    deallocate(work)
        if(allocated(rwork))   deallocate(rwork)
        write(*,*) 'hamiltonian with eigenvector has been given'
      end subroutine