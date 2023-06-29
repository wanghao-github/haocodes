      program hoti_version2

      implicit none

      complex,allocatable:: hops(:,:,:),ham(:,:),hij(:,:,:)

      integer,allocatable:: irvec(:,:),numbers(:)
      integer,allocatable:: nrpts(:),excludeup(:),excludedown(:)

      real,allocatable   :: eigvals(:),btompositions(:,:)
      real,allocatable   :: wannier_center(:,:),atompos(:,:)
      real,allocatable   :: atompositions(:,:),testat(:,:)
      real,allocatable   :: wannier_center_real(:,:)
      real               :: rdum, idum,fermi
      integer            :: ierr, isize, irank, nwannexup,nwannexdown
      integer            :: nslab1,nslab2,num_lines,num_wann,ix,iy,iz
      integer            :: band1,band2,rvecnum,i,j,num_wann_y,preadmax
      integer            :: preadmin,hdim,up,down
      real               :: amat(3,3)
      integer            :: ijmax,numa,numb,h_dim,htotal
      integer,allocatable   :: atomnumbertot(:),atomorbitaltot(:)
      integer,allocatable   :: atomnumbertot2(:),atomorbitaltot2(:)     
      integer,allocatable   :: atomnumber(:),atomorbital(:)
      integer,allocatable   :: layera(:),layerb(:)
      real,allocatable      :: posreal(:,:)

      real,allocatable      :: positions(:,:)
      real           :: t1,t2,hex_a
      integer        ::   kmin,kmax
      logical       :: l_vec  

   
         open(300,file='inp')
         read(300,*)nwannexup,nwannexdown,fermi
         if(nwannexup .eq.0 ) then
           up= 0
           allocate(excludeup(0))  
           read(300,*)     
         else
           up =2
           allocate(excludeup( 2))       
           read(300,*)excludeup(:)
         endif

         if(nwannexdown .eq. 0) then
           down =0
           allocate(excludedown(0))
           read(300,*)       
         else
           down = 2
           allocate(excludedown(2 ))
           read(300,*)excludedown(:)
         endif

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
         read(300,*) numbers(:)
         read(300,*)amat(1,:)
         read(300,*)amat(2,:)
         read(300,*)amat(3,:)
         read(300,*) ,t1,t2,hex_a
         read(300,*) kmin,kmax,l_vec
         close(300)

         write(*,*) 'da,db,d1,d2,d3:',t1,t2,hex_a
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
   
         write(*,*) 'rvecnum',rvecnum
         allocate(nrpts(rvecnum))
         open(14,file='nrpts_inp')
         do j=1,rvecnum/15
            read(14,'(15I5)') (nrpts(15*(j-1)+i) ,i=1,15)
         enddo
         read(14,'(<mod(rvecnum,15)>I5)') (nrpts(15*(rvecnum/15)+i),
     &                                   i=1,mod(rvecnum,15))
         close(14)

         write(*,*) 'size',size(nrpts)
         allocate( wannier_center(3, num_wann))
         wannier_center  = 0.d0
         open(400,file='position')
         do j=1,num_wann/2
            read(400, *) wannier_center(:,j)
         enddo
         close(400)
         do j=1,num_wann/2
            wannier_center(:,j+num_wann/2)=wannier_center(:,j)
         enddo
         numa=nslab1
         numb=nslab2
         ijmax=12
         allocate( wannier_center_real(3, num_wann))
         wannier_center_real=0.0
         do i =1, num_wann
            wannier_center_real(1, num_wann)=
     &        amat(1,1)*wannier_center(1,num_wann)+
     &        amat(2,1)*wannier_center(2,num_wann)+  
     &        amat(3,1)*wannier_center(3,num_wann)
            wannier_center_real(2, num_wann)=
     &        amat(1,2)*wannier_center(1,num_wann)+
     &        amat(2,2)*wannier_center(2,num_wann)+  
     &        amat(3,2)*wannier_center(3,num_wann)
            wannier_center_real(3, num_wann)=
     &        amat(1,3)*wannier_center(1,num_wann)+
     &        amat(2,3)*wannier_center(2,num_wann)+ 
     &        amat(3,3)*wannier_center(3,num_wann)       
         enddo 
    
          write(*,*) 'num_wann= ', num_wann, 'nslab1 = ',nslab1
          write(*,*) 'nwannexup=',nwannexup,'nwannexdown=',nwannexdown
          write(*,*) 'excludeup=',excludeup,'excludedown',excludedown
          write(*,*)'numa,numb=',numa,numb


          htotal = num_wann*numa*numb - nwannexup-nwannexdown
          hdim= num_wann*numa*numb - nwannexup-nwannexdown
          
          allocate(atomnumber(hdim))
          allocate(hij( -4:4,hdim,hdim))
          write(*,*) 'its time for get dim'

         call tri_dim(nwannexup,nwannexdown,num_wann,numa,rvecnum,
     & excludeup,excludedown,numbers,irvec,hops,nrpts,wannier_center,
     & atomnumber,htotal,up,down, numb,h_dim,t1,t2,hex_a)

          write(*,*) 'its time for a tri_cut program'
          allocate(atompositions(3,h_dim))
          allocate(ham(h_dim,h_dim))
          allocate(atomorbital(h_dim))

          call tri_cut(nwannexup,nwannexdown,num_wann,nslab1,rvecnum,
     & excludeup,excludedown,numbers,irvec,hops,nrpts,wannier_center,
     & atompositions,atomnumber,preadmax,preadmin,htotal,up,down,
     & amat,numb,h_dim,ham,atomorbital,t1,t2,hex_a)

          if(allocated(excludeup))       deallocate(excludeup)
          write(*,*) "@1306"
          if(allocated(excludedown))     deallocate(excludedown)
          write(*,*) "@1307"

      
         open(99, file='orbital.dat',recl = 10000)
         do i = 1, h_dim
            write(99,*) i, atompositions(1,i),atompositions(2,i),
     &                             atompositions(3,i)
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
           posreal(1, i)=
     &        amat(1,1)*atompositions(1,i)+
     &        amat(2,1)*atompositions(2,i)+  
     &        amat(3,1)*atompositions(3,i)
           posreal(2, i)=
     &        amat(1,2)*atompositions(1,i)+
     &        amat(2,2)*atompositions(2,i)+  
     &        amat(3,2)*atompositions(3,i)
           posreal(3, i)=
     &        amat(1,3)*atompositions(1,i)+
     &        amat(2,3)*atompositions(2,i)+ 
     &        amat(3,3)*atompositions(3,i)  
            write(199,*) i, posreal(1,i),posreal(2,i),
     &                             posreal(3,i)
         enddo
         close(199)
         write(*,*) "orbitalreal dat yes!"

         if(allocated(atomnumber))      deallocate(atomnumber)
         write(*,*) "@1311"
         if(allocated(atompositions))   deallocate(atompositions)
         write(*,*) "@1312"
         if(allocated(hij))             deallocate(hij)
         write(*,*) "@1313"
         if(allocated(hops))            deallocate(hops)
         write(*,*) "@1314"
         if(allocated(irvec))           deallocate(irvec)
         write(*,*) "@1315"
         if(allocated(nrpts))           deallocate(nrpts)
         if(allocated(wannier_center))  deallocate(wannier_center)
         write(*,*) "@1318"
         if(allocated(btompositions))   deallocate(btompositions)
         if(allocated(posreal))   deallocate(posreal)
 
         if(allocated(atomnumbertot))   deallocate(atomnumbertot)
         if(allocated(atomnumbertot2))   deallocate(atomnumbertot2)
         if(allocated(atomorbitaltot))   deallocate(atomorbitaltot)
         if(allocated(atomorbitaltot2))   deallocate(atomorbitaltot2)
         if(allocated(layera))   deallocate(layera)
         if(allocated(layerb))   deallocate(layerb)
         if(allocated(positions))   deallocate(positions)
 
 
         write(*,*) "@1323"
         
        allocate(eigvals(h_dim))
        write(*,*) 'fermi level is set as ', fermi
        if(l_vec) then
          call solve_ham_vec(h_dim,ham,eigvals,3810,3811,3812,3813,
     &       fermi,kmin,kmax)
        else
          call solve_ham(h_dim,ham,eigvals,fermi)
        endif




        write(*,*) "@1316"
        if(allocated(ham))             deallocate(ham)
        write(*,*) "@1317"

        if(allocated(eigvals))         deallocate(eigvals)
        write(*,*) "@1319"
        if(allocated(excludeup))       deallocate(excludeup)
        write(*,*) "@1320"
        if(allocated(excludedown))     deallocate(excludedown)
        write(*,*) "@1321"
        if(allocated(numbers))         deallocate( numbers)
        write(*,*) "@1322"



      end program hoti_version2

      subroutine tri_cut(nwannexup,nwannexdown,num_wann,numa,rvecnum,
     & excludeup,excludedown,numbers,irvec,hops,nrpts,wanniercenters,
     & atompos,atomnumber,preadmax,preadmin,htotal,up,down,
     & amat,numb,h_dim,ham,atomorbital,t1,t2,hex_a)
        implicit none
        integer,intent(in)  :: nwannexup,nwannexdown,num_wann,numa,numb
        integer,intent(in)    :: rvecnum,up,down
        integer,intent(in)    :: excludeup(up),excludedown(down)
        integer,intent(in)    :: numbers(num_wann),irvec(3,rvecnum) 
        integer,intent(in)    :: nrpts(rvecnum)
        integer           :: htotal

        integer          :: h_dim,i,i1,ii,ix,iy,iz,zvalue,hi,hj,j,i2
        integer               :: spreadmin,spreadmax,preadmax,preadmin
        integer               :: ijmax,m,n
        integer       :: atomnumbertot(htotal),atomorbitaltot(htotal)
        integer       :: atomnumbertot2(htotal),atomorbitaltot2(htotal)
        integer       :: layerdir(htotal),cut_num,layeratot(htotal)
        integer               :: layerbtot(htotal),atomnumber(htotal)
        integer   :: atomorbital(h_dim)
        integer,allocatable   :: layera(:),layerb(:)
        real,intent(in)       :: amat(3,3)
        integer   :: layerx,layery
        

        real                  :: ra,tol
        real,intent(in)       :: wanniercenters(3,num_wann)

        complex,intent(in)    :: hops(num_wann,num_wann,rvecnum)
        complex,allocatable   :: tmp(:,:,:),tmpp(:,:,:,:)
        real,allocatable      :: testat(:,:)
        real                  :: atompos(3,h_dim)
        real,allocatable      :: positions(:,:)
        real       ::positionmed(3,htotal)
        real    :: t1,t2,hex_a

        complex         :: ham(h_dim,h_dim)

        tol=1.0e-8
        cut_num=0
        atomnumbertot2=0
        atomorbitaltot2=0
        layeratot=0
        layerbtot=0
        ijmax=12
        write(*,*) ' in tri_cut htotal = ' , htotal
        allocate(positions(3,htotal))
!        allocate(positions2(300,htotal*100))
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
            positionmed(1,i)=wanniercenters(1,i1)+m-1
            positionmed(2,i)=wanniercenters(2,i1)+n-1
            positionmed(3,i)=wanniercenters(3,i1)

            layerdir(i) = (i+nwannexup-1)/num_wann+1

            if(positionmed(2,i)<t2+hex_a
     & .and. positionmed(2,i)<positionmed(1,i)+t2-t1+hex_a
     & .and. positionmed(2,i)>t2-hex_a 
     & .and. positionmed(2,i)>positionmed(1,i)+t2-t1-hex_a
     & .and. positionmed(1,i)>t1-hex_a
     & .and. positionmed(1,i)<t1+hex_a
     & ) then  

            
            positions(1,i2) = positionmed(1,i)
            positions(2,i2) = positionmed(2,i)
            positions(3,i2) = positionmed(3,i)

            atomnumbertot2(i2)=atomnumbertot(i)
            atomorbitaltot2(i2)=atomorbitaltot(i)
            layeratot(i2)=m
            layerbtot(i2)=n

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
!        allocate(atompos(3,htotal))
        atompos=0.0
        write(*,*) 'allocate array atompos yes!', h_dim
!        allocate(ham(h_dim,h_dim))
        ham=0.0
        write(*,*) 'allocate array ham yes! ' , h_dim

        do i =1,h_dim
          atompos(:,i)=positions(:,i)
        enddo

!        allocate(atomorbital(h_dim))
        atomorbital(:)=atomorbitaltot2(1:h_dim)
        allocate(layera(h_dim))
        layera(:)=layeratot(1:h_dim)
        allocate(layerb(h_dim))
        layerb(:)=layerbtot(1:h_dim)
        deallocate( positions)
        write(*,*) 'atomorbital yes! ' , h_dim

        do i=1, h_dim
!         ! atomorbitaltot(i)=atomorbitaltot(1,i)+50-int(
!      &     atompos(1,i) + 50)

          atomorbital(i)=mod(atomorbital(i),num_wann)
!          write(*,*) 'atomorbital mod i yes! ' , i
   
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
        allocate(tmpp(spreadmin:spreadmax,preadmin:preadmax,
     &                      num_wann,num_wann))

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
          layerx=   layera(j)-layera(i)
          layery=   layerb(j)-layerb(i)
          zvalue =layerdir(j)-layerdir(i)

          if(layerx.ge.spreadmin.and.layerx.le.spreadmax.and.
     &     layery.ge.preadmin.and.layery.le.preadmax)then

       
            i2 =i2+1
            if(atompos(2,i)<t2+hex_a 
     & .and. atompos(2,i)<atompos(1,i)+t2-t1+hex_a
     & .and. atompos(2,i)>t2-hex_a 
     & .and. atompos(2,i)>atompos(1,i)+t2-t1-hex_a
     & .and. atompos(1,i)>t1-hex_a
     & .and. atompos(1,i)< t1+hex_a
     & ) then
   
           if(atompos(2,j)<t2+hex_a 
     & .and. atompos(2,j)<atompos(1,j)+t2-t1+hex_a
     & .and. atompos(2,j)>t2-hex_a 
     & .and. atompos(2,j)>atompos(1,j)+t2-t1-hex_a
     & .and. atompos(1,j)>t1-hex_a
     & .and. atompos(1,j)< t1+hex_a
     & ) then
                
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

        !!!###

        write(*,*) 'all matrix element is ', i2

        write(*,*) 'the Triangle cut is gone'
        
       
        if(allocated(tmp))         deallocate(tmp) 
        write(*,*) "@1303"
        if(allocated(tmpp))        deallocate(tmpp) 
        write(*,*) "@1304"
        deallocate(layera)
        deallocate(layerb)
       end subroutine   

       
       subroutine tri_dim(nwannexup,nwannexdown,num_wann,numa,rvecnum,
     & excludeup,excludedown,numbers,irvec,hops,nrpts,wanniercenters,
     & atomnumber,htotal,up,down, numb,h_dim,t1,t2,hex_a)
        implicit none
        integer,intent(in)  :: nwannexup,nwannexdown,num_wann,numa,numb
        integer,intent(in)    :: rvecnum,up,down
        integer,intent(in)    :: excludeup(up),excludedown(down)
        integer,intent(in)    :: numbers(num_wann),irvec(3,rvecnum) 
        integer,intent(in)    :: nrpts(rvecnum)
        integer           :: htotal

        integer          :: h_dim,i,i1,ii,ix,iy,iz,zvalue,hi,hj,j,i2
        integer               :: ijmax,m,n
        integer       :: atomnumbertot(htotal),atomorbitaltot(htotal)
        integer       :: atomnumbertot2(htotal),atomorbitaltot2(htotal)
        integer       :: layerdir(htotal),cut_num,layeratot(htotal)
        integer               :: layerbtot(htotal),atomnumber(htotal)

        real                  :: tol
        real,intent(in)       :: wanniercenters(3,num_wann)

        complex,intent(in)    :: hops(num_wann,num_wann,rvecnum)
        real,allocatable      :: positions(:,:)
        real       ::positionmed(3,htotal)
        real      :: t1,t2,hex_a


!        htotal = num_wann*numa*numb - nwannexup-nwannexdown
        tol=1.0e-8
        cut_num=0
        atomnumbertot2=0
        atomorbitaltot2=0
        layeratot=0
        layerbtot=0
        ijmax=12
        write(*,*) ' in tri_cut htotal = ' , htotal
        allocate(positions(3,htotal))
!        allocate(positions2(300,htotal*100))
        write(*,*) 'allocate array positions yes! ' 
        positions=0.0
        positionmed=0.0
         i=1
         i2=1

        do m=1,numa
         do n=1,numb
            i1=1

          do i=1,num_wann
            atomnumbertot(i) = numbers(i1)
            atomorbitaltot(i) = i1
            positionmed(1,i)=wanniercenters(1,i1)+m-1
            positionmed(2,i)=wanniercenters(2,i1)+n-1
            positionmed(3,i)=wanniercenters(3,i1)

            layerdir(i) = (i+nwannexup-1)/num_wann+1

            if(positionmed(2,i)<t2+hex_a  
     & .and. positionmed(2,i)<positionmed(1,i)+t2-t1+hex_a
     & .and. positionmed(2,i)>t2-hex_a 
     & .and. positionmed(2,i)>positionmed(1,i)+t2-t1-hex_a 
     & .and. positionmed(1,i)>t1-hex_a 
     & .and. positionmed(1,i)<t1+hex_a
     & ) then
     

            i2=i2+1           
            endif 
            i1=i1+1
          enddo
         enddo
        enddo
        h_dim=i2-1
        write(*,*)  "total_dim:",i2-1

        deallocate(positions)

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
        call zheev('N','U',NN,A,NN,
     &       eigvals,  work,lwork,rwork,info)
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
        if(allocated(work))   deallocate(work)
        write(*,*) "@1"
        if(allocated(rwork))   deallocate(rwork)
        write(*,*) "@2"  
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
         call zheev('V','U',NN,A,NN,
     &       eigvals,  work,lwork,rwork,info)
     
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
         

!write( cTemp,'(i2)' ) k
!open ( 1, file = 'myFile' // trim(adjustl( cTemp )) // '.dat', status =
!'old' )

        do k = kmin,kmax 
         write( tempt,'(i8)' ) k
         open(k,file= trim(adjustl(tempt)) //'.dat',recl = 10000)
         do kp1 = 1,NN
             write(k,*) kp1, abs(A(kp1,k))**2
         enddo
         close(k)
        enddo
  

         
        if(allocated(work))   deallocate(work)
        write(*,*) "@1"
        if(allocated(rwork))   deallocate(rwork)
        write(*,*) "@2"  
        write(*,*) 'hamiltonian with eigenvector has been given'
      end subroutine
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      subroutine to_red_sc(vec_unit,vec_super )
       implicit none
       real,intent(in)     :: vec_unit(2)
       real,intent(out)    :: vec_super(2)
       real                :: trans(2,2)
       real                :: k1,k2
       integer             :: ipiv(2)
       integer             :: info

      
       trans(1,1) =  1
       trans(1,2) =  1
       trans(2,1) = -1
       trans(2,2) =  1
       k1 = vec_unit(1)
       k2 = vec_unit(2) 
       vec_super(1) = k1  
       vec_super(2) = k2
       ipiv = 0
       call dgesv(2,2,transpose(trans),2,ipiv,vec_super,2,info)
    

!       write(*,*) 'new vec_unit is', vec_unit
! it's testified and same as old can't be wrong
      end subroutine



      subroutine make_supercell(rvecnum,num_wann,hoppings,use_sc_red_lat
     & ,self_lat,self_orb, nrpts,irvec,hij, sc_irvec,sc_orbital,num_hop)
        implicit none
        integer               :: self_dim_r,max_R,i,j,t,m,n,hi,hj
        integer               :: num_sc,pair_ind,ijmax,num_hop
        integer,intent(in)    :: num_wann,rvecnum
        integer,intent(in)    :: nrpts(rvecnum),irvec(3,rvecnum)
        integer               :: isc_part(2)
        real                  :: vec(2),tmp_red(2),save_vec(2)
        real,intent(in)       :: use_sc_red_lat(2,2)
        real,intent(in)       :: self_lat(2)
        real                  :: sc_cart_lat(2)
        integer               :: t1,t2
        integer,allocatable   :: sc_cands(:,:)
        real,allocatable      :: sc_vec(:,:)

        real                  :: det,eps_shift,R(2)
        real,intent(in)       :: self_orb(2,num_wann)
        real                  :: sc_part(2),orig_part(2),orb1(2),orb2(2)
        real,allocatable      :: sc_orb(:,:),final_orb(:,:)

        complex,intent(in)    :: hoppings(num_wann,num_wann,rvecnum)
        complex               :: amp
        complex,intent(inout) :: hij(num_hop)
        integer,intent(inout) :: sc_irvec(2,num_hop)
        integer,intent(inout) :: sc_orbital(2,num_hop)
        logical               :: inside

        ijmax = 10
        self_dim_r=2
        write(*,*) use_sc_red_lat
        det = use_sc_red_lat(1,1)*use_sc_red_lat(2,2)-
     &              use_sc_red_lat(1,2)*use_sc_red_lat(2,1)

        if (det .le. 0.0)then
            write(*,*) ' need to form right handed
     &   system.',use_sc_red_lat
            stop
        endif

! For sure , all the high dimensional matrix the firtst is lo
        max_R= 1 *self_dim_r
        allocate(sc_cands(2,(2*max_R+1)*(2*max_R+1)))
        t = 0
        do  i = -max_R,max_R
         do j = -max_R,max_R
            t = t+1
            sc_cands(1,t)= i  
            sc_cands(2,t)= j
         enddo
        enddo
! shift of the grid cant understand unles I make such mistake
        eps_shift= sqrt(2.0)*1.0E-8 
        tmp_red = 0d0      
        num_sc = 0
        do i=1,t
           
           vec(1)=float( sc_cands(1,i))
           vec(2)=float( sc_cands(2,i))
           save_vec = vec
           call to_red_sc(vec,tmp_red,2)
           inside = .True.
           do j =1,2 
            if(tmp_red(j) .le. -1.0*eps_shift .or. tmp_red(j) 
     &            .ge. 1.0-eps_shift) then
               inside = .False.   
            endif  
           enddo           
           if( inside)then
              num_sc = num_sc+1
           endif
        enddo

        write(*,*) 'num_Sc=,', num_sc              
        allocate(sc_vec( 2,num_sc))
        num_sc = 0
        tmp_red = 0d0
        sc_vec  = 0d0
        do i=1,t
           vec(1)=float( sc_cands(1,i))
           vec(2)=float( sc_cands(2,i))
           save_vec = vec
           call to_red_sc(vec,tmp_red)
           inside = .True.
           do j =1,2
            if (tmp_red(j) .le. -1.0*eps_shift .or. tmp_red(j)
     &            .ge. 1.0-eps_shift) then
                  inside = .False.
            endif
           enddo
           if( inside)then
              num_sc = num_sc+1
              sc_vec(1,num_sc) = save_vec(1)
              sc_vec(2,num_sc) = save_vec(2)
           endif
        enddo

        write(*,*) 'num_sc= ',num_sc,' det=',det
        write(*,*) 'value for the supercell is'
        write(*,*) sc_vec

        if ( NINT ( abs( det ) ) .ne. num_sc )then
            write(*,*)'Super-cell generation failed! Wrong
     &            number of super-cell vectors found.'
        endif
        sc_cart_lat=  matmul(use_sc_red_lat,self_lat)

        allocate(sc_orb(2, num_wann*2  ))
        allocate(final_orb(2,num_wann*2))
        final_orb = 0.d0
        sc_orb = 0.d0
        t = 0
        do i=1,num_sc 
           do j =1,num_wann
              t = t+1
              R(1)= self_orb(1,j) + sc_vec(1,i)
              R(2)= self_orb(2,j) + sc_vec(2,i)
              call to_red_sc( R, tmp_red )
              sc_orb(1,t) = tmp_red(1)
              sc_orb(2,t) = tmp_red(2)
           enddo
        enddo
        
        t = 0
        hij = 0.d0
        sc_irvec = 0
        do i =1,num_sc        
          do j =1,rvecnum
            do m =1,num_wann
             do n =1,num_wann
                t = t+1
                t1 = n + (i-1)*num_wann
                t2 = m + (i-1)*num_wann
                orb1=sc_orb(:,t1)
                orb2=sc_orb(:,t2)
                   
                sc_part = 0d0
                amp= hoppings(n,m,j) /float(nrpts(j))
                R(1) = float(irvec(1,j)) + sc_vec(1,i)
                R(2) = float(irvec(2,j)) + sc_vec(2,i)
!    round down
                call to_red_sc( R , sc_part )
                isc_part= int( floor(sc_part))
!    R should be changed so we give it value twice
                R(1) = float(irvec(1,j)) + sc_vec(1,i)
                R(2) = float(irvec(2,j)) + sc_vec(2,i)
                orig_part(1)=R(1) - (isc_part(1) - isc_part(2) )
                orig_part(2)=R(2) - (isc_part(1) + isc_part(2) )
                pair_ind= -1
                if( abs( orig_part(1) - sc_vec(1,1)) < 1e-5  )then
                  if( abs( orig_part(2) - sc_vec(1,2)) < 1e-5  )then
                    pair_ind = 0
                  endif
                endif
!                write(*,*) 'orig_part=',orig_part
                if( abs( orig_part(1) - sc_vec(2,1)) < 1e-5  )then
                  if( abs( orig_part(2) - sc_vec(2,2)) < 1e-5  )then
                     pair_ind = 1
                  endif
                endif

                if( pair_ind .eq. -1)then
                   write(*,*) 'some thing wrong about super cell vector'
                   stop
                endif

                hi= n + (  i-1 )*num_wann
                hj= m + pair_ind*num_wann

                call shift_to_home(isc_part,orb1,orb2)

                if( abs( orb1(1) ) .gt. 1e-5  )then
                  if( abs( orb1(2) ) .gt. 1e-5  )then
                     final_orb(:,t1) = orb1
                  endif
                endif
                sc_orbital(1,t) = hi
                sc_orbital(2,t) = hj
                
                sc_irvec(1,t)=isc_part(1)
                sc_irvec(2,t)=isc_part(2)
                
                hij(t) = hij(t) + amp

              
              enddo
            enddo
           enddo
         enddo
         if(allocated(sc_cands))   deallocate(sc_cands)
             write(*,*) "@1391"
         if(allocated(sc_vec))   deallocate(sc_vec)
             write(*,*) "@1392"
         if(allocated(sc_orb))   deallocate(sc_orb)
             write(*,*) "@1393"
         if(allocated(final_orb))   deallocate(final_orb)
             write(*,*) "@1394"

!         deallocate(sc_cands,sc_vec,sc_orb,final_orb)
         

       
      end subroutine

      subroutine  shift_to_home(R , O1, O2)
        implicit none
        integer,intent(inout) :: R(2)
        real,intent(inout)    :: O1(2),O2(2)
        real                  :: round_orb1(2),round_orb2(2)
        integer               :: dis1(2),dis2(2)

!Shifts all orbital positions to the home unit cell.

        round_orb1(1)=mod( (  O1(1) +1.0E-6) ,1.0)
        round_orb1(2)=mod( (  O1(2) +1.0E-6) ,1.0)

        round_orb2(1)=mod( (  O2(1) +1.0E-6) ,1.0)
        round_orb2(2)=mod( (  O2(2) +1.0E-6) ,1.0)

!find displacement vector needed to bring back to home cell
        dis1 = nint(  O1 - round_orb1 )
        dis2 = nint(  O2 - round_orb2 )

        O1   = O1 - float(dis1)
        O2   = O2 - float(dis2)
        R    = R - dis1 + dis2

      end subroutine

