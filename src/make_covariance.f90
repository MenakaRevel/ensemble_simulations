program make_covariance
    implicit none
    ! designed to make covariance matrix for spatially correlated normal fields 
    ! considering the spatial distance 
    ! r=exp(d/l)
    !$ use omp_lib
    integer                         :: i,j,im,jm,ios,n,ut,l
    integer                         :: ix,iy,iix,iiy
    !$ integer                      :: omp_k
    integer,parameter               :: latpx=720,lonpx=1440
    real,dimension(lonpx,latpx)     :: nextdst,lon,lat
    integer,dimension(lonpx,latpx)  :: nextX,nextY,rivnum 
    real,dimension(lonpx,latpx)     :: WSE
    !real,dimension(365*lonpx*latpx) :: WSE_tmp
    character*256                   :: fname,buf,camadir,rivnumdir,rfile
    character*4                     :: yyyy
    real                            :: c,gaussian_cov,covariance_lag,lag_dist,hubeny_real,sigma,dist
    real                            :: lat1, lon1, lat2, lon2
    integer                         :: ens_num,num
    character*3                     :: numch
    real,dimension(:),allocatable   :: error
    integer,allocatable             :: xlist(:),ylist(:)
    real,allocatable                :: idlist(:),uparea(:),subarea(:)
    real,allocatable                :: cov(:,:)
    
    !$ omp_k=omp_get_max_threads()
    !$ write(*,*)omp_k
    !$ call  omp_set_num_threads(omp_k)
    
    ! ! rivnum 
    ! call getarg(1,buf)
    ! read(buf,*) num
    
    ! pixels
    call getarg(1,buf)
    read(buf,*) l 
    
    ! standard deviation
    call getarg(2,buf)
    read(buf,*) sigma 

    call getarg(3,buf)
    read(buf,*) dist
    
    call getarg(4,buf)
    read(buf,"(A)") camadir

    call getarg(5,buf)
    read(buf,"(A)") rfile
    
    ! read next grid information
    ! read nextX and nextY
    fname=trim(adjustl(camadir))//"map/glb_15min/nextxy.bin"
    open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
    if(ios==0)then
        read(34,rec=1) nextX
        read(34,rec=2) nextY
    else
        write(*,*) "no file nextXY at:",fname    
    end if
    close(34)
    
    ! read distance to next grid
    fname=trim(adjustl(camadir))//"map/glb_15min/nxtdst.bin"
    open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
    if(ios==0)then
        read(34,rec=1) nextdst
    else
        write(*,*) "no file nextdst",fname
    end if
    close(34)
    
    ! read lat lon 
    fname=trim(adjustl(camadir))//"map/glb_15min/lonlat.bin"
    open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
    if(ios==0)then
        read(34,rec=1) lon
        read(34,rec=2) lat
    else
        write(*,*) "no file lonlat",fname
    end if
    close(34)
    
    ! read rivnum.bin 
    !hubeny_realfname=trim(adjustl(rivnumdir))//"rivnum.bin"
    !open(34,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="old",iostat=ios)
    !if(ios==0)then  
    !    read(34,rec=1) rivnum
    !else
    !    write(*,*) "no file at rivnum", fname
    !nd if
    !close(34)
    
    ! read river array
    allocate(idlist(l),xlist(l),ylist(l),uparea(l),subarea(l),cov(l,l))
    ! fname="outsub.txt"
    fname=trim(rfile)
    open(34,file=fname,status='old',access='sequential',form='formatted',action='read')!
    do i=1,l 
      read(34,*) idlist(i),xlist(i),ylist(i),uparea(i),subarea(i)
      !write(*,*) xlist(i),ylist(i)
    end do
    !write(*,23) xlist,ylist,wgt
    close(34)
    !
    cov=0.0d0
    !-- 
    do i=1, l
      ix=xlist(i)
      iy=ylist(i)
      lat1=lat(ix,iy)
      lon1=lon(ix,iy) 
      do j=i, l
        iix=xlist(j)
        iiy=ylist(j)
        !write(*,*) ix,iy,iix,iiy
        !call lag_distance(ix,iy,iix,iiy,lonpx,latpx,nextX,nextY,nextdst,lag_dist)
        !call lag_distance_com(ix,iy,iix,iiy,lonpx,latpx,nextX,nextY,nextdst,lag_dist)
        lat2=lat(iix,iiy)
        lon2=lon(iix,iiy) 
        lag_dist=hubeny_real(lat1, lon1, lat2, lon2)
        if (lag_dist==-9999.0d0) then
          c=0.0d0
        elseif (lag_dist==0.0d0) then
          c=1.0d0
        else
        !   c=gaussian_cov(lag_dist,sigma)
          c=covariance_lag(lag,sigma,dist)
        end if 
        !write(*,*) i,j,lag_dist,c
        cov(i,j)=c
        cov(j,i)=c
      end do
    end do
    ! save cov
    write(*,*) "dimensions of covariance:",l,"x",l
    fname="cov.bin"
    open(69,file=fname,form="unformatted",access="direct",recl=4*l*l,status="replace",iostat=ios)
      if(ios==0)then
         write(69,rec=1) cov 
      end if
    close(69)
    deallocate(idlist,xlist,ylist,uparea,subarea,cov)
    end program make_covariance
    !*****************************************************************
    subroutine lag_distance(i,j,x,y,nx,ny,nextX,nextY,nextdst,lag_dist)
    implicit none 
    !--
    integer                             :: i,j,x,y,nx,ny
    integer,dimension(nx,ny)            :: nextX,nextY
    real,dimension(nx,ny)               :: nextdst
    !--
    real                                :: lag_dist
    integer                             :: ix,iy,iix,iiy,tx,ty,ud
    real                                :: length,rl
    !--
    if (i==x .and. j==y) then
      ud=0
    else
      ud=-1
    end if
    !--
    if (ud==-1) then
      tx=x
      ty=y
      ix=i
      iy=j
      length=0.0
      lag_dist=0.0
      do while (ix/=tx .or. iy/=ty) 
        iix=ix
        iiy=iy 
        ix=nextX(iix,iiy)
        iy=nextY(iix,iiy)
        if (ix==-9 .or. iy==-9) then
          ud=+1
          exit
        end if
        if (ix==-9999 .or. iy==-9999) then
          ud=+1
          exit
        end if
        !-- half of the present grid
        rl=anint((nextdst(ix,iy)/1000.0)*100)/100.0
        length=length+rl!/2.0
      end do
    end if
    !--
    if (ud==+1) then
      tx=i
      ty=j
      ix=x
      iy=y
      length=0.0
      do while (ix/=tx .or. iy/=ty) 
        iix=ix
        iiy=iy
        ix=nextX(iix,iiy)
        iy=nextY(iix,iiy)
        !---
        if (ix==-9 .or. iy==-9) then
          ud=-9999
          exit
        end if
        if (ix==-9999 .or. iy==-9999) then
          ud=-9999
          exit
        end if
        !-- half of the present grid
        rl=anint((nextdst(ix,iy)/1000.0)*100)/100.0
        length=length+rl!/2.0
      end do
    end if
    !-- 
    if (ud==-9999) then
      lag_dist=-9999.0d0
    elseif (ud==0) then
      lag_dist=0.0d0
    else
      lag_dist=length
    end if
    !---
    return
    !---
    end subroutine lag_distance
    !*****************************************************************
    subroutine lag_distance_com(i,j,x,y,nx,ny,nextX,nextY,nextdst,lag_dist)
    implicit none 
    !--
    integer                             :: i,j,x,y,nx,ny
    integer,dimension(nx,ny)            :: nextX,nextY
    real,dimension(nx,ny)               :: nextdst
    !--
    real                                :: lag_dist
    integer                             :: ix,iy,iix,iiy,tx,ty,ud,ii,jj,iii,jjj
    real                                :: length,rl,r1,r2
    !--
    if (i==x .and. j==y) then
      ud=0
    else
      ud=-1
    end if
    !--
    if (ud==-1) then
      tx=x
      ty=y
      ix=i
      iy=j
      length=0.0!anint((nextdst(ix,iy)/1000.0)*100)/100.0
      lag_dist=0.0
      do while (ix/=tx .or. iy/=ty) 
        !rl=anint((nextdst(ix,iy)/1000.0)*100)/100.0
        iix=ix
        iiy=iy 
        ix=nextX(iix,iiy)
        iy=nextY(iix,iiy)
        if (ix==-9 .or. iy==-9) then
          ud=+1
          exit
        end if
        if (ix==-9999 .or. iy==-9999) then
          ud=+1
          exit
        end if
        !-- the present grid nxtdst iix, iiy
        rl=anint((nextdst(iix,iiy)/1000.0)*100)/100.0
        length=length+rl!/2.0
      end do
    end if
    !--
    if (ud==+1) then
      tx=i
      ty=j
      ix=x
      iy=y
      length=0.0
      do while (ix/=tx .or. iy/=ty) 
        !rl=anint((nextdst(ix,iy)/1000.0)*100)/100.0
        iix=ix
        iiy=iy
        ix=nextX(iix,iiy)
        iy=nextY(iix,iiy)
        !---
        if (ix==-9 .or. iy==-9) then
          ud=-9999
          exit
        end if
        if (ix==-9999 .or. iy==-9999) then
          ud=-9999
          exit
        end if
        !-- the present grid nxtdst iix, iiy
        rl=anint((nextdst(iix,iiy)/1000.0)*100)/100.0
        length=length+rl!/2.0
      end do
    end if
    !-- 
    if (ud==-9999) then
      r2=anint((nextdst(i,j)/1000.0)*100)/100.0
      r1=0.0 !anint((nextdst(x,y)/1000.0)*100)/100.0
      ix=x !nextX(x,y)
      iy=y !nextY(x,y)
      ii=nextX(i,j)
      jj=nextY(i,j)
      do while (ix==-9 .or. iy==-9) 
        iix=ix
        iiy=iy
        ix=nextX(iix,iiy)
        iy=nextY(iix,iiy)
        r1=r1+anint((nextdst(iix,iiy)/1000.0)*100)/100.0
        do while ((ii==ix .and. jj==iy) .or. (ii==-9 .or. jj==-9))  
          iii=ii
          jjj=jj
          ii=nextX(iii,jjj)
          jj=nextY(iii,jjj)
          r2=r2++anint((nextdst(ii,jj)/1000.0)*100)/100.0
        end do
      end do
      lag_dist=r1+r2
    elseif (ud==0) then
      lag_dist=0.0d0
    else
      lag_dist=length
    end if
    !---
    return
    !---
    end subroutine lag_distance_com
    !**************************************************
    function gaussian_cov(lag,sigma)
    ! calculate the c=sigma**2*exp(-3*lag/distance)
    implicit none
    real                         :: lag,gaussian_cov,sigma!,L
    !real,parameter               :: sigma=2.5 ! m 
    real,parameter               :: L=100.0d0 ! km
    
    gaussian_cov=(sigma**2)*(exp(-3.0d0*lag/L))
    return 
    end function gaussian_cov
    !*************************************************
    function covariance_lag(lag,sigma,L)
    ! calculate the c=sigma**2*exp(-3*lag/distance)
    implicit none
    real                         :: lag,covariance_lag,sigma,L
    !real,parameter               :: sigma=2.5 ! m 
    ! real,parameter               :: L=100.0d0 ! km
    
    gaussian_cov=(sigma**2)*(exp(-1.0*lag/L))
    return 
    end function gaussian_cov
    !*************************************************
    FUNCTION hubeny_real(lat1, lon1, lat2, lon2)
      implicit none
      !-- for input -----------
      real                                  lat1, lon1, lat2, lon2
      !-- for output-----------
      real                                  hubeny_real  ! (m)
      !-- for calc ------------
      real,parameter                     :: pi = atan(1.0)*4.0
      real,parameter                     :: a  = 6378137
      real,parameter                     :: b  = 6356752.314140
      real,parameter                     :: e2 = 0.00669438002301188
      real,parameter                     :: a_1_e2 = 6335439.32708317
      real                                  M, N, W
      real                                  latrad1, latrad2, lonrad1, lonrad2
      real                                  latave, dlat, dlon
      real                                  dlondeg
      !------------------------
      latrad1   = lat1 * pi / 180.0
      latrad2   = lat2 * pi / 180.0
      lonrad1   = lon1 * pi / 180.0
      lonrad2   = lon2 * pi / 180.0
      !
      latave    = (latrad1 + latrad2)/2.0
      dlat      = latrad2 - latrad1
      dlon      = lonrad2 - lonrad1
      !
      dlondeg   = lon2 - lon1
      if ( abs(dlondeg) .gt. 180.0) then
        dlondeg = 180.0 - mod(abs(dlondeg), 180.0)
        dlon    = dlondeg * pi / 180.0
      end if
      !-------
      W  = sqrt(1.0 - e2 * sin(latave)**2.0 )
      M  =  a_1_e2 / (W**3.0)
      N  =  a / W
      hubeny_real  = sqrt( (dlat * M)**2.0 + (dlon * N * cos(latave))**2.0 )
    RETURN
    END FUNCTION hubeny_real
    !***************************************************