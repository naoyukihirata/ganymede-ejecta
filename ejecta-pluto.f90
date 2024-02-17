

     implicit none
     integer :: i, k,k2,j
     integer, parameter :: num = 10000, num2 = 3600
     real(8) :: c_lon, c_lat,pi,rtem1,rtem3,rtem2
     integer :: loop,loop2,knum
     real(8) :: pp(3),pv(3),dx(3),dp(3),dx_prev(3),dp_prev(3),dt,g,lunch_point
     real(8) :: gm,minbox,distance1,distance2,w, radius,lunch,rotation,c_radius
     real(8) :: c_nz,c_nx,c_ny,Q,c_Q,c_cQ,c_sQ,im_lat,im_lon
     real(8) :: mcx,mcy,mcz,cy,cx,cz,cQ,sQ,nz,ny,nx,polar_lon,polar_lat
     real(8) :: density(-90:90, 0:359) = 0.d0
     real(8) :: old_lunch_point, vol,volaspect
     real(8) :: record(num*num2,6),record2(num*num2,2),clat,clon
     integer*1 :: ini(54), col(1000,-90:90,3),cole
     real(8) :: maxdensity,mindensity,area
     integer :: Target, bina
     real(8) :: colat, coarea,exit_radius 
     real(8) :: c1, h, u, n2, kk, ro, yy, p
     real(8) :: ob_u,ob_c1,ob_kk,ob_c_radius,zeta,impact_angle
     real(8) :: u1,u2,aa,a1,a2
     real(8) :: CENTERLATITUDE, CENTERLONGITUDE,AANGLE
     pi = acos(-1.0d0)

   !!!!!!!!parameter setting
   volaspect= 0. !azimuth of uprange, anti-crockwise from the east
                          ! = 0 means the impactor comes from the exact east 
                          ! = 90 means the impactor comes from the exact north
   Target = 4      ! 1 ~ 8, Housen and Holsapple scaling constants
   impact_angle=70.  != 90 means the vertical impact, = 0 means the horizontal
   c_radius = 400.*1000.   !transient crater radius in meter
   CENTERLATITUDE = 20.    !crater center latitude
   CENTERLONGITUDE = 90.   !crater center east longitude
   AANGLE = 45.            ! ejecta launch angle
   !!!!!!!!parameter setting end


   if (Target .eq. 1) then
     u = 0.55 ; kk = 0.2 ; c1 = 1.5 ; h  = 0.68 ; n2 = 1.5 ; p  = 0.5 ; ro = 1000. ; endif
   if (Target .eq. 2) then 
     u = 0.55 ; kk = 0.3 ; c1 = 1.5 ; h  = 1.1 ; n2 = 1 ; p  = 0.5 ; ro = 3000. ; yy = 30 ; endif
   if (Target .eq. 3) then
     u = 0.46 ; kk = 0.3 ; c1 = 0.18 ; h  = 0.38 ; n2 = 1 ; p  = 0.3 ; ro = 2600. ; yy = 0.45 ; endif
   if ((Target .eq. 4) .or. (Target .eq. 5)) then
     u = 0.41 ; kk = 0.3 ; c1 = 0.55 ; h  = 0.59 ; n2 = 1.3 ; p  = 0.3 ; ro = 1000. ; endif
   if (Target .eq. 6) then
     u = 0.45 ; kk = 0.5 ; c1 = 1. ; h  = 0.8 ; n2 = 1.3 ; p  = 0.3 ; ro = 1500. ; endif
   if (Target .eq. 7) then
     u = 0.4 ; kk = 0.3 ; c1 = 0.55 ; h  = 0.4 ; n2 = 1. ; p  = 0.3 ; ro = 1500. ; yy = 4.e-3 ; endif
   if (Target .eq. 8) then
     u = 0.35 ; kk = 0.32 ; c1 = 0.6 ; h  = 0.81 ; n2 = 1. ; p  = 0.2 ; ro = 1200. ; yy = 2.e-3 ; endif

 open(35,file= 'T1E_vol.txt')
 open(14,file="T1E.bmp", access='direct',recl=1)

    dt = 1.d0   !time step
    radius = 1188.*1000.
    c_lat = CENTERLATITUDE *pi/180.d0 
    c_lon = CENTERLONGITUDE *pi/180.d0 
    c_lon = c_lon -pi

    colat = 1.d0*pi/180.d0
    coarea  = 2.d0*pi*radius*radius*(1.d0-dcos(colat))
    k = 0
    k2=0
    rotation = dble(6.387*24.) ! rotation period
    gm =1.30d22*6.67408d-11  
    g = gm/(radius*radius)
    w = 2.*pi/(rotation*3600.) !rad/s

     u1=gm ! GM of Pluto
     u2=1.90d21*6.67408d-11 ! GM of Charon
     aa=19570.*1000. ! distance between Pluto and Charon
     a1= u2/(u1+u2)
     a2= u1/(u1+u2)
     exit_radius = aa*5.d0  

   do loop2 = 1, num2
      polar_lon = 0.d0*pi/180.d0
      polar_lat = 0.d0*pi/180.d0
      nz = dsin(polar_lat)
      nx = dcos(polar_lon)*dcos(polar_lat)
      ny = dsin(polar_lon)*dcos(polar_lat)
      Q = dble(loop2)*2.d0*pi/dble(num2)
      cQ = dcos(Q)
      sQ = dsin(Q)

   do loop = 1, num-1
    zeta = abs(dble(loop2)*2.d0*pi/dble(num2)-volaspect*pi/180.d0)
    if (zeta .gt. pi ) then;  zeta=2.*pi-zeta ; endif
    ob_u=u*(1.d0+0.5*cos(zeta)*cos(impact_angle*pi/180.d0))
    ob_c1=c1*exp(-5.d0*cos(zeta)*cos(impact_angle*pi/180.d0))
    ob_kk=kk*exp(-0.02*cos(zeta)*cos(impact_angle*pi/180.d0))
    ob_c_radius=c_radius*(1.d0-(90.d0-impact_angle)*cos(zeta)/200.)

    lunch_point = n2*ob_c_radius*dble(loop)/(radius*dble(num)) 
    if ((Target .eq. 1) .or.(Target .eq. 4) .or.(Target .eq. 5) .or.(Target .eq. 6)) then 
    lunch = ob_c1*(h*(4.*pi/3.)**(1./3.))**((-2.-ob_u)/(2.*ob_u))*dsqrt(g*ob_c_radius) &
      *(lunch_point/(ob_c_radius/radius))**(-1./ob_u) &
      *(1.d0 -  lunch_point/(n2*ob_c_radius/radius) )**p
    else
    lunch = ob_c1*(h*(4.*pi/3.)**(1./3.))**(-1./ob_u)*dsqrt(yy/ro) &
      *(lunch_point/(ob_c_radius/radius))**(-1./ob_u) &
      *(1.d0 - lunch_point/(n2*ob_c_radius/radius) )**p
    endif
     pp(2) = radius*dsin(lunch_point) !positive y , 180 degree east longitude　
     pp(1) = radius*dcos(lunch_point) !positive x , 90 degree east longitude
    pp(3) = 0.d0                      !positive z , north pole
    rtem2 = atan2(pp(2), pp(1))
    pv(1) = lunch *cos(rtem2+0.5*pi-AANGLE*pi/180.)
    pv(2) = lunch *sin(rtem2+0.5*pi-AANGLE*pi/180.)
    pv(3) = 0.d0


    !ejecta volume
    old_lunch_point = n2*ob_c_radius*dble(loop-1)/(radius*dble(num))
    vol = abs(ob_kk*( (lunch_point*radius)**3.d0 - (old_lunch_point*radius)**3.d0)/dble(num2))

    cx = pp(1)
    cy = pp(2)
    cz = pp(3)
      mcx = ( cQ+nx*nx*(1.d0-cQ) )*cx      &
           +(nx*ny*(1.d0-cQ)-nz*sQ)*cy      &
           +(nx*nz*(1.d0-cQ)+ny*sQ)*cz
      mcy = (ny*nx*(1.d0-cQ)+nz*sQ)*cx      &
           +(cQ+ny*ny*(1.d0-cQ)  )*cy      &
           +(ny*nz*(1.d0-cQ)-nx*sQ)*cz
      mcz = (nz*nx*(1.d0-cQ)-ny*sQ)*cx      &
           +(nz*ny*(1.d0-cQ)+nx*sQ)*cy      &
           +(cQ+nz*nz*(1.d0-cQ) )*cz
    pp(1) = mcx
    pp(2) = mcy
    pp(3) = mcz

    cx = pv(1)
    cy = pv(2)
    cz = pv(3)
      mcx = ( cQ+nx*nx*(1.d0-cQ) )*cx      &
           +(nx*ny*(1.d0-cQ)-nz*sQ)*cy      &
           +(nx*nz*(1.d0-cQ)+ny*sQ)*cz
      mcy = (ny*nx*(1.d0-cQ)+nz*sQ)*cx      &
           +(cQ+ny*ny*(1.d0-cQ)  )*cy      &
           +(ny*nz*(1.d0-cQ)-nx*sQ)*cz
      mcz = (nz*nx*(1.d0-cQ)-ny*sQ)*cx      &
           +(nz*ny*(1.d0-cQ)+nx*sQ)*cy      &
           +(cQ+nz*nz*(1.d0-cQ) )*cz
    pv(1) = mcx
    pv(2) = mcy
    pv(3) = mcz

      c_nz = 0.d0 
      c_nx = 0.d0 
      c_ny = 1.d0 
      c_Q = -c_lat
      c_cQ = dcos(c_Q)
      c_sQ = dsin(c_Q)
    cx = pp(1)
    cy = pp(2)
    cz = pp(3)
      mcx = ( c_cQ+c_nx*c_nx*(1.d0-c_cQ) )*cx      &
           +(c_nx*c_ny*(1.d0-c_cQ)-c_nz*c_sQ)*cy      &
           +(c_nx*c_nz*(1.d0-c_cQ)+c_ny*c_sQ)*cz
      mcy = (c_ny*c_nx*(1.d0-c_cQ)+c_nz*c_sQ)*cx      &
           +(c_cQ+c_ny*c_ny*(1.d0-c_cQ)  )*cy      &
           +(c_ny*c_nz*(1.d0-c_cQ)-c_nx*c_sQ)*cz
      mcz = (c_nz*c_nx*(1.d0-c_cQ)-c_ny*c_sQ)*cx      &
           +(c_nz*c_ny*(1.d0-c_cQ)+c_nx*c_sQ)*cy      &
           +(c_cQ+c_nz*c_nz*(1.d0-c_cQ) )*cz
    pp(1) = mcx
    pp(2) = mcy
    pp(3) = mcz
    cx = pv(1)
    cy = pv(2)
    cz = pv(3)
      mcx = ( c_cQ+c_nx*c_nx*(1.d0-c_cQ) )*cx      &
           +(c_nx*c_ny*(1.d0-c_cQ)-c_nz*c_sQ)*cy      &
           +(c_nx*c_nz*(1.d0-c_cQ)+c_ny*c_sQ)*cz
      mcy = (c_ny*c_nx*(1.d0-c_cQ)+c_nz*c_sQ)*cx      &
           +(c_cQ+c_ny*c_ny*(1.d0-c_cQ)  )*cy      &
           +(c_ny*c_nz*(1.d0-c_cQ)-c_nx*c_sQ)*cz
      mcz = (c_nz*c_nx*(1.d0-c_cQ)-c_ny*c_sQ)*cx      &
           +(c_nz*c_ny*(1.d0-c_cQ)+c_nx*c_sQ)*cy      &
           +(c_cQ+c_nz*c_nz*(1.d0-c_cQ) )*cz
    pv(1) = mcx
    pv(2) = mcy
    pv(3) = mcz
      c_nz = 1.d0
      c_nx = 0.d0 
      c_ny = 0.d0 
      c_Q = c_lon
      c_cQ = dcos(c_Q)
      c_sQ = dsin(c_Q)
    cx = pp(1)
    cy = pp(2)
    cz = pp(3)
      mcx = ( c_cQ+c_nx*c_nx*(1.d0-c_cQ) )*cx      &
           +(c_nx*c_ny*(1.d0-c_cQ)-c_nz*c_sQ)*cy      &
           +(c_nx*c_nz*(1.d0-c_cQ)+c_ny*c_sQ)*cz
      mcy = (c_ny*c_nx*(1.d0-c_cQ)+c_nz*c_sQ)*cx      &
           +(c_cQ+c_ny*c_ny*(1.d0-c_cQ)  )*cy      &
           +(c_ny*c_nz*(1.d0-c_cQ)-c_nx*c_sQ)*cz
      mcz = (c_nz*c_nx*(1.d0-c_cQ)-c_ny*c_sQ)*cx      &
           +(c_nz*c_ny*(1.d0-c_cQ)+c_nx*c_sQ)*cy      &
           +(c_cQ+c_nz*c_nz*(1.d0-c_cQ) )*cz
    pp(1) = mcx
    pp(2) = mcy
    pp(3) = mcz
    cx = pv(1)
    cy = pv(2)
    cz = pv(3)
      mcx = ( c_cQ+c_nx*c_nx*(1.d0-c_cQ) )*cx      &
           +(c_nx*c_ny*(1.d0-c_cQ)-c_nz*c_sQ)*cy      &
           +(c_nx*c_nz*(1.d0-c_cQ)+c_ny*c_sQ)*cz
      mcy = (c_ny*c_nx*(1.d0-c_cQ)+c_nz*c_sQ)*cx      &
           +(c_cQ+c_ny*c_ny*(1.d0-c_cQ)  )*cy      &
           +(c_ny*c_nz*(1.d0-c_cQ)-c_nx*c_sQ)*cz
      mcz = (c_nz*c_nx*(1.d0-c_cQ)-c_ny*c_sQ)*cx      &
           +(c_nz*c_ny*(1.d0-c_cQ)+c_nx*c_sQ)*cy      &
           +(c_cQ+c_nz*c_nz*(1.d0-c_cQ) )*cz
    pv(1) = mcx
    pv(2) = mcy
    pv(3) = mcz

          im_lat = atan(pp(3) / sqrt(pp(1)**2.d0+pp(2)**2.d0))*180.d0/pi
          im_lon = atan2(pp(2), pp(1))*180.d0/pi
          im_lon = im_lon -180.
          if ( im_lon .le. 0.d0 ) then
           im_lon=im_lon+360.d0
          endif
          k2 = k2 +1
          record2(k2,1) = im_lat
          record2(k2,2) = im_lon
 

    !Pluto position x=-a1 adjustment orbits
     pp(1)=pp(1)-a1


  !Adams-Bashforth Method
      !first step
        !gravity force
           distance1 =sqrt( (pp(1)+a1)**2.+pp(2)**2.+pp(3)**2.)
           distance2 =sqrt( (pp(1)-a2)**2.+pp(2)**2.+pp(3)**2.)
              rtem1 = dt*u1/distance1**3.
              rtem2 = dt*u2/distance2**3.
              dp(1) =   -(pp(1)+a1)*rtem1-(pp(1)-a2)*rtem2
              dp(2) =   -pp(2)*rtem1-pp(2)*rtem2
              dp(3) =   -pp(3)*rtem1-pp(3)*rtem2
         !Coriolis force
             dp(1) = dp(1) +dt*2.d0*w*pv(2)
             dp(2) = dp(2) -dt*2.d0*w*pv(1)
         !centrifugal force
             dp(1) = dp(1) +dt*w*w*pp(1)
             dp(2) = dp(2) +dt*w*w*pp(2)
         dx(1) = pv(1)*dt 
         dx(2) = pv(2)*dt 
         dx(3) = pv(3)*dt 
         do i = 1, 3
         pv(i)=pv(i)+dp(i)
         pp(i)=pp(i)+dx(i)
         dp_prev(i) = dp(i)
         dx_prev(i) = dx(i)
         enddo

      !second step  
      do j = 1, int(1000000./dt)
        !gravity force
       distance1 =sqrt( (pp(1)+a1)**2.+pp(2)**2.+pp(3)**2.)
       distance2 =sqrt( (pp(1)-a2)**2.+pp(2)**2.+pp(3)**2.)
        if (distance1 .gt. exit_radius  ) exit
              rtem1 = dt*u1/distance1**3.
              rtem2 = dt*u2/distance2**3.
              dp(1) =   -(pp(1)+a1)*rtem1-(pp(1)-a2)*rtem2
              dp(2) =   -pp(2)*rtem1-pp(2)*rtem2
              dp(3) =   -pp(3)*rtem1-pp(3)*rtem2
        !Coriolis force
              dp(1) = dp(1) +dt*2.d0*w*pv(2)
              dp(2) = dp(2) -dt*2.d0*w*pv(1)
        !centrifugal force
              dp(1) = dp(1) +dt*w*w*pp(1)
              dp(2) = dp(2) +dt*w*w*pp(2)
         dx(1) = pv(1)*dt
         dx(2) = pv(2)*dt
         dx(3) = pv(3)*dt
         do i = 1, 3
         pv(i)=pv(i)+1.5*dp(i)-0.5*dp_prev(i)
         pp(i)=pp(i)+1.5*dx(i)-0.5*dx_prev(i)
         dp_prev(i) = dp(i)
         dx_prev(i) = dx(i)
         enddo

        !after landing of particles
         if ( distance1 .le. radius ) then
          pp(1)=pp(1)+a1
          im_lat = atan(pp(3) / sqrt(pp(1)**2.d0+pp(2)**2.d0 ) )*180.d0/pi
          im_lon = atan2(pp(2), pp(1))*180.d0/pi
          im_lon = im_lon -180.
            if ( im_lon .le. 0.d0 ) then
            im_lon=im_lon+360.d0
            endif
          k = k +1
          record(k,1) = im_lat
          record(k,2) = im_lon
          record(k,3) = vol
          record(k,4) = pp(1)
          record(k,5) = pp(2)
          record(k,6) = pp(3)
            exit
           endif
     enddo
  enddo !loop
  enddo !loop2


 !calculate ejecta thickness
   knum = k
  rtem3 = pi/180.d0
     do j = -90,90
       clat = dble(j)*rtem3
     do i = 0, 359
       clon = dble(i)*rtem3
      do k = 1, knum
       rtem1 = record(k,1)*rtem3
       if (rtem1 .gt. clat+colat) cycle
       if (rtem1 .le. clat-colat) cycle
       rtem2 = record(k,2)*rtem3
  distance1 = dacos(dsin(rtem1)*dsin(clat)+  &
  dcos(rtem1)*dcos(clat)*dcos(rtem2-clon))
    if ( distance1 .le.  colat) then
     density(j,i) = density(j,i)+ record(k,3)
    endif
       enddo
     density(j,i) = density(j,i)/coarea !num/km2
    enddo
    enddo
    do i = 1, k2 
      density(nint(record2(i,1)),nint(record2(i,2)))=0.d0
    enddo

!     do j = -90,90
!     do i = 0, 360-1
!   write(6,*) i,j,density(j,i)
!    enddo
!    enddo

     do k = 1,knum
       if (density(nint(record(k,1)),nint(record(k,2))) .ne. 0.d0 ) then
          write(35,*) record(k,1),record(k,2),record(k,3),record(k,4),record(k,5),record(k,6)
       endif
     enddo
     close(35)

   maxdensity = maxval( density )
   mindensity = 1000000000.d0
    do j = -90,90
    do k = 0, 360-1
      if (density(j,k) .eq. 0.d0) cycle
     if (mindensity .gt. density(j,k)) then 
     mindensity = density(j,k)
     endif
    enddo
    enddo
   mindensity =  log10(mindensity) 
    do j = -90,90
    do k = 0, 360-1
      if (density(j,k) .eq. 0.d0) then
      density(j,k) =  mindensity 
      else 
      density(j,k) = log10(density(j,k))
      endif
    enddo
    enddo
   maxdensity = dble( int( log10(maxdensity) )+1 )
   mindensity =  dble( int( mindensity )-1 )
   write(6,*) mindensity,maxdensity,rotation,radius,c_radius,c_lat,c_lon


  !color map export
   i=55
   open(12,file="color-bar5.bmp",access='direct',recl=1)
   do k=-90,90
   do j=1,1000
    read(12,rec=i) col(j,k,1)
    i=i+1
    read(12,rec=i) col(j,k,2)
    i=i+1
    read(12,rec=i) col(j,k,3)
    i=i+1
   enddo
   enddo
   close(12)

  ini=0;ini(1)=66;ini(2)=77;ini(3)=18;ini(4)=-4;ini(5)=2;ini(11)=122;ini(15)=108;ini(19)=104;ini(20)=1;ini(23)=-75
  ini(27)=1;ini(29)=24;ini(35)=-104;ini(36)=-5;ini(37)=2;ini(38)=65;ini(39)=92;ini(42)=65;ini(43)=92
  do j=1,54
   write(14,rec=j) ini(j)
  enddo

     i = 123
     do j = -90,90
     do k = 0, 359
     bina= int( (density(j,k)-mindensity)*1000.d0 / (maxdensity-mindensity) )+1
    write(14,rec=i) col(bina,0,1)
    i=i+1
    write(14,rec=i) col(bina,0,2)
    i=i+1
    write(14,rec=i) col(bina,0,3)
    i=i+1
    enddo
    enddo
   close(14)

   end

