!===============================================================================
!     Fortran 90/95 Source Code File
!===============================================================================

!===============================================================================
   subroutine create_synthetics(E,N,Z,dfreq,noise,spol,nop,fast,tlag)
!===============================================================================
      use f90sac  ! use the f90sac module
!-------------------------------------------------------------------------------
      implicit none
      type (SACTrace) :: E,N,Z
      

      integer i,nop,iop
      integer, parameter :: nop_max = 10
      
      real :: fast(nop_max),tlag(nop_max)
      
      real :: noise,delta,dfreq,spol
      
!  ** generate an appropriate sample rate for the required dominant frequency
      delta = 1.0/(dfreq*200.)

!  ** initialise the random noise generator
      call f90sac_init_random()

!  ** make the initial pulse.       
      call source_pulse(350,350,spol,delta,noise,E,N,Z)

!  ** apply the splitting operators
      do iop = 1,nop
         call apply_split_op(N,E,fast(iop),tlag(iop))
      enddo
      
!  ** calculate a sensible analysis window, basically the original wave 
!  ** extended for the sum total of the splitting tlags.
      N%a = -real(150)*delta-sum(tlag(1:nop)) ; 
      N%f =  real(150)*delta+sum(tlag(1:nop)) ; 
      E%a = -real(150)*delta-sum(tlag(1:nop)) ; 
      E%f =  real(150)*delta+sum(tlag(1:nop)) ; 
      Z%a = -real(150)*delta-sum(tlag(1:nop)) ; 
      Z%f =  real(150)*delta+sum(tlag(1:nop)) ; 
         
      N%user0 = N%a ; N%user1 = N%a+delta   
      N%user2 = N%f ; N%user3 = N%f+delta   
      E%user0 = E%a ; E%user1 = E%a+delta   
      E%user2 = E%f ; E%user3 = E%f+delta   
      Z%user0 = Z%a ; Z%user1 = Z%a+delta   
      Z%user2 = Z%f ; Z%user3 = Z%f+delta   
         
   end subroutine create_synthetics
!===============================================================================
   
!===============================================================================
   subroutine apply_split_op(h1,h2,phi,dt)
!===============================================================================
!
!     Apply shear-wave splitting
!
!     parameters:
!  
!     h1,h2 :  (I/O) (SACtrace) input horizontal, orthogonal components
!     dt,phi:  (I)   (real) lag-time, fast direction to correct by
!
!-------------------------------------------------------------------------------
      use f90sac  ! use the f90sac module
!-------------------------------------------------------------------------------
      implicit none
      
      type (SACtrace) :: h1, h2
      real dt , phi
      real cmpazdiff,rotang

!  ** calculate the angle of rotation, based on phi and the azimuth of h1
      rotang = (phi - h1 % cmpaz)
   
!  ** rotate the traces      
      call f90sac_rotate2d(h1,h2,rotang)      

!  ** time shift the slow trace
      call f90sac_tshift(h2,dt)

!  ** rotate back to original orientation    
      call f90sac_rotate2d(h1,h2,-rotang)      

!  * done      
      return
   end
!===============================================================================


!===============================================================================
   subroutine source_pulse(nzpadb,nzpade,spol,delta,noise,E,N,Z)
!===============================================================================
      use f90sac  ! use the f90sac module
!-------------------------------------------------------------------------------
      implicit none

      integer :: nzpadb,nzpade ! number of zeros samples to add to the beginning
                               ! and end of the trace
      
      real :: spol,delta ! sample rate and dominant frequency of the trace                         
      real :: wave(301), noise
      real,parameter :: pi = 3.1415927410125732421875
      type (SACTrace) :: E,N,Z
      
      integer :: nspperiod
                  
!  ** create the North trace and set up the headers
      call f90sac_newtrace(nzpadb+nzpade+301,delta,N)
            
      call f90sac_setdate(N,3001,1,12,00,00,00)
      N%stla = 0.0 ; N%stlo = 0.0
      call givloc(90.0,0.0,80.0,spol,N%evla,N%evlo)
      call stadis2(N%evla,N%evlo,N%stla,N%stlo,N%gcarc,N%az)
      call stadis2(N%stla,N%stlo,N%evla,N%evlo,N%gcarc,N%baz)
      
      N%kstnm = 'SWAV'
      N%evdp = 0.0
      
      N%o = -real(nzpadb+150)*delta
      N%b = -real(nzpadb+150)*delta
      N%e =  real(nzpade+150)*delta

      N%cmpaz=0.0 ; N%cmpinc=90.0; 
      
      call f90sac_clonetrace(N,E)
      E%cmpaz=90.0 ; E%cmpinc=90.0; 

      call f90sac_clonetrace(N,Z)
      Z%cmpaz=0.0 ; Z%cmpinc=0.0; 
      
      call simple_wavelet(wave)
      
      N%trace(nzpadb+1:nzpadb+301) = wave(1:301)*cos(spol*pi/180.0)
      E%trace(nzpadb+1:nzpadb+301) = wave(1:301)*sin(spol*pi/180.0)
      
      call f90sac_addwnoise(N,-noise)
      call f90sac_addwnoise(E,-noise)
      call f90sac_addwnoise(Z,-noise)
      
!      call f90sac_writetrace('splitwave.BHN',N)
!      call f90sac_writetrace('splitwave.BHE',E)
!      call f90sac_writetrace('splitwave.BHZ',Z)
      
         
      return
      end
!===============================================================================

!===============================================================================
   subroutine simple_wavelet(wave)
!===============================================================================
      implicit none
      real :: wave(301),wave_par(301)
      data wave_par / &
       0.000000E+00, -0.409567E-05, -0.929897E-05,  0.131005E-04,  0.809800E-04, &
       0.162307E-03,  0.210853E-03,  0.230160E-03,  0.215660E-03,  0.176950E-03, &
       0.111436E-03,  0.946606E-05, -0.110559E-03, -0.224143E-03, -0.331512E-03, &
       0.433529E-03, -0.508395E-03, -0.538056E-03, -0.405200E-03, -0.112606E-03, &
       0.179249E-03,  0.312222E-03,  0.304611E-03,  0.287822E-03,  0.271034E-03, &
       0.263422E-03,  0.536212E-03,  0.113495E-02,  0.173520E-02,  0.200775E-02, & 
       0.196110E-02,  0.185827E-02,  0.175567E-02,  0.170889E-02,  0.189211E-02, & 
       0.233495E-02,  0.288161E-02,  0.337194E-02,  0.377469E-02,  0.417014E-02, & 
       0.457180E-02,  0.499612E-02,  0.542197E-02,  0.585339E-02,  0.633546E-02, & 
       0.691277E-02,  0.770392E-02,  0.870161E-02,  0.977783E-02,  0.107972E-01, & 
       0.116953E-01,  0.125567E-01,  0.134580E-01,  0.144847E-01,  0.156876E-01, & 
       0.170432E-01,  0.185302E-01,  0.201137E-01,  0.218216E-01,  0.236829E-01, & 
       0.256828E-01,  0.278215E-01,  0.301152E-01,  0.326009E-01,  0.352618E-01, & 
       0.380798E-01,  0.411071E-01,  0.443548E-01,  0.478079E-01,  0.514282E-01, & 
       0.552359E-01,  0.592812E-01,  0.635172E-01,  0.679372E-01,  0.725775E-01, & 
       0.774461E-01,  0.825008E-01,  0.876511E-01,  0.928404E-01,  0.981224E-01, & 
       0.103640E+00,  0.109575E+00,  0.116115E+00,  0.123219E+00,  0.130598E+00, & 
       0.137959E+00,  0.145253E+00,  0.152597E+00,  0.160112E+00,  0.167863E+00, & 
       0.175948E+00,  0.184347E+00,  0.192837E+00,  0.201277E+00,  0.209573E+00, & 
       0.217809E+00,  0.226115E+00,  0.234536E+00,  0.243166E+00,  0.251933E+00, & 
       0.260744E+00,  0.269571E+00,  0.278538E+00,  0.287663E+00,  0.296570E+00, & 
       0.304881E+00,  0.312576E+00,  0.319853E+00,  0.326835E+00,  0.333595E+00, & 
       0.340280E+00,  0.346851E+00,  0.353010E+00,  0.358519E+00,  0.363643E+00, & 
       0.368513E+00,  0.372689E+00,  0.375685E+00,  0.377918E+00,  0.379888E+00, & 
       0.381288E+00,  0.381824E+00,  0.381264E+00,  0.379786E+00,  0.377701E+00, & 
       0.375325E+00,  0.372089E+00,  0.367543E+00,  0.362110E+00,  0.356247E+00, & 
       0.349776E+00,  0.342287E+00,  0.333978E+00,  0.324973E+00,  0.315036E+00, & 
       0.303985E+00,  0.291993E+00,  0.279345E+00,  0.265663E+00,  0.250824E+00, & 
       0.235343E+00,  0.219620E+00,  0.203851E+00,  0.187683E+00,  0.171036E+00, & 
       0.153832E+00,  0.135750E+00,  0.116957E+00,  0.977126E-01,  0.784165E-01, & 
       0.591127E-01,  0.395471E-01,  0.199320E-01,  0.288529E-03, -0.193922E-01, &
      -0.390953E-01, -0.587637E-01, -0.781489E-01, -0.974606E-01, -0.116685E+00, &
      -0.135503E+00, -0.153732E+00, -0.171326E+00, -0.188518E+00, -0.205146E+00, &
      -0.221053E+00, -0.236481E+00, -0.251411E+00, -0.265644E+00, -0.278871E+00, &
      -0.291199E+00, -0.302902E+00, -0.313761E+00, -0.323665E+00, -0.332829E+00, &
      -0.341399E+00, -0.349144E+00, -0.355756E+00, -0.361703E+00, -0.367223E+00, &
      -0.371759E+00, -0.374785E+00, -0.376785E+00, -0.378489E+00, -0.379671E+00, &
      -0.380113E+00, -0.379661E+00, -0.378475E+00, -0.376795E+00, -0.374873E+00, &
      -0.372254E+00, -0.368543E+00, -0.364124E+00, -0.359338E+00, -0.353889E+00, &
      -0.347525E+00, -0.340648E+00, -0.333719E+00, -0.326780E+00, -0.319612E+00, &
      -0.312196E+00, -0.304458E+00, -0.296308E+00, -0.287732E+00, -0.278944E+00, &
      -0.270163E+00, -0.261317E+00, -0.252370E+00, -0.243397E+00, -0.234539E+00, &
      -0.225829E+00, -0.217165E+00, -0.208607E+00, -0.200128E+00, -0.191742E+00, &
      -0.183461E+00, -0.175263E+00, -0.167206E+00, -0.159150E+00, -0.151137E+00, &
      -0.143385E+00, -0.136055E+00, -0.129272E+00, -0.122841E+00, -0.116601E+00, &
      -0.110391E+00, -0.104051E+00, -0.977298E-01, -0.915963E-01, -0.858647E-01, &
      -0.805372E-01, -0.754590E-01, -0.706753E-01, -0.661829E-01, -0.619642E-01, &
      -0.579870E-01, -0.542345E-01, -0.507278E-01, -0.474727E-01, -0.444348E-01, &
      -0.415246E-01, -0.386311E-01, -0.356415E-01, -0.326238E-01, -0.298041E-01, &
      -0.274067E-01, -0.254740E-01, -0.238108E-01, -0.222383E-01, -0.205899E-01, &
      -0.187558E-01, -0.168405E-01, -0.150237E-01, -0.134675E-01, -0.121416E-01, &
      -0.109317E-01, -0.988059E-02, -0.904162E-02, -0.842485E-02, -0.793914E-02, &
      -0.747276E-02, -0.691073E-02, -0.607618E-02, -0.505447E-02, -0.413950E-02, &
      -0.362409E-02, -0.345904E-02, -0.336215E-02, -0.326891E-02, -0.311546E-02, &
      -0.270829E-02, -0.205794E-02, -0.141296E-02, -0.101638E-02, -0.854278E-03, &
      -0.735190E-03, -0.641102E-03, -0.555095E-03, -0.478753E-03, -0.415893E-03, &
      -0.351564E-03, -0.270360E-03, -0.129518E-03,  0.581389E-04,  0.222125E-03, & 
       0.292185E-03,  0.229472E-03,  0.823665E-04, -0.858300E-04, -0.216284E-03, &
      -0.318124E-03, -0.397662E-03, -0.417479E-03, -0.365925E-03, -0.191173E-03, & 
       0.638601E-05,  0.100210E-03,  0.829807E-04,  0.373569E-04,  0.929093E-05, & 
       0.000000E+00 /
         
         wave = wave_par 
         
      return   
   end subroutine simple_wavelet 
!===============================================================================

!===============================================================================
      subroutine STADIS2(qlat,qlon,slat,slon,del,az)
!===============================================================================
! Computes the epicentral distance and azimuth from source to receiver.
! Latitudes are converted to geocentric latitudes prior to performing
! the computations (it is assumed that input latitudes are geographic).
!  Inputs:   qlat  =  quake latitude (degrees)
!            qlon  =  quake longitude (degrees)
!            slat  =  station latitude (degrees)
!            slon  =  station longitude (degrees) 
!  Returns:  del     =  epicentral distance (degrees)
!            az      =  azimuth at quake to receiver, from North (degrees)
!
!  This version calculates doesn't require colats and colons. JW 2007
!
      real co,si,caz,saz
      data rad/57.29578/
      qcolat=90.-qlat
      qcolon=qlon
      if (qcolon.lt.0.) qcolon=qcolon+360.
      scolat=90.-slat
      scolon=slon
      if (scolon.lt.0.) scolon=scolon+360.

      if (qcolat.eq.scolat.and.qcolon.eq.scolon) then
         del=0.
         az=999.
         return
      end if
!      t1deg=scolon
!      p1deg=scolat              !this bug corrected 2/28/90
      t1deg=scolat
      p1deg=scolon
      colat=qcolat
      colon=qcolon
      t1=t1deg/rad
      p1=p1deg/rad
!  first do eq coords.
      t0=geocen(colat/rad)    ! use geocentric e.q.colat }
      p0=colon/rad
      c0=cos(t0)
      s0=sin(t0)
!  now do station coords.
      t2=geocen(t1)           ! use geocentric station colat }
      c1=cos(t2)
      s1=sin(t2)
!  now calculate distance
      dp=p1-p0
      co=c0*c1+s0*s1*cos(dp)
      si=dsqrt(1.d0-co*co)
      del=atan2(si,co)*rad
!  now calculate azimuth
      caz=(c1-c0*co)/(si*s0)
      dp2=-dp
      saz=-s1*sin(dp2)/si
      az=atan2(saz,caz)*rad
      if(az.lt.0.0) az=360.0 + az  !change az to be between 0 and 360}
      return
      end subroutine STADIS2
!===============================================================================

!===============================================================================
      subroutine givloc(colat,colon,del,az,t1,p1)
!===============================================================================
! input:
!   colat,colon = source colatitude and colongitude from sub. source
!   del         = distance in degrees
!   az          = azimuth of station from source
! output:
!   t1 = station latitude  ( + = N, - = S)
!   p1 = station longitude ( + = E, - = W)
!  
!===============================================================================
   implicit real (a-h,o-z)

      integer, parameter :: r4 = selected_real_kind(6,37) ;
      integer, parameter :: r8 = selected_real_kind(15,307) ; 
      integer, parameter :: rs = r8 ; 

      data rad/57.29578/
      delr=del/rad
      azr=az/rad
      t0=geocen(colat/rad)    !convert to geocentric}
      ctheta=sin(delr)*sin(t0)*cos(azr) + cos(t0)*cos(delr)
      t1=acos(ctheta)
      if (t0.eq.0.0) then
        p1=az
      elseif (t1.eq.0.0) then
        p1=0.0
      else
        sphi=sin(delr)*sin(azr)/sin(t1)
        cphi=(cos(delr) - cos(t0)*ctheta)/(sin(t0)*sin(t1))
        p1=colon + atan2(sphi,cphi)*rad
      endif
      t1=90.0-geogrf(t1)*rad      !convert colatitude to geograf. latitude}
      if (p1.gt.360.0) p1 = p1 - 360.0   ! assume p1 never > 720 }
      if (p1.gt.180.0) p1 = p1 - 360.0   !convert colongitude to longitude}
      return
      end subroutine givloc
!===============================================================================


!===============================================================================
      real function geocen(arg)
!===============================================================================
! input:
!   arg    = geographic colatitude (radians)
! output:
!   geocen = geocentric colatitude (radians)
! (n.b. fac=(1-f)**2)
!===============================================================================
      data pi2,fac/1.570796326794895,0.993305621334896/
      geocen=pi2-atan(fac*cos(arg)/amax1(1.e-30,sin(arg)))   
      return
      end function geocen
!===============================================================================

!===============================================================================
      real function geogrf(arg)
!===============================================================================
! input:
!   arg    = geocentric colatitude (radians)
! output:
!   geogrf = geographic colatitude (radians
! (n.b. fac=(1-f)**2)
!
!===============================================================================
      data pi2,fac/1.570796326794895,0.993305621334896/
      geogrf=pi2-atan(cos(arg)/(fac*amax1(1.e-30,sin(arg)))) 
      return                       
      end function geogrf
!===============================================================================
