!=======================================================================
!  S H E B A - Shear-wave Birefringence Analysis
!=======================================================================
!  This software is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!=======================================================================
!
!  James Wookey, School of Earth Sciences, University of Bristol
!  CVS: $Revision: 1.4 $ $Date: 2008/09/19 00:20:08 $
!
!-----------------------------------------------------------------------
!
! SUBROUTINES IN THIS FILE WERE WRITTEN BY N. TEANBY UNIVERSITY OF LEEDS
! BASED ON TEANBY AND KENDALL (2003) (?)

C=======================================================================
      subroutine zsplit(x0,y0,n,wbeg,wend,delta,b,tlag_scale,
     >fast,dfast,tlag,dtlag,spol,dspol,error,error_int,lam1,lam1_int,f,
     >lambda2_min,ndf,snr)
C=======================================================================
c
c      subroutine to perform splitting correction of Silver and Chan 1991
c      on two orthogonal components contained in file1/2.
c
c      this method performs a grid search over a range of fast directions
c      and lag times in order to find the fast direction lag time which 
c      gives the minimum second eigenvalue of the particle motion covaraince
c      matrix
c
c      interpolation in the tlag direction is used to increase resolution
c
c      Method
c      ------
c      -x0 and y0 are two orthogonal components to do the analysis on
c       (ideally orthogonal to ray direction).
c      -The S-wave window, delta, and b are specified in the subroutine call.
c      -Detrend the data (using a least squares fit to datapoints)
c      -Window the data. The data window has 'iwextra' points after the specified
c       end point. this is so that when lagging is done the number of points
c       within the window is constant. Hence, data points outside the window are
c       also used.
c      -Map out the error surface over fast=-90 to 90, 
c       lag=0-40 delta * itlag_step. The error surface is the
c       value of lambda2 (smallest eigenvalue from partical motion
c       covariance matrix).
c      -Interpolate error surface in tlag direction by the factor
c       f=(np2int-1)/(np1-1).
c      -Grid search to find minimum lambda2 of interpolated error surface.
c       This point corresponds to the fast direction and tlag which best
c       correct for the splitting.
c      -Interpolate the time series to same resolution as error surface, so that 
c       the correct value of lag can be used in the calculation of ndf, and
c       source polarisation
c      -Window the interpolated time series.
c      -Rotate, lag, then calc covariance matrix, eigenvalues, and source
c       polarisation
c      -Use source polarisation and fast direction to rotate seismograms into the
c       noise direction and use noise to calc ndf.
c      -Use ndf and lambda2_min to calc the value of the 95% confidence contour
c       and normalise interpolated error surface so 1 = 95% confidence interval.
c       (return  both interpolated and uninterpolated error surfaces because
c       uninterpolated one is used for stacking).
c
c      Variables
c      ---------
c    input:
c      x0(np)/y0(np)      real            two orthogonal seismogram componentstim
c      n                  int            number of points in entire time series
c      wbeg/wend            real            S-wave window to use in analysis
c      delta                  real            sampling interval, read from SAC header
c      b                  real            start time of time series
c      tlag_scale            real            max lag time of error surface
c                                    (reset to an integer multiple of 
c                                    (np2-1)*delta, where integer is =>1)
c      np1                  int            number of grid points for fast direction
c      np2                  int            number of grid points for lag time
c      np2int            int            number of grid points for lag time
c                                     (after interp)
c
c    output:
c      fast                  real            fast direction (deg clock from N)
c      dfast                  real            1.s.d.
c      tlag                  real            lag time in seconds
c      dtlag                  real            1 s.d.
c      spol                  real            s-wave source polarisation(deg clock from N)
c      dspol                  real            Maximum Angular Deviation (MAD)
c      error(np1,np2)      real            uninterpolated error surface
c      error_int(np1,np2int)      real      interpolated error surface
c      f                  int            interpolation factor=(np2int-1)/(np1-1)
c                                     MUST BE and integer for zsplint to work
c      lambda2_min            real            minimum lambda2 (corresponds to solution)
c      ndf                  real            no. deg. freedom
c
c    local:
c      np                  int            array dimension
c                                     (read from SIZE_np.h at compile time)
c      norig                  int            original (uninterpolated) number of points
c                                     in the shear wave window
c      nwindow            int            number of points in the S-wave window
c                                     (including the extra points to allow for 
c                                     the lag)
c      noverlap            int            number of overlapping points after lag has 
c                                     been applied
c      ninterp            int            number of points after interpolation
c
c      iwbeg/iwend            int            indecices of beginning and end of s-wave
c                                     analysis window
c      iwbegx/iwendx      int            indecices of beginning and end of s-wave
c                                     total window
c      iwextra            int            number of extra points to include on the
c                                     end of the window to allow for lagging
c                                     (= max_lag = np2 for uniterpolated case)
c
c      x(np)/y(np)            real            time series after trend removal
c      x/ywindow(np)      real            windowed time series
c      x/yrot(np)            real            rotated time series
c      x/ylag(np)            real            lagged time series
c      x/yinterp(np)      real            interpolated time series
c      x/ynoise(np)      real            seismograms roated into noise (polarisation)
c                                     direction. xnoise is the noise, ynoise is
c                                     the S-wave signal
c      ndf                  int            number of degrees of freedom in noise signal
c                                     (used to calc 95% confidence level)
c      cov(2,2)            real            covariance matrix of particle motion
c      lambda1/2            real            eignevalues of cov (lambda1 is largest)
c      vec1/2            real            eignevectors of cov (vec1 corresponds to
c                                     lambda1)
c      itlag                  int            index of lag corresponding to lambda2_min
c      idtlag            real            1s.d. in lag in units of interp grid spacing
c      itlag_step            int            gridding step on error surface
c
c-----------------------------------------------------------------------
c      N. Teanby      5-8-02      original code
c      N. Teanby      12-5-03      itlag_step introduced so that tlag_scale can
c                              be independent of number of grid points.
c-----------------------------------------------------------------------
      use array_sizes ! use the array_sizes modules
      use sheba_config
c-----------------------------------------------------------------------
      
      implicit none
      integer n,ninterp,norig,noverlap,nwindow
      real wbeg,wend,delta,b,tlag_scale,fast,dfast,tlag,dtlag,spol,dspol
      real error(np1,np2),error_int(np1,np2int),idtlag,idfast
      real lam1(np1,np2),lam1_int(np1,np2int)
      integer f,iwbeg,iwend,ndf,itlag,ifast,itlag_step
      real x(np),y(np),x0(np),y0(np),xinterp(np),yinterp(np)
      real xnoise(np),ynoise(np)
      real xwindow(np),ywindow(np),xrot(np),yrot(np),xlag(np),ylag(np)
      real xlag1(np),ylag1(np),xlag2(np),ylag2(np)
      real lambda1,lambda2,lambda2_min,vec1(2),vec2(2),cov(2,2)
      integer iwextra,iwbegx,iwendx,itlag1,itlag2
      
      integer i
C     ** Hilbert versions

      real hxwindow(np),hywindow(np),hxrot(np),hyrot(np)
      real hxlag1(np),hylag1(np)
      real xtemp(np), ytemp(np)
      
      real snr

c  ** calc itlag_step from tlag_scale **
c  ** itlag_step is the grid spacing in tlag for the grid search **
c  ** it must be an integer greater than 1 **
c  ** use itlag_step to redefine tlag_scale **
      itlag_step = max(1 ,  nint( tlag_scale/(real(np2-1)*delta) )   )
c  ** redefine tlag_scale **
      tlag_scale = real(itlag_step*(np2-1))*delta

c  ** calculate interpolation factor for error surface from np2 and np2int **
c  ** this factor is used to interpolate the error surface in the tlag dirn **
      f = (np2int-1)/(np2-1)
c  ** f needs to be an integer for zsplint to work **
      if (mod(np2int-1,np2-1).ne.0) then
         pause 'ERROR: zsplit: f not a whole number'
      endif

c  ** calculate the index of the S wave window **
      iwbeg = nint((wbeg-b)/delta)+1
      iwend = nint((wend-b)/delta)+1
c      print*,'wbeg,wend',wbeg,wend,iwbeg,iwend
c  ** calc the extra to add to end of the window to allow for lagging 
c      iwextra = maximum lag = np2*itlag_step **
      iwextra = np2*itlag_step
      iwbegx = iwbeg
      iwendx = iwend + iwextra
      if (iwendx.gt.n) then
         print*,iwendx,n
         print*,'wbeg,wend',b,wbeg,wend,itlag_step
         print*,'wbeg,wend',wbeg,wend,iwbeg,iwend
         pause 'ERROR: zsplit: window out of range of data'
      endif
      
      x(:)=0.0
      y(:)=0.0
            
c  ** calc number of points originally in the window **
      norig=iwend-iwbeg+1

c  ** detrend the data **
      call zdetrend(x0,n,np,x)
      call zdetrend(y0,n,np,y)

c  ** window the data **
      call zwindow(x,y,n,np,iwbegx,iwendx,xwindow,ywindow,nwindow)

c  ** map out the error surface **
      call zgrid_lambda2(xwindow,ywindow,nwindow,iwextra,itlag_step,
     >                        lam1,error,delta)

c  ** interpolate error surface in tlag direction **
      call zerror_interp(error,error_int)
      call zerror_interp(lam1,lam1_int)

c  ** find the interpolated minimum position **
      call zerror_min(error_int,np1,np2int,ifast,itlag,lambda2_min)

c  ** convert indices itlag/ifast into tlag and fast **
      tlag  = delta*real(itlag_step*(itlag-1))/real(f)
      fast = -90. + 180.*real(ifast-1)/real(np1-1)

c  ** interpolate the time series to be the same or as high resolution 
c      as error surface
c      because lag index is now in terms of interpolated times **
      call zsplint(x,n,f,xinterp,ninterp)
      call zsplint(y,n,f,yinterp,ninterp)

c  ** error surface resolution = delta *  itlag_step / f **
c  ** time series resolution   = delta *  1 / f **
c  ** itlag is in terms of error surface index **
c  ** convert itlag to time series index for use by zlag**
      itlag = itlag * itlag_step

c  ** window the interpolated data **
      iwbegx=f*(iwbegx-1)+1
      iwendx=f*(iwendx-1)+1
      iwextra=f*iwextra
      call zwindow(xinterp,yinterp,ninterp,np,iwbegx,iwendx,
     >xwindow,ywindow,nwindow)

C  ** Hilbert transforms
      call hilbert(xwindow,nwindow,hxwindow)
      call hilbert(ywindow,nwindow,hywindow)

c  ** rotate, lag, calc the covariance and eigenvalues **
      call zrotate2d(xwindow,ywindow,nwindow,np,fast,xrot,yrot)
      call zrotate2d(hxwindow,hywindow,nwindow,np,fast,hxrot,hyrot)

c  **  split time lag into two stages - JW 2004    
      itlag1 = nint(real(itlag)/2.0)
      itlag2 = itlag - itlag1
      
c  **  first time lag
      call zlag(xrot,yrot,nwindow,np,itlag1,iwextra,
     >               xlag1,ylag1,noverlap)
      call zlag(hxrot,hyrot,nwindow,np,itlag1,iwextra,
     >               hxlag1,hylag1,noverlap)

c  ** ScS correction - JW 2004
c      call scscorr(xlag1,ylag1,hxlag1,hylag1,nwindow,fast,itlag)

c  **  second time lag
      call zlag(xlag1,ylag1,nwindow,np,itlag2,iwextra,
     >               xlag2,ylag2,noverlap)
C  **
      xlag(1:nwindow) = xlag2(1:nwindow) 
      ylag(1:nwindow) = ylag2(1:nwindow)
      
c  ** perform any required post-correction      
      if (config % i_src_corr == 1) then
         print*,'running source correction'
         
c         print*,tlag,itlag,delta 
         
         itlag = nint(config % src_tlag / (delta*1./real(f)) ) 
         
c         print*,config % src_tlag,itlag,delta,config % src_tlag
         
         
         call zrotate2d(xlag,ylag,n,np,
     >                   (config % src_fast-fast),xtemp,ytemp)
         xlag(1:n) = xtemp(1:n) ; ylag(1:n) = ytemp(1:n)
         call zlag(xlag,ylag,n,np,itlag,iwextra,xtemp,ytemp,noverlap)
         xlag(1:n) = xtemp(1:n) ; ylag(1:n) = ytemp(1:n)

         call zrotate2d(xlag,ylag,n,np,
     >                  -(config % src_fast-fast),xtemp,ytemp)
         xlag(1:n) = xtemp(1:n) ; ylag(1:n) = ytemp(1:n)
         
         ! new itlag 
         
c        lag = nint(config % src_lag / delta) 
c        call zrotate2d(xlag,ylag,n,np,
c    >                   config % src_fast,xtemp,ytemp)
c        xlag(1:n) = xtemp(1:n) ; ylag(1:n) = ytemp(1:n)
c        call zlag(xlag,ylag,n,np,lag,iwextra,xtemp,ylag,noverlap)
c        xlag(1:n) = xtemp(1:n) ; ylag(1:n) = ytemp(1:n)
      endif 
      
      
      
C  **       
      call zcovariance(xlag,ylag,noverlap,np,cov)
      call zeigen2x2(cov,lambda1,lambda2,vec1,vec2)
      call zsourcepol(fast,lambda1,lambda2,vec1,vec2,spol,dspol)

c  ** calc the number of degrees of freedom **
c  ** first rotate into spol-fast (so y is signal and x is noise) **
      call zrotate2d(xlag,ylag,noverlap,np,spol-fast,xnoise,ynoise)
      call zndf(xnoise,noverlap,norig,ndf)
c  ** estimate signal to noise ratio *
      call calcsnr(ynoise,xnoise,noverlap,snr)
c  ** normalise error surface and calc errors in fast and lag**
      call zerror95(error_int,ndf,lambda2_min,idfast,idtlag)
      dtlag = delta * idtlag * itlag_step / real(f)
      dfast = 180.  * idfast / real(np1-1)
      
      return
      
      end
C=======================================================================

C=======================================================================
      subroutine zgrid_lambda2(x,y,n,iwextra,itlag_step,
     >                         lambda1grid,lambda2grid,delta)
C=======================================================================
c
c      calculate the second eigenvalue of the particle motion covaraince 
c      matrix over fast direction from -90 - 90 deg (1deg grid spacing) 
c      and lags of 0 to 40 (1 sample point grid spacing).
c
c      variables
c    in:
c      x(np)                  real            time series (local east component)
c      y(np)                  real            time series (local north component)
c      n                  int            number of points
c      iwextra            int            number of extra points included in window
c                                    to allow for the lagging
c      itlag_step            int            gridding step in for lag time
c      [np                  int            array dimension
c                              (read from SIZE_np.h at compile time)]
c      np1/2                  int            array dimension
c    out:
c      lambda2grid(np1,np2)      real      gridded value of lambda2
c    local:
c      fast                  real            fast direction in degrees
c      lag                  int            lag time in sample spacings
c
c     MODIFIED BY J.WOOKEY FOR SHEBA ...
c
c-----------------------------------------------------------------------
c      N. Teanby      4-8-02      Original code
c-----------------------------------------------------------------------
      use sheba_config ! use the sheba_config module
      use array_sizes ! use the array_sizes modules
C-----------------------------------------------------------------------
      implicit none
      integer n,i,j,noverlap,k
      real x(np),y(np),lambda2grid(np1,np2),lambda1grid(np1,np2)
      real fast
      integer lag,iwextra,itlag_step,lag1,lag2
      real cov(2,2),lambda1,lambda2,vec1(2),vec2(2)
      real xrot(np),yrot(np),xlag1(np),ylag1(np)
      real xlag2(np),ylag2(np),xlag(np),ylag(np)
      real xpol(np),ypol(np),spol
c  ** temporary arrays
      real xtemp(np),ytemp(np)
      real delta

C     ** Hilbert transformed versions
      real hx(np),hy(np)
      real hxrot(np),hyrot(np),hxlag1(np),hylag1(np)

c  ** initialise all arrays to zero
      do i=1,np
         xrot(i) = 0.0
         yrot(i) = 0.0
         xlag1(i) = 0.0
         ylag1(i) = 0.0
         xlag2(i) = 0.0
         ylag2(i) = 0.0
         xlag(i) = 0.0
         ylag(i) = 0.0
         xpol(i) = 0.0
         ypol(i) = 0.0
         hx(i) = 0.0
         hy(i) = 0.0
         hxrot(i) = 0.0
         hyrot(i) = 0.0
         hxlag1(i) = 0.0
         hylag1(i) = 0.0
      enddo ! i=1,np

c  ** calculate the Hilbert transforms
      if (config % iscs_corr == 1) call hilbert(x,n,hx)
      if (config % iscs_corr == 1) call hilbert(y,n,hy)

c  ** map out the lambda2 surface **
      do 1 i=1,np1
c         ** set fast direction (range is -90 to 90deg) **
         fast = -90. + 180.*real(i-1)/real(np1-1)

c         ** rotate, lag, calc the covariance and eigenvalues **
         call zrotate2d(x,y,n,np,fast,xrot,yrot)
         if (config % iscs_corr == 1) 
     >       call zrotate2d(hx,hy,n,np,fast,hxrot,hyrot)
         
         do 2 j = 1,np2

c         ** JW do lag in two stages for ScS-correction
            lag = (j - 1)*itlag_step
            lag1 = nint(real(lag)/2.)
            lag2 = lag-lag1
            call zlag(xrot,yrot,n,np,lag1,iwextra,xlag1,ylag1,noverlap)
            if (config % iscs_corr == 1) 
     >         call zlag(hxrot,hyrot,n,np,lag1,iwextra,
     >                   hxlag1,hylag1,noverlap)

C        ** ScS correction
            call scscorr(xlag1,ylag1,hxlag1,hylag1,noverlap,fast,j)
C        ** don't need the Hilbert traces any more ...            

            call zlag(xlag1,ylag1,n,np,lag2,iwextra,
     >                xlag2,ylag2,noverlap)
            
            xlag(1:n) = xlag2(1:n) ; ylag(1:n) = ylag2(1:n)

C        ** apply a source correction (if specified)
            if (config % i_src_corr == 1) then
               lag = nint(config % src_tlag / delta) 
               call zrotate2d(xlag,ylag,n,np,
     >                         config % src_fast-fast,xtemp,ytemp)
               xlag(1:n) = xtemp(1:n) ; ylag(1:n) = ytemp(1:n)
               call zlag(xlag,ylag,n,np,lag,iwextra,
     >                   xtemp,ytemp,noverlap)
               xlag(1:n) = xtemp(1:n) ; ylag(1:n) = ytemp(1:n)
               call zrotate2d(xlag,ylag,n,np,
     >                         -(config % src_fast-fast),xtemp,ytemp)
               xlag(1:n) = xtemp(1:n) ; ylag(1:n) = ytemp(1:n)
            endif 

C           ** calculate the surface point                           
            if (config % imode == 0) then
C           ** TRANSVERSE ENERGY MINIMISATION **
               spol = config % source_pol
             call zrotate2d(xlag,ylag,noverlap,np,(spol-fast),xpol,ypol)
               lambda2grid(i,j)=0.0
               do k=1,n
                  lambda2grid(i,j)=lambda2grid(i,j)+abs(xpol(k))
               enddo                
C           ** do covariance calc. for eigenvalue ratio
               call zcovariance(xlag,ylag,noverlap,np,cov)
               call zeigen2x2(cov,lambda1,lambda2,vec1,vec2)
C           ** SECOND EIGENVALUE MINIMISATION **
            elseif (config % imode == 1) then
               call zcovariance(xlag,ylag,noverlap,np,cov)
               call zeigen2x2(cov,lambda1,lambda2,vec1,vec2)
               lambda2grid(i,j)=lambda2
               lambda1grid(i,j)=lambda1
               
            endif
2         continue ! j = 1,np2
1      continue ! i=1,np1

      return
      end
C=======================================================================
