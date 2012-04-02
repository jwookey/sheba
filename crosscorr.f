!===============================================================================
!  S H E B A - Shear-wave Birefringence Analysis
!===============================================================================
!  This software is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!===============================================================================
!
!  James Wookey, School of Earth Sciences, University of Bristol
!  CVS: $Revision: 1.1 $ $Date: 2011/02/11 16:09:41 $
!

	subroutine crosscorr(s1, s2, n, ss, norm)
c	routine to compute the cross correlation function
c	between two time series of length n.  It is normalized
c	by sqrt(acf1(0)*acf2(0)) if norm  = 1, or not if norm = 0.
c	It is normalized by default.
c	acf1(0) and acf2(0) are stored in common
c  	The result ss is of length 2n + 1.  The zero-lag point
c	is in position n + 1. ss(1) and ss(2n + 1) are set to zero.
c	This is a test of shearsp by doing the cc in time domain
c	Assumes the signal is zero outside the boundaries of the array
c	
	dimension s1(n),s2(n),ss(2*n + 1)
	real acfs1,acfs2
	integer norm
c	data norm/1/
c       *************************************************************
	
	
	n2     = 2*n
	indx   = 1
	ss(1)  = 0.
	ss(n2) = 0.
	
	do lag = -(n-1), n-1, 1
	  indx     = indx + 1
  	  ss(indx) = 0.
c	  find overlap
	  iover = n-iabs(lag)
	  if(lag.le.0) then
	    do i = 1,iover
	      ss(indx)=ss(indx) + s1(i) * s2(i-lag)
	    enddo
	  else
	    do i = 1,iover
	      ss(indx)=ss(indx) + s1(i + lag) * s2(i)
	    enddo
	  endif
	enddo
	
	
	acfs1 = 0.
	acfs2 = 0.
	do i = 1,n	
	  acfs1 = acfs1 + s1(i)**2
	  acfs2 = acfs2 + s2(i)**2
	enddo
	fac  =  sqrt(acfs1 * acfs2)
	
	
	
c	normalize cc
	if (norm.eq.1) then
c	  write(*,*) 'cross correlation normalized'
	  do i = 1,n2	
	     ss(i) = ss(i) / fac 
	  enddo	  
	else
c	  write(*,*) 'cross correlation not normalized'
	endif
	
	
	return
	end


	
c-----------------------------------------------------------------------
      subroutine zerror_interp_xc(error,error_int)
c-----------------------------------------------------------------------
c
c      interpolate the error surface in the tlag direction.
c      the interpolation is done one row at a time
c
c      variables
c    in:
c      error(np1,np2)      real            error surface (i.e. lambda2)
c      np1/2                  int            array dimensions
c      np2int            int            np2 after interpolation
c    out:
c      error_int(np1,np2int)      real      interpolated error surface
c    local:
c      np                  int            array dimension
c                              (read from SIZE_np.h at compile time)
c
c
c-----------------------------------------------------------------------
c      N. Teanby      4-8-02      Original code
c-----------------------------------------------------------------------
      use array_sizes ! use the array_sizes modules
C-----------------------------------------------------------------------
      implicit none
      real error(np1,np2XC),error_int(np1,np2XCint)
      real error_row(np),error_row_int(np)
      integer f,i,j,n_int

c  ** check that np is big enough **
      if (np.lt.2*np2int) then
         pause 'ERROR: zerror_interp: np not big enough'
      endif

c  ** calc the interpolation factor based on np2 and np2int **
      f = (np2XCint-1)/(np2XC-1)
c  ** f needs to be an integer for zsplint to work **
      if (mod(np2XCint-1,np2XC-1).ne.0) then
         pause 'ERROR: zerror_interp: f not a whole number'
      endif

c  ** interpolate error surface in tlag direction **
c  ** do interpolation one row at a time **
      do 1 i=1,np1
c      ** copy ech row to a dummy array **
!         do 11 j=1,2*np2+1
          do 11 j=1,np2XC
            error_row(j)=error(i,j)
11         continue
c      ** interpolate the row data **
         call zsplint(error_row,np2XC,f,error_row_int,n_int)
c     ** check that n_int = np2int
c         (this should not be possible but check anyway)**
         if (np2XCint.ne.n_int) then
            pause 'ERROR: zerror_interp: np2int.ne.n_int'
         endif
c      ** copy interpolated row to output array **
         do 22 j=1,n_int
            error_int(i,j)=error_row_int(j)
22         continue
1      continue

      return
      end
c-----------------------------------------------------------------------
      subroutine zerror_max(error,np1,np2XC,ifast,itlag,xc_max)
c-----------------------------------------------------------------------
c
c      grid search for the maximum value of correlation on the interpolated
c      error surface
c
c      variables
c    in:
c      error(np1,np2)      real      error surface (i.e. lambda2)
c      np1/2                  int      array dimensions
c    out:
c      xc_max            real      maximum value of xc_max
c      itlag                  int      index of lag corresponding to xc_max
c      ifast                  int      index of fast dirn corresponding to xc_max
c
c-----------------------------------------------------------------------
c      N. Teanby      4-8-02      Original code
c-----------------------------------------------------------------------

      implicit none
      integer np1,np2XC
      real error(np1,np2XC),xc_max
      integer i,j,itlag,ifast

c  ** find the minimum lambda2 position **
      xc_max=error(1,1)
      itlag=1
      ifast=1
      do 1 i=1,np1
        do 2 j = 1,np2XC
           if (error(i,j).gt.xc_max) then
              xc_max = error(i,j)
              itlag = j
              ifast = i
           endif
2         continue
1      continue
c       print *,np1,ifast,np2,itlag
      return
      end
	  

c-----------------------------------------------------------------------
      subroutine zerror95XC(error,ndf,lambda2_min,ierror,jerror)
c-----------------------------------------------------------------------
c
c      subroutine to calculate normalise contours so that 1=95% confidence
c      contour. Also calc the errorbars on tlag and fast direction
c      CrossCorrelation needs the Fisher-Transform of the correltion map 
c      see Wuestefeld & Bokelmann BSSA, 2007
c
      use array_sizes ! use the array_sizes modules
C-----------------------------------------------------------------------

      implicit none
      integer ndf,counter
      integer i,j,jmin,jmax,irange,jrange
      integer line(npc)
      real ierror,jerror
      integer k,k1,irange_min,irange_max,istart,line_test(npc)
      real error(np1,np2XCint),lambda2_min,fftable,lambda_crit,Z_crit
	  real Z
      external fftable

c  ** check that npc is big enough **
      if (npc.lt.np1) then
         pause 'ERROR: zerror95: npc < np1'
      endif
      
c     DO THE FISHER TRANSFORM
      Z = atanh(lambda2_min)
c  ** calc value of lambda at 95% confidence limit from tabulated values **
      if (ndf.ge.3) then
		 Z_crit = Z + ( -2.*Z*sign(1.,Z) / (ndf-2) * fftable(ndf-2))
     >                  * sqrt(1./real(ndf-3));
      else
         print*,'WARNING: zerror95: ndf <=2, set to 3)'
         Z_crit= Z*( 1. + 2.*fftable(1))
      endif
c	  BACKTRANSFORM

c      print *, lambda_max,Z,lambda2_min
	  lambda_crit = tanh(Z_crit)
	  
c      print *, lambda_max,lambda2_min
c  ** normalise errors by lambda_max, so that and error of 1 = 95% confidence **
      counter=0
	  ierror=error(1,1)
	  jerror=error(1,1)
	  do i=1,np1
         do j=1,np2XCint 
   		    error(i,j)=error(i,j)/lambda_crit
		    if (error(i,j).le.ierror) then
			   ierror=error(i,j)
			endif
		enddo
      enddo
	  
	  
c  ** find i/j error (half width of 95% confidence contour) **
c  ** find min and max j, simply search the array **
      jmin=np2int
      jmax=1
      do i=1,np1
         do j=1,np2XCint
            if (error(i,j).ge.1.0) then
               jmin = min0(jmin,j)
               jmax = max0(jmax,j)
            endif
	      enddo
      enddo
	  
      jrange=jmax-jmin
c  ** finding min max i is more difficult because of cyclicity of angles **
c  ** sweep a line over all j, set point on line equal to 1 if it falls within the 95% convidence contour for any j. The height of the bounding rectangle is defined by the shortest line which includes all points with a value of line(i)=1. This line is found by searching all line lengths from the minimum = sum_i linr(i) to maximum = np1**
      do i=1,np1
         line(i)=0
      enddo
      do j=1,np2XCint
         do i=1,np1
            if (error(i,j).ge.1.0) then
               line(i)=1
            endif
         enddo
      enddo
c  ** min line length **
      irange_min = 0
      do i=1,np1
         irange_min = irange_min + line(i)
      enddo
c  ** max line length **
      irange_max = np1
c  ** search all line length and starting points to find irange **
      do i = irange_min , irange_max
        do istart = 1,np1
          do k = 1,np1
            line_test(k)=0
          enddo
          do k = istart,istart+i
            if (k.gt.np1) then
               k1 = k - np1
            else
               k1 = k
            endif
            line_test(k1) = 1
          enddo
          do k = 1,np1
            if ((line(k).eq.1).and.(line_test(k).ne.1)) then
              goto 1
            endif
          enddo
          irange = i
          goto 11
1          continue
        enddo
      enddo
11      continue

c  ** one standard deviation = 0.5*95% error half width
c      (so x full width of 95% contour by 0.25 to get 1s.d.)**
      ierror = 0.25*real(irange)
      jerror = 0.25*real(jrange)

      return
      end
  
	  
	  
c-----------------------------------------------------------------------
      real function NULLCRITERION(fastEV,fastXC,tlagEV,tlagXC,
     >     dtlagNORM)
c-----------------------------------------------------------------------	  
c     ** Now doing Quality and Null Check by comparign to methods:
c     ** see Wuestefeld et al, BSSA, 2008
      real Omega, Rat, m1, m2, dis1, dis2
	  
      m1 = max(mod(fastEV+180., 180.), mod(fastXC+180., 180.))
      m2 = min(mod(fastEV+180., 180.), mod(fastXC+180., 180.))
      
      Omega = mod (m1 - m2, 90.)
      if (Omega .gt.45) then
         Omega = 90 - Omega
      endif
	  
c		Now Normalise	  
	  Omega = Omega / 45
	  Rat   = tlagXC / tlagEV
	  
c		normalised Distance:
c       Quality of nulls is the distance from point 0/1
c       Quality of non-nulls is the distance from point 1/0 
	  dis1 = sqrt( Rat**2     + (Omega-1)**2 ) * sqrt(2.)
	  dis2 = sqrt( (Rat-1)**2 + ( Omega )**2 ) * sqrt(2.)
	  


		 
	  if (dis1<dis2) then
		if (dis1>1) then
			dis1=1.
		endif
		NULLCRITERION = -1.*(1.-dis1)
	  else
	  	if (dis2>1) then
	        dis2=1.
	    endif
	    NULLCRITERION =  1.*(1.-dis2)
	  endif
	  
	  
c     ** Now make sure that delay time close to itlag_scale are POOR
       if (Rat.gt.0.9 .AND. tlagEV/dtlagNORM.gt.0.9) then	  
	     NULLCRITERION=NULLCRITERION*(-10*tlagEV/dtlagNORM+10)
       endif
	



	  
	  
	  
c      if (Omega<8.0 .AND. 0.8<=Rat .AND. Rat<=1.1) then
c         if (dfast<8.0 .AND. dtlagNORM<0.15) then
c		     NULLCRITERION = 2
c		 else
c		     NULLCRITERION = 1
c		 endif
c      elseif (Omega<15.0.AND.0.7<=Rat.AND.Rat<=1.2) then 
c    	  if (dfast<15.0 .AND. dtlagNORM<0.33) then
c		     NULLCRITERION = 1
c		 else
c		     NULLCRITERION = 0
c		 endif
c      elseif (32.0<=Omega.AND.Omega<=58.0.AND.Rat<=0.2) then
c         NULLCRITERION = -2
c      elseif (37.0<=Omega.AND.Omega<=53.0.AND.Rat<=0.3) then
c         NULLCRITERION = -1
c      else
c         NULLCRITERION = 0
c      endif	  
      
    
      return
      end
	  
	  
	  
	 
