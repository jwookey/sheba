!-----------------------------------------------------------------------
      subroutine zndf(y,n,norig,ndf)
!-----------------------------------------------------------------------
!
!      subroutine to calculate number of degrees of freedom (ndf) of a time
!      series y, with length n, and physical dimension np. If y has been
!      interpolated then norig is the original number of points in the shear
!      wave window
!
!      calculated according to appendix in Silver and Chan 1991
!
!      For a gaussian white noise ndf should be equal to n.
!      The true ndf should be less than n
!
!      variables
!    in:
!      y(np)                        real      time series
!      n                        int      number of points
!      norig                        int      original (uninterpolated) number of points
!                                    in the shear wave window
!    out:
!      ndf                        int      number of degrees of freedom
!    local:
!      np                        int      array dimension
!                                    (read from SIZE_np.h at compile time)
!      n2                        int      next power of 2 up from n
!      yint(np)                  real      y interpolated onto n2 points
!      yfft(np*2)                   real      fft of interpolated (to n2) y in NR format
!      yfft_amp(np)             real      amplitude spectum
!      nf                        int      number of freq points for n data points
!      yfft_amp                  real      interpolated amp spec on nf points
!      f2,f4                        real      f2 and f4 from Silver and Chan 1991
!                                    used to calc ndf
!
!-----------------------------------------------------------------------
!      N. Teanby      5-7-02   Original code, based on ndf_fun2.f
!                              and ndf_spect.f in DTM
!      J. Wookey      7-4-14   Modified ndf calculation following
!                              Walsh, Arnold and Savage (JGR, 2013)
!-----------------------------------------------------------------------
      use array_sizes ! use the array_sizes modules
!-----------------------------------------------------------------------
      implicit none
      integer ndf
      integer n,n2,nf,norig
      real y(np),yint(np),yfft(np*2),yfft_amp(np)
   	real sce2,sce4
      integer scndf
   	real e2,e4
      integer i

!  -- calc the fft of y --
!  ** calc n2, the next power of 2 from n (or n if n=power of 2) **
      n2 = 2**(int(log(real(n)-0.1)/log(2.0)) + 1)
!  ** interpolate to have n2 points **
      call zlinint(y,n,np,n2,yint)
!  ** perform fft (1 means do forward fft)**
      call zfft(yint,n2,np,1,yfft)
!  ** determine number of frequency points for original time series.
      nf=(norig+2)/2
!  ** calculate amplitude spectrum **
      do i=1,nf
         yfft_amp(i) = sqrt(yfft(2*i-1)**2 + yfft(2*i)**2)
      enddo

   	e2=0.
   	e4=0.
   	do i=1,nf
   !!          yfft_amp(i)**2 is |f_i|^2|g_i|^2
   !!          note that have used i as the index rather than n
   	   e2   = e2 + yfft_amp(i)**2
   !!          yfft_amp(i)**4 is |f_i|^4|g_i|^4
   !!          note that have used i as the index rather than n
   	   e4   = e4 + yfft_amp(i)**4
	      
   	   if(i.eq.1.or.i.eq.nf)then
   !!            coefficients from SC equation A3
   !!            note adding 1/2 is the same as adding 1 then subtracting 1/2	  
   	      sce2 = e2 - 0.5*yfft_amp(i)**2
   !!            coefficients from SC original code (see Walsh eq (33))
   !!            note adding 1/2 is the same as adding 1 then subtracting 1/2	    
   	      sce4 = e4 - 0.5*yfft_amp(i)**4
   	   endif
      enddo
      e2 = sce2


!!    Walsh et al eq 33 shows that the coefficients are
!!    1/3 for n=0,N and 4/3 otherwise
!!    e4 has coefficients 1 for all n so we multiply
!!    e4 by 4/3. The end points (n=0,N) now have coefficients
!!    4/3 so we must subtract off 1 to get 1/3      
      e4 = 4.*e4/3. - yfft_amp(1)**4 - yfft_amp(nf)**4 
   	ndf = nint( 2.0 * (2.0*e2**2/e4 - 1.) )
   	
      ! original version
      scndf = nint( 2.0 * (2.0*sce2**2/sce4 - 1.) )
      
!      print *, "e2,e4,sce2,sce4,e2,e4,ndf,scndf are "
!      print *, e2,e4,sce2,sce4,e2,e4,ndf,scndf   

      return
      end

! MFAST VERSION
! c-----------------------------------------------------------------------
! 	subroutine zndf(y,n,norig,ndf)
! c-----------------------------------------------------------------------
! c
! c	subroutine to calculate number of degrees of freedom (ndf) of a time
! c	series y, with length n, and physical dimension np. If y has been 
! c	interpolated then norig is the original number of points in the shear
! c	wave window
! c
! c	calculated according to appendix in Silver and Chan 1991
! c    Modifed by Walsh, Arnold and Savage 2013. (Silver and Chan Revisited,
! c     submitted to J. Geophysical Research June 18, 2013.
! c
! c	For a gaussian white noise ndf should be equal to n.
! c	The true ndf should be less than n
! c
! c	variables
! c    in:
! c	y(np)				real	time series
! c	n				int	number of points
! c	norig				int	original (uninterpolated) number of points
! c						in the shear wave window
! c    out:
! c	ndf				int	number of degrees of freedom
! c    local:
! c	np				int	array dimension
! c						(read from SIZE_np.h at compile time)
! c	n2				int	next power of 2 up from n
! c	yint(np)			real	y interpolated onto n2 points
! c	yfft(np*2) 			real	fft of interpolated (to n2) y in NR format
! c	yfft_amp(np) 		real	amplitude spectum
! c	nf				int	number of freq points for n data points
! c	yfft_amp			real	interpolated amp spec on nf points
! c	sce2,sce4				real	f2 and f4 from Silver and Chan 1991
! c       scndf,ndf				int	number off degrees of freedom from Silver and Chan--i.e., previous codes
! c						used to calc ndf
! c       e2, e4				real	e2 and e4 from Walsh et al. 2013
! 
! c
! c-----------------------------------------------------------------------
! c	N. Teanby	5-7-02	Original code, based on ndf_fun2.f
! c					and ndf_spect.f in DTM
! c
! c   M. Savage 15 June 2013  Fixing to use modified f2 and f4 from Walsh et al.
! c
! c       Walsh, E., R. Arnold and M. K. Savage (2013) "Silver and Chan 
! c      Revisited), submitted to JGR June 2013.
! c-----------------------------------------------------------------------
! 
! 	implicit none
! 	integer ndf
! 	integer n,np,n2,nf,norig
! 	include "SIZE_np.h"
! 	real y(np),yint(np),yfft(np*2),yfft_amp(np)
! !	real f2,f4
! 	real sce2,sce4
!         integer scndf
! 	real e2,e4
! 	integer i
! 	
! c  -- calc the fft of y --
! c  ** calc n2, the next power of 2 from n (or n if n=power of 2) **
! 	n2 = 2**(int(log(real(n)-0.1)/log(2.0)) + 1)
! c  ** interpolate to have n2 points **
! 	call zlinint(y,n,np,n2,yint)
! c  ** perform fft (1 means do forward fft)**
! 	call zfft(yint,n2,np,1,yfft)
! c  ** determine number of frequency points for original time series.
! 	nf=(norig+2)/2
! c  ** calculate amplitude spectrum **
! 	do 2 i=1,nf
! 	   yfft_amp(i) = sqrt(yfft(2*i-1)**2 + yfft(2*i)**2)
! 2	continue
! 
! !! below is the original code based on the SC codes
! !! they refer to f2 and f4 but these are actually
! !! e2 and e4
! !! SC note that the noise is a convolution E \approx \sum_{n=0}^N a_n x f_n*g_n 
! !! where a is 1/2 for n=0,N and 1 otherwise(SC equation A3)
! !! we cannot access f but we can estimate it by calculating E
! !! from the data
! !!c  -- calc f2 and f4 from Silver and Chan 1991 --
! !! which are actually e2 and e4
! !!	f2=0.
! !!	f4=0.
! !!	do 7 i=1,nf
! !!        yfft_amp(i)**2 is |f_i|^2|g_i|^2
! !!        note that have used i as the index rather than n
! !!	  f2   = f2 + yfft_amp(i)**2
! !!        yfft_amp(i)**4 is |f_i|^4|g_i|^4
! !!        note that have used i as the index rather than n
! !!	  f4   = f4 + yfft_amp(i)**4
! !!cc	  for zero frequency and for Nyquist:
! !!c	  if(i.eq.1.or.i.eq.nf)then
! !!          coefficients from SC equation A3
! !!          note adding 1/2 is the same as adding 1 then subtracting 1/2
! !!c	    f2 = f2 - 0.5*yfft_amp(i)**2
! !!          coefficients from SC original code (see Walsh eq (33))
! !!          note adding 1/2 is the same as adding 1 then subtracting 1/2
! !!c	    f4 = f4 - 0.5*yfft_amp(i)**4
! !!c	  endif
! 
! 
! !!  testing SC coefficients to compare  (equation 33 Walsh et al.)
! !!  relabel coefficients to make it less confusing
! 
! 	e2=0.
! 	e4=0.
! 	do 7 i=1,nf
! !!        yfft_amp(i)**2 is |f_i|^2|g_i|^2
! !!        note that have used i as the index rather than n
! 	  e2   = e2 + yfft_amp(i)**2
! !!        yfft_amp(i)**4 is |f_i|^4|g_i|^4
! !!        note that have used i as the index rather than n
! 	  e4   = e4 + yfft_amp(i)**4
! 	  
! 	  if(i.eq.1.or.i.eq.nf)then
! !!          coefficients from SC equation A3
! !!          note adding 1/2 is the same as adding 1 then subtracting 1/2	  
! 	    sce2 = e2 - 0.5*yfft_amp(i)**2
! !!          coefficients from SC original code (see Walsh eq (33))
! !!          note adding 1/2 is the same as adding 1 then subtracting 1/2	    
! 	    sce4 = e4 - 0.5*yfft_amp(i)**4
! 	  endif
! 7 	continue
! c
! !!  put in Walsh' coefficients
! !!  Replace sce2 and sce4 with e2 and e4 from Walsh et al. (eqn 25 and 27)
! !!  Below eqn 26 it notes that the bn coefficients in eq 25 are equal to the an
! !!  coefficients
! !!  Note that the first equation in section 2.5 (which happens to be un-numbered) 
! !!  mentions that the an coefficients are 1/2 for n=0,N and 1 otherwise
! !!  so the sce2 and e2 are the same
! !!  Walsh et al eq 25
!         e2 = sce2
! !!  below is the main point of different between the two degrees of freedom
! !!  compare Walsh et al (eq 27,28,33)
! !!  note that eq 27 = e4 and eq 33=sce4
! !!  eq 28 which is what is in the SC article has only ever been implemented
! !!  in the original SC codes and later scrapped because eq 33 gave superior
! !!  performance (see Walsh et al section 2.5 and Figure 3)
! !!  you should find that ndf \approx scndf * 3/4 roughly
! 
! !!  Walsh et al eq 33 shows that the coefficients are
! !!  1/3 for n=0,N and 4/3 otherwise
! !!  e4 has coefficients 1 for all n so we multiply
! !!  e4 by 4/3. The end points (n=0,N) now have coefficients
! !!  4/3 so we must subtract off 1 to get 1/3
!         e4 = 4.*e4/3. - yfft_amp(1)**4 - yfft_amp(nf)**4 
! 
! !! estimate the degrees of freedom using (Walsh et al 31 aka SC A12)
! !! based on our defintion of e4 (Walsh et al eq 27)
! 	ndf = nint( 2.0 * (2.0*e2**2/e4 - 1.) )
! !! estimate the degrees of freedom using (Walsh et al 31 aka SC A12)
! !! based on our the original SC code definition of e4 (Walsh et al eq 33)
! 	scndf = nint( 2.0 * (2.0*sce2**2/sce4 - 1.) )
! c  MKS 15 June 2013 replacing with Walsh method (Equation 31)
!         print *, "e2,e4,sce2,sce4,e2,e4,ndf,scndf are "
!         print *, e2,e4,sce2,sce4,e2,e4,ndf,scndf
! 
!         open (unit=59,file="df.txt",action="write"
!      x  ,form="formatted")
! c        open (unit=59,file="df.txt",action="write",status="new"
! c     x  ,form="formatted")
!         write (59,*) e2,e4,sce2,sce4,e2,e4,ndf,scndf
!         close (59)
! c
! 	return
! 	end

      
      
      