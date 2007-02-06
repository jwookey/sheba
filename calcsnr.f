!=======================================================================
!  S H E B A - Shear-wave Birefringence Analysis
!=======================================================================
!  This software is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!=======================================================================
!
!  James Wookey, School of Earth Sciences, University of Bristol
!  CVS: $Revision: 1.2 $ $Date: 2007/02/06 14:03:49 $
!


c-----------------------------------------------------------------------
      subroutine calcsnr(xsig,xnoise,n,snr)
c-----------------------------------------------------------------------
c
c      Estimate signal-to-noise, based on definition of Restivo and 
c      Helffrich (GJI, 1999). Which is the 
c  
c      peak amplitude on 'radial' / 2*std(transverse) 
c
c-----------------------------------------------------------------------
      use array_sizes ! use the array_sizes modules
C-----------------------------------------------------------------------
      implicit none
      integer ndf
      integer norig
      real xsig(np),xnoise(np)
      real snr
      real totsq,signal,noise
      
      integer i,n
      
      totsq=0.0
      signal = 0.0
      
      do i=1,n
         totsq=totsq+xnoise(i)**2.0
         if (abs(xsig(i)).gt.signal) signal = abs(xsig(i))
      enddo
      
      noise = 2.0 * sqrt(totsq/real(n))
      
      snr = signal / noise
      

      return
      end
