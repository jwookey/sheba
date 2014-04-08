!===============================================================================
!  S H E B A - Shear-wave Birefringence Analysis
!===============================================================================
!  This software is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!===============================================================================
!
!  James Wookey, School of Earth Sciences, University of Bristol
!

!=======================================================================
      subroutine desplit(h1,h2,phi,dt)
!=======================================================================
!
!     Remove shear-wave splitting
!
!     parameters:
!  
!     h1,h2 :  (I/O) (SACtrace) input horizontal, orthogonal components
!     dt,phi:  (I)   (real) lag-time, fast direction to correct by
!
!-----------------------------------------------------------------------
      use f90sac  ! use the f90sac module
      use sheba_config ! use the sheba_config module
!-----------------------------------------------------------------------
      implicit none
      
      type (SACtrace) :: h1, h2
      real dt , phi
      real cmpazdiff,rotang

      write(iulog,*) 'Correcting input by (phi,dt)',phi,dt

!  ** calculate the angle of rotation, based on phi and the azimuth of h1
      rotang = (phi - h1 % cmpaz)
   
!  ** rotate the traces      
      call f90sac_rotate2d(h1,h2,rotang)      

!  ** time shift the slow trace
      call f90sac_tshift(h2,-dt)

!  ** rotate back to original orientation    
      call f90sac_rotate2d(h1,h2,-rotang)      

!  * done      
      return
      end
!=======================================================================


