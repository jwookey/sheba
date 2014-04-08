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
      subroutine enz2abc(h1,h2,v,slw)
!=======================================================================
!  
!     Subroutine to rotate the radial and vertical traces into ray-based
!     reference frame based on the slowness. Slowness is converted to
!     angle using IASP91 velocities. Input traces are assumed to be:
!     h1 = N , h2 = E and v = Z
!
!-----------------------------------------------------------------------
      use f90sac  ! use the f90sac module
      use sheba_config ! use the sheba_config module
      use event_info ! use event_info module
!-----------------------------------------------------------------------
      implicit none
      type (SACTrace) :: h1,h2,v ! the re-ordered traces for analysis
      real :: phi, theta
      real :: slw
      real :: svel ! AK135
      
!  ** first rotate horizontals into radial-transverse direction
      phi = h1 % baz - 180.0 ! radial bearing
      
      svel=3.46

      slw = slw / 111.16
      theta = asin(svel*slw) * (180./pi)
      
      call f90sac_unwind(phi)
      call f90sac_rotate2d(h1,h2,phi)
      call f90sac_rotate2d_rz(h1,v,theta)

!  ** re-set the components header
      h1 % cmpinc = 90.0
      v % cmpinc = 0.0

!  ** rotate 'horizontals' back to E-N
      call f90sac_rotate2d(h1,h2,-phi)

      return
      end subroutine enz2abc
!=======================================================================
