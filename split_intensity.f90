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
      subroutine calc_split_intensity(h1,h2,spol, wbeg, wend, SI)
!=======================================================================
!
!     Remove shear-wave splitting
!
!     parameters:
!  
!     h1,h2 :  (I) (SACtrace) input horizontal, orthogonal components
!     wbeg,wend : analysis window
!     SI : split intensity
!
!-----------------------------------------------------------------------
      use f90sac  ! use the f90sac module
      use sheba_config ! use the sheba_config module
!-----------------------------------------------------------------------
      implicit none
      
      type (SACtrace) :: h1, h2, wh1, wh2
      real*4 wbeg,wend,spol,r2
      real*4 SI
      !real*4,allocatable :: R(:,:),RT(:,:),T(:,:)

      integer i,j

      call f90sac_window(h1,wh1,wbeg,wend)
      call f90sac_window(h2,wh2,wbeg,wend)

!  ** rotate the traces to spol, h1 contains radial, h2 transverse      
      call f90sac_rotate2d(wh1,wh2,spol)       

!  ** take the time derivative of the radial trace
      call f90sac_time_derivative(wh1)

!      call f90sac_writetrace('test_gr.sac',wh1)

!  ** calculate the projection
      r2 = 0 
      SI = 0 
      do i = 1,wh1%npts
         r2 = r2 + wh1%trace(i)*wh1%trace(i)
         SI = SI + wh1%trace(i)*wh2%trace(i)
      enddo
!     
      SI = -2.*SI/r2


!  * done      
      return
      end
!=======================================================================


