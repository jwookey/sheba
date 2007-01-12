!===============================================================================
!     S H E B A
!===============================================================================
!     Fortran 90/95 Source Code File
!-------------------------------------------------------------------------------
!
!     PROGRAM : sheba
!     FILE    : misc.f
!     AUTHOR  : James Wookey
!     PLACE   : School of Earth Sciences, University of Leeds
!     DATE    : December 2003
!     PURPOSE : Various subroutines for sheba
!     VERSION : 1.0
!     COMPLETE: No
!     COMMENTS: 
!
!-------------------------------------------------------------------------------
!     This software is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!-------------------------------------------------------------------------------
!
!      SHEBA is a prototype code to implement various algorithms for
!      analysing shear-wave splitting. 
!
!-------------------------------------------------------------------------------
!     Changes log
!-------------------------------------------------------------------------------
!     2003-12-04     * Incept date

!===============================================================================
      subroutine sheba_banner(version)
!===============================================================================
!  
!     Print sheba startup information
!
      implicit none
      real :: version
      
      write(*,'("")')
      write(*,'(a,a)') '========================================', & 
                       '========================================'
      write(*,'(a,a)') '-               S H E B A  -  Shear-wave', & 
                       ' Birefringence Analysis                -'
     
      write(*,'(a,a)') '----------------------------------------', & 
                       '----------------------------------------'
      write(*,'(a,f4.2,a)') '-                       Version ',version, & 
      ' - James Wookey 2004                       -'
      write(*,'(a,a)') '========================================', & 
                       '========================================'
      write(*,'("")')
     
      end subroutine sheba_banner
!===============================================================================

!===============================================================================
      subroutine unwind(angle)
!===============================================================================
!  
!     'Unwind' an angle to between 0 to 360
!
      implicit none
      real :: angle
      integer i
      
      do i = 1 , 1000
         if (angle >= 0.0 .and. angle < 360.0) exit
         if (angle < 0.0) angle = angle + 360.0 
         if (angle >= 360.0) angle = angle - 360.0 
      enddo
    
      end subroutine unwind
!===============================================================================

