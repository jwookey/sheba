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

!===============================================================================
      subroutine sheba_banner()
!===============================================================================
!  
!     Print sheba startup information
!
      implicit none
      
      write(*,'("")')

      write(*,'(a,a)') '=========================================', &
                       '======================================='
      write(*,'(a,a)') '-                                   S H E', & 
                       ' B A                                  -'
      write(*,'(a,a)') '-----------------------------------------', & 
                       '---------------------------------------'
      write(*,'(a,a)') '-                       Shear-wave Birefr', & 
                       'ingence Analysis                      -'
      write(*,'(a,a)') '=========================================', & 
                       '======================================='
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

!===============================================================================
      subroutine unwind_pm90(angle)
!===============================================================================
!  
!     'Unwind' an angle to between -90 to 90
!
      implicit none
      real :: angle
      integer i
      
      do i = 1 , 1000
         if (angle >= -90.0 .and. angle < 90.0) exit
         if (angle < -90.0) angle = angle + 180.0 
         if (angle >= 90.0) angle = angle - 180.0 
      enddo
    
      end subroutine unwind_pm90
!===============================================================================

