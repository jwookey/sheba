!===============================================================================
!  S H E B A - Shear-wave Birefringence Analysis
!===============================================================================
!  This software is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!===============================================================================
!
!  James Wookey, School of Earth Sciences, University of Bristol
!  CVS: $Revision: 1.4 $ $Date: 2011/02/11 16:09:41 $
!

!===============================================================================
   program sheba_main
!===============================================================================
      use f90sac  ! use the f90sac module
      use sheba_config ! use the config module
!-------------------------------------------------------------------------------
      implicit none
         real, parameter :: sheba_version = 1.00

! ** print copyright information
     call sheba_banner(sheba_version) 

!  ** initialise configuration
      call sheba_config_init()

!  ** get input file information
      call read_config()

!  ** run the analysis
      call sheba()

!  ** finish up
      
      stop
   end program sheba_main
!===============================================================================
