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
      module event_info ! Utility module for F90/95 for SAC files
!=======================================================================      
      use array_sizes
      implicit none

!     ** DATA STRUCTURE FOR EVENT INFORMATION
         type event_info_type
            real*4 :: dist,az,baz ! event info
            real*4 :: tlag,dtlag,fast,dfast,spol,dspol ! splitting pars
            real*4 :: tlagXC,dtlagXC,fastXC,dfastXC ! XC splitting pars
            real*4 :: wbeg,wend ! best window
            real*4 :: snr ! signal-to-noise
            real*4 :: eigrat_orig,eigrat_corr
            integer :: ndf
            integer :: ibest ! best window number
            real*4 :: Quality ! AW quality factor (-1 -> 1)
            real*4 :: intensity ! splitting intensity as predicted from
                                ! the best fitting tlag, fast and spol.
            real*4 :: intensity_p ! splitting intensity calculated by
                                  ! projection                    
!        ** other information
            real :: error_grid_tlag_int           

!        ** parameter search space
            real*4 :: fast_vector(np1), tlag_vector(np2int) 

!        ** lam1, lam2 grids before normalisation
            real*4 :: lam1_raw_grid(np1,np2int),lam2_raw_grid(np1,np2int) 

!        ** error grids after f-test normalisation
            real*4 :: lam2_norm_grid(np1,np2int)

!        ** multiwindow information
            integer :: nwindows
            real*4 :: mw_wbeg(npc)
            real*4 :: mw_wend(npc)
            real*4 :: mw_tlag(npc)
            real*4 :: mw_dtlag(npc)
            real*4 :: mw_fast(npc)
            real*4 :: mw_dfast(npc)

!        ** cluster information
            integer :: ncluster
            integer :: kbest
            real*4  :: cluster_xc0(npc)
            real*4  :: cluster_yc0(npc)
            real*4  :: cluster_vxc0(npc)
            real*4  :: cluster_vyc0(npc)

         end type event_info_type     

         type (event_info_type),public :: event
 
!=======================================================================
      end module event_info
!=======================================================================
!     END OF event_info module
!=======================================================================
