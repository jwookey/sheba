!=======================================================================
!     S H E B A - Shear-wave Birefringence Analysis
!=======================================================================
!     Fortran 90/95 Source Code File
!-----------------------------------------------------------------------
!
!     PROGRAM : sheba
!     FILE    : event_info.f
!     AUTHOR  : James Wookey
!     PLACE   : School of Earth Sciences, University of Leeds
!     DATE    : December 2003
!     PURPOSE : 
!     VERSION : 1.0
!     COMPLETE: No
!     COMMENTS: 
!
!-----------------------------------------------------------------------
!     This software is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!-----------------------------------------------------------------------
!
!      SHEBA is a prototype code to implement various algorithms for
!      analysing shear-wave splitting. 
!
!      The basic analysis is after Silver and Chan, 1991.
!      Also incorporated are many routines from Nicholas Teanby,
!      Georg Rumpker, and Chris Chapman. Individual routines are 
!      attributed to their various authors
!
!-----------------------------------------------------------------------
!     Changes log
!-----------------------------------------------------------------------
!     2003-12-04     * Incept date

!=======================================================================
      module event_info ! Utility module for F90/95 for SAC files
!=======================================================================      
      implicit none

!     ** DATA STRUCTURE FOR EVENT INFORMATION
         type event_info_type
            real*4 :: dist,az,baz ! event info
            real*4 :: tlag,dtlag,fast,dfast,spol,dspol ! splitting pars
            real*4 :: wbeg,wend ! best window
            real*4 :: snr ! signal-to-noise
            real*4 :: eigrat_orig,eigrat_corr
            integer :: ndf
            integer :: ibest ! best window number
!        ** other information
            real :: error_grid_tlag_int           
         end type event_info_type     

         type (event_info_type),public :: event
 
!=======================================================================
      end module event_info
!=======================================================================
!     END OF event_info module
!=======================================================================
