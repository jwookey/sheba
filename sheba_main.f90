!=======================================================================
!     S H E B A - Shear-wave Birefringence Analysis
!=======================================================================
!     Fortran 90/95 Source Code File
!-----------------------------------------------------------------------
!
!     PROGRAM : sheba
!     FILE    : sheba.f
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
!      sheba is a prototype code to implement various algorithms for
!      analysing shear-wave splitting. 
!
!-----------------------------------------------------------------------
!     Changes log
!-----------------------------------------------------------------------
!     2003-12-04     * Incept date


!=======================================================================
      program sheba_main
!=======================================================================
      use f90sac  ! use the f90sac module
      use sheba_config ! use the config module
!-----------------------------------------------------------------------
      implicit none
         real, parameter :: sheba_version = 0.99

!        * print copyright information
         call sheba_banner(sheba_version) 

!     * get input file information
         call read_config()

!     * run the analysis
         call sheba()

!     * finish up
      
         stop
      end program sheba_main
!=======================================================================
