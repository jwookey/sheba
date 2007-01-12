!=======================================================================
!     S H E B A - Shear-wave Birefringence Analysis
!=======================================================================
!     Fortran 90/95 Source Code File
!-----------------------------------------------------------------------
!
!     PROGRAM : sheba
!     FILE    : array_sizes.f
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
      module array_sizes ! Utility module for F90/95 for SAC files
!=======================================================================      
      implicit none
 
!        * parameters
         integer, parameter :: np1=181
         integer, parameter :: np2=41
         integer, parameter :: np2int=161
         integer, parameter :: npc=500
         integer, parameter :: np=200000
         integer, parameter :: nsurfmax = 1000 ; ! this one is for stack 
 
!=======================================================================
      end module array_sizes
!=======================================================================
!     END OF sheba_config module
!=======================================================================

