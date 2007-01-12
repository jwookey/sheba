!=======================================================================
!     S H E B A - Shear-wave Birefringence Analysis
!=======================================================================
!     Fortran 90/95 Source Code File
!-----------------------------------------------------------------------
!
!     PROGRAM : sheba
!     FILE    : sheba_config.f
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
      module sheba_config ! Utility module for F90/95 for SAC files
!=======================================================================      
      implicit none
!        * data structure definition 
         type sheba_config_type
!           * filename of input traces
            character (len=50) :: fname_base

!           * filename extensions
            character (len=3) :: comp(3)

!           * cluster analysis specific stuff            
            integer :: nwindbeg, nwindend

!           * shear-wave splitting mode (0=trans. ener. 1=eigenvalue)            
            integer :: imode
            real :: source_pol ! required for transverse energy min.

!           * tlag_scale            
            real :: max_tlag

!           * pre-analysis station correction            
            integer :: iuma_corr
            real :: uma_dt, uma_phi
            
!           * ScS specific stuff
            integer :: iscs_corr ! scs correct flag
            real :: scs_slw ! slowness of ScS phase s/deg

!           * Useful?
            real :: input_h1_pol
            
!        ** post-analysis correction (source-side anisotropy)
            integer :: i_src_corr
            real :: src_tlag, src_fast

!        ** pre-analysis ABC rotation
            integer :: i_rotate_to_ABC
            real :: slw

!           * Cluster analysis parameters, hardwired
            real :: dtlag_max = 2.0
            real :: dfast_max = 20.0
            real :: fast_scale = 90.0
            integer :: min_pts_per_cluster = 5
            integer :: max_no_clusters = 5
 
         end type sheba_config_type     

!        * data structure for various config flags
         type (sheba_config_type),public,save :: config
 
!        * useful stuff
         integer, parameter :: iulog = 8 ! logfile unit number
                  
!        * maths constants and other useful things
         real*8, parameter :: pi = 3.141592653589793238462643d0 ;
         real*8, parameter :: to_rad = 1.74532925199433D-002 ;  
         real*8, parameter :: to_deg = 57.2957795130823d0 ;  
         real*4, parameter :: angtol = 0.001 ; ! angle comparison tol.
         real*8, parameter :: to_km = 111.194926644559 ;
 
!=======================================================================
      end module sheba_config
!=======================================================================
!     END OF sheba_config module
!=======================================================================

