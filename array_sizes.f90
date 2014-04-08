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
      module array_sizes ! Utility module for F90/95 for SAC files
!=======================================================================      
      implicit none
 
!        * parameters
         integer, parameter :: np1=181
         integer, parameter :: np2=41
         integer, parameter :: np2int=161
         integer, parameter :: np2XC=np2
         integer, parameter :: np2XCint=np2int         
         integer, parameter :: npc=500
         integer, parameter :: np=200000
         integer, parameter :: nsurfmax = 1000 ; ! this one is for stack 
 
!=======================================================================
      end module array_sizes
!=======================================================================
!     END OF sheba_config module
!=======================================================================

