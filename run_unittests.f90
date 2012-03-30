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
   program run_unittests
!===============================================================================
   use fruit ! xUnit module.
   use sheba_test ! sheba testing module
!===============================================================================

   call init_fruit()
   call test_split()
   call fruit_summary()

end program run_unittests
