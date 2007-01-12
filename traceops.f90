!=======================================================================
!     S H E B A - Shear-wave Birefringence Analysis
!=======================================================================
!     Fortran 90/95 Source Code File
!-----------------------------------------------------------------------
!
!     PROGRAM : sheba
!     FILE    : traceops.f
!     AUTHOR  : James Wookey
!     PLACE   : School of Earth Sciences, University of Leeds
!     DATE    : December 2003
!     PURPOSE : Subroutines for various trace tasks 
!     VERSION : 1.0
!     COMPLETE: No
!     COMMENTS: 
!
!-----------------------------------------------------------------------
!     This software is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!-----------------------------------------------------------------------

!=======================================================================
      subroutine windowed_coveig(h1,h2,lambda1,lambda2,t1,t2)
!=======================================================================
!
!     Calculate the eigenvalues of the covariance matrix
!
!-----------------------------------------------------------------------
      use f90sac  ! use the f90sac module
!-----------------------------------------------------------------------
      implicit none
      type (SACTrace) :: h1,h2
      type (SACTrace) :: wh1,wh2
      real :: lambda1,lambda2,t1,t2
      
      call f90sac_window(h1,wh1,t1,t2)
      call f90sac_window(h2,wh2,t1,t2)
      
      call coveig(wh1,wh2,lambda1,lambda2) 
      
      return
      end subroutine windowed_coveig
!=======================================================================

!=======================================================================
      subroutine coveig(h1,h2,lambda1,lambda2)
!=======================================================================
!
!     Calculate the eigenvalues of the covariance matrix
!
!-----------------------------------------------------------------------
      use f90sac  ! use the f90sac module
!-----------------------------------------------------------------------
      implicit none
      type (SACTrace) :: h1,h2
      real :: lambda1,lambda2
      real :: cov(2,2),vec1(2),vec2(2)
      integer i

!  ** calc cross correlation matrix **
		cov(1,1) = 0.
		cov(2,2) = 0.
		cov(1,2) = 0.
		cov(2,1) = 0.

		do 1 i=1,h1 % npts
		   cov(1,1) = cov(1,1) + h1 % trace(i)**2
		   cov(2,2) = cov(2,2) + h2 % trace(i)**2
		   cov(1,2) = cov(1,2) + h1 % trace(i)*h2 % trace(i)
1		continue
		cov(2,1) = cov(1,2)

      call zeigen2x2(cov,lambda1,lambda2,vec1,vec2)

      
      return
      end subroutine coveig
!=======================================================================


!=======================================================================
   subroutine check_windows(tr)
!=======================================================================
!
!     Check that the window lengths are sensible given the max_tlag
!     selected
!
!-----------------------------------------------------------------------
   use f90sac  ! use the f90sac module
   use event_info ! use the event info module
   use sheba_config ! use the config module
!-----------------------------------------------------------------------
   implicit none
      type (SACTrace) :: tr
      real :: min_window_length
      character :: ans
!     ** single window     
         if (config%nwindbeg==1 .and. config%nwindend==1) then
            min_window_length = tr%user2 - tr%user0
!     ** multiple windows      
         else
            min_window_length = tr%user2 - tr%user1
         endif

!     ** compare with max_tlag
         if (min_window_length <= config%max_tlag) then
            write(0,'(a)') '---------------------------------------'
            write(0,'(a)') '  WARNING! Smallest window is smaller'
            write(0,'(a)') '  than MAX_TLAG. This is probably not'
            write(0,'(a)') '  a very good idea ...               '
            write(0,'(a)') '---------------------------------------'
            pause ' Hit ^c to abort or return to continue'
         endif
         
         
      return
   end subroutine check_windows
!=======================================================================


!=======================================================================
      subroutine check_trace_headers(tr)
!=======================================================================
!
!     Check that the required headers for the splitting analysis are set
!     These are:
!
!     az,baz,npts,delta,b,e ! Standard SAC variables 
!     user0..3 - window definitions for cluster analysis
!     
!     Also check whether user4..8 are unset. If they are NOT, flag a
!     warning and get permission to continue
!
!-----------------------------------------------------------------------
      use f90sac  ! use the f90sac module
      use event_info ! use the event info module
      use sheba_config ! use the config module
!-----------------------------------------------------------------------
      implicit none
      type (SACTrace) :: tr
      logical miss
      character (len = 80) :: mstr
      
      print*,'> Checking trace header'
      
      miss = .false.
      mstr = ''

      if (tr % az == SAC_rnull) then
         miss = .true.
         mstr = trim(mstr) // 'az '
      endif 

      if (tr % baz == SAC_rnull) then
         miss = .true.
         mstr = trim(mstr) // 'baz '
      endif 

      if (tr % cmpaz == SAC_rnull) then
         miss = .true.
         mstr = trim(mstr) // 'cmpaz '
      endif 

      if (tr % cmpinc == SAC_rnull) then
         miss = .true.
         mstr = trim(mstr) // 'cmpinc '
      endif 

      if (tr % b == SAC_rnull) then
         miss = .true.
         mstr = trim(mstr) // 'b '
      endif 

      if (tr % e == SAC_rnull) then
         miss = .true.
         mstr = trim(mstr) // 'e '
      endif 

      if (tr % delta == SAC_rnull) then
         miss = .true.
         mstr = trim(mstr) // 'delta '
      endif 

      if (tr % npts == SAC_inull) then
         miss = .true.
         mstr = trim(mstr) // 'npts '
      endif 
      
      if (tr % user0 == SAC_rnull) then
         miss = .true.
         mstr = trim(mstr) // 'user0 '
      endif 
      
      if (tr % user1 == SAC_rnull) then
         if (.not.(config%nwindbeg==1 .and. config%nwindend==1)) then
            miss = .true.
            mstr = trim(mstr) // 'user1 '
         endif
      endif 

      if (tr % user2 == SAC_rnull) then
         miss = .true.
         mstr = trim(mstr) // 'user2 '
      endif 
      
      if (tr % user3 == SAC_rnull) then
         if (.not.(config%nwindbeg==1 .and. config%nwindend==1)) then
            miss = .true.
            mstr = trim(mstr) // 'user3 '
         endif   
      endif 

      if (miss) then
         print*,'SHEBA: ERROR! Input trace is missing headers: ', & 
               mstr    
         STOP
      else
!         print*,'> Trace headers OK'
      endif

!  **  check ordering of user0..3
      if (config % nwindbeg == 1  .and. config % nwindend == 1) then
!  **  single window mode
         if (tr % user2 < tr % user0) then ! broken
         print*,'> Bad window positions!'
         print*,tr % user0,  tr % user2
         STOP
         endif
      else
            
         if (tr % user3 > tr % user2 .and. tr % user2 > tr % user1 .and. & 
           tr % user1 > tr % user0) then
!           print*,'> Window positions OK'
         else
            print*,'> Bad window positions!'
            print*,tr % user0, tr % user1, tr % user2, tr % user3
            STOP
         endif
      endif

!  **  check for overwrite
      miss = .false.
      mstr = ''
       
      if (tr % user4 == SAC_rnull) then
         miss = .true.
         mstr = trim(mstr) // 'user4 '
      endif 
      
      if (tr % user5 == SAC_rnull) then
         miss = .true.
         mstr = trim(mstr) // 'user5 '
      endif 

      if (tr % user6 == SAC_rnull) then
         miss = .true.
         mstr = trim(mstr) // 'user6 '
      endif 
      
      if (tr % user7 == SAC_rnull) then
         miss = .true.
         mstr = trim(mstr) // 'user7 '
      endif 

      if (tr % user8 == SAC_rnull) then
         miss = .true.
         mstr = trim(mstr) // 'user8 '
      endif 

!  **  check ordering of user0..3
      if (miss) then
!         print*,'> Overwrite OK'
      else
         print*,'> WARNING! The following positions are not null: ', & 
            mstr
         print*,'> These positions will be overwritten by SHEBA.'  
!         pause
      endif
      
         
      return
      end subroutine check_trace_headers
!=======================================================================

!=======================================================================
      subroutine upload_splitpar(tr)
!=======================================================================
!
!     Upload the splitting parameters (stored in the event structure)
!     into the SAC trace tr. The following storage positions are used
!     
!     user4 = best s-window begin
!     user5 = best s-window end
!     user6 = fast
!     user7 = tlag
!     user8 = spol
!     user9 = ?
!
!-----------------------------------------------------------------------
      use f90sac  ! use the f90sac module
      use event_info ! use the event info module
!-----------------------------------------------------------------------
      implicit none
      type (SACTrace) :: tr
!  **  upload the parameters
      tr % user4 = event % wbeg
      tr % user5 = event % wend
      tr % user6 = event % fast
      tr % user7 = event % tlag
      tr % user8 = event % spol
!  **  done      
      return
      end subroutine upload_splitpar
!=======================================================================
