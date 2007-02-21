!===============================================================================
!-------------------------------------------------------------------------------
!
!  Fortran 90/95 Source Code File
!
!-------------------------------------------------------------------------------
!===============================================================================
!
!  PROGRAM : f90sac_distrib
!  VERSION : 3.2 (Distrib. version)
!  CVS: $Revision: 1.2 $ $Date: 2007/02/21 15:55:32 $
!
!  (C) James Wookey
!  Department of Earth Sciences, University of Bristol
!  Wills Memorial Building, Queen's Road, Bristol, BR8 1RJ, UK
!  j.wookey@bristol.ac.uk
!
!-------------------------------------------------------------------------------
!
!   The module provides data structures and functions for reading,
!   writing and handling SAC files in Fortran 90/95.
!
!   Please report bugs/problems to email address above
!
!   This has been successfully compiled using Solaris f90 and IFC versions
!   6-8.1. NOTE: Check parameter f90sac_32bit_record_length is set to the
!   correct value for your compiler. UPDATE the free compilers g95 and 
!   gfortran also work. Ditto xlf.
!
!   NOTE: This version of the code assumes IO filestream 99 is available 
!         for reading and writing. 
!
!-------------------------------------------------------------------------------
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to: 
!
!   Free Software Foundation, Inc., 59 Temple Place, Suite 330, 
!   Boston, MA  02111-1307  USA
!
!   See http://www.gnu.org/licenses/gpl.txt for full GPL terms.
!
!===============================================================================

   module f90sac ! Utility module for F90/95 for SAC files

!===============================================================================
   implicit none

!  ** DECLARE CONTAINED FUNCTIONS
      public :: f90sac_compare_origin_time
      public :: f90sac_orient2d
      public :: f90sac_rotate2d
      public :: f90sac_rotate2d_rz
      public :: f90sac_unwind
      public :: f90sac_tshift
      public :: f90sac_newtrace
      public :: f90sac_clonetrace
      public :: f90sac_copytraceheader
      public :: f90sac_readtrace
      public :: f90sac_writetrace 
      public :: f90sac_enumhdr
      public :: f90sac_getfhdr
      public :: f90sac_getihdr
      public :: f90sac_getlhdr
      public :: f90sac_getkhdr
      public :: f90sac_setfhdr
      public :: f90sac_setihdr
      public :: f90sac_setlhdr
      public :: f90sac_setkhdr
      public :: f90sac_setdate
      public :: f90sac_ymd2jd
      public :: f90sac_window
      public :: f90sac_destroytrace
      public :: f90sac_setstation
      public :: f90sac_setevent
      public :: f90sac_cattraces
      public :: f90sac_dateseed
                
!  ** define a long (32 bit) integer and 32 bit real      
      integer, parameter, private :: int4 = selected_int_kind(9) ;
      integer, parameter, private :: real4 = selected_real_kind(6,37) ;

!  ** define the record length in a sequential access file for a 32 bit number
!  ** this is compiler dependent: 
!        IFORT/IFC Version >= 8.0 = 1 
!        IFC Version < v8.0 = 4
!        Solaris F90 = 4
!        g95/gfortran = 4
      integer, parameter :: f90sac_32bit_record_length = 4 ;
      
!  ** define the unit number to use for reading and writing (opened and closed
!  ** within each call)
      integer, parameter :: f90sac_iounit = 99 ;


!  ** set whether automatic byteswapping is done
      logical, parameter :: f90sac_force_byteswap = .false. ;

!  ** define the 'active' value for the iTraceActive header
      integer, parameter :: TraceIsActive = 11075 ;
      
!  ** noise generator seed value
      integer, private :: f90sac_random_seed ;      

!=============================================================================== 
!  ** Define a specialised data structure for containing SAC files
!=============================================================================== 
      type SACtrace
!     ** Header floating point part
         real(real4) :: delta,depmin,depmax,scale,odelta,b,e,o,a,internal0
         real(real4) :: t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,f
         real(real4) :: resp0,resp1,resp2,resp3,resp4,resp5
         real(real4) :: resp6,resp7,resp8,resp9
         real(real4) :: stla,stlo,stel,stdp,evla,evlo,evel,evdp,mag
         real(real4) :: user0,user1,user2,user3,user4
         real(real4) :: user5,user6,user7,user8,user9
         real(real4) :: dist,az,baz,gcarc,internal1,internal2,depmen
         real(real4) :: cmpaz,cmpinc
         real(real4) :: xminimum,xmaximum,yminimum,ymaximum
         real(real4) :: unused1,unused2,unused3,unused4
         real(real4) :: unused5,unused6,unused7
!     ** Header integer part
         integer(int4) :: nzyear,nzjday,nzhour,nzmin,nzsec,nzmsec
         integer(int4) :: nvhdr,norid,nevid,npts
         integer(int4) :: internal3,nwfid,nxsize,nysize,unused8
         integer(int4) :: iftype,idep,iztype,unused9
         integer(int4) :: iinst,istreg,ievreg,ievtyp
         integer(int4) :: iqual,isynth,imagtyp,imagsrc
         integer(int4) :: unused10,unused11,unused12,unused13,unused14
         integer(int4) :: unused15,unused16,unused17         
!     ** Header logical part (stored as integers)
         integer(int4) ::  leven,lpspol,lovrok,lcalda,unused18 
!     ** Header character part
         character (len = 16) :: kevnm
         character (len = 8) :: kstnm,khole,ko,ka
         character (len = 8) :: kt0,kt1,kt2,kt3,kt4
         character (len = 8) :: kt5,kt6,kt7,kt8,kt9
         character (len = 8) :: kf,kuser0,kuser1,kuser2
         character (len = 8) :: kcmpnm,knetwk,kdatrd,kinst
!     ** the trace
         real, pointer :: trace(:)
         
!     ** ACTIVE flag, do not alter.
         integer :: iTraceActive
         
      end type SACtrace      

!=============================================================================== 
!  ** Define a data structure for containing SAC xy files
!=============================================================================== 
      type SACxy
!     ** Header floating point part
         real(real4) :: delta,depmin,depmax,scale,odelta,b,e,o,a,internal0
         real(real4) :: t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,f
         real(real4) :: resp0,resp1,resp2,resp3,resp4,resp5
         real(real4) :: resp6,resp7,resp8,resp9
         real(real4) :: stla,stlo,stel,stdp,evla,evlo,evel,evdp,mag
         real(real4) :: user0,user1,user2,user3,user4
         real(real4) :: user5,user6,user7,user8,user9
         real(real4) :: dist,az,baz,gcarc,internal1,internal2,depmen
         real(real4) :: cmpaz,cmpinc
         real(real4) :: xminimum,xmaximum,yminimum,ymaximum
         real(real4) :: unused1,unused2,unused3,unused4
         real(real4) :: unused5,unused6,unused7
!     ** Header integer part
         integer(int4) :: nzyear,nzjday,nzhour,nzmin,nzsec,nzmsec
         integer(int4) :: nvhdr,norid,nevid,npts
         integer(int4) :: internal3,nwfid,nxsize,nysize,unused8
         integer(int4) :: iftype,idep,iztype,unused9
         integer(int4) :: iinst,istreg,ievreg,ievtyp
         integer(int4) :: iqual,isynth,imagtyp,imagsrc
         integer(int4) :: unused10,unused11,unused12,unused13,unused14
         integer(int4) :: unused15,unused16,unused17         
!     ** Header logical part (stored as integers)
         integer(int4) ::  leven,lpspol,lovrok,lcalda,unused18 
!     ** Header character part
         character (len = 16) :: kevnm
         character (len = 8) :: kstnm,khole,ko,ka
         character (len = 8) :: kt0,kt1,kt2,kt3,kt4
         character (len = 8) :: kt5,kt6,kt7,kt8,kt9
         character (len = 8) :: kf,kuser0,kuser1,kuser2
         character (len = 8) :: kcmpnm,knetwk,kdatrd,kinst
!     ** The trace

         real, pointer :: x(:)
         real, pointer :: y(:)
      end type SACxy      

!  ** tolerance for the comparison of angles      
      real, parameter :: f90sac_angle_tolerance = 0.001

!  ** NULL values set in SAC objects      
      real, parameter :: SAC_rnull = -12345.0
      integer, parameter :: SAC_inull = -12345
      integer, parameter :: SAC_lnull = -12345
      character (len = 8) :: SAC_cnull = '-12345'

!  ** OPTIONAL suppression of warnings, set to 1 to supress      
      integer,public :: f90sac_suppress_warnings
      
!===============================================================================
!
!  ** MODULE SUBROUTINES
!
!===============================================================================

   CONTAINS

!===============================================================================
   function f90sac_compare_origin_time(t1,t2)
!===============================================================================
!
!     Compare origin times of two SAC traces (using nzyear,nzjday etc)
!        returns:
!        -1 = t1 is earlier
!         0 = same time
!         1 = t2 is earlier
!
      implicit none
      type (SACtrace) :: t1,t2
      integer f90sac_compare_origin_time
      
!
      f90sac_compare_origin_time = 0
      
      if (t1%nzyear<t2%nzyear) then
         f90sac_compare_origin_time=-1
         return
      elseif (t1%nzyear>t2%nzyear) then
         f90sac_compare_origin_time=1
         return
      else
         if (t1%nzjday<t2%nzjday) then
            f90sac_compare_origin_time=-1
            return
         elseif (t1%nzjday>t2%nzjday) then
            f90sac_compare_origin_time=1
            return
         else
            if (t1%nzhour<t2%nzhour) then
               f90sac_compare_origin_time=-1
               return
            elseif (t1%nzhour>t2%nzhour) then
               f90sac_compare_origin_time=1
               return
            else
               if (t1%nzmin<t2%nzmin) then
                  f90sac_compare_origin_time=-1
                  return
               elseif (t1%nzmin>t2%nzmin) then
                  f90sac_compare_origin_time=1
                  return
               else
                  if (t1%nzsec<t2%nzsec) then
                     f90sac_compare_origin_time=-1
                     return
                  elseif (t1%nzsec>t2%nzsec) then
                     f90sac_compare_origin_time=1
                     return
                  else
                     if (t1%nzmsec<t2%nzmsec) then
                      f90sac_compare_origin_time=-1
                        return
                     elseif (t1%nzmsec>t2%nzmsec) then
                        f90sac_compare_origin_time=1
                        return
                     else
                        f90sac_compare_origin_time=0
                        return
                     endif         
                  endif         
               endif         
            endif   
         endif   
      endif

         
      return
      end function f90sac_compare_origin_time
!===============================================================================


!===============================================================================
   subroutine f90sac_cattraces(t1,t2,tc)
!===============================================================================
!
!     Concatenate two traces. t1 // t2 = tc. Header information is taken from
!     t1, start time is from t1
!
      implicit none
      type (SACtrace) :: t1,t2,tc
      integer :: npts_new

!  ** check delta is the same for both traces. 
      if (t1 % delta /= t2 % delta) then     
         write(0,'(a)') &
         'F90SAC_CATTRACES: Error: Traces have different deltas'
         STOP     
      endif      

      npts_new = t1 % npts + t2 % npts
      
!  ** make a new trace      
      call f90sac_newtrace(npts_new, t1 % delta, tc)

!  ** copy the trace header to the new trace      
      call f90sac_copytraceheader(t1,tc)
      tc % npts = npts_new

!  ** set the new end time (if begin time is not null)
      if (tc%b /= SAC_rnull) then     
         tc % e = tc % b + real(tc % npts) * tc % delta
      endif
      
!  ** copy the traces      
      tc % trace(1:t1 % npts) = t1 % trace(1:t1 % npts)
      tc % trace(t1%npts + 1 : t1%npts + t2%npts) = t2%trace(1:t2 % npts)
      
         
      return
   end subroutine f90sac_cattraces
!===============================================================================

!===============================================================================
   subroutine f90sac_setevent(trace,lat,lon,depth)
!===============================================================================
!
!     Set event information in SAC object trace
!
      implicit none
      type (SACtrace) :: trace
      real :: lat,lon,depth
      
      trace % evla = lat
      trace % evlo = lon
      trace % evdp = depth
      
      return
   end subroutine f90sac_setevent
!===============================================================================

!===============================================================================
   subroutine f90sac_setstation(trace,lat,lon,elevation)
!===============================================================================
!
!     Set station information in SAC object trace
!
      implicit none
      type (SACtrace) :: trace
      real :: lat,lon,elevation
      
      trace % stla = lat
      trace % stlo = lon
      trace % stel = elevation
      
      return
   end subroutine f90sac_setstation
!===============================================================================

!===============================================================================
   subroutine f90sac_window(tr_in,tr_out,t1,t2)
!===============================================================================
!
!     Window the data in SAC object .between t1 and t2 and produce a new
!     trace containing this data
!
      implicit none
      type (SACtrace) :: tr_in,tr_out
      
      integer iwbeg,iwend,new_npts,istatus,i
      real t1,t2
      real new_b,new_e

!  ** calculate indices of window      
      iwbeg = nint((t1-tr_in % b)/tr_in % delta)+1
      iwend = nint((t2-tr_in % b)/tr_in % delta)+1   
      
!  ** Check these values
      if (iwbeg <=0 .or. iwend > tr_in % npts) then
         write(0,'(a)') &
         'F90SAC_WINDOW: Error: Window limits out of range of data.'
         STOP     
      endif      

!  ** Check these values
      if (iwend <= iwbeg) then
         write(0,'(a)') &
         'F90SAC_WINDOW: Error: Resulting trace has no data points.'     
         STOP
      endif      

!  ** calculate new header values
      new_npts = iwend - iwbeg + 1
      new_b = tr_in % b + real(iwbeg-1)*tr_in % delta
      new_e = new_b + real(new_npts-1)*tr_in % delta

!  ** create a new trace
      call f90sac_clonetrace(tr_in,tr_out)

!  ** deallocate it's trace memory, and reallocate the correct amount           
      deallocate (tr_out % trace, stat = istatus)
      allocate (tr_out % trace(new_npts))

!  ** set header values
      tr_out % npts = new_npts
      tr_out % b = new_b
      tr_out % e = new_e

!  ** copy the trace
      do i = 1,new_npts
         tr_out % trace(i) = tr_in % trace( iwbeg + (i-1) )          
      enddo
      return
   end subroutine f90sac_window
!===============================================================================

!===============================================================================
   subroutine f90sac_setdate(trace,iyr,ijd,ihr,imi,ise,ims)
!===============================================================================
!
!     Set the date and time in a SAC object
!
      implicit none
      type (SACtrace) :: trace
      
      integer iyr,ijd,ihr,imi,ise,ims
      
      trace % nzyear = iyr
      trace % nzjday = ijd
      trace % nzhour = ihr
      trace % nzmin = imi
      trace % nzsec = ise
      trace % nzmsec = ims
      
      return
   end subroutine f90sac_setdate
!===============================================================================
     
!===============================================================================
   subroutine f90sac_orient2d(t1,t2)
!===============================================================================
!
!     take 2 traces in arbitrary (orthogonal) orientation and rotate to a 
!     north-east reference frame
!
!     t1 :  (I/O) SAC trace, on output this holds the north component
!     t2 :  (I/O) SAC trace, on output this holds the east component
!
      implicit none
      type (SACtrace) :: t1,t2
      
      real theta,cmpazdiff

!   ** check that cmpaz is set
      if (t1 % cmpaz == SAC_rnull .or. t2 % cmpaz == SAC_rnull) then
         write(0,'(a)') &
         'F90SAC_ORIENT2D: Error: CMPAZ Header is not set'
         STOP
      endif

      
!  ** check for orthogonality
      cmpazdiff = abs(t1 % cmpaz - t2 % cmpaz)
      if ( .not.((abs(cmpazdiff-90.0) <= f90sac_angle_tolerance) &
            .or. (abs(cmpazdiff-270.0) <= f90sac_angle_tolerance))) then
         write(0,'(a)') &
         'F90SAC_ORIENT2D: Error: Input components are not orthogonal'
         STOP
      endif      

!  ** first rotate t1 to the north
      call f90sac_rotate2d(t1,t2,-t1 % cmpaz)

!  ** now check whether t2 needs to be reversed
      if ( (abs(t2 % cmpaz-270.0) <= f90sac_angle_tolerance) ) then
         t2 % trace(1:t2 % npts) = -t2 % trace(1:t2 % npts)
      endif
      
!  ** update component azimuths
      t1 % cmpaz = 0.0   ! north
      t2 % cmpaz = 90.0  ! east

!  ** update component names
      write(t1 % kcmpnm,'(f5.1)') t1 % cmpaz        
      write(t2 % kcmpnm,'(f5.1)') t2 % cmpaz        
      
      return
   end subroutine f90sac_orient2d
!===============================================================================

!===============================================================================
      function f90sac_orth2d(t1,t2)
!===============================================================================
!
!     Check for the (azimuthal) orthogonality between two traces
!
      implicit none
      type (SACtrace) :: t1,t2
      integer f90sac_orth2d
      real cmpazdiff
      
!  ** check for orthogonality
      cmpazdiff = abs(t1 % cmpaz - t2 % cmpaz)
      if ( .not.((abs(cmpazdiff-90.0) <= f90sac_angle_tolerance) &
            .or. (abs(cmpazdiff-270.0) <= f90sac_angle_tolerance))) then
         f90sac_orth2d = 0
      else      
         f90sac_orth2d = 1
      endif
         
      return
      end function f90sac_orth2d
!===============================================================================

!===============================================================================
   subroutine f90sac_rotate2d_rz(t1,t2,theta)
!===============================================================================
!
!     rotate 2 SAC traces in the radial vertical plane
!
!     t1,t2 :  (I/O) SAC traces
!     theta : (I) angle to rotate by (clockwise)
!
      implicit none
      integer :: isamp
      type (SACtrace) :: t1,t2
      
      real theta, rotmat(2,2), sample(2,1), rsample(2,1),temp1,temp2
      real,parameter :: pi = 3.1415927410125732421875
      real :: cmpazdiff

!  ** check CMPINC header is present
      if (t1 % cmpinc == SAC_rnull .or. t2 % cmpinc == SAC_rnull) then
         if (f90sac_suppress_warnings /= 1) then
            write(0,'(a)') &
            'F90SAC_ROTATE2D_RZ: Warning: CMPINC Header is not set,'
            write(0,'(a)') &
            'F90SAC_ROTATE2D_RZ: Warning: setting CMPINC to 0,90'
         endif   
         t1 % cmpinc = 0.0
         t2 % cmpinc = 90.0
      endif

!  ** check for orthogonality
      cmpazdiff = abs(t1 % cmpinc - t2 % cmpinc)
      if ( .not.((abs(cmpazdiff-90.0) <= f90sac_angle_tolerance) &
            .or. (abs(cmpazdiff-270.0) <= f90sac_angle_tolerance))) then
         write(0,'(a)') &
            'F90SAC_ROTATE2D_RZ: Error: Input components are not orthogonal'
         STOP
      endif
            
!  ** check that traces are the same length
      if (t1 % npts /= t2 % npts) then
         write(0,'(a)') &
         'F90SAC_ROTATE2D_RZ: Error: Input components are different lengths'
         STOP
      endif            

!  ** make theta clockwise       
      theta = theta * (-1.0)
      
!  ** build rotation matrix
      rotmat(1,1) = cos(theta*pi/180.)
      rotmat(1,2) = sin(theta*pi/180.)
      rotmat(2,1) = -sin(theta*pi/180.)
      rotmat(2,2) = cos(theta*pi/180.)
      
!  ** rotate traces      
      do isamp = 1 , t1 % npts     
         sample(1,1) = t1 % trace(isamp) ; 
         sample(2,1) = t2 % trace(isamp) ; 
         
         rsample = matmul(rotmat,sample) ;
         
         t1 % trace(isamp) = rsample(1,1) ;
         t2 % trace(isamp) = rsample(2,1) ; 
      enddo

!  ** update component azimuths      
      t1 % cmpinc = t1 % cmpinc + theta
      t2 % cmpinc = t2 % cmpinc + theta

!  ** force to be in the range 0-360 degrees
      do ! forever
         if (t1 % cmpinc >= 0.0 .and. t1 % cmpinc < 360.0) exit
         if (t1 % cmpinc >= 360.0) t1 % cmpinc = t1 % cmpinc - 360.0
         if (t1 % cmpinc < 0.0) t1 % cmpinc = t1 % cmpinc + 360.0      
      enddo      
      
      do ! forever
         if (t2 % cmpinc >= 0.0 .and. t2 % cmpinc < 360.0) exit
         if (t2 % cmpinc >= 360.0) t2 % cmpinc = t2 % cmpinc - 360.0
         if (t2 % cmpinc < 0.0) t2 % cmpinc = t2 % cmpinc + 360.0      
      enddo
      
      return
   end subroutine f90sac_rotate2d_rz
!===============================================================================

!===============================================================================
   subroutine f90sac_rotate2d(t1,t2,theta)
!===============================================================================
!
!     rotate 2 SAC traces in azimuth
!
!     t1,t2 :  (I/O) SAC traces
!     theta : (I) angle to rotate by (clockwise)
!
      implicit none
      integer :: isamp
      type (SACtrace) :: t1,t2
      
      real theta, rotmat(2,2), sample(2,1), rsample(2,1),temp1,temp2
      real,parameter :: pi = 3.1415927410125732421875
      real :: cmpazdiff

!   ** check that cmpaz is set
      if (t1 % cmpaz == SAC_rnull .or. t2 % cmpaz == SAC_rnull) then
        if (f90sac_suppress_warnings /= 1) then
            write(0,'(a)') &
            'F90SAC_ROTATE2D: Warning: CMPAZ Header is not set'
            write(0,'(a)') &
            'F90SAC_ROTATE2D: Warning: Setting CMPAZ to 0,90'
         endif
         t1 % cmpaz = 0.
         t2 % cmpaz = 90.0         
      endif

!  ** check for orthogonality
      cmpazdiff = abs(t1 % cmpaz - t2 % cmpaz)
      if ( .not.((abs(cmpazdiff-90.0) <= f90sac_angle_tolerance) &
            .or. (abs(cmpazdiff-270.0) <= f90sac_angle_tolerance))) then
         write(0,'(a)') &
         'F90SAC_ROTATE2D: Error: Input components are not orthogonal'
         STOP
      endif
            
!  ** check that traces are the same length
      if (t1 % npts /= t2 % npts) then
         write(0,'(a)') &
         'F90SAC_ROTATE2D: Error: Input traces are different lengths'
         STOP
      endif            
      
!  ** build rotation matrix
      rotmat(1,1) = cos(theta*pi/180.)
      rotmat(1,2) = sin(theta*pi/180.)
      rotmat(2,1) = -sin(theta*pi/180.)
      rotmat(2,2) = cos(theta*pi/180.)
      
!  ** rotate traces      
      do isamp = 1 , t1 % npts     
         sample(1,1) = t1 % trace(isamp) ; 
         sample(2,1) = t2 % trace(isamp) ; 
         
         rsample = matmul(rotmat,sample) ;
         
         t1 % trace(isamp) = rsample(1,1) ;
         t2 % trace(isamp) = rsample(2,1) ; 
      enddo

!  ** update component azimuths      
      t1 % cmpaz = t1 % cmpaz + theta
      t2 % cmpaz = t2 % cmpaz + theta

!  ** force to be in the range 0-360 degrees
      do ! forever
         if (t1 % cmpaz >= 0.0 .and. t1 % cmpaz < 360.0) exit
         if (t1 % cmpaz >= 360.0) t1 % cmpaz = t1 % cmpaz - 360.0
         if (t1 % cmpaz < 0.0) t1 % cmpaz = t1 % cmpaz + 360.0      
      enddo      
      
      do ! forever
         if (t2 % cmpaz >= 0.0 .and. t2 % cmpaz < 360.0) exit
         if (t2 % cmpaz >= 360.0) t2 % cmpaz = t2 % cmpaz - 360.0
         if (t2 % cmpaz < 0.0) t2 % cmpaz = t2 % cmpaz + 360.0      
      enddo

!  ** update component names
      write(t1 % kcmpnm,'(f5.1)') t1 % cmpaz        
      write(t2 % kcmpnm,'(f5.1)') t2 % cmpaz        
      
      return
   end subroutine f90sac_rotate2d
!===============================================================================

!===============================================================================
   subroutine f90sac_unwind(angle)
!===============================================================================
!
!  unwind an angle to be in the range 0-360 degrees
!
!     angle :  (I/O) angle to unwind
!
      implicit none
      real :: angle
      real :: dt

      do ! forever
         if (angle >= 0.0 .and. angle < 360.0) exit
         if (angle >= 360.0) angle = angle - 360.0
         if (angle < 0.0) angle = angle + 360.0      
      enddo
      
      return
   end subroutine f90sac_unwind
!===============================================================================


!===============================================================================
   subroutine f90sac_tshift(trace,dt)
!===============================================================================
!
!     time shift a sac trace by dt (to nearest sample), zero additional
!     samples
!
!     trace :  (I/O) SAC trace
!     dt : (I) time to shift by
!
      implicit none
      integer :: isamp,ishift
      real :: dt
      type (SACtrace) :: trace
      
      ishift = nint (dt / (trace % delta) )
            
!  ** check for no shift      
      if (ishift == 0 .and. f90sac_suppress_warnings /= 1) then
         write(0,'(a)') &
            'F90SAC_TSHIFT: Warning: no shift applied, dt too small'
      endif      
         
!  ** shift array      
      trace % trace = cshift((trace % trace),-ishift)
      
!  ** if negative shift, zero last ishift points      
      if (ishift < 0) &
         trace % trace(trace % npts - abs(ishift): trace % npts) = 0.0
!  ** if positive shift, zero first ishift points      
      if (ishift > 0) trace % trace(1:abs(ishift)) = 0.0
      
      return
   end subroutine f90sac_tshift
!===============================================================================


!===============================================================================
   subroutine f90sac_destroytrace(tr)
!===============================================================================
!
!     Destroy an existing SAC trace. i.e. null all its header values and
!     deallocate its trace memory
!
      implicit none
      type (SACtrace) :: tr
      integer :: istatus ! number of samples for trace 

      tr%delta     = SAC_rnull
      tr%depmin    = SAC_rnull
      tr%depmax    = SAC_rnull
      tr%scale     = SAC_rnull
      tr%odelta    = SAC_rnull
      tr%b         = SAC_rnull
      tr%e         = SAC_rnull
      tr%o         = SAC_rnull
      tr%a         = SAC_rnull
      tr%internal0 = SAC_rnull
      tr%t0        = SAC_rnull
      tr%t1        = SAC_rnull
      tr%t2        = SAC_rnull
      tr%t3        = SAC_rnull
      tr%t4        = SAC_rnull
      tr%t5        = SAC_rnull
      tr%t6        = SAC_rnull
      tr%t7        = SAC_rnull
      tr%t8        = SAC_rnull
      tr%t9        = SAC_rnull
      tr%f         = SAC_rnull
      tr%resp0     = SAC_rnull
      tr%resp1     = SAC_rnull
      tr%resp2     = SAC_rnull
      tr%resp3     = SAC_rnull
      tr%resp4     = SAC_rnull
      tr%resp5     = SAC_rnull
      tr%resp6     = SAC_rnull
      tr%resp7     = SAC_rnull
      tr%resp8     = SAC_rnull
      tr%resp9     = SAC_rnull
      tr%stla      = SAC_rnull
      tr%stlo      = SAC_rnull
      tr%stel      = SAC_rnull
      tr%stdp      = SAC_rnull

      tr%evla      = SAC_rnull
      tr%evlo      = SAC_rnull
      tr%evel      = SAC_rnull
      tr%evdp      = SAC_rnull
      tr%mag       = SAC_rnull
      tr%user0     = SAC_rnull
      tr%user1     = SAC_rnull
      tr%user2     = SAC_rnull
      tr%user3     = SAC_rnull
      tr%user4     = SAC_rnull
      tr%user5     = SAC_rnull
      tr%user6     = SAC_rnull
      tr%user7     = SAC_rnull
      tr%user8     = SAC_rnull
      tr%user9     = SAC_rnull
      tr%dist      = SAC_rnull
      tr%az        = SAC_rnull
      tr%baz       = SAC_rnull
      tr%gcarc     = SAC_rnull
      tr%internal1 = SAC_rnull
      tr%internal2 = SAC_rnull
      tr%depmen    = SAC_rnull
      tr%cmpaz     = SAC_rnull
      tr%cmpinc    = SAC_rnull
      tr%xminimum  = SAC_rnull
      tr%xmaximum  = SAC_rnull
      tr%yminimum  = SAC_rnull
      tr%ymaximum  = SAC_rnull
      tr%unused1   = SAC_rnull
      tr%unused2   = SAC_rnull
      tr%unused3   = SAC_rnull
      tr%unused4   = SAC_rnull
      tr%unused5   = SAC_rnull
      tr%unused6   = SAC_rnull
      tr%unused7   = SAC_rnull

      tr%nzyear    = SAC_inull
      tr%nzjday    = SAC_inull
      tr%nzhour    = SAC_inull
      tr%nzmin     = SAC_inull
      tr%nzsec     = SAC_inull
      tr%nzmsec    = SAC_inull
      tr%nvhdr     = SAC_inull
      tr%norid     = SAC_inull
      tr%nevid     = SAC_inull
      tr%npts      = SAC_inull
      tr%internal3 = SAC_inull
      tr%nwfid     = SAC_inull
      tr%nxsize    = SAC_inull
      tr%nysize    = SAC_inull
      tr%unused8   = SAC_inull
      tr%iftype    = SAC_inull
      tr%idep      = SAC_inull
      tr%iztype    = SAC_inull
      tr%unused9   = SAC_inull
      tr%iinst     = SAC_inull
      tr%istreg    = SAC_inull
      tr%ievreg    = SAC_inull
      tr%ievtyp    = SAC_inull
      tr%iqual     = SAC_inull
      tr%isynth    = SAC_inull
      tr%imagtyp   = SAC_inull
      tr%imagsrc   = SAC_inull
      tr%unused10  = SAC_inull
      tr%unused11  = SAC_inull
      tr%unused12  = SAC_inull
      tr%unused13  = SAC_inull
      tr%unused14  = SAC_inull
      tr%unused15  = SAC_inull
      tr%unused16  = SAC_inull
      tr%unused17  = SAC_inull
      tr%leven     = SAC_lnull
      tr%lpspol    = SAC_lnull
      tr%lovrok    = SAC_lnull
      tr%lcalda    = SAC_lnull
      tr%unused18  = SAC_lnull

      tr%kstnm     = SAC_cnull
      tr%kevnm     = SAC_cnull
      tr%khole     = SAC_cnull
      tr%ko        = SAC_cnull
      tr%ka        = SAC_cnull
      tr%kt0       = SAC_cnull
      tr%kt1       = SAC_cnull
      tr%kt2       = SAC_cnull
      tr%kt3       = SAC_cnull
      tr%kt4       = SAC_cnull
      tr%kt5       = SAC_cnull
      tr%kt6       = SAC_cnull
      tr%kt7       = SAC_cnull
      tr%kt8       = SAC_cnull
      tr%kt9       = SAC_cnull
      tr%kf        = SAC_cnull
      tr%kuser0    = SAC_cnull
      tr%kuser1    = SAC_cnull
      tr%kuser2    = SAC_cnull
      tr%kcmpnm    = SAC_cnull
      tr%knetwk    = SAC_cnull
      tr%kdatrd    = SAC_cnull
      tr%kinst     = SAC_cnull
 
      deallocate (tr%trace, stat = istatus)

!  ** flag the trace as inactive      
      tr%iTraceActive = 0 

   end subroutine f90sac_destroytrace
!===============================================================================


!===============================================================================
   subroutine f90sac_newtrace(nsamp,delta,out)
!===============================================================================
      implicit none
      type (SACtrace) :: out
      integer :: nsamp ! number of samples for trace 
      real :: delta

!  ** destroy any pre-existing trace
      if (out%iTraceActive == TraceIsActive) then
         call f90sac_destroytrace(out)
      endif
      
      out%delta     = delta
      out%depmin    = SAC_rnull
      out%depmax    = SAC_rnull
      out%scale     = SAC_rnull
      out%odelta    = SAC_rnull
      out%b         = 0.0
      out%e         = real(nsamp-1)*delta
      out%o         = SAC_rnull
      out%a         = SAC_rnull
      out%internal0 = SAC_rnull
      out%t0        = SAC_rnull
      out%t1        = SAC_rnull
      out%t2        = SAC_rnull
      out%t3        = SAC_rnull
      out%t4        = SAC_rnull
      out%t5        = SAC_rnull
      out%t6        = SAC_rnull
      out%t7        = SAC_rnull
      out%t8        = SAC_rnull
      out%t9        = SAC_rnull
      out%f         = SAC_rnull
      out%resp0     = SAC_rnull
      out%resp1     = SAC_rnull
      out%resp2     = SAC_rnull
      out%resp3     = SAC_rnull
      out%resp4     = SAC_rnull
      out%resp5     = SAC_rnull
      out%resp6     = SAC_rnull
      out%resp7     = SAC_rnull
      out%resp8     = SAC_rnull
      out%resp9     = SAC_rnull
      out%stla      = SAC_rnull
      out%stlo      = SAC_rnull
      out%stel      = SAC_rnull
      out%stdp      = SAC_rnull
      out%evla      = SAC_rnull
      out%evlo      = SAC_rnull
      out%evel      = SAC_rnull
      out%evdp      = SAC_rnull
      out%mag       = SAC_rnull
      out%user0     = SAC_rnull
      out%user1     = SAC_rnull
      out%user2     = SAC_rnull
      out%user3     = SAC_rnull
      out%user4     = SAC_rnull
      out%user5     = SAC_rnull
      out%user6     = SAC_rnull
      out%user7     = SAC_rnull
      out%user8     = SAC_rnull
      out%user9     = SAC_rnull
      out%dist      = SAC_rnull
      out%az        = SAC_rnull
      out%baz       = SAC_rnull
      out%gcarc     = SAC_rnull
      out%internal1 = SAC_rnull
      out%internal2 = SAC_rnull
      out%depmen    = SAC_rnull
      out%cmpaz     = SAC_rnull
      out%cmpinc    = SAC_rnull
      out%xminimum  = SAC_rnull
      out%xmaximum  = SAC_rnull
      out%yminimum  = SAC_rnull
      out%ymaximum  = SAC_rnull
      out%unused1   = SAC_rnull
      out%unused2   = SAC_rnull
      out%unused3   = SAC_rnull
      out%unused4   = SAC_rnull
      out%unused5   = SAC_rnull
      out%unused6   = SAC_rnull
      out%unused7   = SAC_rnull

      out%nzyear    = SAC_inull
      out%nzjday    = SAC_inull
      out%nzhour    = SAC_inull
      out%nzmin     = SAC_inull
      out%nzsec     = SAC_inull
      out%nzmsec    = SAC_inull
      out%nvhdr     = 6 ! default
      out%norid     = SAC_inull
      out%nevid     = SAC_inull
      out%npts      = nsamp ! number of samples
      out%internal3 = SAC_inull
      out%nwfid     = SAC_inull
      out%nxsize    = SAC_inull
      out%nysize    = SAC_inull
      out%unused8   = SAC_inull
      out%iftype    = 1 ! default (time series file)
      out%idep      = 5 ! default
      out%iztype    = 9 ! default
      out%unused9   = SAC_inull
      out%iinst     = SAC_inull
      out%istreg    = SAC_inull
      out%ievreg    = SAC_inull
      out%ievtyp    = 5 ! default
      out%iqual     = SAC_inull
      out%isynth    = SAC_inull
      out%imagtyp   = SAC_inull
      out%imagsrc   = SAC_inull
      out%unused10  = SAC_inull
      out%unused11  = SAC_inull
      out%unused12  = SAC_inull
      out%unused13  = SAC_inull
      out%unused14  = SAC_inull
      out%unused15  = SAC_inull
      out%unused16  = SAC_inull
      out%unused17  = SAC_inull
      out%leven     = 1 ! default
      out%lpspol    = 0
      out%lovrok    = 1
      out%lcalda    = 1
      out%unused18  = 0

      out%kstnm   = SAC_cnull
      out%kevnm   = SAC_cnull
      out%khole   = SAC_cnull
      out%ko      = SAC_cnull
      out%ka      = SAC_cnull
      out%kt0     = SAC_cnull
      out%kt1     = SAC_cnull
      out%kt2     = SAC_cnull
      out%kt3     = SAC_cnull
      out%kt4     = SAC_cnull
      out%kt5     = SAC_cnull
      out%kt6     = SAC_cnull
      out%kt7     = SAC_cnull
      out%kt8     = SAC_cnull
      out%kt9     = SAC_cnull
      out%kf      = SAC_cnull
      out%kuser0  = SAC_cnull
      out%kuser1  = SAC_cnull
      out%kuser2  = SAC_cnull
      out%kcmpnm  = SAC_cnull
      out%knetwk  = SAC_cnull
      out%kdatrd  = SAC_cnull
      out%kinst   = SAC_cnull
 
!  ** allocate memory for the trace *
      allocate (out%trace(out%npts))
   
      out % trace(1:out % npts) = 0.0
      
      out % iTraceActive = TraceIsActive ;

   end subroutine f90sac_newtrace
!===============================================================================

!===============================================================================
   subroutine f90sac_clonetrace(target,clone)
!===============================================================================
!
!     Create a new, identical trace to the target
!
! input: target   SACtrace    target SAC trace object
! output:clone      SACTrace    clone SAC trace object
!
      implicit none
      type (SACtrace) :: target,clone

!  ** make a new trace object     
      call f90sac_newtrace(target % npts, target % delta, clone)
!  ** duplicate target header
      call f90sac_copytraceheader(target,clone) 
!  ** copy data   
      clone % trace(1:clone % npts) = target % trace(1:target % npts)
      
   end subroutine f90sac_clonetrace
!===============================================================================

!===============================================================================
   subroutine f90sac_copytraceheader(source,dest)
!===============================================================================
! input: source   SACtrace    source SAC trace object
! output:dest      SACTrace    destination SAC trace object
!-------------------------------------------------------------------------------
!  modifications:
!     24-03-03    J. Wookey      Modified to use F90 constructs
!-------------------------------------------------------------------------------
      implicit none
      type (SACtrace) :: source,dest

      dest%delta     =  source%delta     
      dest%depmin    =  source%depmin    
      dest%depmax    =  source%depmax    
      dest%scale     =  source%scale     
      dest%odelta    =  source%odelta    
      dest%b         =  source%b         
      dest%e         =  source%e         
      dest%o         =  source%o         
      dest%a         =  source%a         
      dest%internal0 =  source%internal0 
      dest%t0        =  source%t0        
      dest%t1        =  source%t1        
      dest%t2        =  source%t2        
      dest%t3        =  source%t3        
      dest%t4        =  source%t4        
      dest%t5        =  source%t5        
      dest%t6        =  source%t6        
      dest%t7        =  source%t7        
      dest%t8        =  source%t8        
      dest%t9        =  source%t9        
      dest%f         =  source%f         
      dest%resp0     =  source%resp0     
      dest%resp1     =  source%resp1     
      dest%resp2     =  source%resp2     
      dest%resp3     =  source%resp3     
      dest%resp4     =  source%resp4     
      dest%resp5     =  source%resp5     
      dest%resp6     =  source%resp6     
      dest%resp7     =  source%resp7     
      dest%resp8     =  source%resp8     
      dest%resp9     =  source%resp9     
      dest%stla      =  source%stla      
      dest%stlo      =  source%stlo      
      dest%stel      =  source%stel      
      dest%stdp      =  source%stdp      
      dest%evla      =  source%evla      
      dest%evlo      =  source%evlo      
      dest%evel      =  source%evel      
      dest%evdp      =  source%evdp      
      dest%mag       =  source%mag       
      dest%user0     =  source%user0     
      dest%user1     =  source%user1     
      dest%user2     =  source%user2     
      dest%user3     =  source%user3     
      dest%user4     =  source%user4     
      dest%user5     =  source%user5     
      dest%user6     =  source%user6     
      dest%user7     =  source%user7     
      dest%user8     =  source%user8     
      dest%user9     =  source%user9     
      dest%dist      =  source%dist      
      dest%az        =  source%az        
      dest%baz       =  source%baz       
      dest%gcarc     =  source%gcarc     
      dest%internal1 =  source%internal1 
      dest%internal2 =  source%internal2 
      dest%depmen    =  source%depmen    
      dest%cmpaz     =  source%cmpaz     
      dest%cmpinc    =  source%cmpinc    
      dest%xminimum  =  source%xminimum  
      dest%xmaximum  =  source%xmaximum  
      dest%yminimum  =  source%yminimum  
      dest%ymaximum  =  source%ymaximum  
      dest%unused1   =  source%unused1   
      dest%unused2   =  source%unused2   
      dest%unused3   =  source%unused3   
      dest%unused4   =  source%unused4   
      dest%unused5   =  source%unused5   
      dest%unused6   =  source%unused6   
      dest%unused7   =  source%unused7   

      dest%nzyear    =  source%nzyear    
      dest%nzjday    =  source%nzjday    
      dest%nzhour    =  source%nzhour    
      dest%nzmin     =  source%nzmin     
      dest%nzsec     =  source%nzsec     
      dest%nzmsec    =  source%nzmsec    
      dest%nvhdr     =  source%nvhdr     
      dest%norid     =  source%norid     
      dest%nevid     =  source%nevid     
      dest%npts      =  source%npts      
      dest%internal3 =  source%internal3 
      dest%nwfid     =  source%nwfid     
      dest%nxsize    =  source%nxsize    
      dest%nysize    =  source%nysize    
      dest%unused8   =  source%unused8   
      dest%iftype    =  source%iftype    
      dest%idep      =  source%idep      
      dest%iztype    =  source%iztype    
      dest%unused9   =  source%unused9   
      dest%iinst     =  source%iinst     
      dest%istreg    =  source%istreg    
      dest%ievreg    =  source%ievreg    
      dest%ievtyp    =  source%ievtyp    
      dest%iqual     =  source%iqual     
      dest%isynth    =  source%isynth    
      dest%imagtyp   =  source%imagtyp   
      dest%imagsrc   =  source%imagsrc   
      dest%unused10  =  source%unused10  
      dest%unused11  =  source%unused11  
      dest%unused12  =  source%unused12  
      dest%unused13  =  source%unused13  
      dest%unused14  =  source%unused14  
      dest%unused15  =  source%unused15  
      dest%unused16  =  source%unused16  
      dest%unused17  =  source%unused17  
      dest%leven     =  source%leven     
      dest%lpspol    =  source%lpspol    
      dest%lovrok    =  source%lovrok    
      dest%lcalda    =  source%lcalda    
      dest%unused18  =  source%unused18  

      dest%kstnm     =  source%kstnm     
      dest%kevnm     =  source%kevnm     
      dest%khole     =  source%khole     
      dest%ko        =  source%ko        
      dest%ka        =  source%ka        
      dest%kt0       =  source%kt0       
      dest%kt1       =  source%kt1       
      dest%kt2       =  source%kt2       
      dest%kt3       =  source%kt3       
      dest%kt4       =  source%kt4       
      dest%kt5       =  source%kt5       
      dest%kt6       =  source%kt6       
      dest%kt7       =  source%kt7       
      dest%kt8       =  source%kt8       
      dest%kt9       =  source%kt9       
      dest%kf        =  source%kf        
      dest%kuser0    =  source%kuser0    
      dest%kuser1    =  source%kuser1    
      dest%kuser2    =  source%kuser2    
      dest%kcmpnm    =  source%kcmpnm    
      dest%knetwk    =  source%knetwk    
      dest%kdatrd    =  source%kdatrd    
      dest%kinst     =  source%kinst     

   end subroutine f90sac_copytraceheader
!===============================================================================

!===============================================================================
   subroutine f90sac_readtrace(fname,out)
!===============================================================================
!
!  read a SAC time-series object from a file   
!
      implicit none
      integer :: lu
      character (len=4) :: tmp_char1,tmp_char2,tmp_char3,tmp_char4
      character (len=50) :: fname
      integer :: i, istatus
      
      real(real4) :: sacrh(70) ! SAC floating point header
      integer(int4) :: sacih(40) ! SAC floating point header

      type (SACtrace) :: out
      logical UNITOK, UNITOP
            
!  ** destroy any pre-existing trace
      if (out%iTraceActive == TraceIsActive) then
         call f90sac_destroytrace(out)
      endif
            
!  ** Open the file for reading
      open (unit=f90sac_iounit, file=fname,form='unformatted', &
            access='direct',recl=f90sac_32bit_record_length, &
            status='old',err=900)
      
      lu = f90sac_iounit
      
!
!  ** READ IN THE SAC HEADER
!
      do i=1,70
         read(lu,rec=i,err=901) sacrh(i)
      enddo

      if (f90sac_force_byteswap) then
         call f90sac_real32_byteswap(sacrh,70)
      endif

!  ** populate the structure
      out%delta     = sacrh(001)          
      out%depmin    = sacrh(002) 
      out%depmax    = sacrh(003) 
      out%scale     = sacrh(004) 
      out%odelta    = sacrh(005) 
      out%b         = sacrh(006) 
      out%e         = sacrh(007) 
      out%o         = sacrh(008) 
      out%a         = sacrh(009) 
      out%internal0 = sacrh(010) 
      out%t0        = sacrh(011) 
      out%t1        = sacrh(012) 
      out%t2        = sacrh(013) 
      out%t3        = sacrh(014) 
      out%t4        = sacrh(015) 
      out%t5        = sacrh(016) 
      out%t6        = sacrh(017) 
      out%t7        = sacrh(018) 
      out%t8        = sacrh(019) 
      out%t9        = sacrh(020) 
      out%f         = sacrh(021) 
      out%resp0     = sacrh(022) 
      out%resp1     = sacrh(023) 
      out%resp2     = sacrh(024) 
      out%resp3     = sacrh(025) 
      out%resp4     = sacrh(026) 
      out%resp5     = sacrh(027) 
      out%resp6     = sacrh(028) 
      out%resp7     = sacrh(029) 
      out%resp8     = sacrh(030) 
      out%resp9     = sacrh(031) 
      out%stla      = sacrh(032) 
      out%stlo      = sacrh(033) 
      out%stel      = sacrh(034) 
      out%stdp      = sacrh(035) 
      out%evla      = sacrh(036) 
      out%evlo      = sacrh(037) 
      out%evel      = sacrh(038) 
      out%evdp      = sacrh(039) 
      out%mag       = sacrh(040) 
      out%user0     = sacrh(041) 
      out%user1     = sacrh(042) 
      out%user2     = sacrh(043) 
      out%user3     = sacrh(044) 
      out%user4     = sacrh(045) 
      out%user5     = sacrh(046) 
      out%user6     = sacrh(047) 
      out%user7     = sacrh(048) 
      out%user8     = sacrh(049) 
      out%user9     = sacrh(050) 
      out%dist      = sacrh(051) 
      out%az        = sacrh(052) 
      out%baz       = sacrh(053) 
      out%gcarc     = sacrh(054) 
      out%internal1 = sacrh(055) 
      out%internal2 = sacrh(056) 
      out%depmen    = sacrh(057) 
      out%cmpaz     = sacrh(058) 
      out%cmpinc    = sacrh(059) 
      out%xminimum  = sacrh(060) 
      out%xmaximum  = sacrh(061) 
      out%yminimum  = sacrh(062) 
      out%ymaximum  = sacrh(063) 
      out%unused1   = sacrh(064) 
      out%unused2   = sacrh(065) 
      out%unused3   = sacrh(066) 
      out%unused4   = sacrh(067) 
      out%unused5   = sacrh(068) 
      out%unused6   = sacrh(069) 
      out%unused7   = sacrh(070) 
      
      do i=1,40
         read(lu,rec=i+70,err=901) sacih(i)
      enddo

      if (f90sac_force_byteswap) then
         call f90sac_int32_byteswap(sacih,40)
      endif

      out%nzyear    = sacih(001)
      out%nzjday    = sacih(002)
      out%nzhour    = sacih(003)
      out%nzmin     = sacih(004)
      out%nzsec     = sacih(005)
      out%nzmsec    = sacih(006)
      out%nvhdr     = sacih(007)
      out%norid     = sacih(008)
      out%nevid     = sacih(009)
      out%npts      = sacih(010)
      out%internal3 = sacih(011)
      out%nwfid     = sacih(012)
      out%nxsize    = sacih(013)
      out%nysize    = sacih(014)
      out%unused8   = sacih(015)
      out%iftype    = sacih(016)
      out%idep      = sacih(017)
      out%iztype    = sacih(018)
      out%unused9   = sacih(019)
      out%iinst     = sacih(020)
      out%istreg    = sacih(021)
      out%ievreg    = sacih(022)
      out%ievtyp    = sacih(023)
      out%iqual     = sacih(024)
      out%isynth    = sacih(025)
      out%imagtyp   = sacih(026)
      out%imagsrc   = sacih(027)
      out%unused10  = sacih(028)
      out%unused11  = sacih(029)
      out%unused12  = sacih(030)
      out%unused13  = sacih(031)
      out%unused14  = sacih(032)
      out%unused15  = sacih(033)
      out%unused16  = sacih(034)
      out%unused17  = sacih(035)
      out%leven     = sacih(036)
      out%lpspol    = sacih(037)
      out%lovrok    = sacih(038)
      out%lcalda    = sacih(039)
      out%unused18  = sacih(040)
      
      read(lu,rec=111,err=901)tmp_char1
      read(lu,rec=112,err=901)tmp_char2
      out%kstnm = tmp_char1 // tmp_char2
      
      read(lu,rec=113,err=901)tmp_char1
      read(lu,rec=114,err=901)tmp_char2
      read(lu,rec=115,err=901)tmp_char3
      read(lu,rec=116,err=901)tmp_char4
      out%kevnm = tmp_char1//tmp_char2//tmp_char3//tmp_char4
      
      read(lu,rec=117,err=901)tmp_char1
      read(lu,rec=118,err=901)tmp_char2
      out%khole = tmp_char1 // tmp_char2
      
      read(lu,rec=119,err=901)tmp_char1
      read(lu,rec=120,err=901)tmp_char2
      out %ko = tmp_char1 // tmp_char2
      
      read(lu,rec=121,err=901)tmp_char1
      read(lu,rec=122,err=901)tmp_char2
      out %ka = tmp_char1 // tmp_char2
      
      read(lu,rec=123,err=901)tmp_char1
      read(lu,rec=124,err=901)tmp_char2
      out %kt0 = tmp_char1 // tmp_char2
      
      read(lu,rec=125,err=901)tmp_char1
      read(lu,rec=126,err=901)tmp_char2
      out %kt1 = tmp_char1 // tmp_char2
      
      read(lu,rec=127,err=901)tmp_char1
      read(lu,rec=128,err=901)tmp_char2
      out %kt2 = tmp_char1 // tmp_char2
      
      read(lu,rec=129,err=901)tmp_char1
      read(lu,rec=130,err=901)tmp_char2
      out %kt3 = tmp_char1 // tmp_char2
      
      read(lu,rec=131,err=901)tmp_char1
      read(lu,rec=132,err=901)tmp_char2
      out %kt4 = tmp_char1 // tmp_char2
      
      read(lu,rec=133,err=901)tmp_char1
      read(lu,rec=134,err=901)tmp_char2
      out %kt5 = tmp_char1 // tmp_char2
      
      read(lu,rec=135,err=901)tmp_char1
      read(lu,rec=136,err=901)tmp_char2
      out %kt6 = tmp_char1 // tmp_char2
      
      read(lu,rec=137,err=901)tmp_char1
      read(lu,rec=138,err=901)tmp_char2
      out %kt7 = tmp_char1 // tmp_char2
      
      read(lu,rec=139,err=901)tmp_char1
      read(lu,rec=140,err=901)tmp_char2
      out %kt8 = tmp_char1 // tmp_char2
      
      read(lu,rec=141,err=901)tmp_char1
      read(lu,rec=142,err=901)tmp_char2
      out %kt9 = tmp_char1 // tmp_char2
      
      read(lu,rec=143,err=901)tmp_char1
      read(lu,rec=144,err=901)tmp_char2
      out %kf = tmp_char1 // tmp_char2
      
      read(lu,rec=145,err=901)tmp_char1
      read(lu,rec=146,err=901)tmp_char2
      out %kuser0 = tmp_char1 // tmp_char2
      
      read(lu,rec=147,err=901)tmp_char1
      read(lu,rec=148,err=901)tmp_char2
      out %kuser1 = tmp_char1 // tmp_char2
      
      read(lu,rec=149,err=901)tmp_char1
      read(lu,rec=150,err=901)tmp_char2
      out %kuser2 = tmp_char1 // tmp_char2
      
      read(lu,rec=151,err=901)tmp_char1
      read(lu,rec=152,err=901)tmp_char2
      out %kcmpnm = tmp_char1 // tmp_char2
      
      read(lu,rec=153,err=901)tmp_char1
      read(lu,rec=154,err=901)tmp_char2
      out %knetwk = tmp_char1 // tmp_char2
      
      read(lu,rec=155,err=901)tmp_char1
      read(lu,rec=156,err=901)tmp_char2
      out %kdatrd = tmp_char1 // tmp_char2
      
      read(lu,rec=157,err=901)tmp_char1
      read(lu,rec=158,err=901)tmp_char2
      out %kinst = tmp_char1 // tmp_char2
!
!  ** DONE WITH THE HEADER, READ THE TRACE ITSELF
!

!  ** check that the nvhdr header is sensible 
      if (out % nvhdr < 0 .or. out % nvhdr > 10 ) then
         write(0,'(a)') &
          'F90SAC_READTRACE: Error: NVHDR is not sensible, check f90sac_force_byteswap?'
         stop
      endif 

!  ** allocate memory for the trace *
      allocate (out%trace(out%npts))
   
      do i=1,out % npts
            read(lu,rec=158+i) out%trace(i)
      enddo
      close(lu)

!  ** if required, byteswap the trace
      if (f90sac_force_byteswap) then
         call f90sac_real32_byteswap(out%trace,out%npts)
      endif


!  ** set the trace as active      
      out%iTraceActive = TraceIsActive

      
      return
900   write(0,'(a,a)') 'F90SAC_READTRACE: Error: File does not exist: ', &
         trim(fname)
      stop
901   write(0,'(a,a)') 'F90SAC_READTRACE: Error: Bad read in file: ', &
         trim(fname)
      stop     
   end subroutine f90sac_readtrace
!===============================================================================

!===============================================================================
   subroutine f90sac_real32_byteswap(x,n)
!===============================================================================
!
!  byteswap an array of 32 bit floats
!
      integer :: n,i
      real(real4) :: x(n),xx
      integer(int4) :: itmp,itmp2
      
      do i=1,n
         itmp = transfer(x(i),itmp)

         call mvbits( itmp, 24, 8, itmp2, 0  )
         call mvbits( itmp, 16, 8, itmp2, 8  )
         call mvbits( itmp,  8, 8, itmp2, 16 )
         call mvbits( itmp,  0, 8, itmp2, 24 )
         
         x(i) = transfer(itmp2,xx)
      enddo
    
      return     
   end subroutine f90sac_real32_byteswap
!===============================================================================

!===============================================================================
   subroutine f90sac_int32_byteswap(nx,n)
!===============================================================================
!
!  byteswap an array of 32 bit integers
!
      integer :: n,i
      integer(int4) :: nx(n)
      integer(int4) :: itmp,itmp2
      
      do i=1,n
         itmp = nx(i)

         call mvbits( itmp, 24, 8, itmp2, 0  )
         call mvbits( itmp, 16, 8, itmp2, 8  )
         call mvbits( itmp,  8, 8, itmp2, 16 )
         call mvbits( itmp,  0, 8, itmp2, 24 )
         
         nx(i) = itmp2
      enddo
    
      return     
   end subroutine f90sac_int32_byteswap
!===============================================================================

!===============================================================================
   subroutine f90sac_readheader(fname,out)
!===============================================================================
!
!  read a SAC header from a file. This is a trace object but with a single null    
!  value as the trace. 
!
      implicit none
      integer :: lu
      character (len=4) :: tmp_char1,tmp_char2,tmp_char3,tmp_char4
      character (len=50) :: fname
      integer :: i, istatus
      type (SACtrace) :: out
      logical UNITOK, UNITOP
            
      real(real4) :: sacrh(70) ! SAC floating point header
      integer(int4) :: sacih(40) ! SAC floating point header
            
!  ** destroy any pre-existing trace
      if (out%iTraceActive == TraceIsActive) then
         call f90sac_destroytrace(out)
      endif
            
!  ** Open the file for reading
      open (unit=f90sac_iounit, file=fname,form='unformatted', &
            access='direct',recl=f90sac_32bit_record_length, &
            status='old',err=900)
      
      lu = f90sac_iounit
      
!
!  ** READ IN THE SAC HEADER
!
      do i=1,70
         read(lu,rec=i,err=901) sacrh(i)
      enddo

      if (f90sac_force_byteswap) then
         call f90sac_real32_byteswap(sacrh,70)
      endif

!  ** populate the structure
      out%delta     = sacrh(001)          
      out%depmin    = sacrh(002) 
      out%depmax    = sacrh(003) 
      out%scale     = sacrh(004) 
      out%odelta    = sacrh(005) 
      out%b         = sacrh(006) 
      out%e         = sacrh(007) 
      out%o         = sacrh(008) 
      out%a         = sacrh(009) 
      out%internal0 = sacrh(010) 
      out%t0        = sacrh(011) 
      out%t1        = sacrh(012) 
      out%t2        = sacrh(013) 
      out%t3        = sacrh(014) 
      out%t4        = sacrh(015) 
      out%t5        = sacrh(016) 
      out%t6        = sacrh(017) 
      out%t7        = sacrh(018) 
      out%t8        = sacrh(019) 
      out%t9        = sacrh(020) 
      out%f         = sacrh(021) 
      out%resp0     = sacrh(022) 
      out%resp1     = sacrh(023) 
      out%resp2     = sacrh(024) 
      out%resp3     = sacrh(025) 
      out%resp4     = sacrh(026) 
      out%resp5     = sacrh(027) 
      out%resp6     = sacrh(028) 
      out%resp7     = sacrh(029) 
      out%resp8     = sacrh(030) 
      out%resp9     = sacrh(031) 
      out%stla      = sacrh(032) 
      out%stlo      = sacrh(033) 
      out%stel      = sacrh(034) 
      out%stdp      = sacrh(035) 
      out%evla      = sacrh(036) 
      out%evlo      = sacrh(037) 
      out%evel      = sacrh(038) 
      out%evdp      = sacrh(039) 
      out%mag       = sacrh(040) 
      out%user0     = sacrh(041) 
      out%user1     = sacrh(042) 
      out%user2     = sacrh(043) 
      out%user3     = sacrh(044) 
      out%user4     = sacrh(045) 
      out%user5     = sacrh(046) 
      out%user6     = sacrh(047) 
      out%user7     = sacrh(048) 
      out%user8     = sacrh(049) 
      out%user9     = sacrh(050) 
      out%dist      = sacrh(051) 
      out%az        = sacrh(052) 
      out%baz       = sacrh(053) 
      out%gcarc     = sacrh(054) 
      out%internal1 = sacrh(055) 
      out%internal2 = sacrh(056) 
      out%depmen    = sacrh(057) 
      out%cmpaz     = sacrh(058) 
      out%cmpinc    = sacrh(059) 
      out%xminimum  = sacrh(060) 
      out%xmaximum  = sacrh(061) 
      out%yminimum  = sacrh(062) 
      out%ymaximum  = sacrh(063) 
      out%unused1   = sacrh(064) 
      out%unused2   = sacrh(065) 
      out%unused3   = sacrh(066) 
      out%unused4   = sacrh(067) 
      out%unused5   = sacrh(068) 
      out%unused6   = sacrh(069) 
      out%unused7   = sacrh(070) 
      
      do i=1,40
         read(lu,rec=i+70,err=901) sacih(i)
      enddo

      if (f90sac_force_byteswap) then
         call f90sac_int32_byteswap(sacih,40)
      endif

      out%nzyear    = sacih(001)
      out%nzjday    = sacih(002)
      out%nzhour    = sacih(003)
      out%nzmin     = sacih(004)
      out%nzsec     = sacih(005)
      out%nzmsec    = sacih(006)
      out%nvhdr     = sacih(007)
      out%norid     = sacih(008)
      out%nevid     = sacih(009)
      out%npts      = sacih(010)
      out%internal3 = sacih(011)
      out%nwfid     = sacih(012)
      out%nxsize    = sacih(013)
      out%nysize    = sacih(014)
      out%unused8   = sacih(015)
      out%iftype    = sacih(016)
      out%idep      = sacih(017)
      out%iztype    = sacih(018)
      out%unused9   = sacih(019)
      out%iinst     = sacih(020)
      out%istreg    = sacih(021)
      out%ievreg    = sacih(022)
      out%ievtyp    = sacih(023)
      out%iqual     = sacih(024)
      out%isynth    = sacih(025)
      out%imagtyp   = sacih(026)
      out%imagsrc   = sacih(027)
      out%unused10  = sacih(028)
      out%unused11  = sacih(029)
      out%unused12  = sacih(030)
      out%unused13  = sacih(031)
      out%unused14  = sacih(032)
      out%unused15  = sacih(033)
      out%unused16  = sacih(034)
      out%unused17  = sacih(035)
      out%leven     = sacih(036)
      out%lpspol    = sacih(037)
      out%lovrok    = sacih(038)
      out%lcalda    = sacih(039)
      out%unused18  = sacih(040)
      
      read(lu,rec=111,err=901)tmp_char1
      read(lu,rec=112,err=901)tmp_char2
      out%kstnm = tmp_char1 // tmp_char2
      
      read(lu,rec=113,err=901)tmp_char1
      read(lu,rec=114,err=901)tmp_char2
      read(lu,rec=115,err=901)tmp_char3
      read(lu,rec=116,err=901)tmp_char4
      out%kevnm = tmp_char1//tmp_char2//tmp_char3//tmp_char4
      
      read(lu,rec=117,err=901)tmp_char1
      read(lu,rec=118,err=901)tmp_char2
      out%khole = tmp_char1 // tmp_char2
      
      read(lu,rec=119,err=901)tmp_char1
      read(lu,rec=120,err=901)tmp_char2
      out %ko = tmp_char1 // tmp_char2
      
      read(lu,rec=121,err=901)tmp_char1
      read(lu,rec=122,err=901)tmp_char2
      out %ka = tmp_char1 // tmp_char2
      
      read(lu,rec=123,err=901)tmp_char1
      read(lu,rec=124,err=901)tmp_char2
      out %kt0 = tmp_char1 // tmp_char2
      
      read(lu,rec=125,err=901)tmp_char1
      read(lu,rec=126,err=901)tmp_char2
      out %kt1 = tmp_char1 // tmp_char2
      
      read(lu,rec=127,err=901)tmp_char1
      read(lu,rec=128,err=901)tmp_char2
      out %kt2 = tmp_char1 // tmp_char2
      
      read(lu,rec=129,err=901)tmp_char1
      read(lu,rec=130,err=901)tmp_char2
      out %kt3 = tmp_char1 // tmp_char2
      
      read(lu,rec=131,err=901)tmp_char1
      read(lu,rec=132,err=901)tmp_char2
      out %kt4 = tmp_char1 // tmp_char2
      
      read(lu,rec=133,err=901)tmp_char1
      read(lu,rec=134,err=901)tmp_char2
      out %kt5 = tmp_char1 // tmp_char2
      
      read(lu,rec=135,err=901)tmp_char1
      read(lu,rec=136,err=901)tmp_char2
      out %kt6 = tmp_char1 // tmp_char2
      
      read(lu,rec=137,err=901)tmp_char1
      read(lu,rec=138,err=901)tmp_char2
      out %kt7 = tmp_char1 // tmp_char2
      
      read(lu,rec=139,err=901)tmp_char1
      read(lu,rec=140,err=901)tmp_char2
      out %kt8 = tmp_char1 // tmp_char2
      
      read(lu,rec=141,err=901)tmp_char1
      read(lu,rec=142,err=901)tmp_char2
      out %kt9 = tmp_char1 // tmp_char2
      
      read(lu,rec=143,err=901)tmp_char1
      read(lu,rec=144,err=901)tmp_char2
      out %kf = tmp_char1 // tmp_char2
      
      read(lu,rec=145,err=901)tmp_char1
      read(lu,rec=146,err=901)tmp_char2
      out %kuser0 = tmp_char1 // tmp_char2
      
      read(lu,rec=147,err=901)tmp_char1
      read(lu,rec=148,err=901)tmp_char2
      out %kuser1 = tmp_char1 // tmp_char2
      
      read(lu,rec=149,err=901)tmp_char1
      read(lu,rec=150,err=901)tmp_char2
      out %kuser2 = tmp_char1 // tmp_char2
      
      read(lu,rec=151,err=901)tmp_char1
      read(lu,rec=152,err=901)tmp_char2
      out %kcmpnm = tmp_char1 // tmp_char2
      
      read(lu,rec=153,err=901)tmp_char1
      read(lu,rec=154,err=901)tmp_char2
      out %knetwk = tmp_char1 // tmp_char2
      
      read(lu,rec=155,err=901)tmp_char1
      read(lu,rec=156,err=901)tmp_char2
      out %kdatrd = tmp_char1 // tmp_char2
      
      read(lu,rec=157,err=901)tmp_char1
      read(lu,rec=158,err=901)tmp_char2
      out %kinst = tmp_char1 // tmp_char2
!
!  ** DONE WITH THE HEADER, READ THE TRACE ITSELF
!

!  ** check that the nvhdr header is sensible 
      if (out % nvhdr < 0 .or. out % nvhdr > 10 ) then
         write(0,'(a)') &
          'F90SAC_READTRACE: Error: NVHDR is not sensible, byteswap required?'
         stop
      endif 

 
!  ** allocate memory for the trace *
      allocate (out%trace(1))
   
      out%trace(1) = SAC_rnull ;

      close(lu)

!  ** set the trace as active      
      out%iTraceActive = TraceIsActive

      
      return
900   write(0,'(a,a)') 'F90SAC_READTRACE: Error: File does not exist: ', &
         trim(fname)
      stop
901   write(0,'(a,a)') 'F90SAC_READTRACE: Error: Bad read in file: ', &
         trim(fname)
      stop     
   end subroutine f90sac_readheader
!===============================================================================


!===============================================================================
   subroutine f90sac_writetrace(fname,out)
!===============================================================================
!
!  write a SAC time-series object to a file   
!
      implicit none
      integer :: lu
      character (len=4) :: tmp_char1,tmp_char2,tmp_char3,tmp_char4
      character (len=50):: fname
      integer :: i, istatus
      type (SACtrace) :: out
      logical UNITOK, UNITOP

      real(real4) :: sacrh(70) ! SAC floating point header
      integer(int4) :: sacih(40) ! SAC floating point header


!  ** Open the file for writing
      open (unit=f90sac_iounit, file=fname,form='unformatted', &
            access='direct',recl=f90sac_32bit_record_length, &
            status='unknown')

      lu = f90sac_iounit

!
!  ** WRITE OUT THE SAC HEADER
!
!  ** populate the structure
      sacrh(001)  = out%delta             
      sacrh(002)  = out%depmin   
      sacrh(003)  = out%depmax   
      sacrh(004)  = out%scale    
      sacrh(005)  = out%odelta   
      sacrh(006)  = out%b        
      sacrh(007)  = out%e        
      sacrh(008)  = out%o        
      sacrh(009)  = out%a        
      sacrh(010)  = out%internal0
      sacrh(011)  = out%t0       
      sacrh(012)  = out%t1       
      sacrh(013)  = out%t2       
      sacrh(014)  = out%t3       
      sacrh(015)  = out%t4       
      sacrh(016)  = out%t5       
      sacrh(017)  = out%t6       
      sacrh(018)  = out%t7       
      sacrh(019)  = out%t8       
      sacrh(020)  = out%t9       
      sacrh(021)  = out%f        
      sacrh(022)  = out%resp0    
      sacrh(023)  = out%resp1    
      sacrh(024)  = out%resp2    
      sacrh(025)  = out%resp3    
      sacrh(026)  = out%resp4    
      sacrh(027)  = out%resp5    
      sacrh(028)  = out%resp6    
      sacrh(029)  = out%resp7    
      sacrh(030)  = out%resp8    
      sacrh(031)  = out%resp9    
      sacrh(032)  = out%stla     
      sacrh(033)  = out%stlo     
      sacrh(034)  = out%stel     
      sacrh(035)  = out%stdp     
      sacrh(036)  = out%evla     
      sacrh(037)  = out%evlo     
      sacrh(038)  = out%evel     
      sacrh(039)  = out%evdp     
      sacrh(040)  = out%mag      
      sacrh(041)  = out%user0    
      sacrh(042)  = out%user1    
      sacrh(043)  = out%user2    
      sacrh(044)  = out%user3    
      sacrh(045)  = out%user4    
      sacrh(046)  = out%user5    
      sacrh(047)  = out%user6    
      sacrh(048)  = out%user7    
      sacrh(049)  = out%user8    
      sacrh(050)  = out%user9    
      sacrh(051)  = out%dist     
      sacrh(052)  = out%az       
      sacrh(053)  = out%baz      
      sacrh(054)  = out%gcarc    
      sacrh(055)  = out%internal1
      sacrh(056)  = out%internal2
      sacrh(057)  = out%depmen   
      sacrh(058)  = out%cmpaz    
      sacrh(059)  = out%cmpinc   
      sacrh(060)  = out%xminimum 
      sacrh(061)  = out%xmaximum 
      sacrh(062)  = out%yminimum 
      sacrh(063)  = out%ymaximum 
      sacrh(064)  = out%unused1  
      sacrh(065)  = out%unused2  
      sacrh(066)  = out%unused3  
      sacrh(067)  = out%unused4  
      sacrh(068)  = out%unused5  
      sacrh(069)  = out%unused6  
      sacrh(070)  = out%unused7  

      if (f90sac_force_byteswap) then
         call f90sac_real32_byteswap(sacrh,70)
      endif

      do i=1,70
         write(lu,rec=i) sacrh(i)
      enddo

      sacih(001) = out%nzyear   
      sacih(002) = out%nzjday   
      sacih(003) = out%nzhour   
      sacih(004) = out%nzmin    
      sacih(005) = out%nzsec    
      sacih(006) = out%nzmsec   
      sacih(007) = out%nvhdr    
      sacih(008) = out%norid    
      sacih(009) = out%nevid    
      sacih(010) = out%npts     
      sacih(011) = out%internal3
      sacih(012) = out%nwfid    
      sacih(013) = out%nxsize   
      sacih(014) = out%nysize   
      sacih(015) = out%unused8  
      sacih(016) = out%iftype   
      sacih(017) = out%idep     
      sacih(018) = out%iztype   
      sacih(019) = out%unused9  
      sacih(020) = out%iinst    
      sacih(021) = out%istreg   
      sacih(022) = out%ievreg   
      sacih(023) = out%ievtyp   
      sacih(024) = out%iqual    
      sacih(025) = out%isynth   
      sacih(026) = out%imagtyp  
      sacih(027) = out%imagsrc  
      sacih(028) = out%unused10 
      sacih(029) = out%unused11 
      sacih(030) = out%unused12 
      sacih(031) = out%unused13 
      sacih(032) = out%unused14 
      sacih(033) = out%unused15 
      sacih(034) = out%unused16 
      sacih(035) = out%unused17 
      sacih(036) = out%leven    
      sacih(037) = out%lpspol   
      sacih(038) = out%lovrok   
      sacih(039) = out%lcalda   
      sacih(040) = out%unused18 

      if (f90sac_force_byteswap) then
         call f90sac_int32_byteswap(sacih,40)
      endif

      do i=1,40
         write(lu,rec=i+70) sacih(i)
      enddo
      
      write(lu,rec=111) out%kstnm(1:4)
      write(lu,rec=112) out%kstnm(5:8)
      write(lu,rec=113) out%kevnm(1:4)
      write(lu,rec=114) out%kevnm(5:8)
      write(lu,rec=115) out%kevnm(9:12)
      write(lu,rec=116) out%kevnm(13:16)
      write(lu,rec=117) out%khole(1:4)
      write(lu,rec=118) out%khole(5:8)
      write(lu,rec=119) out%ko(1:4)
      write(lu,rec=120) out%ko(5:8)
      write(lu,rec=121) out%ka(1:4)
      write(lu,rec=122) out%ka(5:8)
      write(lu,rec=123) out%kt0(1:4)
      write(lu,rec=124) out%kt0(5:8)
      write(lu,rec=125) out%kt1(1:4)
      write(lu,rec=126) out%kt1(5:8)
      write(lu,rec=127) out%kt2(1:4)
      write(lu,rec=128) out%kt2(5:8)
      write(lu,rec=129) out%kt3(1:4)
      write(lu,rec=130) out%kt3(5:8)
      write(lu,rec=131) out%kt4(1:4)
      write(lu,rec=132) out%kt4(5:8)
      write(lu,rec=133) out%kt5(1:4)
      write(lu,rec=134) out%kt5(5:8)
      write(lu,rec=135) out%kt6(1:4)
      write(lu,rec=136) out%kt6(5:8)
      write(lu,rec=137) out%kt7(1:4)
      write(lu,rec=138) out%kt7(5:8)
      write(lu,rec=139) out%kt8(1:4)
      write(lu,rec=140) out%kt8(5:8)
      write(lu,rec=141) out%kt9(1:4)
      write(lu,rec=142) out%kt9(5:8)
      write(lu,rec=143) out%kf(1:4)
      write(lu,rec=144) out%kf(5:8)
      write(lu,rec=145) out%kuser0(1:4)
      write(lu,rec=146) out%kuser0(5:8)
      write(lu,rec=147) out%kuser1(1:4)
      write(lu,rec=148) out%kuser1(5:8)
      write(lu,rec=149) out%kuser2(1:4)
      write(lu,rec=150) out%kuser2(5:8)
      write(lu,rec=151) out%kcmpnm(1:4)
      write(lu,rec=152) out%kcmpnm(5:8)
      write(lu,rec=153) out%knetwk(1:4)
      write(lu,rec=154) out%knetwk(5:8)
      write(lu,rec=155) out%kdatrd(1:4)
      write(lu,rec=156) out%kdatrd(5:8)
      write(lu,rec=157) out%kinst(1:4)
      write(lu,rec=158) out%kinst(5:8)

!  ** if required, byteswap the trace
      if (f90sac_force_byteswap) then
         call f90sac_real32_byteswap(out%trace,out%npts)
      endif

!  ** Output the trace   
       do i=1,out % npts
            write(lu,rec=158+i) out%trace(i)
      enddo

 
      close(lu)
   end subroutine f90sac_writetrace
!===============================================================================

!===============================================================================
   subroutine f90sac_writeheader(fname,out)
!===============================================================================
!
!  write a SAC time-series object to a file (header only)  
!
      implicit none
      integer :: lu
      character (len=4) :: tmp_char1,tmp_char2,tmp_char3,tmp_char4
      character (len=50):: fname
      integer :: i, istatus
      type (SACtrace) :: out
      logical UNITOK, UNITOP
      real(real4) :: sacrh(70) ! SAC floating point header
      integer(int4) :: sacih(40) ! SAC floating point header


!  ** Open the file for writing
      open (unit=f90sac_iounit, file=fname,form='unformatted', &
            access='direct',recl=f90sac_32bit_record_length, &
            status='unknown')

      lu = f90sac_iounit

!
!  ** WRITE OUT THE SAC HEADER
!
      sacrh(001)  = out%delta             
      sacrh(002)  = out%depmin   
      sacrh(003)  = out%depmax   
      sacrh(004)  = out%scale    
      sacrh(005)  = out%odelta   
      sacrh(006)  = out%b        
      sacrh(007)  = out%e        
      sacrh(008)  = out%o        
      sacrh(009)  = out%a        
      sacrh(010)  = out%internal0
      sacrh(011)  = out%t0       
      sacrh(012)  = out%t1       
      sacrh(013)  = out%t2       
      sacrh(014)  = out%t3       
      sacrh(015)  = out%t4       
      sacrh(016)  = out%t5       
      sacrh(017)  = out%t6       
      sacrh(018)  = out%t7       
      sacrh(019)  = out%t8       
      sacrh(020)  = out%t9       
      sacrh(021)  = out%f        
      sacrh(022)  = out%resp0    
      sacrh(023)  = out%resp1    
      sacrh(024)  = out%resp2    
      sacrh(025)  = out%resp3    
      sacrh(026)  = out%resp4    
      sacrh(027)  = out%resp5    
      sacrh(028)  = out%resp6    
      sacrh(029)  = out%resp7    
      sacrh(030)  = out%resp8    
      sacrh(031)  = out%resp9    
      sacrh(032)  = out%stla     
      sacrh(033)  = out%stlo     
      sacrh(034)  = out%stel     
      sacrh(035)  = out%stdp     
      sacrh(036)  = out%evla     
      sacrh(037)  = out%evlo     
      sacrh(038)  = out%evel     
      sacrh(039)  = out%evdp     
      sacrh(040)  = out%mag      
      sacrh(041)  = out%user0    
      sacrh(042)  = out%user1    
      sacrh(043)  = out%user2    
      sacrh(044)  = out%user3    
      sacrh(045)  = out%user4    
      sacrh(046)  = out%user5    
      sacrh(047)  = out%user6    
      sacrh(048)  = out%user7    
      sacrh(049)  = out%user8    
      sacrh(050)  = out%user9    
      sacrh(051)  = out%dist     
      sacrh(052)  = out%az       
      sacrh(053)  = out%baz      
      sacrh(054)  = out%gcarc    
      sacrh(055)  = out%internal1
      sacrh(056)  = out%internal2
      sacrh(057)  = out%depmen   
      sacrh(058)  = out%cmpaz    
      sacrh(059)  = out%cmpinc   
      sacrh(060)  = out%xminimum 
      sacrh(061)  = out%xmaximum 
      sacrh(062)  = out%yminimum 
      sacrh(063)  = out%ymaximum 
      sacrh(064)  = out%unused1  
      sacrh(065)  = out%unused2  
      sacrh(066)  = out%unused3  
      sacrh(067)  = out%unused4  
      sacrh(068)  = out%unused5  
      sacrh(069)  = out%unused6  
      sacrh(070)  = out%unused7  

      if (f90sac_force_byteswap) then
         call f90sac_real32_byteswap(sacrh,70)
      endif

      do i=1,70
         write(lu,rec=i) sacrh(i)
      enddo

      sacih(001) = out%nzyear   
      sacih(002) = out%nzjday   
      sacih(003) = out%nzhour   
      sacih(004) = out%nzmin    
      sacih(005) = out%nzsec    
      sacih(006) = out%nzmsec   
      sacih(007) = out%nvhdr    
      sacih(008) = out%norid    
      sacih(009) = out%nevid    
      sacih(010) = out%npts     
      sacih(011) = out%internal3
      sacih(012) = out%nwfid    
      sacih(013) = out%nxsize   
      sacih(014) = out%nysize   
      sacih(015) = out%unused8  
      sacih(016) = out%iftype   
      sacih(017) = out%idep     
      sacih(018) = out%iztype   
      sacih(019) = out%unused9  
      sacih(020) = out%iinst    
      sacih(021) = out%istreg   
      sacih(022) = out%ievreg   
      sacih(023) = out%ievtyp   
      sacih(024) = out%iqual    
      sacih(025) = out%isynth   
      sacih(026) = out%imagtyp  
      sacih(027) = out%imagsrc  
      sacih(028) = out%unused10 
      sacih(029) = out%unused11 
      sacih(030) = out%unused12 
      sacih(031) = out%unused13 
      sacih(032) = out%unused14 
      sacih(033) = out%unused15 
      sacih(034) = out%unused16 
      sacih(035) = out%unused17 
      sacih(036) = out%leven    
      sacih(037) = out%lpspol   
      sacih(038) = out%lovrok   
      sacih(039) = out%lcalda   
      sacih(040) = out%unused18 

      if (f90sac_force_byteswap) then
         call f90sac_int32_byteswap(sacih,40)
      endif

      do i=1,40
         write(lu,rec=i+70) sacih(i)
      enddo
      
      write(lu,rec=111) out%kstnm(1:4)
      write(lu,rec=112) out%kstnm(5:8)
      write(lu,rec=113) out%kevnm(1:4)
      write(lu,rec=114) out%kevnm(5:8)
      write(lu,rec=115) out%kevnm(9:12)
      write(lu,rec=116) out%kevnm(13:16)
      write(lu,rec=117) out%khole(1:4)
      write(lu,rec=118) out%khole(5:8)
      write(lu,rec=119) out%ko(1:4)
      write(lu,rec=120) out%ko(5:8)
      write(lu,rec=121) out%ka(1:4)
      write(lu,rec=122) out%ka(5:8)
      write(lu,rec=123) out%kt0(1:4)
      write(lu,rec=124) out%kt0(5:8)
      write(lu,rec=125) out%kt1(1:4)
      write(lu,rec=126) out%kt1(5:8)
      write(lu,rec=127) out%kt2(1:4)
      write(lu,rec=128) out%kt2(5:8)
      write(lu,rec=129) out%kt3(1:4)
      write(lu,rec=130) out%kt3(5:8)
      write(lu,rec=131) out%kt4(1:4)
      write(lu,rec=132) out%kt4(5:8)
      write(lu,rec=133) out%kt5(1:4)
      write(lu,rec=134) out%kt5(5:8)
      write(lu,rec=135) out%kt6(1:4)
      write(lu,rec=136) out%kt6(5:8)
      write(lu,rec=137) out%kt7(1:4)
      write(lu,rec=138) out%kt7(5:8)
      write(lu,rec=139) out%kt8(1:4)
      write(lu,rec=140) out%kt8(5:8)
      write(lu,rec=141) out%kt9(1:4)
      write(lu,rec=142) out%kt9(5:8)
      write(lu,rec=143) out%kf(1:4)
      write(lu,rec=144) out%kf(5:8)
      write(lu,rec=145) out%kuser0(1:4)
      write(lu,rec=146) out%kuser0(5:8)
      write(lu,rec=147) out%kuser1(1:4)
      write(lu,rec=148) out%kuser1(5:8)
      write(lu,rec=149) out%kuser2(1:4)
      write(lu,rec=150) out%kuser2(5:8)
      write(lu,rec=151) out%kcmpnm(1:4)
      write(lu,rec=152) out%kcmpnm(5:8)
      write(lu,rec=153) out%knetwk(1:4)
      write(lu,rec=154) out%knetwk(5:8)
      write(lu,rec=155) out%kdatrd(1:4)
      write(lu,rec=156) out%kdatrd(5:8)
      write(lu,rec=157) out%kinst(1:4)
      write(lu,rec=158) out%kinst(5:8)

 
      close(lu)
   end subroutine f90sac_writeheader
!===============================================================================


!===============================================================================
   subroutine f90sac_getfhdr(tr,id_hdr,val)
!===============================================================================
!
!     Return the value of a floating point header, identified by its ID number
!     (this can be got using f90sac_enumhdr)
!     
      implicit none
      type (SACtrace) :: tr
      integer :: id_hdr,ifound
      real :: val

!  ** check the ID
      if (id_hdr<1 .or. id_hdr>70) then
         write(0,'(a)') &
         'F90SAC_GETFHDR: Error: Bad FP header ID number (not 1-70)'
         STOP     
      endif      

      if (id_hdr == 01) val = tr % delta     
      if (id_hdr == 02) val = tr % depmin    
      if (id_hdr == 03) val = tr % depmax    
      if (id_hdr == 04) val = tr % scale     
      if (id_hdr == 05) val = tr % odelta    
      if (id_hdr == 06) val = tr % b         
      if (id_hdr == 07) val = tr % e         
      if (id_hdr == 08) val = tr % o         
      if (id_hdr == 09) val = tr % a         
      if (id_hdr == 10) val = tr % internal0 
      if (id_hdr == 11) val = tr % t0        
      if (id_hdr == 12) val = tr % t1        
      if (id_hdr == 13) val = tr % t2        
      if (id_hdr == 14) val = tr % t3        
      if (id_hdr == 15) val = tr % t4        
      if (id_hdr == 16) val = tr % t5        
      if (id_hdr == 17) val = tr % t6        
      if (id_hdr == 18) val = tr % t7        
      if (id_hdr == 19) val = tr % t8        
      if (id_hdr == 20) val = tr % t9        
      if (id_hdr == 21) val = tr % f         
      if (id_hdr == 22) val = tr % resp0     
      if (id_hdr == 23) val = tr % resp1     
      if (id_hdr == 24) val = tr % resp2     
      if (id_hdr == 25) val = tr % resp3     
      if (id_hdr == 26) val = tr % resp4     
      if (id_hdr == 27) val = tr % resp5     
      if (id_hdr == 28) val = tr % resp6     
      if (id_hdr == 29) val = tr % resp7     
      if (id_hdr == 30) val = tr % resp8     
      if (id_hdr == 31) val = tr % resp9     
      if (id_hdr == 32) val = tr % stla      
      if (id_hdr == 33) val = tr % stlo      
      if (id_hdr == 34) val = tr % stel      
      if (id_hdr == 35) val = tr % stdp      
      if (id_hdr == 36) val = tr % evla      
      if (id_hdr == 37) val = tr % evlo      
      if (id_hdr == 38) val = tr % evel      
      if (id_hdr == 39) val = tr % evdp      
      if (id_hdr == 40) val = tr % mag       
      if (id_hdr == 41) val = tr % user0     
      if (id_hdr == 42) val = tr % user1     
      if (id_hdr == 43) val = tr % user2     
      if (id_hdr == 44) val = tr % user3     
      if (id_hdr == 45) val = tr % user4     
      if (id_hdr == 46) val = tr % user5     
      if (id_hdr == 47) val = tr % user6     
      if (id_hdr == 48) val = tr % user7     
      if (id_hdr == 49) val = tr % user8     
      if (id_hdr == 50) val = tr % user9     
      if (id_hdr == 51) val = tr % dist      
      if (id_hdr == 52) val = tr % az        
      if (id_hdr == 53) val = tr % baz       
      if (id_hdr == 54) val = tr % gcarc     
      if (id_hdr == 55) val = tr % internal1 
      if (id_hdr == 56) val = tr % internal2 
      if (id_hdr == 57) val = tr % depmen    
      if (id_hdr == 58) val = tr % cmpaz     
      if (id_hdr == 59) val = tr % cmpinc    
      if (id_hdr == 60) val = tr % xminimum  
      if (id_hdr == 61) val = tr % xmaximum  
      if (id_hdr == 62) val = tr % yminimum  
      if (id_hdr == 63) val = tr % ymaximum  
      if (id_hdr == 64) val = tr % unused1   
      if (id_hdr == 65) val = tr % unused2   
      if (id_hdr == 66) val = tr % unused3   
      if (id_hdr == 67) val = tr % unused4   
      if (id_hdr == 68) val = tr % unused5   
      if (id_hdr == 69) val = tr % unused6   
      if (id_hdr == 70) val = tr % unused7         

          
      return
   end subroutine f90sac_getfhdr
!===============================================================================

!===============================================================================
   subroutine f90sac_getihdr(tr,id_hdr,val)
!===============================================================================
!
!     Return the value of an integer header, identified by its ID number
!     (this can be got using f90sac_enumhdr)
!     
      implicit none
      type (SACtrace) :: tr
      integer :: id_hdr,ifound
      integer :: val

!  ** check the ID
      if (id_hdr<71 .or. id_hdr>105) then
         write(0,'(a)') &
         'F90SAC_GETIHDR: Error: Bad Int header ID number (not 71-105)'
         STOP     
      endif      

      if (id_hdr == 071) val = tr % nzyear
      if (id_hdr == 072) val = tr % nzjday
      if (id_hdr == 073) val = tr % nzhour
      if (id_hdr == 074) val = tr % nzmin
      if (id_hdr == 075) val = tr % nzsec
      if (id_hdr == 076) val = tr % nzmsec
      if (id_hdr == 077) val = tr % nvhdr
      if (id_hdr == 078) val = tr % norid
      if (id_hdr == 079) val = tr % nevid
      if (id_hdr == 080) val = tr % npts
      if (id_hdr == 081) val = tr % internal3
      if (id_hdr == 082) val = tr % nwfid
      if (id_hdr == 083) val = tr % nxsize
      if (id_hdr == 084) val = tr % nysize
      if (id_hdr == 085) val = tr % unused8
      if (id_hdr == 086) val = tr % iftype
      if (id_hdr == 087) val = tr % idep
      if (id_hdr == 088) val = tr % iztype
      if (id_hdr == 089) val = tr % unused9
      if (id_hdr == 090) val = tr % iinst
      if (id_hdr == 091) val = tr % istreg
      if (id_hdr == 092) val = tr % ievreg
      if (id_hdr == 093) val = tr % ievtyp
      if (id_hdr == 094) val = tr % iqual
      if (id_hdr == 095) val = tr % isynth
      if (id_hdr == 096) val = tr % imagtyp
      if (id_hdr == 097) val = tr % imagsrc
      if (id_hdr == 098) val = tr % unused10
      if (id_hdr == 099) val = tr % unused11
      if (id_hdr == 100) val = tr % unused12
      if (id_hdr == 101) val = tr % unused13
      if (id_hdr == 102) val = tr % unused14
      if (id_hdr == 103) val = tr % unused15
      if (id_hdr == 104) val = tr % unused16
      if (id_hdr == 105) val = tr % unused17
          
      return
   end subroutine f90sac_getihdr
!===============================================================================

!===============================================================================
   subroutine f90sac_getlhdr(tr,id_hdr,val)
!===============================================================================
!
!     Return the value of a logical header, identified by its ID number
!     (this can be got using f90sac_enumhdr)
!
!     Note that logicals in SAC are integers with the value 0 or 1, not
!     fortran logicals.      
!     
      implicit none
      type (SACtrace) :: tr
      integer :: id_hdr,ifound
      integer :: val

!  ** check the ID
      if (id_hdr<106 .or. id_hdr>110) then
         write(0,'(a)') &
         'F90SAC_GETLHDR: Error: Bad Logical header ID number (not 106-110)'
         STOP     
      endif      

      if (id_hdr == 106) val = tr % leven
      if (id_hdr == 107) val = tr % lpspol
      if (id_hdr == 108) val = tr % lovrok
      if (id_hdr == 109) val = tr % lcalda
      if (id_hdr == 110) val = tr % unused18
          
      return
   end subroutine f90sac_getlhdr
!===============================================================================

!===============================================================================
   subroutine f90sac_getkhdr(tr,id_hdr,val)
!===============================================================================
!
!     Return the value of a character header, identified by its ID number
!     (this can be got using f90sac_enumhdr)
!     
      implicit none
      type (SACtrace) :: tr
      integer :: id_hdr,ifound
      character (len=*) :: val

!  ** check the ID
      if (id_hdr<111 .or. id_hdr>133) then
         write(0,'(a)') &
         'F90SAC_GETKHDR: Error: Bad Character header ID number (not 111-133)'
         STOP     
      endif      

      if (id_hdr == 111) val = tr % kstnm     
      if (id_hdr == 112) val = tr % kevnm     
      if (id_hdr == 113) val = tr % khole     
      if (id_hdr == 114) val = tr % ko        
      if (id_hdr == 115) val = tr % ka        
      if (id_hdr == 116) val = tr % kt0       
      if (id_hdr == 117) val = tr % kt1       
      if (id_hdr == 118) val = tr % kt2       
      if (id_hdr == 119) val = tr % kt3       
      if (id_hdr == 120) val = tr % kt4       
      if (id_hdr == 121) val = tr % kt5       
      if (id_hdr == 122) val = tr % kt6       
      if (id_hdr == 123) val = tr % kt7       
      if (id_hdr == 124) val = tr % kt8       
      if (id_hdr == 125) val = tr % kt9       
      if (id_hdr == 126) val = tr % kf        
      if (id_hdr == 127) val = tr % kuser0    
      if (id_hdr == 128) val = tr % kuser1    
      if (id_hdr == 129) val = tr % kuser2    
      if (id_hdr == 130) val = tr % kcmpnm    
      if (id_hdr == 131) val = tr % knetwk    
      if (id_hdr == 132) val = tr % kdatrd    
      if (id_hdr == 133) val = tr % kinst     
          
      return
   end subroutine f90sac_getkhdr
!===============================================================================

!===============================================================================
   subroutine f90sac_setfhdr(tr,id_hdr,val)
!===============================================================================
!
!     Set the value of a floating point header, identified by its ID number
!     (this can be got using f90sac_enumhdr)
!     
      implicit none
      type (SACtrace) :: tr
      integer :: id_hdr,ifound
      real :: val

!  ** check the ID
      if (id_hdr<1 .or. id_hdr>70) then
         write(0,'(a)') &
         'F90SAC_SETFHDR: Error: Bad FP header ID number (not 1-70)'
         STOP     
      endif      

      if (id_hdr == 01) tr % delta     = val 
      if (id_hdr == 02) tr % depmin    = val 
      if (id_hdr == 03) tr % depmax    = val 
      if (id_hdr == 04) tr % scale     = val 
      if (id_hdr == 05) tr % odelta    = val 
      if (id_hdr == 06) tr % b         = val 
      if (id_hdr == 07) tr % e         = val 
      if (id_hdr == 08) tr % o         = val 
      if (id_hdr == 09) tr % a         = val 
      if (id_hdr == 10) tr % internal0 = val 
      if (id_hdr == 11) tr % t0        = val 
      if (id_hdr == 12) tr % t1        = val 
      if (id_hdr == 13) tr % t2        = val 
      if (id_hdr == 14) tr % t3        = val 
      if (id_hdr == 15) tr % t4        = val 
      if (id_hdr == 16) tr % t5        = val 
      if (id_hdr == 17) tr % t6        = val 
      if (id_hdr == 18) tr % t7        = val 
      if (id_hdr == 19) tr % t8        = val 
      if (id_hdr == 20) tr % t9        = val 
      if (id_hdr == 21) tr % f         = val 
      if (id_hdr == 22) tr % resp0     = val 
      if (id_hdr == 23) tr % resp1     = val 
      if (id_hdr == 24) tr % resp2     = val 
      if (id_hdr == 25) tr % resp3     = val 
      if (id_hdr == 26) tr % resp4     = val 
      if (id_hdr == 27) tr % resp5     = val 
      if (id_hdr == 28) tr % resp6     = val 
      if (id_hdr == 29) tr % resp7     = val 
      if (id_hdr == 30) tr % resp8     = val 
      if (id_hdr == 31) tr % resp9     = val 
      if (id_hdr == 32) tr % stla      = val 
      if (id_hdr == 33) tr % stlo      = val 
      if (id_hdr == 34) tr % stel      = val 
      if (id_hdr == 35) tr % stdp      = val 
      if (id_hdr == 36) tr % evla      = val 
      if (id_hdr == 37) tr % evlo      = val 
      if (id_hdr == 38) tr % evel      = val 
      if (id_hdr == 39) tr % evdp      = val 
      if (id_hdr == 40) tr % mag       = val 
      if (id_hdr == 41) tr % user0     = val 
      if (id_hdr == 42) tr % user1     = val 
      if (id_hdr == 43) tr % user2     = val 
      if (id_hdr == 44) tr % user3     = val 
      if (id_hdr == 45) tr % user4     = val 
      if (id_hdr == 46) tr % user5     = val 
      if (id_hdr == 47) tr % user6     = val 
      if (id_hdr == 48) tr % user7     = val 
      if (id_hdr == 49) tr % user8     = val 
      if (id_hdr == 50) tr % user9     = val 
      if (id_hdr == 51) tr % dist      = val 
      if (id_hdr == 52) tr % az        = val 
      if (id_hdr == 53) tr % baz       = val 
      if (id_hdr == 54) tr % gcarc     = val 
      if (id_hdr == 55) tr % internal1 = val 
      if (id_hdr == 56) tr % internal2 = val 
      if (id_hdr == 57) tr % depmen    = val 
      if (id_hdr == 58) tr % cmpaz     = val 
      if (id_hdr == 59) tr % cmpinc    = val 
      if (id_hdr == 60) tr % xminimum  = val 
      if (id_hdr == 61) tr % xmaximum  = val 
      if (id_hdr == 62) tr % yminimum  = val 
      if (id_hdr == 63) tr % ymaximum  = val 
      if (id_hdr == 64) tr % unused1   = val 
      if (id_hdr == 65) tr % unused2   = val 
      if (id_hdr == 66) tr % unused3   = val 
      if (id_hdr == 67) tr % unused4   = val 
      if (id_hdr == 68) tr % unused5   = val 
      if (id_hdr == 69) tr % unused6   = val 
      if (id_hdr == 70) tr % unused7   = val       

          
      return
   end subroutine f90sac_setfhdr
!===============================================================================

!===============================================================================
   subroutine f90sac_setihdr(tr,id_hdr,val)
!===============================================================================
!
!     Set the value of an integer header, identified by its ID number
!     (this can be got using f90sac_enumhdr)
!     
      implicit none
      type (SACtrace) :: tr
      integer :: id_hdr,ifound
      integer :: val

!  ** check the ID
      if (id_hdr<71 .or. id_hdr>105) then
         write(0,'(a)') &
         'F90SAC_SETIHDR: Error: Bad Int header ID number (not 71-105)'
         STOP     
      endif      

      if (id_hdr == 071)  tr % nzyear    = val
      if (id_hdr == 072)  tr % nzjday    = val
      if (id_hdr == 073)  tr % nzhour    = val
      if (id_hdr == 074)  tr % nzmin     = val
      if (id_hdr == 075)  tr % nzsec     = val
      if (id_hdr == 076)  tr % nzmsec    = val
      if (id_hdr == 077)  tr % nvhdr     = val
      if (id_hdr == 078)  tr % norid     = val
      if (id_hdr == 079)  tr % nevid     = val
      if (id_hdr == 080)  tr % npts      = val
      if (id_hdr == 081)  tr % internal3 = val
      if (id_hdr == 082)  tr % nwfid     = val
      if (id_hdr == 083)  tr % nxsize    = val
      if (id_hdr == 084)  tr % nysize    = val
      if (id_hdr == 085)  tr % unused8   = val
      if (id_hdr == 086)  tr % iftype    = val
      if (id_hdr == 087)  tr % idep      = val
      if (id_hdr == 088)  tr % iztype    = val
      if (id_hdr == 089)  tr % unused9   = val
      if (id_hdr == 090)  tr % iinst     = val
      if (id_hdr == 091)  tr % istreg    = val
      if (id_hdr == 092)  tr % ievreg    = val
      if (id_hdr == 093)  tr % ievtyp    = val
      if (id_hdr == 094)  tr % iqual     = val
      if (id_hdr == 095)  tr % isynth    = val
      if (id_hdr == 096)  tr % imagtyp   = val
      if (id_hdr == 097)  tr % imagsrc   = val
      if (id_hdr == 098)  tr % unused10  = val
      if (id_hdr == 099)  tr % unused11  = val
      if (id_hdr == 100)  tr % unused12  = val
      if (id_hdr == 101)  tr % unused13  = val
      if (id_hdr == 102)  tr % unused14  = val
      if (id_hdr == 103)  tr % unused15  = val
      if (id_hdr == 104)  tr % unused16  = val
      if (id_hdr == 105)  tr % unused17  = val
          
      return
   end subroutine f90sac_setihdr
!===============================================================================

!===============================================================================
   subroutine f90sac_setlhdr(tr,id_hdr,val)
!===============================================================================
!
!     Set the value of a logical header, identified by its ID number
!     (this can be got using f90sac_enumhdr)
!
!     Note that logicals in SAC are integers with the value 0 or 1, not
!     fortran logicals.      
!     
      implicit none
      type (SACtrace) :: tr
      integer :: id_hdr,ifound
      integer :: val

!  ** check the ID
      if (id_hdr<106 .or. id_hdr>110) then
         write(0,'(a)') &
         'F90SAC_SETLHDR: Error: Bad Logical header ID number (not 106-110)'
         STOP     
      endif      

      if (id_hdr == 106)  tr % leven    = val
      if (id_hdr == 107)  tr % lpspol   = val
      if (id_hdr == 108)  tr % lovrok   = val
      if (id_hdr == 109)  tr % lcalda   = val
      if (id_hdr == 110)  tr % unused18 = val
          
      return
   end subroutine f90sac_setlhdr
!===============================================================================

!===============================================================================
   subroutine f90sac_setkhdr(tr,id_hdr,val)
!===============================================================================
!
!     Set the value of a character header, identified by its ID number
!     (this can be got using f90sac_enumhdr)
!     
      implicit none
      type (SACtrace) :: tr
      integer :: id_hdr,ifound
      character (len=*) :: val

!  ** check the ID
      if (id_hdr<111 .or. id_hdr>133) then
         write(0,'(a)') &
         'F90SAC_SETKHDR: Error: Bad Character header ID number (not 111-133)'
         STOP     
      endif      

      if (id_hdr == 111)  tr % kstnm  = val   
      if (id_hdr == 112)  tr % kevnm  = val   
      if (id_hdr == 113)  tr % khole  = val   
      if (id_hdr == 114)  tr % ko     = val   
      if (id_hdr == 115)  tr % ka     = val   
      if (id_hdr == 116)  tr % kt0    = val   
      if (id_hdr == 117)  tr % kt1    = val   
      if (id_hdr == 118)  tr % kt2    = val   
      if (id_hdr == 119)  tr % kt3    = val   
      if (id_hdr == 120)  tr % kt4    = val   
      if (id_hdr == 121)  tr % kt5    = val   
      if (id_hdr == 122)  tr % kt6    = val   
      if (id_hdr == 123)  tr % kt7    = val   
      if (id_hdr == 124)  tr % kt8    = val   
      if (id_hdr == 125)  tr % kt9    = val   
      if (id_hdr == 126)  tr % kf     = val   
      if (id_hdr == 127)  tr % kuser0 = val   
      if (id_hdr == 128)  tr % kuser1 = val   
      if (id_hdr == 129)  tr % kuser2 = val   
      if (id_hdr == 130)  tr % kcmpnm = val   
      if (id_hdr == 131)  tr % knetwk = val   
      if (id_hdr == 132)  tr % kdatrd = val   
      if (id_hdr == 133)  tr % kinst  = val   
          
      return
   end subroutine f90sac_setkhdr
!===============================================================================


!===============================================================================
   subroutine f90sac_enumhdr(hdrstr,id_hdr)
!===============================================================================
!
!     Return the header ID number (just the header's place in the SACfile list)
!     If the header can't be found, -1 is returned
!     
      implicit none
!      type (SACtrace) :: tr
      integer :: id_hdr
      character (len=*) :: hdrstr
      character (len=10) :: headers(133)

!    data headers / &
!    'delta','depmin','depmax','scale','odelta','b','e','o','a','internal0',    &
!    't0','t1','t2','t3','t4','t5','t6','t7','t8','t9','f','resp0','resp1',     &
!    'resp2','resp3','resp4','resp5','resp6','resp7','resp8','resp9','stla',    &
!    'stlo','stel','stdp','evla','evlo','evel','evdp','mag','user0','user1',    &
!    'user2','user3','user4','user5','user6','user7','user8','user9','dist',    &
!    'az','baz','gcarc','internal1','internal2','depmen','cmpaz','cmpinc',      &
!    'xminimum','xmaximum','yminimum','ymaximum','unused1','unused2',           &
!    'unused3','unused4','unused5','unused6','unused7','nzyear','nzjday',       &
!    'nzhour','nzmin','nzsec','nzmsec','nvhdr','norid','nevid','npts',          &
!    'internal3','nwfid','nxsize','nysize','unused8','iftype','idep',           &
!    'iztype','unused9','iinst','istreg','ievreg','ievtyp','iqual','isynth',    &
!    'imagtyp','imagsrc','unused10','unused11','unused12','unused13',           &
!    'unused14','unused15','unused16','unused17','leven','lpspol','lovrok',     &
!    'lcalda','unused18','kstnm','kevnm','khole','ko','ka','kt0','kt1','kt2',   &
!    'kt3','kt4','kt5','kt6','kt7','kt8','kt9','kf','kuser0','kuser1','kuser2', &
!    'kcmpnm','knetwk','kdatrd','kinst'/


!  ** setup the header array
      headers(001)='delta'    ;  headers(068)='unused5'  ;
      headers(002)='depmin'   ;  headers(069)='unused6'  ;
      headers(003)='depmax'   ;  headers(070)='unused7'  ;
      headers(004)='scale'    ;  headers(071)='nzyear'   ;
      headers(005)='odelta'   ;  headers(072)='nzjday'   ;
      headers(006)='b'        ;  headers(073)='nzhour'   ;
      headers(007)='e'        ;  headers(074)='nzmin'    ;
      headers(008)='o'        ;  headers(075)='nzsec'    ;
      headers(009)='a'        ;  headers(076)='nzmsec'   ;
      headers(010)='internal0';  headers(077)='nvhdr'    ;
      headers(011)='t0'       ;  headers(078)='norid'    ;
      headers(012)='t1'       ;  headers(079)='nevid'    ;
      headers(013)='t2'       ;  headers(080)='npts'     ;
      headers(014)='t3'       ;  headers(081)='internal3';
      headers(015)='t4'       ;  headers(082)='nwfid'    ;
      headers(016)='t5'       ;  headers(083)='nxsize'   ;
      headers(017)='t6'       ;  headers(084)='nysize'   ;
      headers(018)='t7'       ;  headers(085)='unused8'  ;
      headers(019)='t8'       ;  headers(086)='iftype'   ;
      headers(020)='t9'       ;  headers(087)='idep'     ;
      headers(021)='f'        ;  headers(088)='iztype'   ;
      headers(022)='resp0'    ;  headers(089)='unused9'  ;
      headers(023)='resp1'    ;  headers(090)='iinst'    ;
      headers(024)='resp2'    ;  headers(091)='istreg'   ;
      headers(025)='resp3'    ;  headers(092)='ievreg'   ;
      headers(026)='resp4'    ;  headers(093)='ievtyp'   ;
      headers(027)='resp5'    ;  headers(094)='iqual'    ;
      headers(028)='resp6'    ;  headers(095)='isynth'   ;
      headers(029)='resp7'    ;  headers(096)='imagtyp'  ;
      headers(030)='resp8'    ;  headers(097)='imagsrc'  ;
      headers(031)='resp9'    ;  headers(098)='unused10' ;
      headers(032)='stla'     ;  headers(099)='unused11' ;
      headers(033)='stlo'     ;  headers(100)='unused12' ;
      headers(034)='stel'     ;  headers(101)='unused13' ;
      headers(035)='stdp'     ;  headers(102)='unused14' ;
      headers(036)='evla'     ;  headers(103)='unused15' ;
      headers(037)='evlo'     ;  headers(104)='unused16' ;
      headers(038)='evel'     ;  headers(105)='unused17' ;
      headers(039)='evdp'     ;  headers(106)='leven'    ;
      headers(040)='mag'      ;  headers(107)='lpspol'   ;
      headers(041)='user0'    ;  headers(108)='lovrok'   ;
      headers(042)='user1'    ;  headers(109)='lcalda'   ;
      headers(043)='user2'    ;  headers(110)='unused18' ;
      headers(044)='user3'    ;  headers(111)='kstnm'    ;
      headers(045)='user4'    ;  headers(112)='kevnm'    ;
      headers(046)='user5'    ;  headers(113)='khole'    ;
      headers(047)='user6'    ;  headers(114)='ko'       ;
      headers(048)='user7'    ;  headers(115)='ka'       ;
      headers(049)='user8'    ;  headers(116)='kt0'      ;
      headers(050)='user9'    ;  headers(117)='kt1'      ;
      headers(051)='dist'     ;  headers(118)='kt2'      ;
      headers(052)='az'       ;  headers(119)='kt3'      ;
      headers(053)='baz'      ;  headers(120)='kt4'      ;
      headers(054)='gcarc'    ;  headers(121)='kt5'      ;
      headers(055)='internal1';  headers(122)='kt6'      ;
      headers(056)='internal2';  headers(123)='kt7'      ;
      headers(057)='depmen'   ;  headers(124)='kt8'      ;
      headers(058)='cmpaz'    ;  headers(125)='kt9'      ;
      headers(059)='cmpinc'   ;  headers(126)='kf'       ;
      headers(060)='xminimum' ;  headers(127)='kuser0'   ;
      headers(061)='xmaximum' ;  headers(128)='kuser1'   ;
      headers(062)='yminimum' ;  headers(129)='kuser2'   ;
      headers(063)='ymaximum' ;  headers(130)='kcmpnm'   ;
      headers(064)='unused1'  ;  headers(131)='knetwk'   ;
      headers(065)='unused2'  ;  headers(132)='kdatrd'   ;
      headers(066)='unused3'  ;  headers(133)='kinst'    ;
      headers(067)='unused4'  ;

!  ** search for the specified header
      do id_hdr = 1,133
         if (trim(hdrstr) == trim(headers(id_hdr))) then
            return ! the do loop
         endif
      enddo   
      
      id_hdr = -1 ! not found condition
          
      return
   end subroutine f90sac_enumhdr
!===============================================================================

!===============================================================================
   function f90sac_dateseed()
!===============================================================================
!
!  Generate a random number seed using the date (day number) and time
!  Updated to use the date_and_time
!
      implicit none
      integer :: dd,mm,yy ! dd,mm,yy
      integer :: dayno,nsec
      integer :: hh,mi,se,ms      

      integer :: f90sac_dateseed ! random number seed 
      character :: dtsdstr*8
      character :: datestr*8,timestr*10

!  ** get the time and date      
      call date_and_time(datestr,timestr)
      read(datestr,'(i4.4,i2.2,i2.2)') yy,mm,dd
      read(timestr,'(i2.2,i2.2,i2.2,1x,i3.3)') hh,mi,se,ms

!  ** get a day (in year) number
!      call f90sac_ymd2jd(yy,mm,dd,dayno)   
      
      f90sac_dateseed = (hh*3600 +mi*60 +se) * 1000 + ms
      
      return
   end function f90sac_dateseed
!===============================================================================

!===============================================================================
   subroutine f90sac_ymd2jd(iyear,imonth,iday,ijd)
!===============================================================================
!     
!  Convert YYYY/MM/DD to DAY number in year
!
      implicit none
      integer iyear,imonth,iday,ijd,i
      integer ndaysin(12)
      character (len=35) :: dstr = '31,28,31,30,31,30,31,31,30,31,30,31'
      logical isleap

!  ** do an internal read to get days-in-month into array      
      read(dstr,*) (ndaysin(i),i=1,12)

!  ** See if it is a leap year
      isleap = .false.
      if (modulo(iyear,4) == 0) isleap = .true.
      if (modulo(iyear,100) == 0) isleap = .false.
      if (modulo(iyear,400) == 0) isleap = .true.
      if (isleap) ndaysin(2) = ndaysin(2) + 1
      
!   ** check the date makes sense
      if (imonth > 12 .or. imonth < 1) then
         write(0,'(a)') 'F90SAC_YMD2JD: Error: Bad date'
         stop
      endif  
 
      if (iday < 0 .or. iday > ndaysin(imonth)) then
         write(0,'(a)') 'F90SAC_YMD2JD: Error: Bad date'
         stop
      endif      
 
      ijd = sum(ndaysin(1:imonth-1)) + iday      
      
      
      return
   end subroutine f90sac_ymd2jd
!===============================================================================

!===============================================================================
   subroutine f90sac_jd2ymd(iyear,ijd,imonth,iday)
!===============================================================================
!     
!  Convert YYYY/MM/DD to DAY number in year
!
      implicit none
      integer iyear,imonth,iday,ijd,i,ijd_temp
      integer ndaysin(12),ndaysinyear
      character (len=35) :: dstr = '31,28,31,30,31,30,31,31,30,31,30,31'
      logical isleap

!  ** do an internal read to get days-in-month into array      
      read(dstr,*) (ndaysin(i),i=1,12)

!  ** See if it is a leap year
      isleap = .false.
      if (modulo(iyear,4) == 0) isleap = .true.
      if (modulo(iyear,100) == 0) isleap = .false.
      if (modulo(iyear,400) == 0) isleap = .true.
      if (isleap) ndaysin(2) = ndaysin(2) + 1

      ndaysinyear = 365 ; if (isleap) ndaysinyear = 366
      
!  ** check the day number makes sense
      if (ijd > ndaysinyear .or. ijd < 1) then
         write(0,'(a)') 'F90SAC_JD2YMD: Error: Bad date'
         stop
      endif  
 
!  ** convert to month/day numbder
      ijd_temp = ijd
      do i=1,12
         if (ijd_temp-ndaysin(i) <= 0) exit
         ijd_temp = ijd_temp - ndaysin(i)
      enddo 
      imonth = i ; iday = ijd_temp      
      
      
      return
   end subroutine f90sac_jd2ymd
!===============================================================================

!===============================================================================
   end module f90sac
!===============================================================================
!  END OF F90SAC module
!===============================================================================

