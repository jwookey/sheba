!===============================================================================
!  S H E B A - Shear-wave Birefringence Analysis
!===============================================================================
!  This software is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!===============================================================================
!
!  James Wookey, School of Earth Sciences, University of Bristol
!  CVS: $Revision: 1.5 $ $Date: 2008/03/20 06:27:11 $
!

!=======================================================================
      subroutine output_result(t1)
!=======================================================================
!  
!     Write a two line file containing event information and result
!     
!  
!-----------------------------------------------------------------------
      use f90sac  ! use the f90sac module
      use sheba_config ! use the sheba_config module
      use event_info ! use the event_info module
!-----------------------------------------------------------------------
      implicit none
      type (SACTrace) :: t1
      character*50 fname

      open(99,file='sheba.result')

      write(99,'(3a)') '%  DATE TIME    EVLA    EVLO    STLA    STLO', &
      '    EVDP    DIST     AZI     BAZ    FAST   DFAST    TLAG   ', &
       'DTLAG    SPOL   DSPOL    WBEG    WEND  STAT'
      write(99,100) t1 % nzyear,t1 % nzjday,  &
                    t1 % nzhour, t1 % nzmin, &
                    t1 % evla, t1 % evlo, &
                    t1 % stla, t1 % stlo, &
                    t1 % evdp, t1 % gcarc, t1 % az, t1 % baz, &
                    event % fast, event % dfast, &
                    event % tlag, event % dtlag, &
                    event % spol, event % dspol, &
                    event % wbeg, event % wend, &
                    t1 % kstnm

       
      close(99) 

      fname = trim(config % fname_base) // '_sheba.result'
      
      open(99,file=fname)
      
      write(99,'(3a)') '%  DATE TIME    EVLA    EVLO    STLA    STLO', &
      '    EVDP    DIST     AZI     BAZ    FAST   DFAST    TLAG   ', &
       'DTLAG    SPOL   DSPOL    WBEG    WEND  STAT'
      write(99,100) t1 % nzyear,t1 % nzjday,  &
                    t1 % nzhour, t1 % nzmin, &
                    t1 % evla, t1 % evlo, &
                    t1 % stla, t1 % stlo, &
                    t1 % evdp, t1 % gcarc, t1 % az, t1 % baz, &
                    event % fast, event % dfast, &
                    event % tlag, event % dtlag, &
                    event % spol, event % dspol, &
                    event % wbeg, event % wend, &
                    t1 % kstnm

       
      close(99) 
       
      fname = trim(config % fname_base) // '_sheba.stats'
      
      open(99,file=fname)
      
      write(99,'(3a)') '%  DATE TIME EIGORIG EIGCORR     SNR   NDF   STAT'
      write(99,101) t1 % nzyear,t1 % nzjday,  &
                    t1 % nzhour, t1 % nzmin, &
                    event % eigrat_orig, event % eigrat_corr, &
                    event % snr, event % ndf, &
                    t1 % kstnm

       
      close(99) 

      return           
100   format(i4.4,i3.3,1x,2i2.2,14f8.2,2f8.2,'  % ',a)      
101   format(i4.4,i3.3,1x,2i2.2,2f8.4,f8.3,x,i5,x,'  % ',a)      
      end subroutine output_result
!=======================================================================

!=======================================================================
      subroutine output_pm_files(t1,t2,fname)
!=======================================================================
!  
!     Calculate and output information needed for GMT plotting
!     This is in the form of a text file which can be sourced
!     in the GMT shell-script
!
!-----------------------------------------------------------------------
      use f90sac  ! use the f90sac module
      use sheba_config ! use the sheba_config module
      use event_info ! use the event_info module
!-----------------------------------------------------------------------
      implicit none
      type (SACTrace) :: t1,t2 ! the traces
      type (SACTrace) :: t1_wind,t2_wind ! windowed traces
      character*50 fname

!     * window the traces
      call f90sac_window(t1,t1_wind,event % wbeg,event % wend)
      call f90sac_window(t2,t2_wind,event % wbeg,event % wend)

!     * create the hybrid file
      call write_xyfile(t1_wind,t2_wind,fname)
      return           
      end subroutine output_pm_files
!=======================================================================


!=======================================================================
      subroutine output_gmt_info()
!=======================================================================
!  
!     Calculate and output information needed for GMT plotting
!     This is in the form of a text file which can be sourced
!     in the GMT shell-script. Also output a short SAC macro
!     to load the event parameters
!
!-----------------------------------------------------------------------
      use sheba_config ! use the sheba_config module
      use event_info ! use the event_info module
      use array_sizes ! use the array sizes module
!-----------------------------------------------------------------------
      implicit none
      character*50 fname
      real tlag_minortick,tlag_majortick
      fname = trim(config % fname_base) // '.gmt'
      
      open(99,file=fname)

!        * output splitting solution parameters
         write(99,'(a,f10.4)') 'set TLAG =', event % tlag 
         write(99,'(a,f10.4)') 'set DTLAG =', event % dtlag 
         write(99,'(a,f10.4)') 'set FAST =', event % fast
         write(99,'(a,f10.4)') 'set DFAST =', event % dfast 
         write(99,'(a,f10.4)') 'set SPOL =', event % spol 
         write(99,'(a,f10.4)') 'set DSPOL =', event % dspol 
!        * generate tick information
         call tickinfo(config%max_tlag,tlag_majortick,tlag_minortick)


!     ** because of the requirement to use integer number of samples
!        need to adjust the maximum slightly for the GMT code to work

!         write(99,'(a,f10.6)') 'set TLAG_SCALE =',config % tlag_scale
         write(99,'(a,f14.10)') 'set TLAG_SCALE =', &
            real(np2int-1) * event % error_grid_tlag_int

         write(99,'(a,f10.6)') 'set TLAG_MAJORTICK =',tlag_majortick
         write(99,'(a,f10.6)') 'set TLAG_MINORTICK =',tlag_minortick
         write(99,'(a,f14.10)') 'set ERROR_GRID_TLAG_INT =', &
                     event % error_grid_tlag_int
         write(99,'(a,f10.4)') 'set WBEG =', event % wbeg 
         write(99,'(a,f10.4)') 'set WEND =', event % wend 
         write(99,'(a,i10)') 'set IBEST =', event % ibest 
           
      close(99)
      fname = trim(config % fname_base) // '.scm'
      
      open(99,file=fname)

!        * output splitting solution parameters
         write(99,'(a)') 'Message "Reading pars"'
         write(99,'(a,f10.4)') 'setbb TLAG ', event % tlag 
         write(99,'(a,f10.4)') 'setbb DTLAG', event % dtlag 
         write(99,'(a,f10.4)') 'setbb FAST ', event % fast
         write(99,'(a,f10.4)') 'setbb DFAST ', event % dfast 
         write(99,'(a,f10.4)') 'setbb SPOL ', event % spol 
         write(99,'(a,f10.4)') 'setbb DSPOL ', event % dspol 
         write(99,'(a,f10.4)') 'setbb WBEG ', event % wbeg 
         write(99,'(a,f10.4)') 'setbb WEND ', event % wend 
           
      close(99)

      return           
      end subroutine output_gmt_info
!=======================================================================


!=======================================================================
      subroutine tickinfo(scale,major,minor)
!=======================================================================
!  
!     Calculate sensible major and minor tick information from the scale
!     (assumes 0->scale). Currently this is simply hardwired for 2,4 or 8
!
!-----------------------------------------------------------------------
      implicit none
      real scale, major,minor

      if (scale <= 0.00001) then
         major=0.0000025 ; minor = 0.0000005
      elseif (scale <= 0.00002) then
         major=0.00005 ; minor = 0.000001
      elseif (scale <= 0.00005) then
         major=0.00001 ; minor = 0.000002


      elseif (scale <= 0.0001) then
         major=0.000025 ; minor = 0.000005
      elseif (scale <= 0.0002) then
         major=0.0005 ; minor = 0.00001
      elseif (scale <= 0.0005) then
         major=0.0001 ; minor = 0.00002

      elseif (scale <= 0.001) then
         major=0.00025 ; minor = 0.00005
      elseif (scale <= 0.002) then
         major=0.0005 ; minor = 0.0001
      elseif (scale <= 0.005) then
         major=0.001 ; minor = 0.0002
 
      elseif (scale <= 0.01) then
         major=0.0025 ; minor = 0.0005
      elseif (scale <= 0.02) then
         major=0.005 ; minor = 0.001
      elseif (scale <= 0.05) then
         major=0.01 ; minor = 0.002
      elseif (scale <= 0.1) then
         major=0.025 ; minor = 0.005            
      elseif (scale <= 0.2) then
         major=0.05 ; minor = 0.01
      elseif (scale <= 0.5) then
         major=0.1 ; minor = 0.02
      elseif (scale <= 1.0) then
         major=0.25 ; minor = 0.05
      elseif (scale <= 2.0) then
         major=0.5 ; minor = 0.10
      elseif (scale <= 5.0) then
         major=1.0 ; minor = 0.2
      elseif (scale <= 10.0) then
         major=2.0 ; minor = 0.4
      else
         major=0.25 ; minor = 0.05      
      endif
      
      return           
      end subroutine tickinfo
!=======================================================================



!=======================================================================
      subroutine output_traces(h1,h2,v,str,iorder)
!=======================================================================
!  
!     Output 3 traces
!
!-----------------------------------------------------------------------
      use f90sac  ! use the f90sac module
      use sheba_config ! use the sheba_config module
!-----------------------------------------------------------------------
      implicit none
      
      type (SACTrace) :: h1,h2,v ! the corrected traces
      character (len=5) :: str ! a descriptive extension to the filename
      character (len=80) :: fn(3) ! filename
      integer :: i,iorder(3) ! index to re-order the traces
!     * make the filenames
      do i = 1,3
         fn(i) = trim(config % fname_base) // trim(str)  &
                     // '.' // config % comp(iorder(i))
      enddo ! do i = 1,3
      call f90sac_writetrace(fn(1),h1)
      call f90sac_writetrace(fn(2),h2)
      call f90sac_writetrace(fn(3),v)
!     * and you're done ...      
      return           
      end subroutine output_traces
!=======================================================================
