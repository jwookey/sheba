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

#ifndef NO_NETCDF

!=======================================================================
   subroutine output_result_nc()
!=======================================================================
!  
!     Output results information to netCDF format
!  
!-----------------------------------------------------------------------
   use sheba_config ! use the sheba_config module
   use event_info ! use the event_info module
   use netcdf, only: nf90_create, nf90_def_dim, nf90_def_var, nf90_enddef, &
      nf90_put_var, nf90_close, NF90_CLOBBER, NF90_FLOAT, nf90_noerr, &
      nf90_strerror, NF90_GLOBAL, nf90_put_att      
!-----------------------------------------------------------------------
      implicit none
      character*50 :: fname
      integer :: j,i
      character*12 :: fmt

      integer :: ncid, tlag_dimid, fast_dimid, window_dimid, cluster_dimid
      integer :: result_dimid ! scalar dimension for single numbers

!  ** variable ids
      integer :: mw_wbeg_varid
      integer :: mw_wend_varid
      integer :: mw_tlag_varid
      integer :: mw_dtlag_varid
      integer :: mw_fast_varid
      integer :: mw_dfast_varid
      integer :: cluster_xc0_varid
      integer :: cluster_yc0_varid
      integer :: cluster_vxc0_varid
      integer :: cluster_vyc0_varid
      integer :: fast_vector_varid     
      integer :: tlag_vector_varid     
      integer :: lam1_raw_grid_varid 
      integer :: lam2_raw_grid_varid 
      integer :: lam2_norm_grid_varid

!  ** open the file
      fname = trim(config % fname_base) // '_sheba_result.nc'
      call check_ncf(nf90_create(trim(fname), NF90_CLOBBER, ncid))

!  ** Make a title
      call check_ncf(nf90_put_att(ncid, NF90_GLOBAL, 'title', &
      'SHEBA splitting result'))

!  ** Save the result details as global attributes
      call check_ncf(nf90_put_att(ncid, NF90_GLOBAL, 'tlag', event % tlag))
      call check_ncf(nf90_put_att(ncid, NF90_GLOBAL, 'dtlag', event % dtlag))
      call check_ncf(nf90_put_att(ncid, NF90_GLOBAL, 'fast', event % fast))
      call check_ncf(nf90_put_att(ncid, NF90_GLOBAL, 'dfast', event % dfast))
      call check_ncf(nf90_put_att(ncid, NF90_GLOBAL, 'spol', event % spol))
      call check_ncf(nf90_put_att(ncid, NF90_GLOBAL, 'dspol', event % dspol))
      call check_ncf(nf90_put_att(ncid, NF90_GLOBAL, 'tlagXC', event % tlagXC))
      call check_ncf(nf90_put_att(ncid, NF90_GLOBAL, 'dtlagXC', event % dtlagXC))
      call check_ncf(nf90_put_att(ncid, NF90_GLOBAL, 'fastXC', event % fastXC))
      call check_ncf(nf90_put_att(ncid, NF90_GLOBAL, 'dfastXC', event % dfastXC))
      call check_ncf(nf90_put_att(ncid, NF90_GLOBAL, 'wbeg', event % wbeg))
      call check_ncf(nf90_put_att(ncid, NF90_GLOBAL, 'wend', event % wend))
      call check_ncf(nf90_put_att(ncid, NF90_GLOBAL, 'snr', event % snr))
      call check_ncf(nf90_put_att(ncid, NF90_GLOBAL, 'eigrat_orig', event % eigrat_orig))
      call check_ncf(nf90_put_att(ncid, NF90_GLOBAL, 'eigrat_corr', event % eigrat_corr))
      call check_ncf(nf90_put_att(ncid, NF90_GLOBAL, 'ndf', event % ndf))
      call check_ncf(nf90_put_att(ncid, NF90_GLOBAL, 'best_window', event % ibest))
      call check_ncf(nf90_put_att(ncid, NF90_GLOBAL, 'best_cluster', event % kbest))
      call check_ncf(nf90_put_att(ncid, NF90_GLOBAL, 'qfactor', event % Quality))
      call check_ncf(nf90_put_att(ncid, NF90_GLOBAL, 'intensity_estimated', event % intensity))
      call check_ncf(nf90_put_att(ncid, NF90_GLOBAL, 'intensity', event % intensity_p))

!  ** Define dimensions
      call check_ncf(nf90_def_dim(ncid, 'result', 1, result_dimid))
      call check_ncf(nf90_def_dim(ncid, 'search_tlag',    np2int, tlag_dimid))
      call check_ncf(nf90_def_dim(ncid, 'search_fast',    np1, fast_dimid))
      call check_ncf(nf90_def_dim(ncid, 'window', event % nwindows, window_dimid))
      call check_ncf(nf90_def_dim(ncid, 'cluster', event % ncluster, cluster_dimid))

!  ** Multiwindow variables
      call check_ncf(nf90_def_var(ncid, 'mw_wbeg', NF90_FLOAT, window_dimid, mw_wbeg_varid))
      call check_ncf(nf90_def_var(ncid, 'mw_wend', NF90_FLOAT, window_dimid, mw_wend_varid))
      call check_ncf(nf90_def_var(ncid, 'mw_tlag', NF90_FLOAT, window_dimid, mw_tlag_varid))
      call check_ncf(nf90_def_var(ncid, 'mw_dtlag', NF90_FLOAT, window_dimid, mw_dtlag_varid))
      call check_ncf(nf90_def_var(ncid, 'mw_fast', NF90_FLOAT, window_dimid, mw_fast_varid))
      call check_ncf(nf90_def_var(ncid, 'mw_dfast', NF90_FLOAT, window_dimid, mw_dfast_varid))

!  ** Cluster variables
      call check_ncf(nf90_def_var(ncid, 'cluster_xc0', NF90_FLOAT, cluster_dimid, cluster_xc0_varid))
      call check_ncf(nf90_def_var(ncid, 'cluster_yc0', NF90_FLOAT, cluster_dimid, cluster_yc0_varid))
      call check_ncf(nf90_def_var(ncid, 'cluster_vxc0', NF90_FLOAT, cluster_dimid, cluster_vxc0_varid))
      call check_ncf(nf90_def_var(ncid, 'cluster_vyc0', NF90_FLOAT, cluster_dimid, cluster_vyc0_varid))

!  ** Grid variables
!call check_ncf(nf90_def_var(ncid, 'amplitude', NF90_FLOAT, &
!      (/x_dimid, y_dimid, twtt_dimid/), amp_varid))
      call check_ncf(nf90_def_var(ncid, 'fast_vector'     , NF90_FLOAT, fast_dimid , fast_vector_varid     ))
      call check_ncf(nf90_def_var(ncid, 'tlag_vector'     , NF90_FLOAT, tlag_dimid , tlag_vector_varid     ))

      call check_ncf(nf90_def_var(ncid, 'lam1_raw_grid' , NF90_FLOAT, (/fast_dimid,tlag_dimid/), lam1_raw_grid_varid ))
      call check_ncf(nf90_def_var(ncid, 'lam2_raw_grid' , NF90_FLOAT, (/fast_dimid,tlag_dimid/), lam2_raw_grid_varid ))
      call check_ncf(nf90_def_var(ncid, 'lam2_norm_grid', NF90_FLOAT, (/fast_dimid,tlag_dimid/), lam2_norm_grid_varid))

!  ** End definition
      call check_ncf(nf90_enddef(ncid))

!  ** Populate variables
      call check_ncf(nf90_put_var(ncid, mw_wbeg_varid       , event % mw_wbeg(1:event % nwindows)))
      call check_ncf(nf90_put_var(ncid, mw_wend_varid       , event % mw_wend(1:event % nwindows)))
      call check_ncf(nf90_put_var(ncid, mw_tlag_varid       , event % mw_tlag(1:event % nwindows)))
      call check_ncf(nf90_put_var(ncid, mw_dtlag_varid      , event % mw_dtlag(1:event % nwindows)))
      call check_ncf(nf90_put_var(ncid, mw_fast_varid       , event % mw_fast(1:event % nwindows)))
      call check_ncf(nf90_put_var(ncid, mw_dfast_varid      , event % mw_dfast(1:event % nwindows)))

      call check_ncf(nf90_put_var(ncid, cluster_xc0_varid   , event % cluster_xc0(1:event % ncluster)))
      call check_ncf(nf90_put_var(ncid, cluster_yc0_varid   , event % cluster_yc0(1:event % ncluster)))
      call check_ncf(nf90_put_var(ncid, cluster_vxc0_varid  , event % cluster_vxc0(1:event % ncluster)))
      call check_ncf(nf90_put_var(ncid, cluster_vyc0_varid  , event % cluster_vyc0(1:event % ncluster)))

      call check_ncf(nf90_put_var(ncid, fast_vector_varid     , event % fast_vector           ))
      call check_ncf(nf90_put_var(ncid, tlag_vector_varid     , event % tlag_vector           ))

      call check_ncf(nf90_put_var(ncid, lam1_raw_grid_varid , event % lam1_raw_grid       ))
      call check_ncf(nf90_put_var(ncid, lam2_raw_grid_varid , event % lam2_raw_grid       ))
      call check_ncf(nf90_put_var(ncid, lam2_norm_grid_varid, event % lam2_norm_grid      ))

!  ** Finalise file
      call check_ncf(nf90_close(ncid))

      return
   end subroutine output_result_nc
!=======================================================================

!===============================================================================
   subroutine check_ncf(status)
!===============================================================================
   use netcdf, only: nf90_noerr, nf90_strerror
      integer, intent(in) :: status
      if (status /= nf90_noerr) then
         write(0,*) 'netCDF error: ' // trim(nf90_strerror(status))
         stop
      endif
   end subroutine check_ncf
!-------------------------------------------------------------------------------

#endif

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
      character*50 :: fname
      integer :: j,i
      character*12 :: fmt

      open(99,file='sheba.result')

      write(99,'(3a)') '%  DATE TIME    EVLA    EVLO    STLA    STLO', &
      '    EVDP    DIST     AZI     BAZ    FAST   DFAST    TLAG   ', &
       'DTLAG    SPOL   DSPOL        WBEG        WEND  STAT'
      write(99,100) t1 % nzyear,t1 % nzjday,  &
                    t1 % nzhour, t1 % nzmin, &
                    t1 % evla, t1 % evlo, &
                    t1 % stla, t1 % stlo, &
                    t1 % evdp, t1 % gcarc, t1 % az, t1 % baz, &
                    event % fast, event % dfast, &
                    event % tlag, event % dtlag, &
                    event % spol, event % dspol, &
                    event % wbeg, event % wend, &
                    trim(t1 % kstnm)

       
      close(99) 

      fname = trim(config % fname_base) // '_sheba.result'
      
      open(99,file=fname)
      
      write(99,'(3a)') '%  DATE TIME    EVLA    EVLO    STLA    STLO', &
      '    EVDP    DIST     AZI     BAZ    FAST   DFAST    TLAG   ', &
       'DTLAG    SPOL   DSPOL        WBEG        WEND  STAT'
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

      fname = trim(config % fname_base) // '_sheba.XCresult'
      
      open(99,file=fname)
      
      write(99,'(3a)') '%  DATE TIME    EVLA    EVLO    STLA    STLO', &
      '    EVDP    DIST     AZI     BAZ    FAST   DFAST    TLAG   ', &
       'DTLAG    SPOL   DSPOL        WBEG        WEND  STAT'
      write(99,100) t1 % nzyear,t1 % nzjday,  &
                    t1 % nzhour, t1 % nzmin, &
                    t1 % evla, t1 % evlo, &
                    t1 % stla, t1 % stlo, &
                    t1 % evdp, t1 % gcarc, t1 % az, t1 % baz, &
                    event % fastXC, event % dfastXC, &
                    event % tlagXC, event % dtlagXC, &
                    event % spol, event % dspol, &
                    event % wbeg, event % wend, &
                    t1 % kstnm

       
      close(99) 
       
      fname = trim(config % fname_base) // '_sheba.stats'
      
      open(99,file=fname)
      
      write(99,'(3a)') '%  DATE TIME EIGORIG EIGCORR      Q     SNR   NDF   SI(Pr)   SI(Pa)  STAT'
      write(99,101) t1 % nzyear,t1 % nzjday,  &
                    t1 % nzhour, t1 % nzmin, &
                    event % eigrat_orig, event % eigrat_corr, &
                    event % Quality, event % snr, event % ndf, &
                    event % intensity_p, &
                    event % intensity, &
                    t1 % kstnm

       
      close(99) 

!  ** Output a chunk of XML for MATISSE
      fname = trim(config % fname_base) // '.mts'
      
      open(99,file=fname)
      
      write(99,'(a)') '         <data>'
!  ** trace names      
      do j=1,3
         write(99,120) j,trim(config % fname_base), &
                         trim(config % comp(config % iorder(j))),j
      enddo
120   format('            <file',i1,'>',a,'.',a,'</file',i1,'>')

!  ** data frame. SHEBA doesn't know this from the data
!     so it is set by the global constant mts_frame
      write(99,'(3a)') '            <frame>',mts_frame,'</frame>'

!  ** window
      write(99,'(a,e12.4,a)') '            <wbeg>', event % wbeg,'</wbeg>'
      write(99,'(a,e12.4,a)') '            <wend>', event % wend,'</wend>'
!  ** SNR
      write(99,'(a,f7.2,a)') '            <snr>', event % snr,'</snr>'
      write(99,'(a,f7.5,a)') &
         '            <l2s>', event % eigrat_corr,'</l2s>'
!  ** NDF      
      write(99,'(a,i4.4,a)') '            <ndf>', event % ndf,'</ndf>'
      
      write(99,'(a)') '         </data>'
      close(99)

!  ** Output a raw lam2/lam1 grid 
      fname = trim(config % fname_base) // '.lamR'
      
      open(99,file=fname)
      write(fmt,'(a1,i5.5,a)') '(',np2int,'f8.5)'
      do i=1,np1
         write(99,fmt) &
            ((event % lam2_raw_grid(i,j)/event % lam2_raw_grid(i,j)),j=1,np2int) 
      enddo

      close(99)

      return            
100   format(i4.4,i3.3,1x,2i2.2,10f8.2,2f8.5,2f8.2,2f12.5,'  % ',a)      
101   format(i4.4,i3.3,1x,2i2.2,2f8.5,f7.3,f8.3,x,i5,&
      x,f8.4,x,f8.4,'  % ',a)      
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
      type (SACTrace) :: pm ! pm file
      character*50 fname

!     * window the traces
      call f90sac_window(t1,t1_wind,event % wbeg,event % wend)
      call f90sac_window(t2,t2_wind,event % wbeg,event % wend)

!     * create the hybrid file
!      call write_xyfile(t1_wind,t2_wind,fname)
      call f90sac_combine_xy(t1_wind,t2_wind,pm)
      call f90sac_writetrace(fname,pm)
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
         write(99,'(a,f10.4)') 'set QUALITY =', event % Quality 

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
         write(99,'(a)') 'Message "Setting result parameters ..."'
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
      subroutine output_traces(h1,h2,v,str)
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
                     // '.' // config % comp(config % iorder(i))
      enddo ! do i = 1,3
      call f90sac_writetrace(fn(1),h1)
      call f90sac_writetrace(fn(2),h2)
      call f90sac_writetrace(fn(3),v)
!     * and you're done ...      
      return           
      end subroutine output_traces
!=======================================================================
