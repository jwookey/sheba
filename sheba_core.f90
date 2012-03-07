!===============================================================================
!  S H E B A - Shear-wave Birefringence Analysis
!===============================================================================
!  This software is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!===============================================================================
!
!  James Wookey, School of Earth Sciences, University of Bristol
!  CVS: $Revision: 1.6 $ $Date: 2011/02/11 16:09:41 $
!

!===============================================================================
      subroutine sheba()
!===============================================================================
!  
!     The main SHEBA soubroutine. This controls the shear-wave 
!     analysis
!
!-------------------------------------------------------------------------------
      use f90sac  ! use the f90sac module
      use sheba_config ! use the sheba_config module
      use event_info ! use event_info module
!-------------------------------------------------------------------------------
      implicit none
      type (SACTrace) :: h1,h2,v ! the re-ordered traces for analysis
      type (SACTrace) :: h1_corr,h2_corr ! the corrected traces
      type (SACTrace) :: h1_proc,h2_proc ! traces to run in the analysis
      
      real plat
      real lam1,lam2,lam1_corr,lam2_corr,eigrat,eigrat_corr
      character*50 fname

      open(iulog,file='sheba.log')

!  ** read the input files for analysis
      call get_traces(h1,h2,v)      

!  ** check the trace headers to make sure there is 
!  ** sufficient information to perform the splitting analysis
      call check_trace_headers(h1)      
      call check_trace_headers(h2)      
      call check_trace_headers(v)      

!  ** check window length
      call check_windows(h1)

!  ** save the input polarisation
      config % input_h1_pol = h1 % cmpaz

!  ** orient the traces to NE reference frame
      call f90sac_orient2d(h1,h2)

     
!  ** copy the traces before UMA correction
      call f90sac_clonetrace(h1,h1_proc)
      call f90sac_clonetrace(h2,h2_proc)

!  **  if required pre-correct for UMA
      if (config % iuma_corr == 1) then
         call desplit(h1_proc,h2_proc,config % uma_phi,config % uma_dt)
      endif

!  **  run the Silver and Chan analysis
      call cluster_split(h1_proc,h2_proc)

!  ** output information for plotting
      call output_gmt_info()      

!  ** upload the splitting parameters into the trace headers
      call upload_splitpar(h1)
      call upload_splitpar(h2)
      call upload_splitpar(v)


!  ** create clone traces to correct
      call f90sac_clonetrace(h1,h1_corr)
      call f90sac_clonetrace(h2,h2_corr)


!  **  if required pre-correct for UMA
      if (config % iuma_corr == 1) then
         call desplit(h1,h2,config % uma_phi,config % uma_dt)
      endif

!  **  if required pre-correct for UMA
      if (config % iuma_corr == 1) then
         call desplit(h1_corr,h2_corr,config % uma_phi,config % uma_dt)
      endif
      
!  ** correct the traces for the best splitting anisotropy
      call desplit(h1_corr,h2_corr,event % fast, event % tlag)

!  ** source anisotropy pre-correction
      if (config % i_src_corr == 1) then
         call desplit(h1_corr,h2_corr,  & 
        config % src_fast, config % src_tlag)
      endif
      
!  ** compare before / after eigenvalue ratios
      call windowed_coveig(h1,h2,lam1,lam2,event % wbeg,event % wend)
      call windowed_coveig(h1_corr,h2_corr, & 
                     lam1_corr,lam2_corr,event % wbeg,event % wend)

      eigrat = lam2/lam1 * 100.0
      eigrat_corr = lam2_corr/lam1_corr * 100.0

!  ** save in the event data structure
      event % eigrat_orig = lam2/lam1
      event % eigrat_corr = lam2_corr/lam1_corr

      write(*,'(a,f6.3,a,f6.2,a,f6.2,a,a,f6.2)') &
         ' Q = ', event%Quality, &
         '  L2/L1(%) =  ', & 
         eigrat,' (pre-corr.) ', & 
         eigrat_corr,' (post-corr.)', &
         '  SNR = ',event%snr                          
      write(*,"(80('='))")

!  ** output result information
      call output_result(h1) 
      
!  ** output extra information
      
      

!  ** output particle motion files
      fname = trim(config % fname_base) // '.xy1'
      call output_pm_files(h1,h2,fname)

      fname = trim(config % fname_base) // '.xy2'
      call output_pm_files(h1_corr,h2_corr, fname)

!  ** write out the corrected traces
      call output_traces(h1_corr,h2_corr,v,'_corr')      

!  ** destroy working traces
      call f90sac_deletetrace(h1)
      call f90sac_deletetrace(h2)
      call f90sac_deletetrace(h1_corr)
      call f90sac_deletetrace(h2_corr)
      call f90sac_deletetrace(h1_proc)
      call f90sac_deletetrace(h2_proc)

      return
      end subroutine sheba
!===============================================================================
