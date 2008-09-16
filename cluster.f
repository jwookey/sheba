!=======================================================================
!  S H E B A - Shear-wave Birefringence Analysis
!=======================================================================
!  This software is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!=======================================================================
!
!  James Wookey, School of Earth Sciences, University of Bristol
!  CVS: $Revision: 1.3 $ $Date: 2008/09/16 17:05:28 $
!
!-----------------------------------------------------------------------
!
! SUBROUTINES IN THIS FILE WERE WRITTEN BY N. TEANBY UNIVERSITY OF LEEDS
! BASED ON TEANBY AND KENDALL (2003) (?)
! 
c     subroutine cluster_split(event,trace1,trace2,
c    >nwbeg,nwend,dtlag_max,dfast_max,tlag_scale,fast_scale,
c    >max_no_clusters,nmin,
c    >OPT_verbose,OPT_outfiles,
c    >wbeg_best,wend_best,ierr)
c-----------------------------------------------------------------------
		subroutine cluster_split(trace1,trace2)
c-----------------------------------------------------------------------
c		AUTOMATIC SEISMIC SPLITTING (ASS)
c
c		subroutine to automatically find the optimum shear wave window for
c		use in shear wave splitting.
c
c		the method is based on a grid search of a range of windows. for each 
c		trial window shear wave splitting is performed. the resulting set of  
c		tlag and fast directions is subjected to cluster analysis in order to  
c		find stable regions of the window space. the optimum window is the  
c		one used in the best measurement from the best cluster.
c
c		windowing:
c		||-------|------------||--------|---------------------------||
c		wbeg spick-t_off_beg  spick    spick + t_off_end            wend
c		
c		-The s-wave window is from wbeg to wend (in seconds)
c		-The position of the window is based on the S-wave onset time (spick)
c		which is read from the SAC header (T5)
c		-The analysis is cycled over nwbeg window beginnings and nwend window
c		endings (nwindows=nwbeg*nwend in total). The end of the window is most
c		critical so nwend>nwbeg normally.
c		-steps wbeg by dt_beg, and wend by dt_end seconds
c
c		-splitting analysis is done for each of the windows
c
c		-cluster analysis used to cluster the data, find optimum no. of clusters
c		and the best measurement within the best cluster
c
c
c     This version has been modified for use with SHEBA. It should no longer
c     be considered compatible with the original code. It also now uses
c     the SAC library f90sac. 
c
c		variables
c		---------
c    in:
c		event				char*50		name of event to analyse
c		ext1/2		char*50		extention of the two components to do analysis on
c		lu				int				logical unit to open files on
c
c     h1,h2       SACTrace       SAC traces to operate on
c
c		nwbeg				int				number of start positions for S-wave window
c		nwend				int				number of end positions for S-wave window
c		dt_beg		real				increment for start of window (in seconds)
c		dt_end		real				increment for end of window (in seconds)
c		t_off_beg		real				maximum window beginning (relative to spick)
c		t_off_end		real				minimum window end (relative to spick)
c		dtlag_max		real				max allowable error in lag time for inclusion 
c										in clustering
c		dfast_max		real				" fast direction "
c		tlag_scale		real				range of tlag scale in seconds
c		fast_scale		real				range of fast direction scale in degrees
c		max_no_clusters		int		max. number of clusters
c		nmin				int				minimum number of points in an acceptable cluster
c		OPT_verbose		logical		true for verbose output
c		OPT_outfiles logical		true if write outfiles for gmt plots
c    out:
c		wbeg_best		real				optimum start of S-wave window
c		wend_best		real				optimum end of S-wave window
c		ierr				int				=0 if there is an acceptable solution
c										=1 if no acceptable solutions
c    other:
c		x0/y0(np)		real				x and y seismic data
c		n				int				number of data points
c		ppick				real				p-wave pick read from sac header
c		spick				real				s-wave pick read from sac header
c		delta				real				sampling interval (s) read from sac header
c		as				real				beginning of s-wave window (hand pick, not used)
c		fs				real				end of s-wave window (hand pick, not used)
c		phi				real				rotation angle (clock from north) read from header
c		theta				real				rotation angle from vertical read from header
c				phi/theta define the frame relative to the ENZ reference frame
c		fast(npc)		real				fast direction for ith window
c		dfast(npc)		real				s.d. of fast direction for ith window
c		tlag(npc)		real				tlag direction for ith window
c		dtlag(npc)		real				s.d. of tlag direction for ith window
c
c-----------------------------------------------------------------------
c		n.teanby		20-8-02		original code
c     j.wookey    02-02-04    modified for SHEBA
c-----------------------------------------------------------------------
      use f90sac  ! use the f90sac module
      use sheba_config ! use the sheba_config module
      use event_info ! use the event_info module
      use array_sizes ! use the array_sizes modules
C-----------------------------------------------------------------------
		implicit none
		integer n,f,ndf,nwindows,nwbeg,nwend,ierr,lu
		logical OPT_verbose,OPT_outfiles
		real wbeg_best,wend_best,wbeg(npc),wend(npc)
		real delta,b,as,fs,ppick,spick
		real fast(npc),dfast(npc),tlag(npc),dtlag(npc)
		real fast_best,dfast_best,tlag_best,dtlag_best
		real spol_best,dspol_best
		real spol,dspol
c      real theta,ev_x,ev_y,ev_z,ev_time,phi
		real lambda2_min,tlag_min,fast_min,dtlag_max,dfast_max
		real lam2m(npc), min_lam2m ! required for second multiple windowing
		integer imin_lam2m
		real error(np1,np2),error_int(np1,np2int)
		real x0(np),y0(np),tlag_scale,fast_scale
		real xc0(npc),yc0(npc),vxc0(npc),vyc0(npc)
		real fastSTR,dfastSTR,fastDIP,dfastDIP
		real spolSTR,dspolSTR,spolDIP,dspolDIP
		character*50 event_name
      type (SACTrace) trace1,trace2
		integer i,j,l,ibest,kbest,k
		integer nmin,max_no_clusters
		real dt_beg,dt_end
		integer cluster(npc,npc)
		real t_off_beg,t_off_end
      real temp_wbeg,temp_wend
		character*50 file_clustxy,file_clusters,file_soln,file_log,ext
      character*50 file_error
      character*12 fmt
c
      integer itlag_step
      real test_tlag,test_fast
      real snr ! signal to noise ratio

      write(iulog,*) 'Running multiple window shear-wave splitting' 
      write(*,*) '> Running multiple window shear-wave splitting' 

c     ** hardwire some parameters **
      event_name = trim(config % fname_base)
      nwbeg = config % nwindbeg
      nwend = config % nwindend
      dtlag_max=config % dtlag_max
      dfast_max=config % dfast_max
      nmin = config % min_pts_per_cluster
      max_no_clusters = config % max_no_clusters
      lu = 99
 		tlag_scale = config % max_tlag
      fast_scale = config % fast_scale
      
      OPT_verbose = .true.
      
c  ** calc the grid spacing **
		tlag_min = tlag_scale/real(np2int)
		fast_min = fast_scale/real(np1)

c  ** cluster analysis parameters
c     subroutine cluster_split(event,trace1,trace2,
c    >nwbeg,nwend,dtlag_max,dfast_max,tlag_scale,fast_scale,
c    >max_no_clusters,nmin,
c    >OPT_verbose,OPT_outfiles,
c    >wbeg_best,wend_best,ierr)
           

C  NEED TO RECODE THIS SECTION !! Populate 		
c  ** read in data **
c		file1=fstrcat(event,ext1)
c		file2=fstrcat(event,ext2)
c		call zreaddataSAC(file1,file2,lu,.true.,
c     >x0,y0,n,np,b,delta,as,fs,ppick,spick,phi,theta,
c     >ev_x,ev_y,ev_z,ev_time)

C     NEED TO SET x0,y0,n,b,delta,as,fs,ppick spick 
c		dt_beg		real				increment for start of window (in seconds)
c		dt_end		real				increment for end of window (in seconds)
c		t_off_beg		real				maximum window beginning (relative to spick)
c		t_off_end		real				minimum window end (relative to spick)

C
C  The beginning and end window postions are required to be stored in 
C  user0 (beg_beg) user1(beg_end) user2(end_beg) user3(end_end) 
C      
      n = trace1 % npts 
      do i = 1,n
         x0(i) = trace2 % trace(i) ! EAST
         y0(i) = trace1 % trace(i) ! NORTH
      enddo
      b = trace1 % b
      delta = trace1 % delta

C     * calculate window dt
      if (nwbeg==1) then
         dt_beg = 0.0
      else
         dt_beg = ((trace1 % user1)-(trace1 % user0)) / real(nwbeg-1)
      endif   

      if (nwend==1) then
          dt_end = 0.0
      else
         dt_end = ((trace1 % user3)-(trace1 % user2)) / real(nwend-1)
      endif   

      t_off_beg = trace1 % user0
      t_off_end = trace1 % user2


      
c  ** calc number of windows and check it's not too many **
		nwindows=nwend*nwbeg
		print*,'no. windows to cycle = ',nwindows
		if (nwindows.gt.npc) then
		   pause 'ERROR: zass: number of windows to search is too large'
		endif

C     ** run the cluster analysis UNLESS only one window was specified
      if (nwbeg==1 .and. nwend==1) then

         kbest=1
         ibest=1
         wbeg(1) = t_off_beg
         wend(1) = t_off_end
         
      else

C  ** run the cluster analysis 
c  ** loop over window end points **
		do i=1,nwend
		   do j=1,nwbeg
		      l = (i-1)*nwbeg + j
C        ** refine window to be relative to 4 header values
            wbeg(l) = t_off_beg + real(j-1)*dt_beg
 		      wend(l) = t_off_end + real(i-1)*dt_end
c            print*,l,wbeg(l),wend(l),temp_wbeg,temp_wend
		      call zsplit(x0,y0,n,wbeg(l),wend(l),delta,b,tlag_scale,
     >                  fast(l),dfast(l),tlag(l),dtlag(l),
     >                  spol,dspol,error,error_int,f,
     >                  lam2m(l),ndf,snr)
         
            if (l==1) then
               min_lam2m = lam2m(l)
               imin_lam2m = 1
            else
               if (lam2m(l)<min_lam2m) then
                  min_lam2m = lam2m(l)
                  imin_lam2m = l
               endif   
            endif      
		      write(iulog,99) l,wbeg(l),wend(l),tlag(l),
     >         dtlag(l),fast(l),dfast(l)

		      write(*,99) l,wbeg(l),wend(l),tlag(l),
     >         dtlag(l),fast(l),dfast(l)

99		      format(i3,2x,f8.3,' -',f8.3,'s  tlag =',f8.5,' +/-',f8.5,
     >      '  fast =',f8.3,' +/-',f8.3)
            
		   enddo
		enddo

c-----------------------------------------------------------------------
      if (config % WindowSelectionMode .eq. 1) then
c
c     Select window based on cluster analysis of splitting results 
c     (ASS method)
c
		   print*,'Done all windows, running cluster analysis'
         
         
		   call zpackresults(wbeg,wend,tlag,dtlag,fast,dfast,nwindows,
     >   npc,dtlag_max,dfast_max)
         
c  **    do cluster analysis and find number of clusters **
		   call zcluster(tlag,dtlag,fast,dfast,nwindows,
     >   		tlag_scale,fast_scale,tlag_min,fast_min,max_no_clusters,
     >   		xc0,yc0,vxc0,vyc0,cluster,k)
         
         
c  **    find best cluster **
		   call zselect_cluster(dtlag,dfast,vxc0,vyc0,nwindows,
     >   		tlag_scale,fast_scale,cluster,nmin,k,
     >   		kbest)
         
c  **    find best measurement **
		   call zselect_measurement(dtlag,dfast,nwindows,
     >   		wbeg,wend,delta,spick,tlag_scale,fast_scale,cluster,k,
     >   		kbest,ibest)
         
c  **    write out measurements to file **
         file_clustxy = trim(config % fname_base) // '.clustxy'
		   call zwrite_clustxy(lu,file_clustxy,nwindows,npc,
     >   		wbeg,wend,fast,dfast,tlag,dtlag)
         
c  **    write out clusters to file **
         file_clusters = trim(config % fname_base) // '.clusters'
		   call zwrite_clusters(lu,file_clusters,k,npc,
     >   		xc0,yc0,vxc0,vyc0)

      elseif (config % WindowSelectionMode .eq. 2) then
c
c     Select window based on best minimisation of lambda2
c
c     NOTE, DEVELOPMENT OF THIS OPTION HAS BEEN SUSPENDED. IT MAY BE USEFUL
c     UNDER CIRCUMSTANCES, BUT NOT THE ONES IT WAS ORIGINALLY INTENDED. THE
c     FRAMEWORK CODE HAS BEEN LEFT IN IN CASE OF FUTURE DEVELOPMENT OF 
c     DIFFERENT WINDOW SELECTION ALGORITHMS. 
c
c     Therefore WindowSelectionMode is hardwired to 1 in sheba_config.
c
		   print*,'Done all windows, selecting on minimum lambda2'
		   print*,' '
		   print*,'!!! Warning !!!'
		   print*,' '
		   print*,'This algorithm is under development, and should be'
		   print*,'used only with all due caution.'
		   
         ibest = imin_lam2m
         
      endif ! of MULTIPLE WINDOW MODE IF
c-----------------------------------------------------------------------------

      endif ! end of IF MULTIPLE WINDOWS

c  ** do splitting and write log file **		
		if (kbest.eq.0 .and. config % WindowSelectionMode.eq.1) then
		   ierr = 1
		   wbeg_best=0.
		   wend_best=0.
		   print*,'NO CLUSTERS MEET THE CRITERIA'
		   print*,'NO SOLUTION FOUND FOR: ',event_name
         STOP
		else
      
! to implement stacking, need to use k (optimum number of clusters), kbest (best
! cluster (based on what?)), cluster (assignment of windows to particular 
! clusters). So windows where cluster(l,k) == kbest should be included
!
! could also use the stacked result to select cluster? no ...
!     
      
		   ierr=0
		   wbeg_best=wbeg(ibest)
		   wend_best=wend(ibest)
c  		** do splitting analysis **
		   call zsplit(x0,y0,n,wbeg_best,wend_best,
     >   		delta,b,tlag_scale,
     >   		fast_best,dfast_best,tlag_best,dtlag_best,
     >   		spol_best,dspol_best,error,error_int,f,
     >   		lambda2_min,ndf,snr)

      if (nwbeg==1 .and. nwend==1) then
 
c  ** write out empty cluster analysis files ONE WINDOW ONLY
         file_clustxy = trim(config % fname_base) // '.clustxy'
         file_clusters = trim(config % fname_base) // '.clusters'
      
         open(99,file=file_clusters)
         write(99,*) tlag_best,fast_best,0.,0.
         close(99)

         open(99,file=file_clustxy)
		   write(99,100) 1,wbeg_best,wend_best,fast_best,dfast_best,
     >                tlag_best,dtlag_best
100	   format(i6,2f12.4,f8.3,f7.3,2f10.6)	

         close(99)
      endif
      
c  ** upload splitting parameters to event_info modules **
      event % tlag = tlag_best
      event % dtlag = dtlag_best
      event % fast = fast_best
      event % dfast = dfast_best
      event % spol = spol_best
      event % dspol = dspol_best

      event % wbeg = wbeg_best
      event % wend = wend_best
      
      event % snr = snr
      event % ndf = ndf

      event % ibest = ibest
      
c  ** calc itlag_step from tlag_scale **
c  ** itlag_step is the grid spacing in tlag for the grid search **
c  ** it must be an integer greater than 1 **
		itlag_step = nint( tlag_scale/(real(np2-1)*delta) )
      
      print*,f,itlag_step
      
      file_error = trim(config % fname_base) // '.error'
		open(lu,file=file_error)
		do i=1,np1
		   do j=1,np2int
		     test_tlag = delta*real((j-1)*itlag_step)/real(f)
		     test_fast = -90. + 180.*real(i-1)/real(np1-1)
		     write(lu,*) test_tlag,test_fast,error_int(i,j)
		   enddo
		enddo
		close(lu)
      event % error_grid_tlag_int = real(itlag_step)*delta/real(f)
      file_error = trim(config % fname_base) // '.lam2'
		open(lu,file=file_error)
      write(lu,'(2i5,a)') np1,np2int,'   % NPfast,NPtlag'
      write(lu,'(f12.4,a)') event % error_grid_tlag_int,'   % dtlag'
      write(lu,'(i5,a)') ndf,'   % NDF'
      write(lu,'(f12.4,a)') snr,'   % SNR'
      
      write(fmt,'(a1,i5.5,a)') '(',np2int,'f12.4)'
      print*,fmt
		do i=1,np1
		     write(lu,fmt) (error_int(i,j),j=1,np2int) 
		enddo
		close(lu)



c		** print output message **		   
       write(*,199)
		  print*,' RESULT for:',event_name
       write(*,199)
		  print*,'s-wave window =',wbeg_best,' - ',wend_best,' seconds'
        if (config % i_rotate_to_ABC == 1) print*,'[Rotated frame]'
		  print*,'lag      =',tlag_best,'+/-',dtlag_best,' seconds'
		  print*,'fast     =',fast_best,'+/-',dfast_best,' degrees'
		  print*,'spol     =',spol_best,'+/-',dspol_best,' degrees'
        write(*,198)

199     format(80('='))
198     format(80('-'))
		endif


C     


		return
		end
