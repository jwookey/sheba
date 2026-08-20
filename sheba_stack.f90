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
!-------------------------------------------------------------------------------
!
!      SHEBA_STACK stacks single event error surfaces, after 
!      Wolfe & Silver (1998). Simplified method after 
!      Restivo and Helffrich (GJI, 1999)
!
!      Update: ignore commented lines in input file
!
!      Update: added some command line options:
!              -wgt [snr|filename|one] snr = SNR ala Restivo and Helffrich '99
!                                 filename = file with weights
!                                 one = equal weighting (the default)  
!              -bs N : bootstrapped solution (N times).
!
!===============================================================================
   program sheba_stack
!===============================================================================
   use array_sizes
   use statistics
!-------------------------------------------------------------------------------
   implicit none
      real :: error_in(np1,np2int,nsurfmax), error_stack(np1,np2int)
      integer :: n1,n2,nsurf,ndf_stack,ndf_in(nsurfmax),ndf
      real :: dtlag_step,lam2min,dfast,dtlag,fast,tlag
      real :: wgt(nsurfmax)
      integer :: ind_surf(nsurfmax) ! array containing ordering of surfaces 
      integer ifast,itlag,i,j,iskip
      integer :: iargc
      real idfast,idtlag
      character (len = 12) :: fmt
      character (len = 80) :: arg,wfname ! weight file name

      integer :: iwmode ! 1=even weighting, 2=snr weighting, 3=file of weights

!  ** Bootstrap mode variables
      integer :: nbootstrap ! Number of stacks to bootstrap
      integer :: ibootstrap
      real :: bs_pdf(np1,np2int)
      real(dp) :: bs_fast(nbsmax), bs_tlag(nbsmax)
      real(dp) :: bs_fast_ave,bs_fast_std,bs_tlag_ave,bs_tlag_std


!  ** parse command line options
      iwmode = 1
      nbootstrap = 1
      iskip = 0      
      do 5 i=1,iargc()
         if (i .le. iskip) go to 5
         call getarg(i,arg)
         if (arg(1:4) == '-wgt') then
            call getarg(i+1,arg)
            if (arg(1:3) == 'one') then
               iwmode = 1
            elseif (arg(1:3) == 'snr') then
               iwmode = 2
            else
               iwmode = 3
               wfname = arg
            endif
            iskip = i+1
         elseif (arg(1:3) == '-bs') then  
            call getarg(i+1,arg)
            read(arg,*,err=900) nbootstrap
            iskip = i+1 
         else
            write(0,*) '**Unrecognized: ',arg(1:index(arg,' '))
         endif
5     continue

!  ** load error surfaces from the file
      write(0,'(a)') 'Loading error surfaces ...'
      call load_error_surfaces(iwmode,wfname,nsurf,n1,n2, &
                               dtlag_step,ndf_in,wgt,error_in)

!-------------------------------------------------------------------------------
      if (nbootstrap==1) then ! Normal mode. 
!-------------------------------------------------------------------------------
         write(0,'(a)') 'Stacking ...'

         ind_surf = [(i, i = 1, nsurf)]

!     ** stack the error surfaces
         call stack_error_surfaces(nsurf,n1,n2,error_stack,&
            ind_surf,error_in,ndf_in,wgt,ndf)

!     ** find the minimum position
         call zerror_min(error_stack,np1,np2int,ifast,itlag,lam2min)
         print*,'Minimum at:', -90.0+real(ifast-1),0+real(itlag-1)*dtlag_step
         fast = -90.0+real(ifast-1)
         tlag = 0.0+real(itlag-1)*dtlag_step

!     ** calculate errors from the 95% confidence interval
         call zerror95(error_stack,ndf,lam2min,idfast,idtlag)
!         print*,idfast,idtlag,ndf
         dtlag = idtlag * dtlag_step
         dfast = real(idfast) 
   
         print*,'Minimum at:'
         print*,'   FAST',fast,' +/- ',dfast
         print*,'   TLAG',tlag,' +/- ',dtlag
         
!     ** output results to read in MATLAB
   
         open(31,file='sheba_stack.err')
         open(32,file='sheba_stack.sol')     
   
         write(fmt,'(a1,i5.5,a)') '(',np2int,'f12.4)'
         do i=1,np1
              write(31,fmt) (error_stack(i,j),j=1,np2int) 
         enddo
         write(32,'(a)') '% FAST DFAST TLAG DTLAG, NSTACKED, TLAG_STEP'
         write(32,'(4f12.4,i5,f12.4)') fast,dfast,tlag,dtlag,nsurf,dtlag_step
         close(31)
         close(32)
!-------------------------------------------------------------------------------
      else ! Bootstrap mode
!-------------------------------------------------------------------------------
!     ** init random number generator, date/time based seed        
         call random_init(repeatable=.false.,image_distinct=.true.)
         write(0,'(a)') 'Bootstrapping stacking ...'

!     ** zero PPDF
         bs_pdf(:,:) = 0.0

!     ** loop over the number of BS
         do ibootstrap = 1,nbootstrap

!        ** generate a random index
            call generate_stack_index(nsurf,ind_surf)

            !print*,ind_surf(1:nsurf)
!        ** stack the surfaces
            call stack_error_surfaces(nsurf,n1,n2,error_stack,&
               ind_surf,error_in,ndf_in,wgt,ndf)

!        ** find the minimum position
            call zerror_min(error_stack,np1,np2int,ifast,itlag,lam2min)

!        ** record            
            bs_fast(ibootstrap) = -90.0+real(ifast-1)
            bs_tlag(ibootstrap) = 0+real(itlag-1)*dtlag_step
            bs_pdf(ifast,itlag) = bs_pdf(ifast,itlag)+1.0

         enddo

         write(0,'(a)') 'Bootstrapping complete.'

!     ** calculate mean and std of solutions, normalise PDF
         call axial_mean_std(bs_fast(1:nbootstrap), bs_fast_ave, bs_fast_std)
         call linear_mean_std(bs_tlag(1:nbootstrap), bs_tlag_ave, bs_tlag_std)

         bs_pdf(:,:) = bs_pdf(:,:) / real(nbootstrap)

         print*,'Minimum at:'
         print*,'   FAST',bs_fast_ave,' +/- ',2.0_dp*bs_fast_std
         print*,'   TLAG',bs_tlag_ave,' +/- ',2.0_dp*bs_tlag_std

!     ** output solution, statistics and PDF
         open(31,file='sheba_stack.bs_pdf')
         open(32,file='sheba_stack.bs_sol')     
   
         write(fmt,'(a1,i5.5,a)') '(',np2int,'f12.4)'
         do i=1,np1
              write(31,fmt) (bs_pdf(i,j),j=1,np2int) 
         enddo
         write(32,'(a)') &
            '% FAST DFAST TLAG DTLAG, NSTACKED, TLAG_STEP, NBOOTSTRAP'
         write(32,'(4f12.4,i5,f12.4,i8)') &
            bs_fast_ave, 2.0_dp * bs_fast_std, &
            bs_tlag_ave, 2.0_dp * bs_tlag_std,nsurf,dtlag_step, nbootstrap
         close(31)
         close(32)
      endif  
     
      stop
900   write(0,'(a)') 'Error in command line argument.'      
   end program sheba_stack
!===============================================================================

!===============================================================================
   subroutine generate_stack_index(nsurf,ind_surf)
!===============================================================================
!
!   Generate a randomised stack index for bootstrapping. 
!
   use array_sizes
!-------------------------------------------------------------------------------
   implicit none
      integer :: ind_surf(nsurfmax) ! array containing ordering of surfaces
      integer :: nsurf
      real :: rnum(nsurfmax) 

      call random_number(rnum)
      ind_surf(1:nsurf) = 1 + (rnum(1:nsurf)*nsurf) 

      return
   end subroutine generate_stack_index

!===============================================================================
   subroutine stack_error_surfaces(nsurf,nfast,ntlag,error_stack, &
                                   ind_surf,error_in,ndf_in,wgt_in,ndf_stack)
!===============================================================================
!
!   Stack error surfaces in the error_in 3D array. 
!
   use array_sizes
!-------------------------------------------------------------------------------
   implicit none
      real :: error_in(np1,np2int,nsurfmax), error_stack(np1,np2int)
      real :: error_local(np1,np2int,nsurfmax) ! local copy
      integer :: nfast,ntlag,nsurf,ndf_stack,ndf(nsurfmax),ndf_in(nsurfmax)
      integer :: ind_surf(nsurfmax) ! array containing ordering of surfaces 
      real :: rndf
      integer :: isurf,ifast,itlag

      real :: lam2max,lam2min ! maximum lamba in surface      
      real :: wgt(nsurfmax),wgt_in(nsurfmax)
      
      real :: sumwgt ! sum of weights applied

!  ** create local surface array based on index
      do isurf = 1,nsurf
         error_local(:,:,isurf) = error_in(:,:,ind_surf(isurf))
         ndf(isurf) = ndf_in(ind_surf(isurf))
         wgt(isurf) = wgt_in(ind_surf(isurf))
      enddo
      
!  ** stack error surfaces
      error_stack(:,:) = 0.0
      rndf = 0.0
      sumwgt = 0.0
      
      do isurf = 1 , nsurf
!     ** get lam2max,lam2min for this surface
         lam2max = 0.0  
         lam2min = error_local(1,1,isurf)  
         do ifast = 1,nfast
            do itlag = 1,ntlag
               if (error_local(ifast,itlag,isurf) > lam2max) &
                  lam2max = error_local(ifast,itlag,isurf)
               if (error_local(ifast,itlag,isurf) < lam2min) &
                  lam2min = error_local(ifast,itlag,isurf)
                  
            enddo ! ifast = 1,nfast
         enddo ! itlag = 1,ntlag

!     ** normalize error surface lam2min to 1   
         error_local(:,:,isurf) = error_local(:,:,isurf)/lam2min

!     ** apply weighting (e.g.,Restivo & Helffrich, GJI, 1999)
         error_local(:,:,isurf) = error_local(:,:,isurf)*wgt(isurf)
         
!     ** add to the stack  
         error_stack(:,:) = error_stack(:,:) + error_local(:,:,isurf)
         rndf = rndf + real(ndf(isurf))*wgt(isurf)
         sumwgt = sumwgt + wgt(isurf)
            
      enddo ! isurf = 1 , nsurf
            
!  ** average error surface
      error_stack(:,:) = error_stack(:,:) / sumwgt       
      
!  ** calculate weighted NDF
      rndf = rndf/sumwgt*real(nsurf)
      ndf_stack = nint(rndf)

      return
   end subroutine stack_error_surfaces
!===============================================================================

!===============================================================================
   function funwgt(x)
!===============================================================================
!  Function funwth defines weight to give to samples with a determined S/N rat.
!  It is defined so that samples below S/N* = 1.0 are given a weight that from
!  0.01 tends asynthotically to 0.0, while samples over S/N* = 21.0 are given a
!  weight that from 0.99 tends to 1.0. Values of S/N* in between return a weight
!  which steadily increases with S/N.
!  S/N* indicates the S/N ratio calculated by the shear program (biased higher).
!
!  Restivo + Helffrich, GJI, 1999
!-------------------------------------------------------------------------------
      real k,mu
   
      epstop = 21.0
      epsbot = 1.0
      width = epstop-epsbot
      mu = width/2.
      funetp = 0.99
      funebt = 0.01

      k = (log(funebt**2./funetp**2.))/width

      funwgt = 1./(exp(k*(x-mu))+1.)
      return
   end function funwgt
!===============================================================================


!===============================================================================
   subroutine load_error_surfaces(iwmode,wfname,nsurf,n1,n2,dtlag_step, &
                                   ndf_in,wgt,error_in)
!===============================================================================
!
!    load error surfaces from files listed in sheba_stack.in
!
    use array_sizes
!-------------------------------------------------------------------------------
   implicit none
      character (len = 80) :: filename
      character (len = 80) :: wfname ! weight file name
      real :: error_in(np1,np2int,nsurfmax),snr,wgt(nsurfmax)
      integer :: n1,n2,nsurf,ndf_in(nsurfmax)
      integer :: n1first,n2first,i,j
      real :: dtlag_step,dtlagfirst
      real :: funwgt ! weighting function
      integer :: iwmode ! 1=even weighting, 2=snr weighting, 3=file of weights
         
!  ** loop through file list, and get file names
!  ** for each file read in error surface and NDF 
      open(10,file='sheba_stack.in') ! open error surface file list

!  ** if necessary open the weight file
      if (iwmode==3) then
         open(30,file=wfname)
      endif   

      nsurf = 0
      do ! forever
50       read(10,*,end=100) filename
         nsurf = nsurf + 1
         if (filename(1:1)=='%' .or. filename(1:1)=='!') then ! comment line
            nsurf = nsurf - 1
            goto 50 ! sorry
         endif
         
!     ** check whether NSURF has been exceeded. 
         if (nsurf>nsurfmax) then
            write(0,'(a)') &
            'Maximum number of surfaces exceeded. Check the value of NSURFMAX.'
            stop
         endif
         open(20,file=filename)

         read(20,*) n1, n2
         read(20,*) dtlag_step
!     ** check that the new grid is the same size as the first
         if (nsurf == 1) then
            n1first = n1 ; n2first = n2; dtlagfirst = dtlag_step
         else
            if (n1/=n1first .or. n2/=n2first .or. dtlag_step/=dtlagfirst) then
               print*,'Warning! Grid ',nsurf,' has a different number of nodes'
               print*,'         or a different DTLAG'
               print*,'Ignoring ...'
               nsurf = nsurf - 1
               goto 50 ! sorry
            endif   
         endif
!     ** if we get to here then we can go ahead and read the file
         read(20,*) ndf_in(nsurf)
         read(20,*) snr
!     ** set the weight based on options         
         if (iwmode == 1) then
            wgt(nsurf) = 1.0
         elseif (iwmode == 2) then
            wgt(nsurf) = funwgt(snr)
            print*,wgt(nsurf),snr
         else
!     ** read in the weight
            read(30,*) wgt(nsurf)
         endif
         
         print*,'File ',nsurf, ' SNR ',snr, ' IWMODE ',iwmode
         print*,'N1 N2 DTLAG NDF WGT',n1,n2,dtlag_step,ndf_in(nsurf),wgt(nsurf)
         do i=1,n1
            read(20,*) (error_in(i,j,nsurf),j=1,n2)
         enddo
         
         close(20)   
      enddo   
100   continue      
!  ** loop is finished
      
      print*,'Read ',nsurf, ' files'

!  ** produce the error stack

      return
   end subroutine load_error_surfaces
!===============================================================================
