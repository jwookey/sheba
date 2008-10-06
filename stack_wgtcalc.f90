!===============================================================================
!  S H E B A - Shear-wave Birefringence Analysis
!===============================================================================
!  This software is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!===============================================================================
!
!  James Wookey, School of Earth Sciences, University of Bristol
!  CVS: $Revision: 1.2 $ $Date: 2008/10/06 11:04:40 $
!
!-------------------------------------------------------------------------------
!
!   STACK_WGTCALC - This calculates stacking weights for splitting error 
!   surfaces based on a list of parameters supplied on the standard input.
!   Command line options define how the parameters are interpreted. The 
!   output is a list of normalised aggregate weights. 
!
!===============================================================================
   program stack_wgtcalc
!===============================================================================
   implicit none
      integer i,iskip,j,ip
      integer :: iargc
      real idfast,idtlag
      character (len = 12) :: fmt
      character (len = 80) :: arg,wfname ! weight file name
      
      integer,parameter :: nfmax = 1000
      integer,parameter :: npmax = 1000
      integer :: np,nf,ibin
      integer :: iwmode(npmax) ! 0 = raw, 1 = SNR, 2 = SPOL
      real :: wgts(nfmax,npmax),awgt(nfmax)
      real :: params(nfmax,npmax)
      
      real :: funwgt,ang0_180 ! SNR weighting function, angle unwinder
      integer :: ninbin(6) ! 30 degree angular bins
      real :: srcpol
      
      if (iargc()==0) then
         write(0,'(a)') &
'STACK_WGTCALC - This calculates stacking weights for splitting error  ',&
'surfaces based on a list of parameters supplied on the standard input.',&
'Command line options define how the parameters are interpreted. The   ',&
'output is a list of normalised aggregate weights.',&
'',&
'Usage: stack_wgtcalc [parameter flags]',&
'Available parameter flags:',&
'     -raw  : Use this parameter without modification (except normalisation',&
'             based on the maxium value',&
'     -snr  : Restivo and Helffrich (1999) SNR weight',&
'     -spol : Weight based on binning of source polarisation. This',&
'             provides a mechanism for correcting for oversampling',&
''
          stop
      endif
      
!  ** parse command line options
      np = 0
      iskip = 0      
      do 5 i=1,iargc()
         if (i .le. iskip) go to 5
         call getarg(i,arg)
         if (arg(1:4) == '-raw') then
            np = np + 1
            iwmode(np) = 0
            iskip = i
         elseif (arg(1:4) == '-snr') then
            np = np + 1
            iwmode(np) = 1
            iskip = i
         elseif (arg(1:5) == '-spol') then
            np = np + 1
            iwmode(np) = 2
            iskip = i
         else
            write(0,*) '**Unrecognized: ',arg(1:index(arg,' '))
            stop   
         endif
5     continue

!  ** now read the parameter table
      nf = 0
      do ! forever
         nf = nf + 1
         read(*,*,end=100,err=900) (params(nf,j),j=1,np)
      enddo   
100   nf = nf - 1         
      
!==============================================================================
!  ** process the list, working over the columns
      do ip=1,np
         if (iwmode(ip)==0) then
!        ** (normalised) raw weights
            wgts(1:nf,ip) = params(1:nf,ip) / maxval(params(1:nf,ip))
!            wgts(1:nf,ip) = (wgts(1:nf,ip)-minval(params(1:nf,ip))) / &
!                      (maxval(params(1:nf,ip))-minval(params(1:nf,ip)))

!==============================================================================
         elseif (iwmode(ip)==1) then
!        ** Restivo and Helffrich SNR weight
            do i=1,nf
               wgts(i,ip) = funwgt(params(i,ip)) 
            enddo          
!==============================================================================
         elseif (iwmode(ip)==2) then
!        ** Weight-based on the source polarisation
!           first fill the bins
            ninbin(:) = 0
            do i=1,nf
               srcpol = ang0_180(params(i,ip)) 
               ibin = 1+int(srcpol/30.0)
               ninbin(ibin) = ninbin(ibin) + 1 
            enddo          
!           then generate the weights
            do i=1,nf
               srcpol = ang0_180(params(i,ip)) 
               ibin = 1+int(srcpol/30.0) 
               wgts(i,ip) = 1.0 / real(ninbin(ibin)) 
            enddo          
!==============================================================================
         else
            write(0,'(a)') 'Unsupported weight mode'
            stop            
         endif      
!==============================================================================
      enddo   

!  ** combine the weights in the columns (product)
      awgt(1:nf) = 1.0
      do ip=1,np
         awgt(1:nf) = awgt(1:nf)*wgts(1:nf,ip)    
      enddo   

!  ** normalise the sum of the weights to 1.0
      awgt(1:nf) = awgt(1:nf)/sum(awgt(1:nf))
      
!  ** output         
      do i=1,nf
!         write(*,*) (wgts(i,ip),ip=1,np)
         write(*,'(f10.8)') awgt(i)
      enddo   
      
      stop
900   write(0,*) 'Error reading parameter list.'      
   end program stack_wgtcalc
!===============================================================================

!===============================================================================
   function ang0_180(x)
!===============================================================================
!  Put angle into 0-180 range (assuming symmetry)
!-------------------------------------------------------------------------------
      real x
      do 
         if (x>=0. .and. x<180.) exit ! done
         if (x>360.0) x=x-360.
         if (x<0.) x=x+360.
         if (x>=180.0) x = x - 180.
      enddo
      ang0_180 = x   
      return
   end function ang0_180
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
