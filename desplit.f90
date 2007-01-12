!=======================================================================
!     S H E B A - Shear-wave Birefringence Analysis
!=======================================================================
!     Fortran 90/95 Source Code File
!-----------------------------------------------------------------------
!
!     PROGRAM : sheba
!     FILE    : desplit.f
!     AUTHOR  : James Wookey
!     PLACE   : School of Earth Sciences, University of Leeds
!     DATE    : December 2003
!     PURPOSE : Various subroutines for sheba
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
      subroutine desplit(h1,h2,phi,dt)
!=======================================================================
!
!     Remove shear-wave splitting
!
!     parameters:
!  
!     h1,h2 :  (I/O) (SACtrace) input horizontal, orthogonal components
!     dt,phi:  (I)   (real) lag-time, fast direction to correct by
!
!-----------------------------------------------------------------------
      use f90sac  ! use the f90sac module
      use sheba_config ! use the sheba_config module
!-----------------------------------------------------------------------
      implicit none
      
      type (SACtrace) :: h1, h2
      real dt , phi
      real cmpazdiff,rotang

      write(iulog,*) 'Correcting input by (phi,dt)',phi,dt

!  ** calculate the angle of rotation, based on phi and the azimuth of h1
      rotang = (phi - h1 % cmpaz)
   
!  ** rotate the traces      
      call f90sac_rotate2d(h1,h2,rotang)      

!  ** time shift the slow trace
      call f90sac_tshift(h2,-dt)

!  ** rotate back to original orientation    
      call f90sac_rotate2d(h1,h2,-rotang)      

!  * done      
      return
      end
!=======================================================================

!=======================================================================
      subroutine desplit_scs(h1,h2,phi,dt,plat)
!=======================================================================
!
!     Remove shear-wave splitting - special version for ScS phase
!
!     parameters:
!  
!     h1,h2 :  (I/O) (SACtrace) input horizontal, orthogonal components
!     dt,phi:  (I)   (real) lag-time, fast direction to correct by
!     plat: (I) (real) horizontal slowness in s/km
!
!-----------------------------------------------------------------------
      use f90sac  ! use the f90sac module
      use sheba_config ! use the sheba_config module
      use array_sizes ! use the array_sizes modules
!-----------------------------------------------------------------------
      implicit none
      type (SACtrace) :: h1, h2
      real dt , phi, plat
      real cmpazdiff,rotang

      real h1r(np),h2r(np),h1h(np),h2h(np),h1c(np),h2c(np)
      integer n,i

!  ** reflection matrix variables
      complex rsv, rsh, r11, r12, r21, r22, temp1, temp2
      complex rr(2,2),rri(2,2) ! complex reflection matrix
      real rri_real(2,2),rri_imag(2,2)

      write(iulog,*) 'Correcting input by (phi,dt) with ScS corr',phi,dt
!      write(*,*) '> Correcting input by (phi,dt)',phi,dt

!  ** calculate the angle of rotation, based on phi and the azimuth of h1
      rotang = (phi - h1 % cmpaz)
      
!  ** rotate the traces      
      call f90sac_rotate2d(h1,h2,rotang)
      
!  ** time shift the slow trace (half the time)
      call f90sac_tshift(h2,-dt/2.0)

!  ** calculate Hilbert transforms            
      n = h1 % npts

      do i=1,n
         h1r(i) = h1 % trace(i)
         h2r(i) = h2 % trace(i)
      enddo ! i=1,n 

      call hilbert(h1r,n,h1h)
      call hilbert(h2r,n,h2h)

!  ** dump out to a file
      open(98,file='corr.dat') 
      
      do i=1,n
         write(98,'(4f14.4)') h1r(i),h1h(i),h2r(i),h2h(i)
      enddo ! i=1,n 


!  ** calculate SV / SH reflection coefficients
      call refsub(plat,rsv,rsh)      
!  ** calculate reflection matrix for fast-slow reference frame
      call dppref(rotang,rsv,rsh,r11,r12,r21,r22)

      rr(1,1)=r11
      rr(1,2)=r12
      rr(2,1)=r21
      rr(2,2)=r22

!  ** calculate inverse matrix 
      call xinv2(rr,rri)

!  ** split matrix into real and imaginary parts
      rri_real(1:2,1:2) = real(rri(1:2,1:2))
      rri_imag(1:2,1:2) = aimag(rri(1:2,1:2))    

!  ** output rri matrix to log file
      write(iulog,*) 'ScS reflection matrix...'
      write(iulog,199) rri_real(1,1),rri_imag(1,1), & 
                       rri_real(1,2),rri_imag(1,2)
      write(iulog,199) rri_real(2,1),rri_imag(2,1), & 
                       rri_real(2,2),rri_imag(2,2)
      write(iulog,'(a,f5.2)') ' ScS slowness used =',plat * to_km
      
      if (rri_imag(1,1)==0.0 .and. rri_imag(1,2)==0.0 .and. & 
          rri_imag(2,1)==0.0 .and. rri_imag(2,2)==0.0) then
         write(*,'(a)') ' ScS      =   no phase shift'
      else
         write(*,'(a)') ' ScS      =   phase shift'
      endif   
      
199   format(f7.4,' + ',f7.4,'i  ,   ',f7.4,'+',f7.4,'i')      
      
!  ** apply the correction to the traces using the analytic signals         
      do i=1,n
            h1c(i) = h1r(i)*rri_real(1,1) +  & 
                  h2r(i)*rri_real(1,2) +  & 
                  h1h(i)*rri_imag(1,1) +  & 
                  h2h(i)*rri_imag(1,2)
            h2c(i) = h1r(i)*rri_real(2,1) +  & 
                  h2r(i)*rri_real(2,2) +  & 
                  h1h(i)*rri_imag(2,1) +  & 
                  h2h(i)*rri_imag(2,2)
      enddo

!  ** load the corrected traces back into the SAC traces
      do i=1,n
         h1 % trace(i) = h1c(i)
         h2 % trace(i) = h2c(i)
      enddo    

!  ** final time shift
      call f90sac_tshift(h2,-dt/2.0)
    
!  ** rotate back to original orientation    
      call f90sac_rotate2d(h1,h2,-rotang)      

!  ** done      
      return
      end
!=======================================================================

