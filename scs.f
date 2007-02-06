!=======================================================================
!  S H E B A - Shear-wave Birefringence Analysis
!=======================================================================
!  This software is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!=======================================================================
!
!  James Wookey, School of Earth Sciences, University of Bristol
!  CVS: $Revision: 1.2 $ $Date: 2007/02/06 14:03:49 $
!
     
      
C=======================================================================
      subroutine scscorr(x,y,xh,yh,n,fast,ilag)
C=======================================================================
C
C     Correct for the effects of the ScS reflection. If iscs_corr is set
C     to zero, do nothing     
C
C     parameters:
C  
C     x,y : two time series to operate on, with n points
C     xh,yh : the Hilbert transforms of x and y
C     fast : fast direction azimuth
C 
C-----------------------------------------------------------------------
      use f90sac  ! use the f90sac module
      use sheba_config ! use the sheba_config module
      use event_info ! use the event info module
      use array_sizes ! use the array_sizes modules
C-----------------------------------------------------------------------
      implicit none
      integer n
C  * coefficients for ScS reflection calculation
      complex rsv, rsh, r11, r12, r21, r22, temp1, temp2
      complex rr(2,2),rri(2,2) ! complex reflection matrix
      real rri_real(2,2),rri_imag(2,2)
      integer i,j
      real x(np),y(np),xc,yc,xh(np),yh(np)
      real p_lat,fast ! phase slowness (s/km) and fast-direction
      integer ilag
      real radial_dir
      real orientation ! the angle between the radial and phi
      
      if (config % iscs_corr == 1) then     
      
c     ** calculate p_lat (s/km) from slowness (s/deg)     
         p_lat = config % scs_slw / to_km
c     ** calculate angle between fast-direction and radial
         radial_dir = event % baz + 180.0
         call unwind(radial_dir)
         orientation = fast - radial_dir
c         print*,orientation,n,p_lat,ilag

c     ** calculate SV / SH reflection coefficient
         call refsub(p_lat,rsv,rsh)      

c     ** calculate reflection matrix for fast-slow reference frame
         call dppref(orientation,rsv,rsh,r11,r12,r21,r22)

         rr(1,1)=r11
         rr(1,2)=r12
         rr(2,1)=r21
         rr(2,2)=r22

C     ** calculate inverse matrix
         call xinv2(rr,rri)
         
         rri_real(1:2,1:2) = real(rri(1:2,1:2))
         rri_imag(1:2,1:2) = aimag(rri(1:2,1:2))    

C     ** apply correction
         do i=1,n
            yc = y(i)*rri_real(1,1) +
     >              x(i)*rri_real(1,2) + 
     >             yh(i)*rri_imag(1,1) + 
     >             xh(i)*rri_imag(1,2)
            xc = y(i)*rri_real(2,1) +
     >              x(i)*rri_real(2,2) + 
     >             yh(i)*rri_imag(2,1) + 
     >             xh(i)*rri_imag(2,2)
            y(i)=yc              
            x(i)=xc
         enddo  

      endif
c  ** done      

      return
      end
