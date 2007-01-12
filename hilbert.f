!=======================================================================
!     S H E B A - Shear-wave Birefringence Analysis
!=======================================================================
!     Fortran 90/95 Source Code File
!-----------------------------------------------------------------------
!
!     PROGRAM : sheba
!     FILE    : scs.f
!     AUTHOR  : James Wookey
!     PLACE   : School of Earth Sciences, University of Leeds
!     DATE    : December 2003
!     PURPOSE : Subroutine for correcting ScS for the reflection 
!     VERSION : 1.0
!     COMPLETE: No
!     COMMENTS: 
!
!-----------------------------------------------------------------------
!     This software is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!-----------------------------------------------------------------------
      
      
C=======================================================================
      subroutine hilbert(x,n,xh)
C=======================================================================
C
C     Perform a Hilbert transform on a real trace and return the
C     imaginary part of the analytical signal
C  
C     Uses: four1.f    
C  
C-----------------------------------------------------------------------
      use array_sizes ! use the array_sizes modules
C-----------------------------------------------------------------------
      implicit none
      integer n,n2
      integer i
      real x(np),xh(np)
      complex xfd(np)
      
      
      
C     * get a power of 2 for n2
      n2 = 1
      do i = 1,100
         if (n2>np) then
            print*,'HILBERT: np too small',n,n2,np
            stop
         endif
C              
         if (n2>(n*2)) exit
         n2 = n2 * 2
      enddo ! i = 1,100          

C     * initialize arrays      
c      print*,n2,n
      do i=1,np
         if (i<n) then
            xfd(i) = cmplx(x(i),0.0)
            xh(i) = 0.0
         else
            xfd(i) = cmplx(0.0,0.0)    
         endif      
      enddo ! i=1,np

C     * FFT
      call four1(xfd,n2,1)

C     * Hilbert transform, w>0.0  * -i , w=0.0  * 0.0 , w<0.0  * +i      
      xfd(1) = cmplx(0.0,0.0)
 
      do i=2,n2/2-1
         xfd(i) = xfd(i) * cmplx(0.,-1.) 
      enddo
 
      do i=n2/2,n2
         xfd(i) = xfd(i) * cmplx(0.,1.) 
      enddo

C     * IFFT
      call four1(xfd,n2,-1)
      
      do i=1,n
         xh(i)=real(xfd(i))/real(-n2)
      enddo
      
      return
      end
C=======================================================================

C=======================================================================
      SUBROUTINE FOUR1(DATA,NN,ISIGN)
C=======================================================================
c
c numerical recipes fft subroutine
c nn=2**n !!!!, devide by nn on inverse transform
c
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
      DIMENSION DATA(*)
c     had to disable this ...
c      if (mod(log(float(nn)),log(2.)).ne.0.0) then
c         print*,'nn=',nn,mod(log(float(nn)),log(2.))
c        stop'nn not a power of 2 in four1'
c      end if
      N=2*NN
      J=1
      DO 11 I=1,N,2
        IF(J.GT.I)THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
1       IF ((M.GE.2).AND.(J.GT.M)) THEN
          J=J-M
          M=M/2
        GO TO 1
        ENDIF
        J=J+M
11    CONTINUE
      MMAX=2
2     IF (N.GT.MMAX) THEN
        ISTEP=2*MMAX
        THETA=6.28318530717959D0/(ISIGN*MMAX)
        WPR=-2.D0*DSIN(0.5D0*THETA)**2
        WPI=DSIN(THETA)
        WR=1.D0
        WI=0.D0
        DO 13 M=1,MMAX,2
          DO 12 I=M,N,ISTEP
            J=I+MMAX
            TEMPR=SNGL(WR)*DATA(J)-SNGL(WI)*DATA(J+1)
            TEMPI=SNGL(WR)*DATA(J+1)+SNGL(WI)*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
12        CONTINUE
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
13      CONTINUE
        MMAX=ISTEP
      GO TO 2
      ENDIF
      RETURN
      END
C=======================================================================

