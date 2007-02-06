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
!-----------------------------------------------------------------------
!
!     SUBROUTINES IN THIS FILE WERE WRITTEN BY G.RUEMPKER, GFZ POTSDAM
!
C=======================================================================
      subroutine refsub(p_lat,rsv,rsh)
C=======================================================================
c forward calculation of reflection coefficient at the CMB
c for given lateral slowness component in (s/km)
c
c modified by James Wookey, 2004, hardwired AK135 velocities, instead of
c calling PREMPS
c
      implicit none
      real*4 p_lat
      real*8 rad
      integer i
      integer layp
      real*8 vp,vs,rho,xx

      integer idrt(10)
      real rt11, rt12, rt13 
      real rt21, rt22, rt23 
      complex coeff, refampsv, refampsh
      complex rsv, rsh
     
      logical logsh    

c premps: double precision
c coeff: single precision

C RADII only needed in version for layered 1D background model
c        data radii/0.0d0,1221.5d0,3480.0d0,
c     .             3630.0d0,5600.0d0,5701.0d0,
c     .             5771.0d0,5971.0d0,6151.0d0,
c     .             6291.0d0,6346.6d0,6356.0d0,
c     .             6371.0d0/

c p_lat denotes the horizontal slowness in units of s/km
c IASPEI travel time tables contain p_lat in s/deg
           
           do i=1,10
           idrt(i) = 0
           enddo
           idrt(4) = 1

c AK135 velocities just above CMB
           vp = 13.6602d0
           vs = 7.2811d0
           rho = 5.5515d0

c       ** altered velocities 
c          vp = 13.3d0
c          vs = 7.3d0
c          rho = 5.5d0
           
           rt11 = 1.0/sngl(vp)
           rt12 = 1.0/sngl(vs)
           rt13 = rho

c AK135 velocities just below CMB
           vp = 8.d0
           vs = 0.0000d0
           rho = 9.9145d0

           rt21 = 1.0/sngl(vp)
           rt22 = 0.0
           rt23 = rho


c SH reflection coefficient
           logsh=.true.             !SH refl/trans used in coeff
          refampsh = coeff(logsh,p_lat,rt11,rt12,rt13,
     .                       rt21,rt22,rt23,idrt)

c SV reflection coefficient
           logsh=.false.            !SH refl/trans used in coeff
          refampsv = coeff(logsh,p_lat,rt11,rt12,rt13,
     .                       rt21,rt22,rt23,idrt)
           rsh = refampsh
           rsv = refampsv

      return
      end
C=======================================================================

C=======================================================================
      subroutine xinv2(a,ai)
C=======================================================================
      implicit none
c reflection coefficients at D"
c in radial and transverse coordinates
c inversion of a 2x2 (complex) matrix
      complex a(2,2),ai(2,2), det
      det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
      if (cabs(det).lt.1.e-15) stop'det too small in xinv2'
      ai(1,1)=a(2,2)/det
      ai(2,2)=a(1,1)/det
      ai(1,2)=-a(1,2)/det
      ai(2,1)=-a(2,1)/det
      return
      end
C=======================================================================


C=======================================================================
      subroutine dppref(phii,rsv,rsh,r11,r12,r21,r22)
C=======================================================================
      implicit none
c reflection coefficients at D"
c in radial and transverse coordinates
      complex rsv,rsh,r11,r12,r21,r22
      real pi,phii,phi,cos2,sin2,cossin
      pi=4.*atan(1.)

      phi=phii*pi/180.
      cos2=cos(phi)*cos(phi)
      sin2=sin(phi)*sin(phi)
      cossin=cos(phi)*sin(phi)
      r11=rsv*cos2+rsh*sin2
      r12=(rsh-rsv)*cossin
      r21=(rsh-rsv)*cossin
      r22=rsv*sin2+rsh*cos2
      return
      end
C=======================================================================

C=====================================================================72
      COMPLEX FUNCTION COEFF (logsh,P,A1,B1,D1,A2,B2,D2,ID)             
C=====================================================================72
C RETURNS PRODUCT OF COMPLEX COEFFICIENTS AT INTERFACE .
C THIS VERSION HAS BEEN WRITTEN TO REDUCE THE AMOUNT OF UNNECESSARY
C COMPLEX ARITHMETIC .
C ID(1) IS MULTIPLICITY OF P1P1, ID(2) IS P1P2 PLUS P2P1, ID(3) IS P2P2
C ID(4) IS S1S1, ID(5) IS S1S2 PLUS S2S1, ID(6) IS S2S2,
C ID(7) IS P1S1 PLUS S1P1, ID(8) IS P1S2 PLUS S2P1,
C ID(9) IS P2S1 PLUS S1P2, AND ID(10) IS P2S2 PLUS S2P2.
C NOTE COEFFICIENTS ARE TAKEN ON UPPER RIEMANN SHEET WITH IMAG(Q).GE.0.
C COMPLEX FUNCTION TO RETURN PRODUCT OF ENERGY FLUX COEFFICIENTS
C NOTE THE S COEFFICIENTS ARE FOR SV WAVES . FOR SH WAVES A DIFFERENT
C ( MUCH SIMPLER ) SUBROUTINE MUST BE SUBSTITUTED .

      logical logsh
c
      REAL DUM(22)
      INTEGER ID(10)
      COMPLEX DELTA,AELP,AELM,BELP,BELM,ELEB1P,ELEB1M,ELEB2P,ELEB2M
      EQUIVALENCE (DUM(1),EL1R),(DUM(2),EL1I),(DUM(3),EL2R),
     .(DUM(4),EL2I),(DUM(5),EB1R),(DUM(6),EB1I),(DUM(7),EB2R),
     .(DUM(8),EB2I),(DUM(9),AEL1R),(DUM(10),AEL1I),(DUM(11),AEL2R),
     .(DUM(12),AEL2I),(DUM(13),BEL1R),(DUM(14),BEL1I),
     .(DUM(15),BEL2R),(DUM(16),BEL2I),(DUM(17),ELEB1R),
     .(DUM(18),ELEB1I),(DUM(19),ELEB2R),(DUM(20),ELEB2I),
     .(DUM(21),AMU1),(DUM(22),AMU2)
      COEFF=CMPLX(1.,0.)
      DO 101 I=1,10
      IF (ID(I).NE.0) GO TO 102
101   CONTINUE
      RETURN
C NOTE NO CHECK IS MADE THAT COEFFICIENTS AND MODEL ARE COMPATIBLE
C P1P1 AND S1S1 MAY BE REPLACED BY WKBJ COEFFICIENT (0.,-1.) BUT NOT
C P2P2 OR S2S2
C INCIDENT AND GENERATED WAVES ARE ASSUMED TO BE REAL
102   PP=P*P
      PP4=4.*PP
      DO 103 I=1,22
103   DUM(I)=0.
C SET UP REAL OR IMAGINARY WAVESLOWNESSES
      X2=(A2-P)*(A2+P)
      IF (X2) 104,105,106
104   EL2I=SQRT(-X2)
      AEL2I=D1*EL2I
      GO TO 105
106   EL2R=SQRT(X2)
      AEL2R=D1*EL2R
105   IF (D1.EQ.0.) GO TO 1
      X1=(A1-P)*(A1+P)
      IF (X1) 107,108,109
107   EL1I=SQRT(-X1)
      AEL1I=D2*EL1I
      GO TO 108
109   EL1R=SQRT(X1)
      AEL1R=D2*EL1R
108   AELP=CMPLX(AEL1R+AEL2R,AEL1I+AEL2I)
      AELM=CMPLX(AEL1R-AEL2R,AEL1I-AEL2I)
      IF (B1.EQ.0..AND.B2.EQ.0.) GO TO 6
      IF (B1.EQ.0.) GO TO 111
      AMU1=D1/B1**2
      X=(B1-P)*(B1+P)
      IF (X) 110,111,112
110   EB1I=SQRT(-X)
      BEL1I=D2*EB1I
      IF (X1) 113,111,114
113   ELEB1R=-EL1I*EB1I
      GO TO 111
114   ELEB1I=EL1R*EB1I
      GO TO 111
112   EB1R=SQRT(X)
      BEL1R=D2*EB1R
      IF (X1) 115,111,116
115   ELEB1I=EB1R*EL1I
      GO TO 111
116   ELEB1R=EB1R*EL1R
111   IF (B2.NE.0.) AMU2=D2/B2**2
      DMU=AMU1-AMU2
      X=PP4*DMU*ELEB1I
      ELEB1P=CMPLX(PP4*(DMU*(ELEB1R+PP)-D1),X)
      ELEB1M=CMPLX(PP4*(DMU*(ELEB1R-PP)+D1),X)
      IF (B2.EQ.0.) GO TO 5
      X=(B2-P)*(B2+P)
      IF (X) 117,118,119
117   EB2I=SQRT(-X)
      BEL2I=D1*EB2I
      IF (X2) 120,118,121
120   ELEB2R=-EL2I*EB2I
      GO TO 118
121   ELEB2I=EL2R*EB2I
      GO TO 118
119   EB2R=SQRT(X)
      BEL2R=D1*EB2R
      IF (X2) 122,118,123
122   ELEB2I=EB2R*EL2I
      GO TO 118
123   ELEB2R=EB2R*EL2R
118   X=DMU*ELEB2I
      ELEB2P=CMPLX(DMU*(ELEB2R+PP)+D2,X)
      ELEB2M=CMPLX(DMU*(ELEB2R-PP)-D2,X)
      IF (B1.EQ.0.) GO TO 4
      DD=PP*(D1+D2)**2
      DMU2=2.*PP*DMU
      IF(LOGSH) GO TO 1000                            !SH MOD}
      DM1=D1-DMU2
      DM2=D2+DMU2
      BELP=CMPLX(BEL1R+BEL2R,BEL1I+BEL2I)
      BELM=CMPLX(BEL1R-BEL2R,BEL1I-BEL2I)
      DELTA=AELP*BELP+ELEB1P*ELEB2P+CMPLX(DD,0.)
C  P1P1 REFLECTION COEFF.
      K=ID(1)
      IF (K.EQ.0) GO TO 18
      IF (P.GE.A1) GO TO 30
      COEFF=COEFF*((AELM*BELP+ELEB1M*ELEB2P-CMPLX(DD,0.))/DELTA)**K
C  P1P2+P2P1 REFLECTION COEFF.
18    K=ID(2)
      IF (K.EQ.0) GO TO 19
      X=2.*SQRT(AEL1R*AEL2R)
      COEFF=COEFF*(CMPLX(X*(DM2*EB1R+DM1*EB2R),X*(DM2*EB1I+DM1*EB2I))/
     .DELTA)**K
C  P2P2 REFLECTION COEFF.
19    K=ID(3)
      IF (K.EQ.0) GO TO 20
      COEFF=COEFF*((ELEB1P*ELEB2M-AELM*BELP-CMPLX(DD,0.))/DELTA)**K
C  S1S1 REFLECTION COEFF
20    K=ID(4)
      IF (K.EQ.0) GO TO 21
      IF (P.GE.B1) GO TO 31
      COEFF=COEFF*((CMPLX(DD,0.)-AELP*BELM-ELEB1M*ELEB2P)/DELTA)**K
C  S1S2+S2S1 REFLECTION COEFF.
21    K=ID(5)
      IF (K.EQ.0) GO TO 22
      X=2.*SQRT(BEL1R*BEL2R)
      COEFF=COEFF*(CMPLX(X*(DM2*EL1R+DM1*EL2R),X*(DM2*EL1I+DM1*EL2I))/
     .DELTA)**K
C  S2S2 REFLECTION COEFF.
22    K=ID(6)
      IF (K.EQ.0) GO TO 23
      COEFF=COEFF*((AELP*BELM-ELEB1P*ELEB2M+CMPLX(DD,0.))/DELTA)**K
C  P1S1+S1P1 REFLECTION COEFF.
23    K=ID(7)
      IF (K.EQ.0) GO TO 24
      Y=2.*DMU*DM1
      X=SQRT(PP4*ELEB1R)
      COEFF=COEFF*(CMPLX(X*(Y*ELEB2R+(DM1-D2)*DM2),X*Y*ELEB2I)/
     .DELTA)**K
C P1S2+S2P1 REFLECTION COEFF.
24    K=ID(8)
      IF (K.EQ.0) GO TO 25
      Y=DMU+DMU
      X=SQRT(PP4*AEL1R*BEL2R)
      COEFF=COEFF*(CMPLX(X*(Y*(EB1R*EL2R-EB1I*EL2I)+DM2-D1),
     .X*Y*(EB1I*EL2R+EB1R*EL2I))/DELTA)**K
C  P2S1+S1P2 REFLECTION COEFF.
25    K=ID(9)
      IF (K.EQ.0) GO TO 26
      Y=DMU+DMU
      X=-SQRT(PP4*AEL2R*BEL1R)
      COEFF=COEFF*(CMPLX(X*(Y*(EL1R*EB2R-EL1I*EB2I)+DM2-D1),
     .X*Y*(EL1R*EB2I+EL1I*EB2R))/DELTA)**K
C  P2S2+S2P2 REFLECTION COEFF.
26    K=ID(10)
      IF (K.EQ.0) GO TO 100
      Y=2.*DM2*DMU
      X=-SQRT(PP4*ELEB2R)
      COEFF=COEFF*(CMPLX(X*(Y*ELEB1R+(DM1-D2)*DM1),X*Y*ELEB1I)/DELTA)**K
      GO TO 100
30    COEFF=COEFF*CMPLX(0.,-1.)**K
      GO TO 18
C SV TURNING RAY CHANGES SIGN OF HORIZONTAL COMPONENT SO COEFFICIENT
C IS -1.*(0.,-1.) ( P TURNING RAY IS THE WKBJ COEFFICIENT (0.,-1.) )
31    COEFF=COEFF*CMPLX(0.,1.)**K
      GO TO 21
C
C FLUID ABOVE
C
4     X=PP4/B2**2
      ELEB1P=-CMPLX(X*EL1R,X*EL1I)
      DELTA=AELP+ELEB1P*ELEB2P
C  P1P1 REFLECTION COEFF.
      K=ID(1)
      IF (K.EQ.0) GO TO 7
      IF (P.GE.A1) GO TO 34
      COEFF=COEFF*((DELTA-CMPLX(AEL2R+AEL2R,AEL2I+AEL2I))/DELTA)**K
C  P1P2+P2P1 REFLECTION COEFF.
7     K=ID(2)
      IF (K.EQ.0) GO TO 8
      COEFF=COEFF*(CMPLX((2.-X)*SQRT(AEL1R*AEL2R),0.)/DELTA)**K
C  P2P2 REFLECTION COEFF.
8     K=ID(3)
      IF (K.EQ.0) GO TO 9
      COEFF=COEFF*((ELEB1P*ELEB2M-AELM)/DELTA)**K
C  S2S2 REFLECTION COEFF.
9     K=ID(6)
      IF (K.EQ.0) GO TO 10
      IF(.NOT.LOGSH) COEFF=COEFF*((AELP-ELEB1P*ELEB2M)/DELTA)**K
C     PREVIOUS STATEMENT CHANGED IN SH MOD+++++++++++++++++++++++
      IF(LOGSH) COEFF=COEFF*CMPLX(1.,0.)              !SH MOD}
C P1S2+S2P1 REFLECTION COEFF.
10    K=ID(8)
      IF (K.EQ.0) GO TO 11
      Y=-X*SQRT(D1*AEL1R*EB2R)/P
      COEFF=COEFF*(CMPLX(Y*EL2R,Y*EL2I)/DELTA)**K
C  P2S2+S2P2 REFLECTION COEFF.
11    K=ID(10)
      IF (K.EQ.0) GO TO 100
      Y=X*(D2-2.*PP*AMU2)*SQRT(ELEB2R)/P
      COEFF=COEFF*(CMPLX(Y*EL1R,Y*EL1I)/DELTA)**K
      GO TO 100
34    COEFF=COEFF*CMPLX(0.,-1.)**K
      GO TO 7
C
C FLUID BELOW
C
5     X=AMU1/D1
      ELEB2P=CMPLX(X*EL2R,X*EL2I)
      X=PP4*X
      DELTA=AELP+ELEB1P*ELEB2P
C P1P1 REFLECTION COEFF.
      K=ID(1)
      IF (K.EQ.0) GO TO 12
      IF (P.GE.A1) GO TO 36
      COEFF=COEFF*((AELM+ELEB2P*ELEB1M)/DELTA)**K
C  P1P2+P2P1 REFLECTION COEFF.
12    K=ID(2)
      IF (K.EQ.0) GO TO 13
      COEFF=COEFF*(CMPLX((2.-X)*SQRT(AEL1R*AEL2R),0.)/DELTA)**K
C P2P2 REFLECTION COEFF.
13    K=ID(3)
      IF (K.EQ.0) GO TO 14
      COEFF=COEFF*((DELTA-CMPLX(AEL1R+AEL1R,AEL1I+AEL1I))/DELTA)**K
C S1S1 REFLECTION COEFF.
14    K=ID(4)
      IF (K.EQ.0) GO TO 15
      IF (P.GE.B1) GO TO 37
      IF(.NOT.LOGSH) COEFF=COEFF*((AELP-ELEB2P*ELEB1M)/DELTA)**K
C     PREVIOUS STATEMENT CHANGED IN SH MOD+++++++++++++++++++++++++
      IF(LOGSH) COEFF=CMPLX(1.,0.)*COEFF              !SH MOD}
C  P1S1+S1P1 REFLECTION COEFFICIENT
15    K=ID(7)
      IF (K.EQ.0) GO TO 16
      Y=X*(D1-2.*PP*AMU1)*SQRT(ELEB1R)/P
      COEFF=COEFF*(CMPLX(Y*EL2R,Y*EL2I)/DELTA)**K
C  P2S1+S1P2 REFLECTION COEFF.
16    K=ID(9)
      IF (K.EQ.0) GO TO 100
      Y=-X*SQRT(D2*AEL2R*EB1R)/P
      COEFF=COEFF*(CMPLX(Y*EL1R,Y*EL1I)/DELTA)**K
      GO TO 100
36    COEFF=COEFF*CMPLX(0.,-1.)**K
      GO TO 12
C SV TURNING RAY CHANGES SIGN OF HORIZONTAL COMPONENT SO COEFFICIENT
C IS -1.*(0.,-1.) ( P TURNING RAY IS THE WKBJ COEFFICIENT (0.,-1.) )
37    IF(.NOT.LOGSH) COEFF=COEFF*CMPLX(0.,1.)**K
      IF(LOGSH) COEFF=COEFF*CMPLX(0.,-1.)**K
      IF(LOGSH) COEFF=-COEFF                               !SHMOD}
C  SH TURNING RAY DOES NOT SHIFT SIGN                      !SH MOD}
      GO TO 15
C
C FLUID-FLUID
C
6     K=ID(1)
      IF (K.EQ.0) GO TO 17
      IF (P.GE.A1) GO TO 40
      COEFF=COEFF*(AELM/AELP)**K
17    K=ID(3)
      IF (K.EQ.0) GO TO 42
      COEFF=COEFF*(-AELM/AELP)**K
42    K=ID(2)
      IF (K.EQ.0) GO TO 100
      COEFF=COEFF*(CMPLX(2.*SQRT(AEL1R*AEL2R),0.)/AELP)**K
      GO TO 100
40    COEFF=COEFF*CMPLX(0.,-1.)**K
      GO TO 17
C
C FREE SURFACE
C
1     IF (B2.EQ.0.) GO TO 2
      IF(LOGSH) COEFF=COEFF*CMPLX(1.,0.)              !SH MOD}
      IF(LOGSH)  RETURN                               !SH MOD}
      X=(B2-P)*(B2+P)
      Y=X-PP
      X2=Y**2
      X=PP4*SQRT(X)
      ELEB2R=X*EL2R
      ELEB2I=X*EL2I
      DELTA=CMPLX(X2+ELEB2R,ELEB2I)
      K=ID(10)
      IF (K.EQ.0) GO TO 3
      COEFF=COEFF*(CMPLX(2.*Y*SQRT(ELEB2R),0.)/DELTA)**K
3     K=ID(3)+ID(6)
      IF (K.EQ.0) GO TO 100
      COEFF=COEFF*(CMPLX(X2-ELEB2R,-ELEB2I)/DELTA)**K
2     IF (MOD(ID(3),2).EQ.1) COEFF=-COEFF
100   RETURN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  THE LAYER INTERFACE MULTIPLIES BELOW WERE ADDED IN THE SH
C  MODIFICATION.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
1000  BELP=AMU1*CMPLX(EB1R,EB1I)
      BELM=AMU2*CMPLX(EB2R,EB2I)
      DELTA=BELP+BELM
      K=ID(4)
      IF(K.EQ.0) GO TO 1001
      IF(P.GE.B1) GO TO 1030
      COEFF=COEFF*((BELP-BELM)/DELTA)**K
1001  K=ID(5)
      IF(K.EQ.0) GO TO 1002
      COEFF=COEFF*(2.*CSQRT(BELP)*CSQRT(BELM)/DELTA)**K
1002  K=ID(6)
      IF(K.EQ.0) GO TO 1003
      COEFF=COEFF*((BELM-BELP)/DELTA)**K
1003  RETURN
1030  COEFF=COEFF*CMPLX(0.,1.)**K
      RETURN
C  END OF THE SH MODIFICATIONS TO LAYER MULTIPLIES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      END                                                               

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine xinver(icmb,nt,dt,phip,urad,utra,plat,
     & alag,angle,uradi,utrai)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c inverse operator

      complex urad(*), utra(*)
      complex uradi(*), utrai(*)

      complex  a, b, as, bs
      complex ccm(2,2), cicm(2,2)
      complex rsv,rsh,r11,r12,r21,r22
      complex rr(2,2), ci1(2,2), ci2(2,2),rri(2,2)

      pi=4.*atan(1.)

      delta=alag
      phi=angle

      if (icmb.eq.1) then 
c adjusted for CMB reflection
      call refsub(plat,rsv,rsh)   
      call dppref(phip,rsv,rsh,r11,r12,r21,r22)
c define reflection matrix
      rr(1,1)=r11
      rr(1,2)=r12
      rr(2,1)=r21
      rr(2,2)=r22
      else
      rr(1,1)=cmplx(1.,0.)
      rr(1,2)=cmplx(0.,0.)
      rr(2,1)=cmplx(0.,0.)
      rr(2,2)=cmplx(1.,0.)
      end if
      call xinv2(rr,rri)
         print*,rri(1,1),rri(1,2)
         print*,rri(2,1),rri(2,2)
      print*,phip

      do iomega=1,nt/2+1
      omega=2.*pi/(float(nt)*dt) * float(iomega-1)

c set up inverse operator 
        theta=omega*delta/2.
        alpha=2.*(phi-phip)*pi/180.
        ar=cos(theta)
        ai=-sin(theta)*cos(alpha)
        br=0.
        bi=-sin(theta)*sin(alpha)
        a=cmplx(ar,ai)
        b=cmplx(br,bi)
        as=cmplx(ar,-ai)
        bs=cmplx(-br,bi)

        ccm(1,1)=a
        ccm(1,2)=b
        ccm(2,1)=bs
        ccm(2,2)=as

       call xmult2(rr,ccm,ci1)
       call xmult2(ccm,ci1,ci2)
       call xinv2(ci2,cicm)   

       uradi(iomega)=cicm(1,1)*urad(iomega)+cicm(1,2)*utra(iomega)
       utrai(iomega)=cicm(2,1)*urad(iomega)+cicm(2,2)*utra(iomega)

      end do

      return
      end

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine xinver2(icmb,nt,dt,phip,urad,utra,plat,
     & alag,angle,uradi,utrai)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c inverse operator

      complex urad(*), utra(*)
      complex uradi(*), utrai(*)

      complex rsv,rsh,r11,r12,r21,r22
      complex rr(2,2), rri(2,2)

      pi=4.*atan(1.)

      delta=alag
      phi=angle

      if (icmb.eq.1) then 
c adjusted for CMB reflection
      call refsub(plat,rsv,rsh)   
      call dppref(phip,rsv,rsh,r11,r12,r21,r22)
c define reflection matrix
      rr(1,1)=r11
      rr(1,2)=r12
      rr(2,1)=r21
      rr(2,2)=r22
      else
      rr(1,1)=cmplx(1.,0.)
      rr(1,2)=cmplx(0.,0.)
      rr(2,1)=cmplx(0.,0.)
      rr(2,2)=cmplx(1.,0.)
      end if
      call xinv2(rr,rri)
         print*,rri(1,1),rri(1,2)
         print*,rri(2,1),rri(2,2)
      print*,phip

      do iomega=1,nt/2+1
      omega=2.*pi/(float(nt)*dt) * float(iomega-1)


       uradi(iomega)=rri(1,1)*urad(iomega)+rri(1,2)*utra(iomega)
       utrai(iomega)=rri(2,1)*urad(iomega)+rri(2,2)*utra(iomega)

      end do

      return
      end



      subroutine xmult2(a,b,c)
      complex a(2,2),b(2,2),c(2,2)
c multiplication of 2X2 (complex) matrix
      c(1,1)=a(1,1)*b(1,1)+a(1,2)*b(2,1)
      c(1,2)=a(1,1)*b(1,2)+a(1,2)*b(2,2)
      c(2,1)=a(2,1)*b(1,1)+a(2,2)*b(2,1)
      c(2,2)=a(2,1)*b(1,2)+a(2,2)*b(2,2)
      return
      end
