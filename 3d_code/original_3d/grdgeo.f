*
*      file grdgeo.f
************************************************************************
      SUBROUTINE GRDGEO(XP,YP,ZP,XE,YN,ZT,DIEP,DJNP,DKTP,
     C                  DISE,DISN,DIST,AREP,ARNP,ARTP,VOLP,
     C                  IBW,IEW,JBW,JEW,KBW,KEW,
     C                  XNET,YNET,ZNET,IB,IE,JB,JE,KB,KE,ID,JD,KD) 
*
*     Subroutine to calculate the geometrical properties for a simple
*     orthogonal grid. All properties are based on control volumes 
*     obtained by joining the corner points with straight lines.
*
*     XP(ID) x location of P node for each c.v.; output
*     YP(JD) y location of P node for each c.v.; output
*     ZP(KD) z location of P node for each c.v.; output
*     XE(ID) x location of e integration point;  output
*     YN(JD) y location of n integration point;  output
*     ZT(KD) z location of t integration point;  output
*
*     DIEP(ID) distance from P to E through e; output
*     DJNP(JD) distance from P to N through n; output
*     DKTP(KD) distance from P to T through t; output
*     
*     DISE(ID) distance in i from P to e face; output
*     DISN(JD) distance in j from P to n face; output
*     DIST(KD) distance in k from P to t face; output
*
*     AREP(JD,KD) c.v. area of face at e point; output
*     ARNP(ID,KD) c.v. area of face at n point; output
*     ARTP(ID,JD) c.v. area of face at t point; output
*     VOLP(ID,JD,KD) volume of c.v. around P node; output
*
*     XNET(ID) x location of net corner points; input
*     YNET(JD) y location of net corner points; input
*     ZNET(KD) z location of net corner points; input
*
*     INTEGER IB,IE first and last interior indices in i; input
*     INTEGER JB,JE first and last interior indices in j; input
*     INTEGER KB,KE first and last interior indices in k; input
*     INTEGER ID,JD,KD array dimensions; input
*
*     Notes: 1) This is setup for an orthogonal grid.
*
*            2) The locations of the e, n, and t integration points 
*               should be calculated as the face midpoints and stored
*               in arrays.
*
************************************************************************
*
      IMPLICIT NONE
      REAL*8 XP(ID),YP(JD),ZP(KD)
      REAL*8 XE(ID),YN(JD),ZT(KD)
      REAL*8 DIEP(ID),DJNP(JD),DKTP(KD)
      REAL*8 DISE(ID),DISN(JD),DIST(KD)
      REAL*8 AREP(JD,KD),ARNP(ID,KD),ARTP(ID,JD),VOLP(ID,JD,KD)
      REAL*8 XNET(ID),YNET(JD),ZNET(KD)
      INTEGER IBW,IEW,JBW,JEW,KBW,KEW
      INTEGER IB,IE,JB,JE,KB,KE,IBM1,IEP1,JBM1,JEP1,KBM1,KEP1
      INTEGER I,J,K,ID,JD,KD
*
      IBM1=IB-1
      IEP1=IE+1
      JBM1=JB-1
      JEP1=JE+1
      KBM1=KB-1
      KEP1=KE+1
*
*     Calculate locations of P nodes for interior and boundary c.v.
*
      XP(IBM1)=XNET(IBM1)
      DO 10 I=IB,IE
        XP(I)=0.5*(XNET(I)+XNET(I-1))
   10 CONTINUE
      XP(IEP1)=XNET(IE)
*
      YP(JBM1)=YNET(JBM1)
      DO 20 J=JB,JE
        YP(J)=0.5*(YNET(J)+YNET(J-1))
   20 CONTINUE
      YP(JEP1)=YNET(JE)
*
      ZP(KBM1)=ZNET(KBM1)
      DO 30 K=KB,KE
        ZP(K)=0.5*(ZNET(K)+ZNET(K-1))
   30 CONTINUE
      ZP(KEP1)=ZNET(KE)
*
*  Compute grid geometry parameters
*
      DO 50 I=IBM1,IE
        DIEP(I)=XP(I+1)-XP(I)
 50   CONTINUE
*
      DO 60 J=JBM1,JE
        DJNP(J)=YP(J+1)-YP(J)
 60   CONTINUE
*
      DO 65 K=KBM1,KE
        DKTP(K)=ZP(K+1)-ZP(K)
 65   CONTINUE
*
*  Calculate midface locations
*     
      DO 70 I=IBM1,IE
        XE(I)=XNET(I)
 70   CONTINUE
      XE(IEP1)=XNET(IE)
*
      DO 80 J=JBM1,JE
        YN(J)=YNET(J)
 80   CONTINUE
      YN(JEP1)=YNET(JE)
*
      DO 90 K=KBM1,KE
        ZT(K)=ZNET(K)
 90   CONTINUE
      ZT(KEP1)=ZNET(KE)
*
      DO 100 I=IBM1,IEP1
        DISE(I)=XE(I)-XP(I)
 100  CONTINUE
*
      DO 120 J=JBM1,JEP1
        DISN(J)=YN(J)-YP(J)
 120  CONTINUE
*
      DO 130 K=KBM1,KEP1
        DIST(K)=ZT(K)-ZP(K)
 130  CONTINUE
*
*  Calculate planar areas and volumes for Cartesian grid
*
      DO 201 I=IB,IE
       DO 200 K=KB,KE
        ARNP(I,K) = 4.0*DISE(I)*DIST(K)
 200   CONTINUE
 201  CONTINUE
      DO 211 J=JB,JE
       DO 210 K=KB,KE
        AREP(J,K) = 4.0*DISN(J)*DIST(K)
 210   CONTINUE
 211  CONTINUE
      DO 216 I=IB,IE
       DO 215 J=JB,JE
        ARTP(I,J) = 4.0*DISE(I)*DISN(J)
 215   CONTINUE
 216  CONTINUE
*     
      DO 240 K=KB,KE
       DO 230 J=JB,JE
        DO 220 I=IB,IE
         VOLP(I,J,K) = 2.0*ARTP(I,J)*DIST(K)
 220    CONTINUE
 230   CONTINUE
 240  CONTINUE
*     
      RETURN
      END
