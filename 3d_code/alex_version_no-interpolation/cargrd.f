*           file cargrd.f
***********************************************************************
*
      SUBROUTINE CARGRD(XNET,YNET,ZNET, GRDX1,GRDX2,GRDX3,GRDX,
     C            GRDY1,GRDY2,GRDY3,GRDY,GRDZ1,GRDZ2,GRDZ3,GRDZ,
     C            IRATX1,IRATX2,IRATX3,IRATX4,
     C            IRATY1,IRATY2,IRATY3,IRATY4,
     C            IRATZ1,IRATZ2,IRATZ3,IRATZ4,
     C            IBW,IEW,JBW,JEW,KBW,KEW,X_IB,X_IE,Y_JB,Y_JE,
     C            Z_KB,Z_KE,XP,YP,ZP,XE,YN,ZT,DIEP,DJNP,DKTP,
     C            DISE,DISN,DIST,AREP,ARNP,ARTP,VOLP,LVLGEO,
     C            POROUS,IPB,IPE,JPB,JPE,KPB,KPE,
     C            IPB2,IPE2,JPB2,JPE2,KPB2,KPE2,
     C            SOLID,ISB,ISE,JSB,JSE,KSB,KSE,
     C            ISB2,ISE2,JSB2,JSE2,KSB2,KSE2,
     C            CVTYPE,IB,IE1,IE2,IE3,IE,JB,JE1,JE2,JE3,JE,
     C            KB,KE1,KE2,KE3,KE,ID,JD,KD,NNB)
*
*  Routine to organize the creation of a cartesian grid.
*
***********************************************************************
*
      IMPLICIT NONE
      REAL*8 XNET(ID),YNET(JD),ZNET(KD)
      REAL*8 GRDX1,GRDX2,GRDX3,GRDX,GRDY1,GRDY2,GRDY3,GRDY
      REAL*8 GRDZ1,GRDZ2,GRDZ3,GRDZ
      REAL*8 IRATX1,IRATX2,IRATX3,IRATX4
      REAL*8 IRATY1,IRATY2,IRATY3,IRATY4
      REAL*8 IRATZ1,IRATZ2,IRATZ3,IRATZ4            
      REAL*8 XP(ID),YP(JD),ZP(KD),XE(ID),YN(JD),ZT(KD)
      REAL*8 DIEP(ID),DJNP(JD),DKTP(KD),DISE(ID),DISN(JD),DIST(KD)
      REAL*8 AREP(JD,KD),ARNP(ID,KD),ARTP(ID,JD),VOLP(ID,JD,KD)
      REAL*8 X_IB,X_IE,Y_JB,Y_JE,Z_KB,Z_KE
*
      INTEGER IBW,IEW,JBW,JEW,KBW,KEW,LVLGEO
      INTEGER POROUS,IPB,IPE,JPB,JPE,KPB,KPE
      INTEGER IPB2,IPE2,JPB2,JPE2,KPB2,KPE2
      INTEGER SOLID,ISB,ISE,JSB,JSE,KSB,KSE
      INTEGER ISB2,ISE2,JSB2,JSE2,KSB2,KSE2
      INTEGER IB,IE1,IE2,IE3,IE,JB,JE1,JE2,JE3,JE,KB,KE1,KE2,KE3,KE
      INTEGER CVTYPE(ID,JD,KD,NNB+1),ID,JD,KD,NNB
*
      IBW = 0
      IEW = 0
      JBW = 0
      JEW = 0
      KBW = 0
      KEW = 0
*
*  Define integration points and c.v. corners; generate grid
*
      CALL MAKGRD(XNET,YNET,ZNET, GRDX1,GRDX2,GRDX3,GRDX,
     C            GRDY1,GRDY2,GRDY3,GRDY,GRDZ1,GRDZ2,GRDZ3,GRDZ,
     C            IRATX1,IRATX2,IRATX3,IRATX4,
     C            IRATY1,IRATY2,IRATY3,IRATY4,
     C            IRATZ1,IRATZ2,IRATZ3,IRATZ4,
     C            IBW,IEW,JBW,JEW,KBW,KEW,X_IB,X_IE,Y_JB,Y_JE,
     C            Z_KB,Z_KE,IB,IE1,IE2,IE3,IE,JB,JE1,JE2,JE3,JE,
     C            KB,KE1,KE2,KE3,KE,ID,JD,KD)
      CALL GRDGEO(XP,YP,ZP,XE,YN,ZT,DIEP,DJNP,DKTP,DISE,DISN,DIST,
     C            AREP,ARNP,ARTP,VOLP,
     C            IBW,IEW,JBW,JEW,KBW,KEW,
     C            XNET,YNET,ZNET,IB,IE,JB,JE,KB,KE,ID,JD,KD) 
      CALL TYPE(CVTYPE, POROUS,IPB,IPE,JPB,JPE,KPB,KPE,
     C          IPB2,IPE2,JPB2,JPE2,KPB2,KPE2,
     C          SOLID,ISB,ISE,JSB,JSE,KSB,KSE,
     C          ISB2,ISE2,JSB2,JSE2,KSB2,KSE2,
     C          IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)     
*
      RETURN
      END
