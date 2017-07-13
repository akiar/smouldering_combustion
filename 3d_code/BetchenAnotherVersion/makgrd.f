*
*           file makgrd.f
********************************************************************
      SUBROUTINE MAKGRD(XNET,YNET,ZNET, GRDX1,GRDX2,GRDX3,GRDX,
     C            GRDY1,GRDY2,GRDY3,GRDY,GRDZ1,GRDZ2,GRDZ3,GRDZ,
     C            IRATX1,IRATX2,IRATX3,IRATX4,
     C            IRATY1,IRATY2,IRATY3,IRATY4,
     C            IRATZ1,IRATZ2,IRATZ3,IRATZ4,
     C            IBW,IEW,JBW,JEW,KBW,KEW,X_IB,X_IE,Y_JB,Y_JE,
     C            Z_KB,Z_KE,IB,IE1,IE2,IE3,IE,JB,JE1,JE2,JE3,JE,
     C            KB,KE1,KE2,KE3,KE,ID,JD,KD)
*
*     Subroutine to calculate the corner locations of asimple 
*     orthogonal grid. The grid is formed from a pin-jointed
*     structure with non-uniform properties.
*
*     XNET(ID)   x location of northeast top corners; output
*     YNET(JD)   y location of northeast top corners; output
*     ZNET(KD)   z location of northeast top corners; output
*
*     GRDX,GRDY,GRDZ       grid lengths; input
*     GRDAB    position of interface; input
*     IRATAB   ratio of 1st to last c.v. in part B of direction A; input
*     INTEGER IB    index of first interior volume in i direction; input
*     INTEGER IE    index of last interior volume in i direction; input
*     INTEGER JB    index of first interior volume in j direction; input
*     INTEGER JE    index of last interior volume in j direction; input
*     INTEGER KB    index of first interior volume in k direction; input
*     INTEGER KE    index of last interior volume in k direction; input
*     INTEGER ID,JD,KD dimensions of arrays
*
*     Notes: 1) Corner locations must be calculated for volumes with 
*               indices running from IB-1 to IE, JB-1 to JE, and
*               KB-1 to KE
*
********************************************************************
*
      IMPLICIT NONE
      REAL*8 XNET(ID),YNET(JD),ZNET(KD)
      REAL*8 GRDX1,GRDX2,GRDX3,GRDX,GRDY1,GRDY2,GRDY3,GRDY
      REAL*8 GRDZ1,GRDZ2,GRDZ3,GRDZ
      REAL*8 IRATX1,IRATX2,IRATX3,IRATX4
      REAL*8 IRATY1,IRATY2,IRATY3,IRATY4
      REAL*8 IRATZ1,IRATZ2,IRATZ3,IRATZ4
      REAL*8 SUMDX,SUMDY,SUMDZ,DX1,DY1,DZ1,GRADEX,GRADEY,GRADEZ
      REAL*8 X_IB,X_IE,Y_JB,Y_JE,Z_KB,Z_KE
      INTEGER IBW,IEW,JBW,JEW,KBW,KEW
      INTEGER IB,IE1,IE2,IE3,IE,JB,JE1,JE2,JE3,JE
      INTEGER KB,KE1,KE2,KE3,KE,ID,JD,KD
      INTEGER NCVX,NCVY,NCVZ,I,J,K
*
*  Check for input errors
*
      IF(IE3.GT.IE.OR.IE2.GT.IE3.OR.IE1.GT.IE2.OR.JE3.GT.JE
     C   .OR.JE2.GT.JE3.OR.JE1.GT.JE2.OR.KE3.GT.KE.OR.KE2.GT.KE3
     C   .OR.KE1.GT.KE2) GOTO 111
*
      NCVX= IE-IB+1
      NCVY= JE-JB+1
      NCVZ= KE-KB+1
*
*  Set position of lower left corner
*
      XNET(IB-1)= 0.0
      YNET(JB-1)= 0.0
      ZNET(KB-1)= 0.0
*
*  Construct 'x' direction
*
*  (a) from IB to IE1
*
      IF(IBW.EQ.1) THEN
        XNET(IB)= X_IB
      ENDIF
      NCVX= IE1-IB+1-IBW
      IF(NCVX.GT.1) THEN
        SUMDX= 0.0
        DO 1 I=1,NCVX
          SUMDX= SUMDX+(1.0-(I-1)*(1.0-1.0/IRATX1)/(NCVX-1))
  1     CONTINUE
        DX1= (GRDX1-X_IB)/SUMDX
        GRADEX= DX1/(NCVX-1)*(1.0-1.0/IRATX1)
        XNET(IB+IBW)= XNET(IB+IBW-1)+DX1
        DO 2 I=2,NCVX
          XNET(I+IB+IBW-1)= XNET(I+IB+IBW-2)+(DX1-(I-1)*GRADEX)
  2     CONTINUE
      ELSE IF(NCVX.EQ.1) THEN
        XNET(IB)= GRDX
      ELSE
        GOTO 111
      ENDIF
*
*  (b) from IE1 to IE2
*
      NCVX= IE2-IE1
      IF(NCVX.LT.1) GOTO 3
        SUMDX= 0.0
        DO 4 I=1,NCVX
          SUMDX= SUMDX+(1.0-(I-1)*(1.0-1.0/IRATX2)/(NCVX-1))
  4     CONTINUE
        DX1= GRDX2/SUMDX
        GRADEX= DX1/(NCVX-1)*(1.0-1.0/IRATX2)
        XNET(IE1+1)= XNET(IE1)+DX1
        DO 5 I=2,NCVX
          XNET(IE1+I)= XNET(IE1+I-1)+(DX1-(I-1)*GRADEX)
  5     CONTINUE
*
*  (c) from IE2 to IE3
*
  3   NCVX= IE3-IE2
      IF(NCVX.LT.1) GOTO 6
        SUMDX= 0.0
        DO 7 I=1,NCVX
          SUMDX= SUMDX+(1.0-(I-1)*(1.0-1.0/IRATX3)/(NCVX-1))
  7     CONTINUE
        DX1= GRDX3/SUMDX
        GRADEX= DX1/(NCVX-1)*(1.0-1.0/IRATX3)
        XNET(IE2+1)= XNET(IE2)+DX1
        DO 8 I=2,NCVX
          XNET(IE2+I)= XNET(IE2+I-1)+(DX1-(I-1)*GRADEX)
  8     CONTINUE  
*
*  (d) from IE3 to IE
*
  6   NCVX= IE-IE3-IEW
      IF(NCVX.LT.1) GOTO 9
        SUMDX= 0.0
        DO 10 I=1,NCVX
          SUMDX= SUMDX+(1.0-(I-1)*(1.0-1.0/IRATX4)/(NCVX-1))
  10     CONTINUE
        DX1= (GRDX-GRDX3-GRDX1-GRDX2-X_IE)/SUMDX
        GRADEX= DX1/(NCVX-1)*(1.0-1.0/IRATX4)
        XNET(IE3+1)= XNET(IE3)+DX1
        DO 11 I=2,NCVX
          XNET(IE3+I)= XNET(IE3+I-1)+(DX1-(I-1)*GRADEX)
  11     CONTINUE
        IF(IEW.EQ.1) THEN
          XNET(IE)= GRDX
        ENDIF
  9     XNET(IE+1)= XNET(IE)
*
*  Construct 'y' direction
*
*  (a) from JB to JE1
*
      IF(JBW.EQ.1) THEN
        YNET(JB)= Y_JB
      ENDIF
      NCVY= JE1-JB+1-JBW
      IF(NCVY.GT.1) THEN
        SUMDY= 0.0
        DO 16 J=1,NCVY
          SUMDY= SUMDY+(1.0-(J-1)*(1.0-1.0/IRATY1)/(NCVY-1))
  16    CONTINUE
        DY1= (GRDY1-Y_JB)/SUMDY
        GRADEY= DY1/(NCVY-1)*(1.0-1.0/IRATY1)
        YNET(JB+JBW)= YNET(JB+JBW-1)+DY1
        DO 17 J=2,NCVY
          YNET(J+JB+JBW-1)= YNET(J+JB+JBW-2)+(DY1-(J-1)*GRADEY)
  17    CONTINUE
      ELSE IF(NCVY.EQ.1) THEN
        YNET(JB)= GRDY
      ELSE
        GOTO 111
      ENDIF
*
*  (b) from JE1 to JE2
*
      NCVY= JE2-JE1
      IF(NCVY.LT.1) GOTO 18
        SUMDY= 0.0
        DO 19 J=1,NCVY
          SUMDY= SUMDY+(1.0-(J-1)*(1.0-1.0/IRATY2)/(NCVY-1))
  19    CONTINUE
        DY1= GRDY2/SUMDY
        GRADEY= DY1/(NCVY-1)*(1.0-1.0/IRATY2)
        YNET(JE1+1)= YNET(JE1)+DY1
        DO 20 J=2,NCVY
          YNET(JE1+J)= YNET(JE1+J-1)+(DY1-(J-1)*GRADEY)
  20    CONTINUE
*
*  (c) from JE2 to JE3
*
  18  NCVY= JE3-JE2
      IF(NCVY.LT.1) GOTO 21
        SUMDY= 0.0
        DO 22 J=1,NCVY
          SUMDY= SUMDY+(1.0-(J-1)*(1.0-1.0/IRATY3)/(NCVY-1))
  22    CONTINUE
        DY1= GRDY3/SUMDY
        GRADEY= DY1/(NCVY-1)*(1.0-1.0/IRATY3)
        YNET(JE2+1)= YNET(JE2)+DY1
        DO 23 J=2,NCVY
          YNET(JE2+J)= YNET(JE2+J-1)+(DY1-(J-1)*GRADEY)
  23    CONTINUE  
*
*  (d) from JE3 to JE
*
  21  NCVY= JE-JE3-JEW
      IF(NCVY.LT.1) GOTO 24
        SUMDY= 0.0
        DO 25 J=1,NCVY
          SUMDY= SUMDY+(1.0-(J-1)*(1.0-1.0/IRATY4)/(NCVY-1))
  25    CONTINUE
        DY1= (GRDY-GRDY3-GRDY1-GRDY2-Y_JE)/SUMDY
        GRADEY= DY1/(NCVY-1)*(1.0-1.0/IRATY4)
        YNET(JE3+1)= YNET(JE3)+DY1
        DO 26 J=2,NCVY
          YNET(JE3+J)= YNET(JE3+J-1)+(DY1-(J-1)*GRADEY)
  26    CONTINUE
        IF(JEW.EQ.1) THEN
          YNET(JE)= GRDY
        ENDIF
  24    YNET(JE+1)= YNET(JE)
*
*  Construct 'z' direction
*
*  (a) from KB to KE1
*
      IF(KBW.EQ.1) THEN
        ZNET(KB)= Z_KB
      ENDIF
      NCVZ= KE1-KB+1-KBW
      IF(NCVZ.GT.1) THEN
        SUMDZ= 0.0
        DO 27 K=1,NCVZ
          SUMDZ= SUMDZ+(1.0-(K-1)*(1.0-1.0/IRATZ1)/(NCVZ-1))
  27    CONTINUE
        DZ1= (GRDZ1-Z_KB)/SUMDZ
        GRADEZ= DZ1/(NCVZ-1)*(1.0-1.0/IRATZ1)
        ZNET(KB+KBW)= ZNET(KB+KBW-1)+DZ1
        DO 28 K=2,NCVZ
          ZNET(K+KB+KBW-1)= ZNET(K+KB+KBW-2)+(DZ1-(K-1)*GRADEZ)
  28    CONTINUE
      ELSE IF(NCVZ.EQ.1) THEN
        ZNET(KB)= GRDZ
      ELSE
        GOTO 111
      ENDIF
*
*  (b) from KE1 to KE2
*
      NCVZ= KE2-KE1
      IF(NCVZ.LT.1) GOTO 29
        SUMDZ= 0.0
        DO 30 K=1,NCVZ
          SUMDZ= SUMDZ+(1.0-(K-1)*(1.0-1.0/IRATZ2)/(NCVZ-1))
  30    CONTINUE
        DZ1= GRDZ2/SUMDZ
        GRADEZ= DZ1/(NCVZ-1)*(1.0-1.0/IRATZ2)
        ZNET(KE1+1)= ZNET(KE1)+DZ1
        DO 31 K=2,NCVZ
          ZNET(KE1+K)= ZNET(KE1+K-1)+(DZ1-(K-1)*GRADEZ)
  31    CONTINUE
*
*  (c) from KE2 to KE3
*
  29  NCVZ= KE3-KE2
      IF(NCVZ.LT.1) GOTO 32
        SUMDZ= 0.0
        DO 33 K=1,NCVZ
          SUMDZ= SUMDZ+(1.0-(K-1)*(1.0-1.0/IRATZ3)/(NCVZ-1))
  33    CONTINUE
        DZ1= GRDZ3/SUMDZ
        GRADEZ= DZ1/(NCVZ-1)*(1.0-1.0/IRATZ3)
        ZNET(KE2+1)= ZNET(KE2)+DZ1
        DO 34 K=2,NCVZ
          ZNET(KE2+K)= ZNET(KE2+K-1)+(DZ1-(K-1)*GRADEZ)
  34    CONTINUE  
*
*  (d) from KE3 to KE
*
  32  NCVZ= KE-KE3-KEW
      IF(NCVZ.LT.1) GOTO 35
        SUMDZ= 0.0
        DO 36 K=1,NCVZ
          SUMDZ= SUMDZ+(1.0-(K-1)*(1.0-1.0/IRATZ4)/(NCVZ-1))
  36    CONTINUE
        DZ1= (GRDZ-GRDZ3-GRDZ1-GRDZ2-Z_KE)/SUMDZ
        GRADEZ= DZ1/(NCVZ-1)*(1.0-1.0/IRATZ4)
        ZNET(KE3+1)= ZNET(KE3)+DZ1
        DO 37 K=2,NCVZ
          ZNET(KE3+K)= ZNET(KE3+K-1)+(DZ1-(K-1)*GRADEZ)
  37    CONTINUE
        IF(KEW.EQ.1) THEN
          ZNET(KE)= GRDZ
        ENDIF
  35    ZNET(KE+1)= ZNET(KE)
*
      RETURN
*
 111  WRITE(*,*) 'ERROR IN makgrd.f'
      STOP
      END

