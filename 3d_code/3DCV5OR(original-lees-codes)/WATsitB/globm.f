*
*     file globm.f
*********************************************************************
*
      SUBROUTINE GLOBMI(GLOBM,RHS,IA,JA,KA,KVLST,NEQA,
     C                  ATP,ATW,ATE,ATS,ATN,ATB,ATT,BT,
     C                  IB,IE,JB,JE,KB,KE,ID,JD,KD,N,
     C                  MG,MR,NELEMX,NJAMAX)
*
*     Subroutine to construct the coefficient array in the structure 
*     required by WATsit1.
*
*     GLOBM(MG)      global array in IA(),JA() format; output
*     RHS(MR)        right-hand-side vector; output
*     
*     AT#(ID,JD,KD)  active coefficient arrays
*     N              number of unknowns in block (=1 for scalar)
*     IB,IE          starting and ending indices in 'x' direction
*     JB,JE          starting and ending indices in 'y' direction
*     KB,KE          starting and ending indices in 'z' direction
*     ID,JD,KD       array dimensions of block coefficients
*
*     NOTE:  The structure of the matrix is set up such that the
*            diagonal of each row is the first element in the packed
*            global matrix.
*
*********************************************************************
*
      IMPLICIT NONE
      DOUBLE PRECISION GLOBM(MG),RHS(MR)
      REAL*8 ATP(ID,JD,KD),ATW(ID,JD,KD),ATE(ID,JD,KD)
      REAL*8 ATS(ID,JD,KD),ATN(ID,JD,KD),BT(ID,JD,KD)
      REAL*8 ATB(ID,JD,KD),ATT(ID,JD,KD)
      INTEGER IA(NELEMX+1),JA(NJAMAX),KA(NJAMAX+1)
      INTEGER KVLST(NELEMX),NEQA(NELEMX)
      INTEGER IB,IE,JB,JE,KB,KE,E,N,ID,JD,KD
      INTEGER MG,MR,NELEMX,NJAMAX,ROW
      INTEGER W,Z,Y,I,J,K,NODE,NNODES
*
      W=IE-IB+1
      Z=JE-JB+1
*
*-------------------------------------------------------
*  [1] Form the IA(), JA() arrays
*-------------------------------------------------------
*
*  Initialize indices and counters
*
      I=IB
      J=JB
      K=KB
      ROW=1
      E=1
*
*  Bottom face control-volumes
*
      IA(ROW)=E
      JA(E)=ROW
      E=E+1
      IF(IE .GT. IB) THEN
        JA(E)=ROW+1
        E=E+1
      ENDIF
      IF(JE .GT. JB) THEN
        JA(E)=ROW+W
        E=E+1
      ENDIF
      IF(KE .GT. KB) THEN
        JA(E)=ROW+W*Z
        E=E+1
      ENDIF
*
      DO 5 I=IB+1,IE-1
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
	  JA(E+1)=ROW-1
	  JA(E+2)=ROW+1
 	  E=E+3
        IF(JE .GT. JB) THEN
          JA(E)=ROW+W
          E=E+1
        ENDIF
        IF(KE .GT. KB) THEN
          JA(E)=ROW+W*Z
          E=E+1
        ENDIF
  5   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
	  JA(E+1)=ROW-1
	  E=E+2
        IF(JE .GT. JB) THEN
          JA(E)=ROW+W
          E=E+1
        ENDIF
        IF(KE .GT. KB) THEN
          JA(E)=ROW+W*Z
          E=E+1
        ENDIF
      ENDIF
*
      DO 7 J=JB+1,JE-1
*
      I=IB
      ROW=ROW+1
      IA(ROW)=E
      JA(E)=ROW
      JA(E+1)=ROW-W
      E=E+2
      IF(IE .GT. IB) THEN
        JA(E)=ROW+1
        E=E+1
      ENDIF
      JA(E)=ROW+W
      E=E+1
      IF(KE .GT. KB) THEN
        JA(E)=ROW+W*Z
        E=E+1
      ENDIF
*
      DO 6 I=IB+1,IE-1
        ROW=ROW+1
	  IA(ROW)=E
        JA(E  )=ROW
        JA(E+1)=ROW-W
        JA(E+2)=ROW-1
	  JA(E+3)=ROW+1
        JA(E+4)=ROW+W
 	  E=E+5
        IF(KE .GT. KB) THEN
	    JA(E)=ROW+W*Z
	    E=E+1
        ENDIF
  6   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W
	  JA(E+2)=ROW-1
        JA(E+3)=ROW+W
	  E=E+4
        IF(KE .GT. KB) THEN
	    JA(E)=ROW+W*Z
	    E=E+1
        ENDIF
      ENDIF
*
  7   CONTINUE
*      
      IF(JE .GT. JB) THEN
*
      J=JE
      I=IB
      ROW=ROW+1
      IA(ROW)=E
      JA(E)=ROW
      JA(E+1)=ROW-W
      E=E+2
      IF(IE .GT. IB) THEN
	  JA(E)=ROW+1
	  E=E+1
      ENDIF
      IF(KE .GT. KB) THEN
	  JA(E)=ROW+W*Z
	  E=E+1
      ENDIF
*
      DO 8 I=IB+1,IE-1
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W
	  JA(E+2)=ROW-1
	  JA(E+3)=ROW+1
 	  E=E+4
        IF(KE .GT. KB) THEN
	    JA(E)=ROW+W*Z
	    E=E+1
        ENDIF
  8   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W
	  JA(E+2)=ROW-1
	  E=E+3
        IF(KE .GT. KB) THEN
	    JA(E)=ROW+W*Z
	    E=E+1
        ENDIF
      ENDIF
*
      ENDIF
*
*  Interior control-volumes
*
      DO 15 K=KB+1,KE-1
*
      I=IB
      J=JB
      ROW=ROW+1
      IA(ROW)=E
      JA(E)=ROW
      JA(E+1)=ROW-W*Z
      E=E+2
      IF(IE .GT. IB) THEN
	  JA(E)=ROW+1
	  E=E+1
      ENDIF
      IF(JE .GT. JB) THEN
	  JA(E)=ROW+W
	  E=E+1
      ENDIF
      JA(E)=ROW+W*Z
      E=E+1

*
      DO 9 I=IB+1,IE-1
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W*Z
	  JA(E+2)=ROW-1
	  JA(E+3)=ROW+1
 	  E=E+4
        IF(JE .GT. JB) THEN
	    JA(E)=ROW+W
	    E=E+1
        ENDIF
        JA(E)=ROW+W*Z
        E=E+1
  9   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W*Z
	  JA(E+2)=ROW-1
	  E=E+3
        IF(JE .GT. JB) THEN
	    JA(E)=ROW+W
	    E=E+1
        ENDIF
        JA(E)=ROW+W*Z
        E=E+1
      ENDIF
*
      DO 11 J=JB+1,JE-1
*
      I=IB
      ROW=ROW+1
      IA(ROW)=E
      JA(E)=ROW
      JA(E+1)=ROW-W*Z
      JA(E+2)=ROW-W
      E=E+3
      IF(IE .GT. IB) THEN
	  JA(E)=ROW+1
	  E=E+1
      ENDIF
      JA(E)=ROW+W
      JA(E+1)=ROW+W*Z
      E=E+2
*
      DO 10 I=IB+1,IE-1
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W*Z
        JA(E+2)=ROW-W
	  JA(E+3)=ROW-1
	  JA(E+4)=ROW+1
        JA(E+5)=ROW+W
        JA(E+6)=ROW+W*Z
        E=E+7
 10   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W*Z
        JA(E+2)=ROW-W
	  JA(E+3)=ROW-1
        JA(E+4)=ROW+W
        JA(E+5)=ROW+W*Z
        E=E+6
      ENDIF
*
 11   CONTINUE
*      
      IF(JE .GT. JB) THEN
*
      J=JE
      I=IB
      ROW=ROW+1
      IA(ROW)=E
      JA(E)=ROW
      JA(E+1)=ROW-W*Z
      JA(E+2)=ROW-W
      E=E+3
      IF(IE .GT. IB) THEN
        JA(E)=ROW+1
        E=E+1
      ENDIF
      JA(E)=ROW+W*Z
      E=E+1
*
      DO 12 I=IB+1,IE-1
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W*Z
        JA(E+2)=ROW-W
	  JA(E+3)=ROW-1
	  JA(E+4)=ROW+1
        JA(E+5)=ROW+W*Z
        E=E+6
 12   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W*Z
	  JA(E+2)=ROW-W
	  JA(E+3)=ROW-1
        JA(E+4)=ROW+W*Z
        E=E+5
      ENDIF
*
      ENDIF
*
 15   CONTINUE
*
*  Top face control-volumes
*
      IF(KE .GT. KB) THEN
*
      I=IB
      J=JB
      K=KE
      ROW=ROW+1
      IA(ROW)=E
      JA(E)=ROW
      JA(E+1)=ROW-W*Z
      E=E+2
      IF(IE .GT. IB) THEN
        JA(E)=ROW+1
        E=E+1
      ENDIF
      IF(JE .GT. JB) THEN
        JA(E)=ROW+W
        E=E+1
      ENDIF
*
      DO 16 I=IB+1,IE-1
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W*Z
	  JA(E+2)=ROW-1
	  JA(E+3)=ROW+1
 	  E=E+4
        IF(JE .GT. JB) THEN
          JA(E)=ROW+W
          E=E+1
        ENDIF
 16   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W*Z
        JA(E+2)=ROW-1
        E=E+3
        IF(JE .GT. JB) THEN
          JA(E)=ROW+W
          E=E+1
        ENDIF
      ENDIF
*
      DO 18 J=JB+1,JE-1
*
      I=IB
      ROW=ROW+1
      IA(ROW)=E
      JA(E)=ROW
      JA(E+1)=ROW-W*Z
      JA(E+2)=ROW-W
      E=E+3
      IF(IE .GT. IB) THEN
        JA(E)=ROW+1
        E=E+1
      ENDIF
      JA(E)=ROW+W
      E=E+1
*
      DO 17 I=IB+1,IE-1
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W*Z
        JA(E+2)=ROW-W
	  JA(E+3)=ROW-1
	  JA(E+4)=ROW+1
        JA(E+5)=ROW+W
        E=E+6
 17   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W*Z
        JA(E+2)=ROW-W
	  JA(E+3)=ROW-1
        JA(E+4)=ROW+W
        E=E+5
      ENDIF
*
 18   CONTINUE
*      
      IF(JE .GT. JB) THEN
*
      J=JE
      I=IB
      ROW=ROW+1
      IA(ROW)=E
      JA(E)=ROW
      JA(E+1)=ROW-W*Z
      JA(E+2)=ROW-W
      E=E+3
      IF(IE .GT. IB) THEN
	  JA(E)=ROW+1
	  E=E+1
      ENDIF
*
      DO 19 I=IB+1,IE-1
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W*Z
        JA(E+2)=ROW-W
	  JA(E+3)=ROW-1
	  JA(E+4)=ROW+1
 	  E=E+5
 19   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W*Z
        JA(E+2)=ROW-W
	  JA(E+3)=ROW-1
        E=E+4
      ENDIF
*
      ENDIF
*
      ENDIF
*
*  Indicate the position after the last entry of JA()
*
      IA(ROW+1)=E
      E=E-1
*
*-------------------------------------------------------
*  [2] Form the NEQA(), KVLST() and KA() arrays
*-------------------------------------------------------
*
      NNODES= (IE-IB+1)*(JE-JB+1)*(KE-KB+1)
      DO 25 NODE=1,NNODES
	 NEQA(NODE)= N
  25  CONTINUE
*
      KVLST(1)= 1
      DO 30 NODE=1,NNODES
	 KVLST(NODE+1)= KVLST(NODE) + NEQA(NODE)
  30  CONTINUE
*
      KA(1)= 1
      DO 40 NODE=1,NNODES
	 DO 35 Y=IA(NODE),IA(NODE+1)-1
	  KA(Y+1)= KA(Y) + NEQA(NODE)*NEQA(JA(Y))
  35   CONTINUE
  40  CONTINUE
*
*-------------------------------------------------------
*  [3] Fill in the non-zero elements in matrix GLOBM()
*-------------------------------------------------------
*
*  Initialize indices and counters
*
      I=IB
      J=JB
      K=KB
      ROW=1
      E=1
*
*  Bottom face control-volumes
*
      GLOBM(KA(E))= DBLE(ATP(I,J,K))
      E=E+1
      IF(IE .GT. IB) THEN
        GLOBM(KA(E))= DBLE(-ATE(I,J,K))
	E=E+1
      ENDIF
      IF(JE .GT. JB) THEN
        GLOBM(KA(E))= DBLE(-ATN(I,J,K))
	E=E+1
      ENDIF
      IF(KE .GT. KB) THEN
  	  GLOBM(KA(E))= DBLE(-ATT(I,J,K))
	  E=E+1
      ENDIF
      RHS(KVLST(ROW))= DBLE(BT(I,J,K))
*
      DO 45 I=IB+1,IE-1
        ROW=ROW+1
        GLOBM(KA(E))= DBLE(ATP(I,J,K))
        GLOBM(KA(E+1))= DBLE(-ATW(I,J,K))
        GLOBM(KA(E+2))= DBLE(-ATE(I,J,K))
        E=E+3
        IF(JE .GT. JB) THEN
	    GLOBM(KA(E))= DBLE(-ATN(I,J,K))
	    E=E+1
        ENDIF
        IF(KE .GT. KB) THEN
	    GLOBM(KA(E))= DBLE(-ATT(I,J,K))
	    E=E+1
        ENDIF
        RHS(KVLST(ROW))= DBLE(BT(I,J,K))
 45   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
	  GLOBM(KA(E))= DBLE(ATP(I,J,K))
        GLOBM(KA(E+1))= DBLE(-ATW(I,J,K))
	  E=E+2
        IF(JE .GT. JB) THEN
	    GLOBM(KA(E))= DBLE(-ATN(I,J,K))
	    E=E+1
        ENDIF
        IF(KE .GT. KB) THEN
	    GLOBM(KA(E))= DBLE(-ATT(I,J,K))
	    E=E+1
        ENDIF
        RHS(KVLST(ROW))= DBLE(BT(I,J,K))
      ENDIF
*
      DO 47 J=JB+1,JE-1
*
      I=IB
      ROW=ROW+1
      GLOBM(KA(E))= DBLE(ATP(I,J,K))
      GLOBM(KA(E+1))= DBLE(-ATS(I,J,K))
      E=E+2
      IF(IE .GT. IB) THEN
        GLOBM(KA(E))= DBLE(-ATE(I,J,K))
	  E=E+1
      ENDIF
      GLOBM(KA(E))= DBLE(-ATN(I,J,K))
	E=E+1
      IF(KE .GT. KB) THEN
	  GLOBM(KA(E))= DBLE(-ATT(I,J,K))
	  E=E+1
      ENDIF
      RHS(KVLST(ROW))= DBLE(BT(I,J,K))
*
      DO 46 I=IB+1,IE-1
        ROW=ROW+1
	  GLOBM(KA(E))= DBLE(ATP(I,J,K))
        GLOBM(KA(E+1))= DBLE(-ATS(I,J,K))
	  GLOBM(KA(E+2))= DBLE(-ATW(I,J,K))
	  GLOBM(KA(E+3))= DBLE(-ATE(I,J,K))
        GLOBM(KA(E+4))= DBLE(-ATN(I,J,K))
        E=E+5
        IF(KE .GT. KB) THEN
	    GLOBM(KA(E))= DBLE(-ATT(I,J,K))
	    E=E+1
        ENDIF
        RHS(KVLST(ROW))= DBLE(BT(I,J,K))
 46   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
	  GLOBM(KA(E))= DBLE(ATP(I,J,K))
        GLOBM(KA(E+1))= DBLE(-ATS(I,J,K))
	  GLOBM(KA(E+2))= DBLE(-ATW(I,J,K))
        GLOBM(KA(E+3))= DBLE(-ATN(I,J,K))
	  E=E+4
        IF(KE .GT. KB) THEN
	    GLOBM(KA(E))= DBLE(-ATT(I,J,K))
	    E=E+1
        ENDIF
        RHS(KVLST(ROW))= DBLE(BT(I,J,K))
      ENDIF
*
 47   CONTINUE
*      
      IF(JE .GT. JB) THEN
*
      J=JE
      I=IB
      ROW=ROW+1
      GLOBM(KA(E))= DBLE(ATP(I,J,K))
      GLOBM(KA(E+1))= DBLE(-ATS(I,J,K))
      E=E+2
      IF(IE .GT. IB) THEN
	  GLOBM(KA(E))= DBLE(-ATE(I,J,K))
	  E=E+1
      ENDIF
      IF(KE .GT. KB) THEN
	  GLOBM(KA(E))= DBLE(-ATT(I,J,K))
	  E=E+1
      ENDIF
      RHS(KVLST(ROW))= DBLE(BT(I,J,K))
*
      DO 48 I=IB+1,IE-1
        ROW=ROW+1
	  GLOBM(KA(E))= DBLE(ATP(I,J,K))
        GLOBM(KA(E+1))= DBLE(-ATS(I,J,K))
	  GLOBM(KA(E+2))= DBLE(-ATW(I,J,K))
	  GLOBM(KA(E+3))= DBLE(-ATE(I,J,K))
 	  E=E+4
        IF(KE .GT. KB) THEN
	    GLOBM(KA(E))= DBLE(-ATT(I,J,K))
	    E=E+1
        ENDIF
        RHS(KVLST(ROW))= DBLE(BT(I,J,K))
 48   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
	  GLOBM(KA(E))= DBLE(ATP(I,J,K))
        GLOBM(KA(E+1))= DBLE(-ATS(I,J,K))
	  GLOBM(KA(E+2))= DBLE(-ATW(I,J,K))
	  E=E+3
        IF(KE .GT. KB) THEN
	    GLOBM(KA(E))= DBLE(-ATT(I,J,K))
	    E=E+1
        ENDIF
        RHS(KVLST(ROW))= DBLE(BT(I,J,K))
      ENDIF
*
      ENDIF
*
*  Interior control-volumes
*
      DO 55 K=KB+1,KE-1
*
      I=IB
      J=JB
      ROW=ROW+1
      GLOBM(KA(E))= DBLE(ATP(I,J,K))
      GLOBM(KA(E+1))= DBLE(-ATB(I,J,K))
      E=E+2
      IF(IE .GT. IB) THEN
	  GLOBM(KA(E))= DBLE(-ATE(I,J,K))
	  E=E+1
      ENDIF
      IF(JE .GT. JB) THEN
	  GLOBM(KA(E))= DBLE(-ATN(I,J,K))
	  E=E+1
      ENDIF
      GLOBM(KA(E))= DBLE(-ATT(I,J,K))
      E=E+1
      RHS(KVLST(ROW))= DBLE(BT(I,J,K))
*
      DO 49 I=IB+1,IE-1
        ROW=ROW+1
	  GLOBM(KA(E))= DBLE(ATP(I,J,K))
        GLOBM(KA(E+1))= DBLE(-ATB(I,J,K))
	  GLOBM(KA(E+2))= DBLE(-ATW(I,J,K))
	  GLOBM(KA(E+3))= DBLE(-ATE(I,J,K))
 	  E=E+4
        IF(JE .GT. JB) THEN
	    GLOBM(KA(E))= DBLE(-ATN(I,J,K))
	    E=E+1
        ENDIF
        GLOBM(KA(E))= DBLE(-ATT(I,J,K))
        E=E+1
        RHS(KVLST(ROW))= DBLE(BT(I,J,K))
 49   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
	  GLOBM(KA(E))= DBLE(ATP(I,J,K))
        GLOBM(KA(E+1))= DBLE(-ATB(I,J,K))
	  GLOBM(KA(E+2))= DBLE(-ATW(I,J,K))
	  E=E+3
        IF(JE .GT. JB) THEN
	    GLOBM(KA(E))= DBLE(-ATN(I,J,K))
	    E=E+1
        ENDIF
        GLOBM(KA(E))= DBLE(-ATT(I,J,K))
        E=E+1
        RHS(KVLST(ROW))= DBLE(BT(I,J,K))
      ENDIF
*
      DO 51 J=JB+1,JE-1
*
      I=IB
      ROW=ROW+1
      GLOBM(KA(E))= DBLE(ATP(I,J,K))
      GLOBM(KA(E+1))= DBLE(-ATB(I,J,K))
      GLOBM(KA(E+2))= DBLE(-ATS(I,J,K))
      E=E+3
      IF(IE .GT. IB) THEN
	  GLOBM(KA(E))= DBLE(-ATE(I,J,K))
	  E=E+1
      ENDIF
      GLOBM(KA(E))= DBLE(-ATN(I,J,K))
      GLOBM(KA(E+1))= DBLE(-ATT(I,J,K))
      E=E+2
      RHS(KVLST(ROW))= DBLE(BT(I,J,K))
*
      DO 50 I=IB+1,IE-1
        ROW=ROW+1
	  GLOBM(KA(E))= DBLE(ATP(I,J,K))
        GLOBM(KA(E+1))= DBLE(-ATB(I,J,K))
        GLOBM(KA(E+2))= DBLE(-ATS(I,J,K))
	  GLOBM(KA(E+3))= DBLE(-ATW(I,J,K))
	  GLOBM(KA(E+4))= DBLE(-ATE(I,J,K))
        GLOBM(KA(E+5))= DBLE(-ATN(I,J,K))
        GLOBM(KA(E+6))= DBLE(-ATT(I,J,K))
        E=E+7
        RHS(KVLST(ROW))= DBLE(BT(I,J,K))
 50   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
	  GLOBM(KA(E))= DBLE(ATP(I,J,K))
        GLOBM(KA(E+1))= DBLE(-ATB(I,J,K))
        GLOBM(KA(E+2))= DBLE(-ATS(I,J,K))
	  GLOBM(KA(E+3))= DBLE(-ATW(I,J,K))
        GLOBM(KA(E+4))= DBLE(-ATN(I,J,K))
        GLOBM(KA(E+5))= DBLE(-ATT(I,J,K))
        E=E+6
        RHS(KVLST(ROW))= DBLE(BT(I,J,K))
      ENDIF
*
 51   CONTINUE
*      
      IF(JE .GT. JB) THEN
*
      J=JE
      I=IB
      ROW=ROW+1
      GLOBM(KA(E))= DBLE(ATP(I,J,K))
      GLOBM(KA(E+1))= DBLE(-ATB(I,J,K))
      GLOBM(KA(E+2))= DBLE(-ATS(I,J,K))
      E=E+3
      IF(IE .GT. IB) THEN
	  GLOBM(KA(E))= DBLE(-ATE(I,J,K))
	  E=E+1
      ENDIF
      GLOBM(KA(E))= DBLE(-ATT(I,J,K))
      E=E+1
      RHS(KVLST(ROW))= DBLE(BT(I,J,K))
*
      DO 52 I=IB+1,IE-1
        ROW=ROW+1
	  GLOBM(KA(E))= DBLE(ATP(I,J,K))
        GLOBM(KA(E+1))= DBLE(-ATB(I,J,K))
        GLOBM(KA(E+2))= DBLE(-ATS(I,J,K))
	  GLOBM(KA(E+3))= DBLE(-ATW(I,J,K))
	  GLOBM(KA(E+4))= DBLE(-ATE(I,J,K))
        GLOBM(KA(E+5))= DBLE(-ATT(I,J,K))
 	  E=E+6
        RHS(KVLST(ROW))= DBLE(BT(I,J,K))
 52   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
	  GLOBM(KA(E))= DBLE(ATP(I,J,K))
        GLOBM(KA(E+1))= DBLE(-ATB(I,J,K))
        GLOBM(KA(E+2))= DBLE(-ATS(I,J,K))
	  GLOBM(KA(E+3))= DBLE(-ATW(I,J,K))
        GLOBM(KA(E+4))= DBLE(-ATT(I,J,K))
	  E=E+5
        RHS(KVLST(ROW))= DBLE(BT(I,J,K))
      ENDIF
*
      ENDIF
*
 55   CONTINUE
*
*  Top face control-volumes
*
      IF(KE .GT. KB) THEN
*
      I=IB
      J=JB
      K=KE
      ROW=ROW+1
      GLOBM(KA(E))= DBLE(ATP(I,J,K))
      GLOBM(KA(E+1))= DBLE(-ATB(I,J,K))
      E=E+2
      IF(IE .GT. IB) THEN
	  GLOBM(KA(E))= DBLE(-ATE(I,J,K))
	  E=E+1
      ENDIF
      IF(JE .GT. JB) THEN
	  GLOBM(KA(E))= DBLE(-ATN(I,J,K))
	  E=E+1
      ENDIF
      RHS(KVLST(ROW))= DBLE(BT(I,J,K))
*
      DO 56 I=IB+1,IE-1
        ROW=ROW+1
        GLOBM(KA(E))= DBLE(ATP(I,J,K))
        GLOBM(KA(E+1))= DBLE(-ATB(I,J,K))
	  GLOBM(KA(E+2))= DBLE(-ATW(I,J,K))
	  GLOBM(KA(E+3))= DBLE(-ATE(I,J,K))
 	  E=E+4
        IF(JE .GT. JB) THEN
	    GLOBM(KA(E))= DBLE(-ATN(I,J,K))
	    E=E+1
        ENDIF
        RHS(KVLST(ROW))= DBLE(BT(I,J,K))
 56   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
	  GLOBM(KA(E))= DBLE(ATP(I,J,K))
        GLOBM(KA(E+1))= DBLE(-ATB(I,J,K))
	  GLOBM(KA(E+2))= DBLE(-ATW(I,J,K))
	  E=E+3
        IF(JE .GT. JB) THEN
	    GLOBM(KA(E))= DBLE(-ATN(I,J,K))
	    E=E+1
        ENDIF
        RHS(KVLST(ROW))= DBLE(BT(I,J,K))
      ENDIF
*
      DO 58 J=JB+1,JE-1
*
      I=IB
      ROW=ROW+1
      GLOBM(KA(E))= DBLE(ATP(I,J,K))
      GLOBM(KA(E+1))= DBLE(-ATB(I,J,K))
      GLOBM(KA(E+2))= DBLE(-ATS(I,J,K))
      E=E+3
      IF(IE .GT. IB) THEN
	  GLOBM(KA(E))= DBLE(-ATE(I,J,K))
	  E=E+1
      ENDIF
      GLOBM(KA(E))= DBLE(-ATN(I,J,K))
      E=E+1
      RHS(KVLST(ROW))= DBLE(BT(I,J,K))
*
      DO 57 I=IB+1,IE-1
        ROW=ROW+1
	  GLOBM(KA(E))= DBLE(ATP(I,J,K))
        GLOBM(KA(E+1))= DBLE(-ATB(I,J,K))
        GLOBM(KA(E+2))= DBLE(-ATS(I,J,K))
	  GLOBM(KA(E+3))= DBLE(-ATW(I,J,K))
	  GLOBM(KA(E+4))= DBLE(-ATE(I,J,K))
        GLOBM(KA(E+5))= DBLE(-ATN(I,J,K))
 	  E=E+6
        RHS(KVLST(ROW))= DBLE(BT(I,J,K))
 57   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
	  GLOBM(KA(E))= DBLE(ATP(I,J,K))
        GLOBM(KA(E+1))= DBLE(-ATB(I,J,K))
        GLOBM(KA(E+2))= DBLE(-ATS(I,J,K))
	  GLOBM(KA(E+3))= DBLE(-ATW(I,J,K))
        GLOBM(KA(E+4))= DBLE(-ATN(I,J,K))
	  E=E+5
        RHS(KVLST(ROW))= DBLE(BT(I,J,K))
      ENDIF
*
 58   CONTINUE
*      
      IF(JE .GT. JB) THEN
*
      J=JE
      I=IB
      ROW=ROW+1
      GLOBM(KA(E))= DBLE(ATP(I,J,K))
      GLOBM(KA(E+1))= DBLE(-ATB(I,J,K))
      GLOBM(KA(E+2))= DBLE(-ATS(I,J,K))
      E=E+3
      IF(IE .GT. IB) THEN
	  GLOBM(KA(E))= DBLE(-ATE(I,J,K))
	  E=E+1
      ENDIF
      RHS(KVLST(ROW))= DBLE(BT(I,J,K))
*
      DO 59 I=IB+1,IE-1
        ROW=ROW+1
	  GLOBM(KA(E))= DBLE(ATP(I,J,K))
        GLOBM(KA(E+1))= DBLE(-ATB(I,J,K))
        GLOBM(KA(E+2))= DBLE(-ATS(I,J,K))
	  GLOBM(KA(E+3))= DBLE(-ATW(I,J,K))
	  GLOBM(KA(E+4))= DBLE(-ATE(I,J,K))
 	  E=E+5
        RHS(KVLST(ROW))= DBLE(BT(I,J,K))
 59   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
	  GLOBM(KA(E))= DBLE(ATP(I,J,K))
        GLOBM(KA(E+1))= DBLE(-ATB(I,J,K))
        GLOBM(KA(E+2))= DBLE(-ATS(I,J,K))
	  GLOBM(KA(E+3))= DBLE(-ATW(I,J,K))
	  E=E+4
        RHS(KVLST(ROW))= DBLE(BT(I,J,K))
      ENDIF
*
      ENDIF
*
      ENDIF
*
      RETURN
      END
*
*********************************************************************
*
      SUBROUTINE GLOBMS(GLOBM,RHS,IA,JA,KA,KVLST,NEQA,
     C                  ATP,ATW,ATE,ATS,ATN,ATB,ATT,BT,
     C                  IB,IE,JB,JE,KB,KE,ID,JD,KD,N,
     C                  MG,MR,NELEMX,NJAMAX)
*
*     Subroutine to construct the coefficient array in the structure 
*     required by WATsitB.
*
*     GLOBM(MG)      global array in IA(),JA(),KA() format; output
*     RHS(MR)        right-hand-side vector; output
*     
*     AT#(N,N,ID,JD,KD) NxN active coefficient block arrays
*     N              number of unknowns in block (=1 for scalar)
*     IB,IE          starting and ending indices in 'x' direction
*     JB,JE          starting and ending indices in 'y' direction
*     KB,KE          starting and ending indices in 'z' direction
*     ID,JD,KD          array dimensions of block coefficients
*
*     NOTE:  The structure of the matrix is set up such that the
*            diagonal of each row is the first element in the packed
*            global matrix.
*
*********************************************************************
*
      IMPLICIT NONE
      DOUBLE PRECISION GLOBM(MG),RHS(MR)
      REAL*8 ATP(N,N,ID,JD,KD),ATW(N,N,ID,JD,KD)
      REAL*8 ATE(N,N,ID,JD,KD)
      REAL*8 ATS(N,N,ID,JD,KD),ATN(N,N,ID,JD,KD),BT(N,ID,JD,KD)
      REAL*8 ATB(N,N,ID,JD,KD),ATT(N,N,ID,JD,KD)
      INTEGER IA(NELEMX+1),JA(NJAMAX),KA(NJAMAX+1)
      INTEGER KVLST(NELEMX),NEQA(NELEMX)
      INTEGER IB,IE,JB,JE,KB,KE,E,N,ID,JD,KD
      INTEGER MG,MR,NELEMX,NJAMAX,ROW
      INTEGER Y,W,Z,I,J,K,NODE,NNODES
*
      W=IE-IB+1
      Z=JE-JB+1
*
*-------------------------------------------------------
*  [1] Form the IA(), JA() arrays
*-------------------------------------------------------
*
*  Initialize indices and counters
*
      I=IB
      J=JB
      K=KB
      ROW=1
      E=1
*
*  Bottom face control-volumes
*
      IA(ROW)=E
      JA(E)=ROW
      E=E+1
      IF(IE .GT. IB) THEN
        JA(E)=ROW+1
        E=E+1
      ENDIF
      IF(JE .GT. JB) THEN
        JA(E)=ROW+W
        E=E+1
      ENDIF
      IF(KE .GT. KB) THEN
        JA(E)=ROW+W*Z
        E=E+1
      ENDIF
*
      DO 5 I=IB+1,IE-1
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
	  JA(E+1)=ROW-1
	  JA(E+2)=ROW+1
 	  E=E+3
        IF(JE .GT. JB) THEN
          JA(E)=ROW+W
          E=E+1
        ENDIF
        IF(KE .GT. KB) THEN
          JA(E)=ROW+W*Z
          E=E+1
        ENDIF
  5   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
	  JA(E+1)=ROW-1
	  E=E+2
        IF(JE .GT. JB) THEN
          JA(E)=ROW+W
          E=E+1
        ENDIF
        IF(KE .GT. KB) THEN
          JA(E)=ROW+W*Z
          E=E+1
        ENDIF
      ENDIF
*
      DO 7 J=JB+1,JE-1
*
      I=IB
      ROW=ROW+1
      IA(ROW)=E
      JA(E)=ROW
      JA(E+1)=ROW-W
      E=E+2
      IF(IE .GT. IB) THEN
        JA(E)=ROW+1
        E=E+1
      ENDIF
      JA(E)=ROW+W
      E=E+1
      IF(KE .GT. KB) THEN
        JA(E)=ROW+W*Z
        E=E+1
      ENDIF
*
      DO 6 I=IB+1,IE-1
        ROW=ROW+1
	  IA(ROW)=E
        JA(E  )=ROW
        JA(E+1)=ROW-W
        JA(E+2)=ROW-1
	  JA(E+3)=ROW+1
        JA(E+4)=ROW+W
 	  E=E+5
        IF(KE .GT. KB) THEN
	    JA(E)=ROW+W*Z
	    E=E+1
        ENDIF
  6   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W
	  JA(E+2)=ROW-1
        JA(E+3)=ROW+W
	  E=E+4
        IF(KE .GT. KB) THEN
	    JA(E)=ROW+W*Z
	    E=E+1
        ENDIF
      ENDIF
*
  7   CONTINUE
*      
      IF(JE .GT. JB) THEN
*
      J=JE
      I=IB
      ROW=ROW+1
      IA(ROW)=E
      JA(E)=ROW
      JA(E+1)=ROW-W
      E=E+2
      IF(IE .GT. IB) THEN
	  JA(E)=ROW+1
	  E=E+1
      ENDIF
      IF(KE .GT. KB) THEN
	  JA(E)=ROW+W*Z
	  E=E+1
      ENDIF
*
      DO 8 I=IB+1,IE-1
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W
	  JA(E+2)=ROW-1
	  JA(E+3)=ROW+1
 	  E=E+4
        IF(KE .GT. KB) THEN
	    JA(E)=ROW+W*Z
	    E=E+1
        ENDIF
  8   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W
	  JA(E+2)=ROW-1
	  E=E+3
        IF(KE .GT. KB) THEN
	    JA(E)=ROW+W*Z
	    E=E+1
        ENDIF
      ENDIF
*
      ENDIF
*
*  Interior control-volumes
*
      DO 15 K=KB+1,KE-1
*
      I=IB
      J=JB
      ROW=ROW+1
      IA(ROW)=E
      JA(E)=ROW
      JA(E+1)=ROW-W*Z
      E=E+2
      IF(IE .GT. IB) THEN
	  JA(E)=ROW+1
	  E=E+1
      ENDIF
      IF(JE .GT. JB) THEN
	  JA(E)=ROW+W
	  E=E+1
      ENDIF
      JA(E)=ROW+W*Z
      E=E+1

*
      DO 9 I=IB+1,IE-1
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W*Z
	  JA(E+2)=ROW-1
	  JA(E+3)=ROW+1
 	  E=E+4
        IF(JE .GT. JB) THEN
	    JA(E)=ROW+W
	    E=E+1
        ENDIF
        JA(E)=ROW+W*Z
        E=E+1
  9   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W*Z
	  JA(E+2)=ROW-1
	  E=E+3
        IF(JE .GT. JB) THEN
	    JA(E)=ROW+W
	    E=E+1
        ENDIF
        JA(E)=ROW+W*Z
        E=E+1
      ENDIF
*
      DO 11 J=JB+1,JE-1
*
      I=IB
      ROW=ROW+1
      IA(ROW)=E
      JA(E)=ROW
      JA(E+1)=ROW-W*Z
      JA(E+2)=ROW-W
      E=E+3
      IF(IE .GT. IB) THEN
	  JA(E)=ROW+1
	  E=E+1
      ENDIF
      JA(E)=ROW+W
      JA(E+1)=ROW+W*Z
      E=E+2
*
      DO 10 I=IB+1,IE-1
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W*Z
        JA(E+2)=ROW-W
	  JA(E+3)=ROW-1
	  JA(E+4)=ROW+1
        JA(E+5)=ROW+W
        JA(E+6)=ROW+W*Z
        E=E+7
 10   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W*Z
        JA(E+2)=ROW-W
	  JA(E+3)=ROW-1
        JA(E+4)=ROW+W
        JA(E+5)=ROW+W*Z
        E=E+6
      ENDIF
*
 11   CONTINUE
*      
      IF(JE .GT. JB) THEN
*
      J=JE
      I=IB
      ROW=ROW+1
      IA(ROW)=E
      JA(E)=ROW
      JA(E+1)=ROW-W*Z
      JA(E+2)=ROW-W
      E=E+3
      IF(IE .GT. IB) THEN
        JA(E)=ROW+1
        E=E+1
      ENDIF
      JA(E)=ROW+W*Z
      E=E+1
*
      DO 12 I=IB+1,IE-1
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W*Z
        JA(E+2)=ROW-W
	  JA(E+3)=ROW-1
	  JA(E+4)=ROW+1
        JA(E+5)=ROW+W*Z
        E=E+6
 12   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W*Z
	  JA(E+2)=ROW-W
	  JA(E+3)=ROW-1
        JA(E+4)=ROW+W*Z
        E=E+5
      ENDIF
*
      ENDIF
*
 15   CONTINUE
*
*  Top face control-volumes
*
      IF(KE .GT. KB) THEN
*
      I=IB
      J=JB
      K=KE
      ROW=ROW+1
      IA(ROW)=E
      JA(E)=ROW
      JA(E+1)=ROW-W*Z
      E=E+2
      IF(IE .GT. IB) THEN
        JA(E)=ROW+1
        E=E+1
      ENDIF
      IF(JE .GT. JB) THEN
        JA(E)=ROW+W
        E=E+1
      ENDIF
*
      DO 16 I=IB+1,IE-1
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W*Z
	  JA(E+2)=ROW-1
	  JA(E+3)=ROW+1
 	  E=E+4
        IF(JE .GT. JB) THEN
          JA(E)=ROW+W
          E=E+1
        ENDIF
 16   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W*Z
        JA(E+2)=ROW-1
        E=E+3
        IF(JE .GT. JB) THEN
          JA(E)=ROW+W
          E=E+1
        ENDIF
      ENDIF
*
      DO 18 J=JB+1,JE-1
*
      I=IB
      ROW=ROW+1
      IA(ROW)=E
      JA(E)=ROW
      JA(E+1)=ROW-W*Z
      JA(E+2)=ROW-W
      E=E+3
      IF(IE .GT. IB) THEN
        JA(E)=ROW+1
        E=E+1
      ENDIF
      JA(E)=ROW+W
      E=E+1
*
      DO 17 I=IB+1,IE-1
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W*Z
        JA(E+2)=ROW-W
	  JA(E+3)=ROW-1
	  JA(E+4)=ROW+1
        JA(E+5)=ROW+W
        E=E+6
 17   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W*Z
        JA(E+2)=ROW-W
	  JA(E+3)=ROW-1
        JA(E+4)=ROW+W
        E=E+5
      ENDIF
*
 18   CONTINUE
*      
      IF(JE .GT. JB) THEN
*
      J=JE
      I=IB
      ROW=ROW+1
      IA(ROW)=E
      JA(E)=ROW
      JA(E+1)=ROW-W*Z
      JA(E+2)=ROW-W
      E=E+3
      IF(IE .GT. IB) THEN
	  JA(E)=ROW+1
	  E=E+1
      ENDIF
*
      DO 19 I=IB+1,IE-1
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W*Z
        JA(E+2)=ROW-W
	  JA(E+3)=ROW-1
	  JA(E+4)=ROW+1
 	  E=E+5
 19   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
	  IA(ROW)=E
	  JA(E  )=ROW
        JA(E+1)=ROW-W*Z
        JA(E+2)=ROW-W
	  JA(E+3)=ROW-1
        E=E+4
      ENDIF
*
      ENDIF
*
      ENDIF
*
*  Indicate the position after the last entry of JA()
*
      IA(ROW+1)=E
      E=E-1
*
*-------------------------------------------------------
*  [2] Form the NEQA(), KVLST() and KA() arrays
*-------------------------------------------------------
*
      NNODES= (IE-IB+1)*(JE-JB+1)*(KE-KB+1)
      DO 25 NODE=1,NNODES
	 NEQA(NODE)= N
  25  CONTINUE
*
      KVLST(1)= 1
      DO 30 NODE=1,NNODES
	 KVLST(NODE+1)= KVLST(NODE) + NEQA(NODE)
  30  CONTINUE
*
      KA(1)= 1
      DO 40 NODE=1,NNODES
	 DO 35 Y=IA(NODE),IA(NODE+1)-1
	  KA(Y+1)= KA(Y) + NEQA(NODE)*NEQA(JA(Y))
  35   CONTINUE
  40  CONTINUE
*
*-------------------------------------------------------
*  [3] Fill in the non-zero elements in matrix GLOBM()
*-------------------------------------------------------
*
*  Initialize indices and counters
*
      I=IB
      J=JB
      K=KB
      ROW=1
      E=1
*
*  Bottom face control-volumes
*
      CALL GMIJKL(GLOBM, KA(E),ATP,I,J,K,1.0,ID,JD,KD,MG,N)
      E=E+1
      IF(IE .GT. IB) THEN
         CALL GMIJKL(GLOBM, KA(E),ATE,I,J,K,-1.0,ID,JD,KD,MG,N)
	E=E+1
      ENDIF
      IF(JE .GT. JB) THEN
         CALL GMIJKL(GLOBM, KA(E),ATN,I,J,K,-1.0,ID,JD,KD,MG,N)
	E=E+1
      ENDIF
      IF(KE .GT. KB) THEN
         CALL GMIJKL(GLOBM, KA(E),ATT,I,J,K,-1.0,ID,JD,KD,MG,N)
	E=E+1
      ENDIF
      CALL RHSIKL(RHS, KVLST(ROW),BT,I,J,K,1.0,ID,JD,KD,MR,N)
*
      DO 45 I=IB+1,IE-1
        ROW=ROW+1
        CALL GMIJKL(GLOBM, KA(E),ATP,I,J,K,1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+1),ATW,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+2),ATE,I,J,K,-1.0,ID,JD,KD,MG,N)
        E=E+3
        IF(JE .GT. JB) THEN
          CALL GMIJKL(GLOBM, KA(E),ATN,I,J,K,-1.0,ID,JD,KD,MG,N)
	    E=E+1
        ENDIF
        IF(KE .GT. KB) THEN
          CALL GMIJKL(GLOBM, KA(E),ATT,I,J,K,-1.0,ID,JD,KD,MG,N)
	    E=E+1
        ENDIF
        CALL RHSIKL(RHS, KVLST(ROW),BT,I,J,K,1.0,ID,JD,KD,MR,N)
 45   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
        CALL GMIJKL(GLOBM, KA(E),ATP,I,J,K,1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+1),ATW,I,J,K,-1.0,ID,JD,KD,MG,N)
	  E=E+2
        IF(JE .GT. JB) THEN
          CALL GMIJKL(GLOBM, KA(E),ATN,I,J,K,-1.0,ID,JD,KD,MG,N)
	    E=E+1
        ENDIF
        IF(KE .GT. KB) THEN
          CALL GMIJKL(GLOBM, KA(E),ATT,I,J,K,-1.0,ID,JD,KD,MG,N)
	    E=E+1
        ENDIF
        CALL RHSIKL(RHS, KVLST(ROW),BT,I,J,K,1.0,ID,JD,KD,MR,N)
      ENDIF
*
      DO 47 J=JB+1,JE-1
*
      I=IB
      ROW=ROW+1
      CALL GMIJKL(GLOBM, KA(E),ATP,I,J,K,1.0,ID,JD,KD,MG,N)
      CALL GMIJKL(GLOBM, KA(E+1),ATS,I,J,K,-1.0,ID,JD,KD,MG,N)
      E=E+2
      IF(IE .GT. IB) THEN
        CALL GMIJKL(GLOBM, KA(E),ATE,I,J,K,-1.0,ID,JD,KD,MG,N)
	  E=E+1
      ENDIF
      CALL GMIJKL(GLOBM, KA(E),ATN,I,J,K,-1.0,ID,JD,KD,MG,N)
      E=E+1
      IF(KE .GT. KB) THEN
        CALL GMIJKL(GLOBM, KA(E),ATT,I,J,K,-1.0,ID,JD,KD,MG,N)
	  E=E+1
      ENDIF
      CALL RHSIKL(RHS, KVLST(ROW),BT,I,J,K,1.0,ID,JD,KD,MR,N)
*
      DO 46 I=IB+1,IE-1
        ROW=ROW+1
        CALL GMIJKL(GLOBM, KA(E),ATP,I,J,K,1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+1),ATS,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+2),ATW,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+3),ATE,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+4),ATN,I,J,K,-1.0,ID,JD,KD,MG,N)
 	  E=E+5
        IF(KE .GT. KB) THEN
          CALL GMIJKL(GLOBM, KA(E),ATT,I,J,K,-1.0,ID,JD,KD,MG,N)
	    E=E+1
        ENDIF
        CALL RHSIKL(RHS, KVLST(ROW),BT,I,J,K,1.0,ID,JD,KD,MR,N)
 46   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
        CALL GMIJKL(GLOBM, KA(E),ATP,I,J,K,1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+1),ATS,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+2),ATW,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+3),ATN,I,J,K,-1.0,ID,JD,KD,MG,N)
        E=E+4
        IF(KE .GT. KB) THEN
          CALL GMIJKL(GLOBM, KA(E),ATT,I,J,K,-1.0,ID,JD,KD,MG,N)
	    E=E+1
        ENDIF
        CALL RHSIKL(RHS, KVLST(ROW),BT,I,J,K,1.0,ID,JD,KD,MR,N)
      ENDIF
*
 47   CONTINUE
*      
      IF(JE .GT. JB) THEN
*
      J=JE
      I=IB
      ROW=ROW+1
      CALL GMIJKL(GLOBM, KA(E),ATP,I,J,K,1.0,ID,JD,KD,MG,N)
      CALL GMIJKL(GLOBM, KA(E+1),ATS,I,J,K,-1.0,ID,JD,KD,MG,N)
      E=E+2
      IF(IE .GT. IB) THEN
        CALL GMIJKL(GLOBM, KA(E),ATE,I,J,K,-1.0,ID,JD,KD,MG,N)
        E=E+1
      ENDIF
      IF(KE .GT. KB) THEN
        CALL GMIJKL(GLOBM, KA(E),ATT,I,J,K,-1.0,ID,JD,KD,MG,N)
	  E=E+1
      ENDIF
      CALL RHSIKL(RHS, KVLST(ROW),BT,I,J,K,1.0,ID,JD,KD,MR,N)
*
      DO 48 I=IB+1,IE-1
        ROW=ROW+1
        CALL GMIJKL(GLOBM, KA(E),ATP,I,J,K,1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+1),ATS,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+2),ATW,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+3),ATE,I,J,K,-1.0,ID,JD,KD,MG,N)
 	  E=E+4
        IF(KE .GT. KB) THEN
          CALL GMIJKL(GLOBM, KA(E),ATT,I,J,K,-1.0,ID,JD,KD,MG,N)
	    E=E+1
        ENDIF
        CALL RHSIKL(RHS, KVLST(ROW),BT,I,J,K,1.0,ID,JD,KD,MR,N)
 48   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
        CALL GMIJKL(GLOBM, KA(E),ATP,I,J,K,1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+1),ATS,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+2),ATW,I,J,K,-1.0,ID,JD,KD,MG,N)
	  E=E+3
        IF(KE .GT. KB) THEN
          CALL GMIJKL(GLOBM, KA(E),ATT,I,J,K,-1.0,ID,JD,KD,MG,N)
	    E=E+1
        ENDIF
        CALL RHSIKL(RHS, KVLST(ROW),BT,I,J,K,1.0,ID,JD,KD,MR,N)
      ENDIF
*
      ENDIF
*
*  Interior control-volumes
*
      DO 55 K=KB+1,KE-1
*
      I=IB
      J=JB
      ROW=ROW+1
      CALL GMIJKL(GLOBM, KA(E),ATP,I,J,K,1.0,ID,JD,KD,MG,N)
      CALL GMIJKL(GLOBM, KA(E+1),ATB,I,J,K,-1.0,ID,JD,KD,MG,N)
      E=E+2
      IF(IE .GT. IB) THEN
        CALL GMIJKL(GLOBM, KA(E),ATE,I,J,K,-1.0,ID,JD,KD,MG,N)
	  E=E+1
      ENDIF
      IF(JE .GT. JB) THEN
        CALL GMIJKL(GLOBM, KA(E),ATN,I,J,K,-1.0,ID,JD,KD,MG,N)
	  E=E+1
      ENDIF
      CALL GMIJKL(GLOBM, KA(E),ATT,I,J,K,-1.0,ID,JD,KD,MG,N)
      E=E+1
      CALL RHSIKL(RHS, KVLST(ROW),BT,I,J,K,1.0,ID,JD,KD,MR,N)
*
      DO 49 I=IB+1,IE-1
        ROW=ROW+1
        CALL GMIJKL(GLOBM, KA(E),ATP,I,J,K,1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+1),ATB,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+2),ATW,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+3),ATE,I,J,K,-1.0,ID,JD,KD,MG,N)
 	  E=E+4
        IF(JE .GT. JB) THEN
          CALL GMIJKL(GLOBM, KA(E),ATN,I,J,K,-1.0,ID,JD,KD,MG,N)
	    E=E+1
        ENDIF
        CALL GMIJKL(GLOBM, KA(E),ATT,I,J,K,-1.0,ID,JD,KD,MG,N)
        E=E+1
        CALL RHSIKL(RHS, KVLST(ROW),BT,I,J,K,1.0,ID,JD,KD,MR,N)
 49   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
        CALL GMIJKL(GLOBM, KA(E),ATP,I,J,K,1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+1),ATB,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+2),ATW,I,J,K,-1.0,ID,JD,KD,MG,N)
	  E=E+3
        IF(JE .GT. JB) THEN
          CALL GMIJKL(GLOBM, KA(E),ATN,I,J,K,-1.0,ID,JD,KD,MG,N)
	    E=E+1
        ENDIF
        CALL GMIJKL(GLOBM, KA(E),ATT,I,J,K,-1.0,ID,JD,KD,MG,N)
        E=E+1
        CALL RHSIKL(RHS, KVLST(ROW),BT,I,J,K,1.0,ID,JD,KD,MR,N)
      ENDIF
*
      DO 51 J=JB+1,JE-1
*
      I=IB
      ROW=ROW+1
      CALL GMIJKL(GLOBM, KA(E),ATP,I,J,K,1.0,ID,JD,KD,MG,N)
      CALL GMIJKL(GLOBM, KA(E+1),ATB,I,J,K,-1.0,ID,JD,KD,MG,N)
      CALL GMIJKL(GLOBM, KA(E+2),ATS,I,J,K,-1.0,ID,JD,KD,MG,N)
      E=E+3
      IF(IE .GT. IB) THEN
        CALL GMIJKL(GLOBM, KA(E),ATE,I,J,K,-1.0,ID,JD,KD,MG,N)
	  E=E+1
      ENDIF
      CALL GMIJKL(GLOBM, KA(E),ATN,I,J,K,-1.0,ID,JD,KD,MG,N)
      CALL GMIJKL(GLOBM, KA(E+1),ATT,I,J,K,-1.0,ID,JD,KD,MG,N)
      E=E+2
      CALL RHSIKL(RHS, KVLST(ROW),BT,I,J,K,1.0,ID,JD,KD,MR,N)
*
      DO 50 I=IB+1,IE-1
        ROW=ROW+1
        CALL GMIJKL(GLOBM, KA(E),ATP,I,J,K,1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+1),ATB,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+2),ATS,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+3),ATW,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+4),ATE,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+5),ATN,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+6),ATT,I,J,K,-1.0,ID,JD,KD,MG,N)
        E=E+7
        CALL RHSIKL(RHS, KVLST(ROW),BT,I,J,K,1.0,ID,JD,KD,MR,N)
 50   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
        CALL GMIJKL(GLOBM, KA(E),ATP,I,J,K,1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+1),ATB,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+2),ATS,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+3),ATW,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+4),ATN,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+5),ATT,I,J,K,-1.0,ID,JD,KD,MG,N)
        E=E+6
        CALL RHSIKL(RHS, KVLST(ROW),BT,I,J,K,1.0,ID,JD,KD,MR,N)
      ENDIF
*
 51   CONTINUE
*      
      IF(JE .GT. JB) THEN
*
      J=JE
      I=IB
      ROW=ROW+1
      CALL GMIJKL(GLOBM, KA(E),ATP,I,J,K,1.0,ID,JD,KD,MG,N)
      CALL GMIJKL(GLOBM, KA(E+1),ATB,I,J,K,-1.0,ID,JD,KD,MG,N)
      CALL GMIJKL(GLOBM, KA(E+2),ATS,I,J,K,-1.0,ID,JD,KD,MG,N)
      E=E+3
      IF(IE .GT. IB) THEN
        CALL GMIJKL(GLOBM, KA(E),ATE,I,J,K,-1.0,ID,JD,KD,MG,N)
	  E=E+1
      ENDIF
      CALL GMIJKL(GLOBM, KA(E),ATT,I,J,K,-1.0,ID,JD,KD,MG,N)
      E=E+1
      CALL RHSIKL(RHS, KVLST(ROW),BT,I,J,K,1.0,ID,JD,KD,MR,N)
*
      DO 52 I=IB+1,IE-1
        ROW=ROW+1
        CALL GMIJKL(GLOBM, KA(E),ATP,I,J,K,1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+1),ATB,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+2),ATS,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+3),ATW,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+4),ATE,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+5),ATT,I,J,K,-1.0,ID,JD,KD,MG,N)
        E=E+6
        CALL RHSIKL(RHS, KVLST(ROW),BT,I,J,K,1.0,ID,JD,KD,MR,N)
 52   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
        CALL GMIJKL(GLOBM, KA(E),ATP,I,J,K,1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+1),ATB,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+2),ATS,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+3),ATW,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+4),ATT,I,J,K,-1.0,ID,JD,KD,MG,N)
	  E=E+5
        CALL RHSIKL(RHS, KVLST(ROW),BT,I,J,K,1.0,ID,JD,KD,MR,N)
      ENDIF
*
      ENDIF
*
 55   CONTINUE
*
*  Top face control-volumes
*
      IF(KE .GT. KB) THEN
*
      I=IB
      J=JB
      K=KE
      ROW=ROW+1
      CALL GMIJKL(GLOBM, KA(E),ATP,I,J,K,1.0,ID,JD,KD,MG,N)
      CALL GMIJKL(GLOBM, KA(E+1),ATB,I,J,K,-1.0,ID,JD,KD,MG,N)
      E=E+2
      IF(IE .GT. IB) THEN
        CALL GMIJKL(GLOBM, KA(E),ATE,I,J,K,-1.0,ID,JD,KD,MG,N)
	  E=E+1
      ENDIF
      IF(JE .GT. JB) THEN
        CALL GMIJKL(GLOBM, KA(E),ATN,I,J,K,-1.0,ID,JD,KD,MG,N)
	  E=E+1
      ENDIF
      CALL RHSIKL(RHS, KVLST(ROW),BT,I,J,K,1.0,ID,JD,KD,MR,N)
*
      DO 56 I=IB+1,IE-1
        ROW=ROW+1
        CALL GMIJKL(GLOBM, KA(E),ATP,I,J,K,1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+1),ATB,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+2),ATW,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+3),ATE,I,J,K,-1.0,ID,JD,KD,MG,N)
 	  E=E+4
        IF(JE .GT. JB) THEN
          CALL GMIJKL(GLOBM, KA(E),ATN,I,J,K,-1.0,ID,JD,KD,MG,N)
	    E=E+1
        ENDIF
        CALL RHSIKL(RHS, KVLST(ROW),BT,I,J,K,1.0,ID,JD,KD,MR,N)
 56   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
        CALL GMIJKL(GLOBM, KA(E),ATP,I,J,K,1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+1),ATB,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+2),ATW,I,J,K,-1.0,ID,JD,KD,MG,N)
	  E=E+3
        IF(JE .GT. JB) THEN
          CALL GMIJKL(GLOBM, KA(E),ATN,I,J,K,-1.0,ID,JD,KD,MG,N)
	    E=E+1
        ENDIF
        CALL RHSIKL(RHS, KVLST(ROW),BT,I,J,K,1.0,ID,JD,KD,MR,N)
      ENDIF
*
      DO 58 J=JB+1,JE-1
*
      I=IB
      ROW=ROW+1
      CALL GMIJKL(GLOBM, KA(E),ATP,I,J,K,1.0,ID,JD,KD,MG,N)
      CALL GMIJKL(GLOBM, KA(E+1),ATB,I,J,K,-1.0,ID,JD,KD,MG,N)
      CALL GMIJKL(GLOBM, KA(E+2),ATS,I,J,K,-1.0,ID,JD,KD,MG,N)
      E=E+3
      IF(IE .GT. IB) THEN
        CALL GMIJKL(GLOBM, KA(E),ATE,I,J,K,-1.0,ID,JD,KD,MG,N)
	  E=E+1
      ENDIF
      CALL GMIJKL(GLOBM, KA(E),ATN,I,J,K,-1.0,ID,JD,KD,MG,N)
      E=E+1
      CALL RHSIKL(RHS, KVLST(ROW),BT,I,J,K,1.0,ID,JD,KD,MR,N)
*
      DO 57 I=IB+1,IE-1
        ROW=ROW+1
        CALL GMIJKL(GLOBM, KA(E),ATP,I,J,K,1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+1),ATB,I,J,K,-1.0,ID,JD,KD,MG,N) 
        CALL GMIJKL(GLOBM, KA(E+2),ATS,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+3),ATW,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+4),ATE,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+5),ATN,I,J,K,-1.0,ID,JD,KD,MG,N)
        E=E+6
        CALL RHSIKL(RHS, KVLST(ROW),BT,I,J,K,1.0,ID,JD,KD,MR,N)
 57   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
        CALL GMIJKL(GLOBM, KA(E),ATP,I,J,K,1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+1),ATB,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+2),ATS,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+3),ATW,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+4),ATN,I,J,K,-1.0,ID,JD,KD,MG,N)
        E=E+5
        CALL RHSIKL(RHS, KVLST(ROW),BT,I,J,K,1.0,ID,JD,KD,MR,N)
      ENDIF
*
 58   CONTINUE
*      
      IF(JE .GT. JB) THEN
*
      J=JE
      I=IB
      ROW=ROW+1
      CALL GMIJKL(GLOBM, KA(E),ATP,I,J,K,1.0,ID,JD,KD,MG,N)
      CALL GMIJKL(GLOBM, KA(E+1),ATB,I,J,K,-1.0,ID,JD,KD,MG,N)
      CALL GMIJKL(GLOBM, KA(E+2),ATS,I,J,K,-1.0,ID,JD,KD,MG,N)
      E=E+3
      IF(IE .GT. IB) THEN
        CALL GMIJKL(GLOBM, KA(E),ATE,I,J,K,-1.0,ID,JD,KD,MG,N)
	  E=E+1
      ENDIF
      CALL RHSIKL(RHS, KVLST(ROW),BT,I,J,K,1.0,ID,JD,KD,MR,N)

*
      DO 59 I=IB+1,IE-1
        ROW=ROW+1
        CALL GMIJKL(GLOBM, KA(E),ATP,I,J,K,1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+1),ATB,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+2),ATS,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+3),ATW,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+4),ATE,I,J,K,-1.0,ID,JD,KD,MG,N)
        E=E+5
        CALL RHSIKL(RHS, KVLST(ROW),BT,I,J,K,1.0,ID,JD,KD,MR,N)
 59   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        ROW=ROW+1
        CALL GMIJKL(GLOBM, KA(E),ATP,I,J,K,1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+1),ATB,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+2),ATS,I,J,K,-1.0,ID,JD,KD,MG,N)
        CALL GMIJKL(GLOBM, KA(E+3),ATW,I,J,K,-1.0,ID,JD,KD,MG,N)
        E=E+4
        CALL RHSIKL(RHS, KVLST(ROW),BT,I,J,K,1.0,ID,JD,KD,MR,N)
      ENDIF
*
      ENDIF
*
      ENDIF
*
      RETURN
      END
*
*
*  subroutine GMIJKL
*********************************************************************
*
      SUBROUTINE GMIJKL(GLOBM, KLSTRT,ABLK,I,J,K,SIGN,
     C                  ID,JD,KD,M,N)
*
*     Subroutine GMIJKL inserts the entries of ABLK into the global
*     stiffness matrix starting at KLSTRT.
*
*********************************************************************
*
      IMPLICIT NONE
      DOUBLE PRECISION GLOBM(M)
      REAL*8 ABLK(N,N,ID,JD,KD),DSIGN
      REAL SIGN
      INTEGER KLSTRT,H,I,J,K,L,ID,JD,KD,M,N,CNT
*
      DSIGN = DBLE(SIGN)
      CNT= KLSTRT
      DO 20 H=1,N
       DO 10 L=1,N
          GLOBM(CNT)= DSIGN*ABLK(H,L,I,J,K)
          CNT=CNT+1
  10   CONTINUE
  20  CONTINUE
*
      RETURN
      END
*
*  subroutine RHSIKL
*********************************************************************
*
      SUBROUTINE RHSIKL(RHS, KLSTRT,BVEC,I,J,K,SIGN,ID,JD,KD,M,N)
*
*     Subroutine RHSIKL inserts the entries of BVEC into the right-
*     hand-side vector starting at KLSTRT.
*
*********************************************************************
*
      IMPLICIT NONE
      DOUBLE PRECISION RHS(M)
      REAL*8 BVEC(N,ID,JD,KD),DSIGN
      REAL SIGN
      INTEGER KLSTRT,I,J,K,L,ID,JD,KD,M,N,CNT
*
      DSIGN = DBLE(SIGN)
      CNT= KLSTRT
      DO 10 L=1,N
        RHS(CNT)= DSIGN*BVEC(L,I,J,K)
        CNT=CNT+1
  10  CONTINUE
*
      RETURN
      END
