*      This file contains two subroutines: INITAL and SAVE
*
*************************************************************************
*
      SUBROUTINE INITAL(TF,TS,P,U,V,W,UHE,VHN,WHT,
     C                  T0,P0,U0,V0,W0,
     C                  CVTYPE,IRSTRT,IRSI,
     C                  IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB,
     C                  RHO,COND,VISC,BETA,CONDFE,          !
     C                  RHO0,COND0,VISC0,BETA0,CONDFE0)     ! 
*
*     Subroutine to set initial fields by either fixing values at 
*     0 (IRSTRT=0) or by reading in values from a restart file
*     (IRSTRT=1).
*
*     U0   initial guess for uniform velocity distribution; input
*     V0   initial guess for uniform velocity distribution; input
*     W0   initial guess for uniform velocity distribution; input
*     P0   initial guess for uniform pressure distribution; input
*     T0   initial guess for uniform temperature distribution; input
*
*     INTEGER IRSTRT    restart parameter =0 uniform =1 readin; input
*     INTEGER IRSI      file number for restart data file; input
*     INTEGER IB,IE     first and last interior indices in i; input
*     INTEGER JB,JE     first and last interior indices in j; input
*     INTEGER KB,KE     first and last interior indices in k; input
*     INTEGER ID,JD,KD  array dimensions; input
*
*************************************************************************
*
      IMPLICIT NONE
      REAL*8 U(ID,JD,KD),V(ID,JD,KD),W(ID,JD,KD),P(ID,JD,KD)
      REAL*8 TF(ID,JD,KD),TS(ID,JD,KD)
      REAL*8 UHE(ID,JD,KD),VHN(ID,JD,KD),WHT(ID,JD,KD)
      REAL*8 U0,V0,W0,P0,T0
      INTEGER IRSTRT,IRSI,CVTYPE(ID,JD,KD,NNB+1),IB,IE,JB,JE,KB,KE
      INTEGER I,J,K,ID,JD,KD,NNB
      INTEGER IBM1,IEP1,JBM1,JEP1,KBM1,KEP1
*
      REAL*8 RHO(ID,JD,KD),COND(ID,JD,KD),VISC(ID,JD,KD),BETA(ID,JD,KD)
      REAL*8 CONDFE(ID,JD,KD)
      REAL*8 RHO0,COND0,VISC0,BETA0,CONDFE0
*     
*
*     Initialize fluid air properties arrays to inputs in in.dat file
*
      DO 9500 I=IB-1,IE+1
      DO 9400 J=JB-1,JE+1
      DO 9300 K=KB-1,KE+1
       RHO(I,J,K)    = RHO0
       COND(I,J,K)   = COND0
       VISC(I,J,K)   = VISC0
       BETA(I,J,K)   = BETA0
       CONDFE(I,J,K) = CONDFE0
 9300 CONTINUE
 9400 CONTINUE
 9500 CONTINUE
*
      IBM1=IB-1
      IEP1=IE+1
      JBM1=JB-1
      JEP1=JE+1
      KBM1=KB-1
      KEP1=KE+1
*
      IF (IRSTRT .EQ. 1) GO TO 100
*
      DO 50 K=KB,KE
       DO 40 J=JB,JE
        DO 30 I=IB,IE
*        
*     Conditions in solid CVs
*
         IF(CVTYPE(I,J,K,1).EQ.3) THEN
          TS(I,J,K)=  T0
          TF(I,J,K)= 0.0
          P(I,J,K)=  0.0
          U(I,J,K)=  0.0
          V(I,J,K)=  0.0
          W(I,J,K)=  0.0          
          UHE(I,J,K)=0.0
          VHN(I,J,K)=0.0
          WHT(I,J,K)=0.0          
*
*     Conditions in clear fluid and porous CVs
*
         ELSE
          IF(CVTYPE(I,J,K,1).EQ.0) THEN
           TS(I,J,K)= 0.0
          ELSE
           TS(I,J,K)= T0
          ENDIF
          TF(I,J,K)= T0
          P(I,J,K)=  P0
          U(I,J,K)=  U0
          V(I,J,K)=  V0
          W(I,J,K)=  W0          
          IF(CVTYPE(I,J,K,3).EQ.3) THEN
           UHE(I,J,K)=0.0
          ELSE
           UHE(I,J,K)=U0
          ENDIF
          IF(CVTYPE(I,J,K,5).EQ.3) THEN
           VHN(I,J,K)=0.0
          ELSE
           VHN(I,J,K)=V0
          ENDIF
          IF(CVTYPE(I,J,K,7).EQ.3) THEN
           WHT(I,J,K)=0.0
          ELSE
           WHT(I,J,K)=W0
          ENDIF
         ENDIF
*
 30     CONTINUE
 40    CONTINUE
 50   CONTINUE
*
*     Conditions at south and north boundary nodes
*
      DO 60 K=KBM1,KEP1
       DO 55 I=IBM1,IEP1
        J=JBM1
        IF(CVTYPE(I,J,K,5).EQ.3) THEN
         TS(I,J,K)= T0       
         TF(I,J,K)= 0.0
         P(I,J,K)=  0.0
         U(I,J,K)=  0.0
         V(I,J,K)=  0.0
         W(I,J,K)=  0.0         
         UHE(I,J,K)=0.0
         VHN(I,J,K)=0.0
         WHT(I,J,K)=0.0         
        ELSE
         IF(CVTYPE(I,J,K,5).EQ.0) THEN
          TS(I,J,K)= 0.0
         ELSE
          TS(I,J,K)= T0
         ENDIF
         TF(I,J,K)= T0
         P(I,J,K)=  P0
         U(I,J,K)=  U0
         V(I,J,K)=  V0
         W(I,J,K)=  W0         
         UHE(I,J,K)=U0
         VHN(I,J,K)=V0
         WHT(I,J,K)=W0
        ENDIF
* 
        J=JEP1
        IF(CVTYPE(I,J,K,4).EQ.3) THEN
         TS(I,J,K)= T0
         TF(I,J,K)= 0.0
         P(I,J,K)=  0.0
         U(I,J,K)=  0.0
         V(I,J,K)=  0.0
         W(I,J,K)=  0.0         
         UHE(I,J,K)=0.0
         VHN(I,J,K)=0.0
         WHT(I,J,K)=0.0
        ELSE
         IF(CVTYPE(I,J,K,4).EQ.0) THEN
          TS(I,J,K)= 0.0
         ELSE
          TS(I,J,K)= T0
         ENDIF
         TF(I,J,K)= T0
         P(I,J,K)=  P0
         U(I,J,K)=  U0
         V(I,J,K)=  V0
         W(I,J,K)=  W0         
         UHE(I,J,K)=U0
         VHN(I,J,K)=V0
         WHT(I,J,K)=W0
        ENDIF
 55    CONTINUE       
 60   CONTINUE
*
*     Conditions at west and east boundary nodes
*
      DO 70 K=KBM1,KEP1
       DO 65 J=JBM1,JEP1
        I=IBM1
        IF(CVTYPE(I,J,K,3).EQ.3) THEN
         TS(I,J,K)= T0       
         TF(I,J,K)= 0.0
         P(I,J,K)=  0.0
         U(I,J,K)=  0.0
         V(I,J,K)=  0.0
         W(I,J,K)=  0.0         
         UHE(I,J,K)=0.0
         VHN(I,J,K)=0.0
         WHT(I,J,K)=0.0         
        ELSE
         IF(CVTYPE(I,J,K,3).EQ.0) THEN
          TS(I,J,K)= 0.0
         ELSE
          TS(I,J,K)= T0
         ENDIF
         TF(I,J,K)= T0
         P(I,J,K)=  P0
         U(I,J,K)=  U0
         V(I,J,K)=  V0
         W(I,J,K)=  W0         
         UHE(I,J,K)=U0
         VHN(I,J,K)=V0
         WHT(I,J,K)=W0
        ENDIF
* 
        I=IEP1
        IF(CVTYPE(I,J,K,2).EQ.3) THEN
         TS(I,J,K)= T0
         TF(I,J,K)= 0.0
         P(I,J,K)=  0.0
         U(I,J,K)=  0.0
         V(I,J,K)=  0.0
         W(I,J,K)=  0.0         
         UHE(I,J,K)=0.0
         VHN(I,J,K)=0.0
         WHT(I,J,K)=0.0
        ELSE
         IF(CVTYPE(I,J,K,2).EQ.0) THEN
          TS(I,J,K)= 0.0
         ELSE
          TS(I,J,K)= T0
         ENDIF
         TF(I,J,K)= T0
         P(I,J,K)=  P0
         U(I,J,K)=  U0
         V(I,J,K)=  V0
         W(I,J,K)=  W0         
         UHE(I,J,K)=U0
         VHN(I,J,K)=V0
         WHT(I,J,K)=W0
        ENDIF
 65    CONTINUE       
 70   CONTINUE
*
*     Conditions at bottom and top boundary nodes
*
      DO 80 J=JBM1,JEP1
       DO 75 I=IBM1,IEP1
        K=KBM1
        IF(CVTYPE(I,J,K,7).EQ.3) THEN
         TS(I,J,K)= T0       
         TF(I,J,K)= 0.0
         P(I,J,K)=  0.0
         U(I,J,K)=  0.0
         V(I,J,K)=  0.0
         W(I,J,K)=  0.0         
         UHE(I,J,K)=0.0
         VHN(I,J,K)=0.0
         WHT(I,J,K)=0.0         
        ELSE
         IF(CVTYPE(I,J,K,7).EQ.0) THEN
          TS(I,J,K)= 0.0
         ELSE
          TS(I,J,K)= T0
         ENDIF
         TF(I,J,K)= T0
         P(I,J,K)=  P0
         U(I,J,K)=  U0
         V(I,J,K)=  V0
         W(I,J,K)=  W0         
         UHE(I,J,K)=U0
         VHN(I,J,K)=V0
         WHT(I,J,K)=W0
        ENDIF
* 
        K=KEP1
        IF(CVTYPE(I,J,K,6).EQ.3) THEN
         TS(I,J,K)= T0
         TF(I,J,K)= 0.0
         P(I,J,K)=  0.0
         U(I,J,K)=  0.0
         V(I,J,K)=  0.0
         W(I,J,K)=  0.0         
         UHE(I,J,K)=0.0
         VHN(I,J,K)=0.0
         WHT(I,J,K)=0.0
        ELSE
         IF(CVTYPE(I,J,K,6).EQ.0) THEN
          TS(I,J,K)= 0.0
         ELSE
          TS(I,J,K)= T0
         ENDIF
         TF(I,J,K)= T0
         P(I,J,K)=  P0
         U(I,J,K)=  U0
         V(I,J,K)=  V0
         W(I,J,K)=  W0         
         UHE(I,J,K)=U0
         VHN(I,J,K)=V0
         WHT(I,J,K)=W0
        ENDIF
 75    CONTINUE       
 80   CONTINUE
      GOTO 200
*
 100  CONTINUE
      READ(IRSI) (((TF(I,J,K),I=IBM1,IEP1),J=JBM1,JEP1),K=KBM1,KEP1)
      READ(IRSI) (((TS(I,J,K),I=IBM1,IEP1),J=JBM1,JEP1),K=KBM1,KEP1)      
      READ(IRSI) (((P(I,J,K),I=IBM1,IEP1),J=JBM1,JEP1),K=KBM1,KEP1)
      READ(IRSI) (((U(I,J,K),I=IBM1,IEP1),J=JBM1,JEP1),K=KBM1,KEP1)
      READ(IRSI) (((V(I,J,K),I=IBM1,IEP1),J=JBM1,JEP1),K=KBM1,KEP1)
      READ(IRSI) (((W(I,J,K),I=IBM1,IEP1),J=JBM1,JEP1),K=KBM1,KEP1)
      READ(IRSI) (((UHE(I,J,K),I=IBM1,IE),J=JB,JE),K=KB,KE)
      READ(IRSI) (((VHN(I,J,K),I=IB,IE),J=JBM1,JE),K=KB,KE)
      READ(IRSI) (((WHT(I,J,K),I=IB,IE),J=JB,JE),K=KBM1,KE)
*
 200  CONTINUE
      RETURN
      END
*
*
************************************************************************
*
      SUBROUTINE SAVE(TF,TS,P,U,V,W,UHE,VHN,WHT,XP,YP,ZP,
     C                SGENT,SGEN,SGENF,SGENS,RATIOT,RATIOS,RATIO,
     C                AREP,ME,CP,DIEP,DISE,
     C                DTF,DTS,DNF,DNS,DEF,DES,DTMX,
     C                FNAPP,IRSO,IB,IE,JB,JE,KB,KE,ID,JD,KD,
     C                CVTYPE,KNTOUT,INDIC,IPB,IPE,NNB)
*
*     Routine to save solution arrays.
*
*************************************************************************
*
      IMPLICIT NONE
      REAL*8 TF(ID,JD,KD),TS(ID,JD,KD)
      REAL*8 P(ID,JD,KD),U(ID,JD,KD),V(ID,JD,KD),W(ID,JD,KD)
      REAL*8 ULID,UHE(ID,JD,KD),VHN(ID,JD,KD),WHT(ID,JD,KD)
      REAL*8 XP(ID),YP(JD),ZP(KD)
      REAL*8 SGENT,SGEN(ID,JD,KD),SGENF(ID,JD,KD),SGENS(ID,JD,KD)
      REAL*8 RATIOT,RATIOS,RATIO(ID,JD,KD)
      REAL*8 AREP(JD,KD),ME(ID,JD,KD)
      REAL*8 DIEP(ID),DISE(ID),CP,DTF(ID,JD,KD),DTS(ID,JD,KD)
      REAL*8 DNF(ID,JD,KD),DNS(ID,JD,KD)
      REAL*8 DEF(ID,JD,KD),DES(ID,JD,KD)
      REAL*8 PFIN,PFOUT,AIN,AOUT,PIN,POUT,DELTAP
      REAL*8 EIN,EOUT,QTOT,MTOT,TIN,TOUT,LTNE,DTMX,QBACK,QF,QFBYQT
      PARAMETER(ULID = 1.0)
      INTEGER IRSO,IB,IE,JB,JE,KB,KE,ID,JD,KD,FNUM,KNUM
      INTEGER KNTOUT,INDIC,IPB,IPE,NNB,CVTYPE(ID,JD,KD,NNB+1)
      INTEGER IBM1,IEP1,JBM1,JEP1,KBM1,KEP1,I,J,K
      CHARACTER*5 CNAME,FNAPP
      CHARACTER*6 COL
      CHARACTER*8 CFNUM
      CHARACTER*17 FNAME
*
      IBM1 = IB - 1
      IEP1 = IE + 1
      JBM1 = JB - 1
      JEP1 = JE + 1
      KBM1 = KB - 1
      KEP1 = KE + 1
*
      IF (INDIC .EQ. 1) THEN
*
       OPEN(UNIT=20,FILE='tec'//FNAPP//'.dat')
       OPEN(UNIT=21,FILE='matlab'//FNAPP//'.dat')
       OPEN(UNIT=22,FILE='xymidu'//FNAPP//'.dat')
       OPEN(UNIT=23,FILE='yzmidw'//FNAPP//'.dat')   
       OPEN(UNIT=24,FILE='midplane'//FNAPP//'.dat')           
       OPEN(UNIT=25,FILE='overall'//FNAPP//'.dat')       
*     
       WRITE(IRSO) (((TF(I,J,K),I=IBM1,IEP1),J=JBM1,JEP1),K=KBM1,KEP1)
       WRITE(IRSO) (((TS(I,J,K),I=IBM1,IEP1),J=JBM1,JEP1),K=KBM1,KEP1)       
       WRITE(IRSO) (((P(I,J,K),I=IBM1,IEP1),J=JBM1,JEP1),K=KBM1,KEP1)
       WRITE(IRSO) (((U(I,J,K),I=IBM1,IEP1),J=JBM1,JEP1),K=KBM1,KEP1)
       WRITE(IRSO) (((V(I,J,K),I=IBM1,IEP1),J=JBM1,JEP1),K=KBM1,KEP1)
       WRITE(IRSO) (((W(I,J,K),I=IBM1,IEP1),J=JBM1,JEP1),K=KBM1,KEP1)
       WRITE(IRSO) (((UHE(I,J,K),I=IBM1,IE),J=JB,JE),K=KB,KE)
       WRITE(IRSO) (((VHN(I,J,K),I=IB,IE),J=JBM1,JE),K=KB,KE)
       WRITE(IRSO) (((WHT(I,J,K),I=IB,IE),J=JB,JE),K=KBM1,KE)
*
*  Write output in TECplot format
*
       WRITE(20,320)
       WRITE(20,322) 0,0,1,(IE-IB+3),(JE-JB+3),(KE-KB+3)
       DO 30 K=KBM1,KEP1
        DO 25 J=JBM1,JEP1
         DO 10 I=IBM1,IEP1
*
*          Calculate percent local thermal non-equilibrium.
*         
           IF(CVTYPE(I,J,K,1).EQ.2) THEN
            LTNE = 100*ABS(TF(I,J,K)-TS(I,J,K))/DTMX
           ELSEIF(CVTYPE(I,J,K,1).EQ.1 .AND. (CVTYPE(I,J,K,2).EQ.2 
     C      .OR. CVTYPE(I,J,K,3).EQ.2 .OR. CVTYPE(I,J,K,4).EQ.2
     C      .OR. CVTYPE(I,J,K,5).EQ.2 .OR. CVTYPE(I,J,K,6).EQ.2
     C      .OR. CVTYPE(I,J,K,7).EQ.2)) THEN
            LTNE = 100*ABS(TF(I,J,K)-TS(I,J,K))/DTMX            
           ELSE
            LTNE = 0
           ENDIF
           WRITE(20,325) XP(I),YP(J),ZP(K),P(I,J,K),U(I,J,K)/ULID,
     C                   V(I,J,K)/ULID,W(I,J,K)/ULID,
     C                   TF(I,J,K),TS(I,J,K),LTNE,SGENF(I,J,K),
     C                   SGENS(I,J,K),SGEN(I,J,K),RATIO(I,J,K)
  10     CONTINUE
  25    CONTINUE
  30   CONTINUE
*
*  Write output for Matlab postprocessing
*
       DO 32 K=KBM1,KEP1
        KNUM = 100 + K
        WRITE(COL,*) KNUM
        CNAME = 'u' // COL(4:6) // '='
        WRITE(21,*) CNAME
        DO 31 J = JBM1,JEP1
         WRITE(21,326) U(IEP1,J,K)
  31    CONTINUE
  32   CONTINUE
*
*  Write output for Excel postprocessing
*
       DO 33 K=KBM1,KEP1
         WRITE(22,327) ZP(K),U(IE/2+1,JE/2+1,K)/ULID
  33   CONTINUE
       DO 34 I=IBM1,IEP1
         WRITE(23,327) XP(I),W(I,JE/2+1,KE/2+1)/ULID
  34   CONTINUE
       DO 36 I=IBM1,IEP1
        DO 35 K=KBM1,KEP1
         J = JB
         WRITE(24,328) XP(I),ZP(K),U(I,J,K),W(I,J,K)
  35    CONTINUE
  36   CONTINUE  
*
*  Calculate and output overall heat transfer and pressure drop,
*  as well as back conduction and fraction of heat transfer directly
*  into fluid constituent
*
       I=IBM1
       AIN = 0.0
       PFIN = 0.0
       DO 38 J=JB,JE
        DO 37 K=KB,KE
         IF(CVTYPE(I,J,K,3).NE.3) THEN
          AIN = AIN + AREP(J,K)
          PFIN = PFIN + P(I,J,K)*AREP(J,K)
         ENDIF
  37    CONTINUE
  38   CONTINUE
*
       PIN = PFIN/AIN   
       DELTAP = PIN - 0.0        
*       
       K = KBM1
       QTOT = 0.0
       QF = 0.0
       DO 40 J=JB,JE
        DO 39 I=IB,IE
         QTOT = QTOT + DTF(I,J,K)*(TF(I,J,K)-TF(I,J,K+1))
     C          + DTS(I,J,K)*(TS(I,J,K)-TS(I,J,K+1))
         QF = QF + DTF(I,J,K)*(TF(I,J,K)-TF(I,J,K+1))
  39    CONTINUE
  40   CONTINUE   
       QFBYQT = QF/QTOT    
*     
       I = IBM1
       QBACK = 0.0
       DO 42 J=JB,JE
        DO 41 K=KB,KE
         QBACK = QBACK + DEF(I,J,K)*(TF(I+1,J,K)-TF(I,J,K))
     C          + DES(I,J,K)*(TS(I+1,J,K)-TS(I,J,K))
  41    CONTINUE
  42   CONTINUE       
*     
       WRITE(25,331) DELTAP
       WRITE(25,335) QTOT                     
       WRITE(25,336) QBACK       
       WRITE(25,337) QFBYQT  
       WRITE(25,338) SGENT       
       WRITE(25,339) RATIOT 
       WRITE(25,340) RATIOS
*
      ELSE
*
       FNUM = KNTOUT+10000
       WRITE(CFNUM,*) FNUM
       FNAME = 'tec'//CFNUM(4:8)//FNAPP//'.dat'
       OPEN(UNIT=FNUM,FILE=FNAME)
*
       WRITE(FNUM,320)
       WRITE(FNUM,322) 0,0,1,(IE-IB+3),(JE-JB+3),(KE-KB+3)
       DO 45 K=KBM1,KEP1
        DO 44 J=JBM1,JEP1
         DO 43 I=IBM1,IEP1
           IF(CVTYPE(I,J,K,1).EQ.2) THEN
            LTNE = 100*ABS(TF(I,J,K)-TS(I,J,K))/DTMX
           ELSEIF(CVTYPE(I,J,K,1).EQ.1 .AND. (CVTYPE(I,J,K,2).EQ.2
     C      .OR. CVTYPE(I,J,K,3).EQ.2 .OR. CVTYPE(I,J,K,4).EQ.2
     C      .OR. CVTYPE(I,J,K,5).EQ.2 .OR. CVTYPE(I,J,K,6).EQ.2
     C      .OR. CVTYPE(I,J,K,7).EQ.2)) THEN
            LTNE = 100*ABS(TF(I,J,K)-TS(I,J,K))/DTMX            
           ELSE
            LTNE = 0
           ENDIF         
           WRITE(FNUM,325) XP(I),YP(J),ZP(K),P(I,J,K),U(I,J,K)/ULID,
     C                     V(I,J,K)/ULID,W(I,J,K)/ULID,
     C                     TF(I,J,K),TS(I,J,K),LTNE,SGENF(I,J,K),
     C                     SGENS(I,J,K),SGEN(I,J,K),RATIO(I,J,K)
  43     CONTINUE
  44    CONTINUE
  45   CONTINUE
      ENDIF
*       
 320  FORMAT('VARIABLES="XP","YP","ZP","P","U","V","W","TF","TS","LTNE"
     C,"SGENF","SGENS","SGEN","RATIO"')
 322  FORMAT('ZONE T="Zone ',I1,I1,I1,'" I=',I6,' J=',I6,' K=',I6,'
     C F=POINT')
 325  FORMAT(14(2X,1PE12.5))
 326  FORMAT(1(2X,1PE12.5))
 327  FORMAT(2(2X,1PE12.5))
 328  FORMAT(4(2X,1PE12.5)) 
 329  FORMAT('Pin = ',1PD13.5,' Pa')
 330  FORMAT('Pout = ',1PD13.5,' Pa') 
 331  FORMAT('DeltaP = ',1PD13.5,' Pa') 
 332  FORMAT('Tin = ',1PD13.5,' C')
 333  FORMAT('Tout = ',1PD13.5,' C') 
 334  FORMAT('Mtot = ',1PD13.5,' kg/s') 
 335  FORMAT('Qtot = ',1PD13.5,' W') 
 336  FORMAT('Qback = ',1PD13.5,' W')
 337  FORMAT('Qf/Qtot = ',1PD13.5)  
 338  FORMAT('Sgen = ',1PD13.5,' W/K') 
 339  FORMAT('Be = ',1PD13.5)
 340  FORMAT('Sgen,s/Sgen = ',1PD13.5)
*
       CALL FLUSH(20)
       CLOSE(UNIT=20)
       CALL FLUSH(FNUM)
       CLOSE(UNIT=FNUM)
      RETURN
      END
