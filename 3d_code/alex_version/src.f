*
*     This file contains 5 subroutines: SRCTF,SCRTS,SRCU,SRCV and SRCW.
*
************************************************************************
      SUBROUTINE SRCTF(QT,RT, DCCE,DCCN,DCCT,HSF,SPECSA,VOLP,
     C                 CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB,TTIME,
     C                 HEATERTIME,HEATER,
     C                 DTIME,GRDY)
*
*     Subroutine to calculate the net source of T in each interior
*     control volume for the entire volume. Fluid-phase.
*     Net Source = Q + R*T
*
*     Q(ID,JD,KD) fixed source coefficient; output
*     R(ID,JD,KD) linearized source coefficient; output
*
************************************************************************
*
      IMPLICIT NONE
      REAL*8 QT(ID,JD,KD),RT(ID,JD,KD)
      REAL*8 DCCE(ID,JD,KD),DCCN(ID,JD,KD),DCCT(ID,JD,KD)
      REAL*8 HSF(ID,JD,KD),SPECSA(ID,JD,KD),VOLP(ID,JD,KD)
      REAL*8 TTIME
      INTEGER IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB,CVTYPE(ID,JD,KD,NNB+1)
      INTEGER I,J,K,HEATER
      REAL*8 INTGEN,HEATERTIME,DTIME,GRDY
*
      DO 30 K=KB,KE
       DO 20 J=JB,JE
        DO 10 I=IB,IE
*
*        Set where the internal source acts - CHANGE WITH NUMBER OF Y CONTROL VOLUMES
*
         IF ((J == HEATER).AND.(TTIME<=HEATERTIME)) THEN
          !INTGEN IS PER CONTROL VOLUME - divided by total I and K CVs
          INTGEN = 25000*DTIME/(GRDY/(JE-JB+2))/(IE-IB+2)/(KE-KB+2) ! Zanoni et al Table 6
         ELSE 
          INTGEN = 0.0
         END IF
*        
*        Deferred corrections
*
         QT(I,J,K) = -1.0*DCCE(I,J,K)+DCCE(I-1,J,K)
     C               -DCCN(I,J,K)+DCCN(I,J-1,K)
     C               -DCCT(I,J,K)+DCCT(I,J,K-1)
     C               +INTGEN*VOLP(I,J,K)
*
*        Interfacial exchange
* 
         IF(CVTYPE(I,J,K,1).EQ.2) THEN
           RT(I,J,K) = -1.0*HSF(I,J,K)*SPECSA(I,J,K)*VOLP(I,J,K)
         ELSE
           RT(I,J,K) = 0.0
         ENDIF
 10     CONTINUE
 20    CONTINUE
 30   CONTINUE
      RETURN
      END
*
************************************************************************
      SUBROUTINE SRCTS(QT,RT, HSF,SPECSA,VOLP,
     C                 CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB,TTIME,
     C                 HEATERTIME,HEATER,
     C                 DTIME,GRDY)
*
*     Subroutine to calculate the net source of T in each interior
*     control volume for the entire volume.Solid-phase.
*     Net Source = Q + R*T 
*
*     Q(ID,JD,KD) fixed source coefficient; output
*     R(ID,JD,KD) linearized source coefficient; output
*
************************************************************************
*
      IMPLICIT NONE
      REAL*8 QT(ID,JD,KD),RT(ID,JD,KD)
      REAL*8 HSF(ID,JD,KD),SPECSA(ID,JD,KD),VOLP(ID,JD,KD)
      REAL*8 TTIME
      INTEGER IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB,CVTYPE(ID,JD,KD,NNB+1)
      INTEGER I,J,K,HEATER
      REAL*8 INTGEN,HEATERTIME,DTIME,GRDY
*
      DO 30 K=KB,KE
       DO 20 J=JB,JE
        DO 10 I=IB,IE
*
*        Set where the internal source acts - CHANGE WITH NUMBER OF Y CONTROL VOLUMES
*
         IF ((J == HEATER).AND.(TTIME<=HEATERTIME)) THEN
          INTGEN = 25000*DTIME/(GRDY/(JE-JB+2))/(IE-IB+2)/(KE-KB+2) !750000.0 !39603.95 ! Zanoni et al Table 6
         ELSE 
          INTGEN = 0.0
         END IF
*
         QT(I,J,K) = 0.0 + INTGEN*VOLP(I,J,K)
*
*        Interfacial exchange
*
         IF(CVTYPE(I,J,K,1).EQ.2) THEN
           RT(I,J,K) = -1.0*HSF(I,J,K)*SPECSA(I,J,K)*VOLP(I,J,K)
         ELSE
           RT(I,J,K) = 0.0
         ENDIF
 10     CONTINUE
 20    CONTINUE
 30   CONTINUE
      RETURN
      END
*
************************************************************************
      SUBROUTINE SRCU(QU,RU, DCCE,DCCN,DCCT,VOLP,
     C                VISC,KPERM,RHO,CFORCH,PRSTY,U,V,W,
     C                CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
*
*     Subroutine to calculate the net source of U in each interior
*     control volume for the entire volume.
*     Net Source = Q + R*U 
*
*     Q(ID,JD,KD) fixed source coefficient; output
*     R(ID,JD,KD) linearized source coefficient; output
*
************************************************************************
*
      IMPLICIT NONE
      REAL*8 QU(ID,JD,KD),RU(ID,JD,KD)
      REAL*8 DCCE(ID,JD,KD),DCCN(ID,JD,KD),DCCT(ID,JD,KD)
      REAL*8 VOLP(ID,JD,KD),U(ID,JD,KD),V(ID,JD,KD),W(ID,JD,KD)
      REAL*8 KPERM(ID,JD,KD),CFORCH(ID,JD,KD),PRSTY(ID,JD,KD)
      REAL*8 VEL
      INTEGER IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB,CVTYPE(ID,JD,KD,NNB+1)
      INTEGER I,J,K
*
      REAL*8 VISC(ID,JD,KD),RHO(ID,JD,KD)
*   
      DO 30 K=KB,KE  
       DO 20 J=JB,JE
        DO 10 I=IB,IE
*
*        Deferred corrections
*                
         QU(I,J,K) = -1.0*DCCE(I,J,K)+DCCE(I-1,J,K)
     C               -DCCN(I,J,K)+DCCN(I,J-1,K)
     C               -DCCT(I,J,K)+DCCT(I,J,K-1)     
*
*        Darcy and Forchheimer terms
*                  
         IF(CVTYPE(I,J,K,1).EQ.2) THEN
          VEL = (U(I,J,K)**2+V(I,J,K)**2+W(I,J,K)**2)**0.5
          RU(I,J,K) = PRSTY(I,J,K)*VOLP(I,J,K)
     C                *(-1.0*VISC(I,J,K)/KPERM(I,J,K)-RHO(I,J,K)
     C                *CFORCH(I,J,K)*VEL/((KPERM(I,J,K))**0.5))
         ELSE
          RU(I,J,K) = 0.0
         ENDIF
 10     CONTINUE
 20    CONTINUE
 30   CONTINUE
      RETURN
      END
*
************************************************************************
      SUBROUTINE SRCV(QV,RV, DCCE,DCCN,DCCT,VOLP,
     C                VISC,KPERM,RHO,CFORCH,PRSTY,U,V,W,
     C                CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB,
     C                TF,BETA)
*
*     Subroutine to calculate the net source of V in each interior
*     control volume for the entire volume.
*     Net Source = Q + R*V 
*
*     Q(ID,JD,KD) fixed source coefficient; output
*     R(ID,JD,KD) linearized source coefficient; output
*
************************************************************************
*
      IMPLICIT NONE
      REAL*8 QV(ID,JD,KD),RV(ID,JD,KD)
      REAL*8 DCCE(ID,JD,KD),DCCN(ID,JD,KD),DCCT(ID,JD,KD)
      REAL*8 VOLP(ID,JD,KD),U(ID,JD,KD),V(ID,JD,KD),W(ID,JD,KD)
      REAL*8 KPERM(ID,JD,KD),CFORCH(ID,JD,KD),PRSTY(ID,JD,KD)
      REAL*8 VEL      
      INTEGER IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB,CVTYPE(ID,JD,KD,NNB+1)
      INTEGER I,J,K
*
      REAL*8 VISC(ID,JD,KD),RHO(ID,JD,KD),BETA(ID,JD,KD)
      REAL*8 TF(ID,JD,KD),TREF,GEY
      PARAMETER(TREF=293.0)
      PARAMETER(GEY=-9.81)
*
      DO 30 K=KB,KE  
       DO 20 J=JB,JE
        DO 10 I=IB,IE
*
*        Deferred corrections
*
         QV(I,J,K) = -1.0*DCCE(I,J,K)+DCCE(I-1,J,K)
     C               -DCCN(I,J,K)+DCCN(I,J-1,K)
     C               -DCCT(I,J,K)+DCCT(I,J,K-1)
     C          -RHO(I,J,K)*BETA(I,J,K)*GEY*(TF(I,J,K)-TREF)*VOLP(I,J,K)   !Buoyant source term
*
*        Darcy and Forchheimer terms
*
         IF(CVTYPE(I,J,K,1).EQ.2) THEN
          VEL = (U(I,J,K)**2+V(I,J,K)**2+W(I,J,K)**2)**0.5
          RV(I,J,K) = PRSTY(I,J,K)*VOLP(I,J,K)
     C                *(-1.0*VISC(I,J,K)/KPERM(I,J,K)-RHO(I,J,K)
     C                *CFORCH(I,J,K)*VEL/((KPERM(I,J,K))**0.5))
         ELSE
          RV(I,J,K) = 0.0
         ENDIF
 10     CONTINUE
 20    CONTINUE
 30   CONTINUE
      RETURN
      END
*
************************************************************************
      SUBROUTINE SRCW(QW,RW, DCCE,DCCN,DCCT,
     C                TF,RHO,VOLP,GEE,BETA,
     C                VISC,KPERM,CFORCH,PRSTY,U,V,W,
     C                CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
*
*     Subroutine to calculate the net source of W in each interior
*     control volume for the entire volume.
*     Net Source = Q + R*W 
*
*     Q(ID,JD,KD) fixed source coefficient; output
*     R(ID,JD,KD) linearized source coefficient; output
*
************************************************************************
*
      IMPLICIT NONE
      REAL*8 QW(ID,JD,KD),RW(ID,JD,KD)
      REAL*8 DCCE(ID,JD,KD),DCCN(ID,JD,KD),DCCT(ID,JD,KD)
      REAL*8 TF(ID,JD,KD),VOLP(ID,JD,KD),GEE,BETA
      REAL*8 U(ID,JD,KD),V(ID,JD,KD),W(ID,JD,KD)
      REAL*8 KPERM(ID,JD,KD),CFORCH(ID,JD,KD),PRSTY(ID,JD,KD)
      REAL*8 VEL
      INTEGER IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB,CVTYPE(ID,JD,KD,NNB+1)
      INTEGER I,J,K
*
      REAL*8 RHO(ID,JD,KD),VISC(ID,JD,KD)
*
      DO 30 K=KB,KE  
       DO 20 J=JB,JE
        DO 10 I=IB,IE
*
*        Deferred corrections
*        
         QW(I,J,K) = -1.0*DCCE(I,J,K)+DCCE(I-1,J,K)
     C               -DCCN(I,J,K)+DCCN(I,J-1,K)
     C               -DCCT(I,J,K)+DCCT(I,J,K-1)
*
*        Darcy and Forchheimer terms
*
         IF(CVTYPE(I,J,K,1).EQ.2) THEN
          VEL = (U(I,J,K)**2+V(I,J,K)**2+W(I,J,K)**2)**0.5
          RW(I,J,K) = PRSTY(I,J,K)*VOLP(I,J,K)
     C                *(-1.0*VISC(I,J,K)/KPERM(I,J,K)-RHO(I,J,K)
     C                *CFORCH(I,J,K)*VEL/((KPERM(I,J,K))**0.5))      
         ELSE
          RW(I,J,K) = 0.0
         ENDIF
 10     CONTINUE
 20    CONTINUE
 30   CONTINUE
      RETURN
      END