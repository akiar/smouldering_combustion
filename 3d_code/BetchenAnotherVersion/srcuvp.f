*
*         file srcuvp.f
*********************************************************************
*
      SUBROUTINE SRCUVP(AP,AW,AE,AS,AN,AB,AT,B,
     C                  P,PEAST,PNORTH,PTOP,PRSTY,   
     C                  VOLP,DIEP,DISE,DJNP,DISN,DKTP,DIST,
     C                  CVTYPE,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB)
*
*     Subroutine to add the pressure source term in the U & V
*     momentum equations.  The order of the solver is P,U,V.
*
*     Note: As per thesis, at porous/pure fluid interfaces standard
*           CDS estimate for gradient is used implicitly, with
*           interface pressure estimate accounted for by deferred 
*           correction.
*
*********************************************************************
*
      IMPLICIT NONE
      REAL*8 AP(N,N,ID,JD,KD),AW(N,N,ID,JD,KD),AE(N,N,ID,JD,KD)
      REAL*8 AS(N,N,ID,JD,KD),AN(N,N,ID,JD,KD),AB(N,N,ID,JD,KD)
      REAL*8 AT(N,N,ID,JD,KD),B(N,ID,JD,KD)
      REAL*8 P(ID,JD,KD),PRSTY(ID,JD,KD),VOLP(ID,JD,KD)
      REAL*8 PEAST(ID,JD,KD),PNORTH(ID,JD,KD),PTOP(ID,JD,KD)
      REAL*8 DIEP(ID),DISE(ID),DJNP(JD),DISN(JD),DKTP(KD),DIST(KD)
      INTEGER IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB
      INTEGER CVTYPE(ID,JD,KD,NNB+1),I,J,K
*  
      DO 20 I=IB,IE
       DO 15 J=JB,JE
        DO 10 K=KB,KE
*        
        IF(CVTYPE(I,J,K,1).NE.3) THEN
*
*       Cases for dP/dx:
*
*        Clear fluid/porous media interface at east face
*
         IF((CVTYPE(I,J,K,1).EQ.0.AND.CVTYPE(I,J,K,3).EQ.2)
     C   .OR.(CVTYPE(I,J,K,1).EQ.2.AND.CVTYPE(I,J,K,3).EQ.0)) THEN
          AP(2,1,I,J,K)= 0.5*PRSTY(I,J,K)*VOLP(I,J,K)
     C                   *(1.0/DIEP(I-1)-1.0/DIEP(I))
          AW(2,1,I,J,K)= 0.5*PRSTY(I,J,K)*VOLP(I,J,K)/DIEP(I-1)
          AE(2,1,I,J,K)=-0.5*PRSTY(I,J,K)*VOLP(I,J,K)/DIEP(I)
          AS(2,1,I,J,K)= 0.0
          AN(2,1,I,J,K)= 0.0
          AB(2,1,I,J,K)= 0.0
          AT(2,1,I,J,K)= 0.0          
          B(2,I,J,K)= B(2,I,J,K)-0.5*PRSTY(I,J,K)*VOLP(I,J,K)
     C                *((PEAST(I,J,K)-P(I,J,K))/DISE(I)
     C                -(P(I+1,J,K)-P(I,J,K))/DIEP(I))
*
*        Clear fluid/porous media interface at west face
*
         ELSEIF((CVTYPE(I,J,K,1).EQ.0.AND.CVTYPE(I,J,K,2).EQ.2) 
     C   .OR.(CVTYPE(I,J,K,1).EQ.2.AND.CVTYPE(I,J,K,2).EQ.0)) THEN
          AP(2,1,I,J,K)= 0.5*PRSTY(I,J,K)*VOLP(I,J,K)
     C                   *(1.0/DIEP(I-1)-1.0/DIEP(I))
          AW(2,1,I,J,K)= 0.5*PRSTY(I,J,K)*VOLP(I,J,K)/DIEP(I-1)
          AE(2,1,I,J,K)=-0.5*PRSTY(I,J,K)*VOLP(I,J,K)/DIEP(I)
          AS(2,1,I,J,K)= 0.0
          AN(2,1,I,J,K)= 0.0
          AB(2,1,I,J,K)= 0.0
          AT(2,1,I,J,K)= 0.0          
          B(2,I,J,K)= B(2,I,J,K)-0.5*PRSTY(I,J,K)*VOLP(I,J,K)
     C                *((P(I,J,K)-PEAST(I-1,J,K))/DISE(I)
     C                -(P(I,J,K)-P(I-1,J,K))/DIEP(I-1))
*
*        Solid interface at east face
*
         ELSEIF(CVTYPE(I,J,K,3).EQ.3) THEN
          AP(2,1,I,J,K)= PRSTY(I,J,K)*VOLP(I,J,K)/DIEP(I-1)
          AW(2,1,I,J,K)= PRSTY(I,J,K)*VOLP(I,J,K)/DIEP(I-1)
          AE(2,1,I,J,K)= 0.0
          AS(2,1,I,J,K)= 0.0
          AN(2,1,I,J,K)= 0.0
          AB(2,1,I,J,K)= 0.0
          AT(2,1,I,J,K)= 0.0          
*
*        Solid interface at west face
*
         ELSEIF(CVTYPE(I,J,K,2).EQ.3) THEN
          AP(2,1,I,J,K)= PRSTY(I,J,K)*VOLP(I,J,K)*(-1.0/DIEP(I))
          AW(2,1,I,J,K)= 0.0
          AE(2,1,I,J,K)= PRSTY(I,J,K)*VOLP(I,J,K)*(-1.0/DIEP(I))
          AS(2,1,I,J,K)= 0.0
          AN(2,1,I,J,K)= 0.0
          AB(2,1,I,J,K)= 0.0
          AT(2,1,I,J,K)= 0.0
*
*        No interface
*
         ELSE
          AP(2,1,I,J,K)= 0.5*PRSTY(I,J,K)*VOLP(I,J,K)
     C                   *(1.0/DIEP(I-1)-1.0/DIEP(I))
          AW(2,1,I,J,K)= 0.5*PRSTY(I,J,K)*VOLP(I,J,K)/DIEP(I-1)
          AE(2,1,I,J,K)=-0.5*PRSTY(I,J,K)*VOLP(I,J,K)/DIEP(I)
          AS(2,1,I,J,K)= 0.0
          AN(2,1,I,J,K)= 0.0
          AB(2,1,I,J,K)= 0.0
          AT(2,1,I,J,K)= 0.0          
         ENDIF
*
*       Cases for dP/dy:
*
*        Clear fluid/porous media interface at north face
*
         IF((CVTYPE(I,J,K,1).EQ.0.AND.CVTYPE(I,J,K,5).EQ.2)
     C   .OR.(CVTYPE(I,J,K,1).EQ.2.AND.CVTYPE(I,J,K,5).EQ.0)) THEN
          AP(3,1,I,J,K)= 0.5*PRSTY(I,J,K)*VOLP(I,J,K)
     C                   *(1.0/DJNP(J-1)-1.0/DJNP(J))
          AW(3,1,I,J,K)= 0.0
          AE(3,1,I,J,K)= 0.0
          AS(3,1,I,J,K)= 0.5*PRSTY(I,J,K)*VOLP(I,J,K)/DJNP(J-1)
          AN(3,1,I,J,K)=-0.5*PRSTY(I,J,K)*VOLP(I,J,K)/DJNP(J)
          AB(3,1,I,J,K)= 0.0
          AT(3,1,I,J,K)= 0.0          
          B(3,I,J,K)= B(3,I,J,K)-0.5*PRSTY(I,J,K)*VOLP(I,J,K)
     C                *((PNORTH(I,J,K)-P(I,J,K))/DISN(J)
     C                -(P(I,J+1,K)-P(I,J,K))/DJNP(J))
*
*        Clear fluid/porous media interface at south face
*
         ELSEIF((CVTYPE(I,J,K,1).EQ.0.AND.CVTYPE(I,J,K,4).EQ.2)
     C   .OR.(CVTYPE(I,J,K,1).EQ.2.AND.CVTYPE(I,J,K,4).EQ.0)) THEN
          AP(3,1,I,J,K)= 0.5*PRSTY(I,J,K)*VOLP(I,J,K)
     C                   *(1.0/DJNP(J-1)-1.0/DJNP(J))
          AW(3,1,I,J,K)= 0.0
          AE(3,1,I,J,K)=0.0
          AS(3,1,I,J,K)= 0.5*PRSTY(I,J,K)*VOLP(I,J,K)/DJNP(J-1)
          AN(3,1,I,J,K)=-0.5*PRSTY(I,J,K)*VOLP(I,J,K)/DJNP(J)
          AB(3,1,I,J,K)= 0.0
          AT(3,1,I,J,K)= 0.0          
          B(3,I,J,K)= B(3,I,J,K)-0.5*PRSTY(I,J,K)*VOLP(I,J,K)
     C                *((P(I,J,K)-PNORTH(I,J-1,K))/DISN(J)
     C                -(P(I,J,K)-P(I,J-1,K))/DJNP(J-1))
*
*
*        Solid interface at north face
*
         ELSEIF(CVTYPE(I,J,K,5).EQ.3) THEN
          AP(3,1,I,J,K)= PRSTY(I,J,K)*VOLP(I,J,K)/DJNP(J-1)
          AW(3,1,I,J,K)= 0.0
          AE(3,1,I,J,K)= 0.0
          AS(3,1,I,J,K)= PRSTY(I,J,K)*VOLP(I,J,K)/DJNP(J-1)
          AN(3,1,I,J,K)= 0.0
          AB(3,1,I,J,K)= 0.0
          AT(3,1,I,J,K)= 0.0          
*
*        Solid interface at south face
*
         ELSEIF(CVTYPE(I,J,K,4).EQ.3) THEN
          AP(3,1,I,J,K)= PRSTY(I,J,K)*VOLP(I,J,K)*(-1.0/DJNP(J))
          AW(3,1,I,J,K)= 0.0
          AE(3,1,I,J,K)= 0.0
          AS(3,1,I,J,K)= 0.0
          AN(3,1,I,J,K)= PRSTY(I,J,K)*VOLP(I,J,K)*(-1.0/DJNP(J))
          AB(3,1,I,J,K)= 0.0
          AT(3,1,I,J,K)= 0.0          
*
*        No interface
*
         ELSE
          AP(3,1,I,J,K)= 0.5*PRSTY(I,J,K)*VOLP(I,J,K)
     C                   *(1.0/DJNP(J-1)-1.0/DJNP(J))
          AW(3,1,I,J,K)= 0.0
          AE(3,1,I,J,K)= 0.0
          AS(3,1,I,J,K)= 0.5*PRSTY(I,J,K)*VOLP(I,J,K)/DJNP(J-1)
          AN(3,1,I,J,K)=-0.5*PRSTY(I,J,K)*VOLP(I,J,K)/DJNP(J) 
          AB(3,1,I,J,K)= 0.0
          AT(3,1,I,J,K)= 0.0          
         ENDIF
*
*       Cases for dP/dz:
*
*        Clear fluid/porous media interface at top face
*
         IF((CVTYPE(I,J,K,1).EQ.0.AND.CVTYPE(I,J,K,7).EQ.2)
     C   .OR.(CVTYPE(I,J,K,1).EQ.2.AND.CVTYPE(I,J,K,7).EQ.0)) THEN
          AP(4,1,I,J,K)= 0.5*PRSTY(I,J,K)*VOLP(I,J,K)
     C                   *(1.0/DKTP(K-1)-1.0/DKTP(K))
          AW(4,1,I,J,K)= 0.0
          AE(4,1,I,J,K)= 0.0
          AS(4,1,I,J,K)= 0.0
          AN(4,1,I,J,K)= 0.0
          AB(4,1,I,J,K)= 0.5*PRSTY(I,J,K)*VOLP(I,J,K)/DKTP(K-1)
          AT(4,1,I,J,K)=-0.5*PRSTY(I,J,K)*VOLP(I,J,K)/DKTP(K)          
          B(4,I,J,K)= B(4,I,J,K)-0.5*PRSTY(I,J,K)*VOLP(I,J,K)
     C                *((PTOP(I,J,K)-P(I,J,K))/DIST(K)
     C                -(P(I,J,K+1)-P(I,J,K))/DKTP(K))
*
*        Clear fluid/porous media interface at bottom face
*
         ELSEIF((CVTYPE(I,J,K,1).EQ.0.AND.CVTYPE(I,J,K,6).EQ.2)
     C   .OR.(CVTYPE(I,J,K,1).EQ.2.AND.CVTYPE(I,J,K,6).EQ.0)) THEN
          AP(4,1,I,J,K)= 0.5*PRSTY(I,J,K)*VOLP(I,J,K)
     C                   *(1.0/DKTP(K-1)-1.0/DKTP(K))
          AW(4,1,I,J,K)= 0.0
          AE(4,1,I,J,K)=0.0
          AS(4,1,I,J,K)= 0.0
          AN(4,1,I,J,K)= 0.0
          AB(4,1,I,J,K)= 0.5*PRSTY(I,J,K)*VOLP(I,J,K)/DKTP(K-1)
          AT(4,1,I,J,K)=-0.5*PRSTY(I,J,K)*VOLP(I,J,K)/DKTP(K)          
          B(4,I,J,K)= B(4,I,J,K)-0.5*PRSTY(I,J,K)*VOLP(I,J,K)
     C                *((P(I,J,K)-PTOP(I,J,K-1))/DIST(K)
     C                -(P(I,J,K)-P(I,J,K-1))/DKTP(K-1))
*
*
*        Solid interface at top face
*
         ELSEIF(CVTYPE(I,J,K,7).EQ.3) THEN
          AP(4,1,I,J,K)= PRSTY(I,J,K)*VOLP(I,J,K)/DKTP(K-1)
          AW(4,1,I,J,K)= 0.0
          AE(4,1,I,J,K)= 0.0
          AS(4,1,I,J,K)= 0.0
          AN(4,1,I,J,K)= 0.0
          AB(4,1,I,J,K)= PRSTY(I,J,K)*VOLP(I,J,K)/DKTP(K-1)
          AT(4,1,I,J,K)= 0.0          
*
*        Solid interface at bottom face
*
         ELSEIF(CVTYPE(I,J,K,6).EQ.3) THEN
          AP(4,1,I,J,K)= PRSTY(I,J,K)*VOLP(I,J,K)*(-1.0/DKTP(K))
          AW(4,1,I,J,K)= 0.0
          AE(4,1,I,J,K)= 0.0
          AS(4,1,I,J,K)= 0.0
          AN(4,1,I,J,K)= 0.0
          AB(4,1,I,J,K)= 0.0
          AT(4,1,I,J,K)= PRSTY(I,J,K)*VOLP(I,J,K)*(-1.0/DKTP(K))          
*
*        No interface
*
         ELSE
          AP(4,1,I,J,K)= 0.5*PRSTY(I,J,K)*VOLP(I,J,K)
     C                   *(1.0/DKTP(K-1)-1.0/DKTP(K))
          AW(4,1,I,J,K)= 0.0
          AE(4,1,I,J,K)= 0.0
          AS(4,1,I,J,K)= 0.0
          AN(4,1,I,J,K)= 0.0 
          AB(4,1,I,J,K)= 0.5*PRSTY(I,J,K)*VOLP(I,J,K)/DKTP(K-1)
          AT(4,1,I,J,K)=-0.5*PRSTY(I,J,K)*VOLP(I,J,K)/DKTP(K)          
         ENDIF
        ENDIF        
*        
 10     CONTINUE
 15    CONTINUE
 20   CONTINUE
*
      RETURN
      END