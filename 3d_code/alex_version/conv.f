*
*     file conv.f
************************************************************************
      SUBROUTINE CONV(HSF, U,V,W,PRSTY,PD,DEQ,
     C                RHO,CP,COND,VISC,
     C                PCONV,CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB,
     C                COND0,RHO0,VISC0,CP0)
*
*     Subroutine to calculate the interface exchange coefficient for
*     the non-equilibrium heat transfer model in porous volumes using 
*     a relevant Nusselt number correlation. 
*
*     HSF(ID,JD,KD)  interfacial exchange coefficient at porous 
*                    volumes; output
*     U,V,W  fluid velocity components at each volume; input
*     PRSTY(ID,JD,KD)  porosity array; input
*     PD  pore diameter of porous medium; input
*     DEQ(ID,JD,KD)  equivalent solid particle diameter array; input
*     RHO  fluid density; input
*     CP  fluid specific heat capacity; input
*     COND  fluid thermal conductivity; input
*     VISC  fluid viscosity; input
*     PCONV sets convection correlation; input
*     CVTYPE(ID,JD,KD,NNB+1)  indicates type of each CV; input
*     IB,IE  first and last interior indices in i; input
*     JB,JE  first and last interior indices in j; input
*     KB,KE  first and last interior indices in k; input
*     ID,JD,KD  array dimensions; input
*     NNB  number of neighbors for interior nodes; input
*
************************************************************************
*
      IMPLICIT NONE
      REAL*8 HSF(ID,JD,KD)
      REAL*8 U(ID,JD,KD),V(ID,JD,KD),W(ID,JD,KD)
      REAL*8 PRSTY(ID,JD,KD),PD,DEQ(ID,JD,KD)
      REAL*8 PR,RE,CT,VEL,MEXP,NEXP
      INTEGER IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB,CVTYPE(ID,JD,KD,NNB+1)
      INTEGER PCONV,I,J,K
*
      REAL*8 RHO(ID,JD,KD),CP(ID,JD,KD),COND(ID,JD,KD),VISC(ID,JD,KD)
      REAL*8 COND0,RHO0,VISC0,CP0
*
*     Set coefficients for specified correlation
*
      IF(PCONV.EQ.0) THEN
       CT = 0.52
       MEXP = 0.5
       NEXP = 0.37
      ELSEIF(PCONV.EQ.1) THEN
       CT = 0.018
       MEXP = 0.27
       NEXP = 1.0/3.0
      ELSEIF(PCONV.EQ.2) THEN
       CT = 0.0075
       MEXP = 0.36
       NEXP = 1.0/3.0
      ELSEIF(PCONV.EQ.3) THEN
       CT = 1000.0
       MEXP = 0.0
       NEXP = 0.0
      ELSEIF(PCONV.EQ.4) THEN
       CT = 0.97
       MEXP = 0.6
       NEXP = 0.0
      ELSEIF(PCONV.EQ.5) THEN
       CT = 0.001
       MEXP = 1.97
       NEXP = 1.0/3.0
      ENDIF
*
      DO 30 K=KB,KE
       DO 20 J=JB,JE
        DO 10 I=IB,IE
         IF(CVTYPE(I,J,K,1).EQ.2) THEN
           VEL = (U(I,J,K)**2+V(I,J,K)**2+W(I,J,K)**2)**0.5
           PR = CP(I,J,K)*VISC(I,J,K)/COND(I,J,K)                         !
           IF(PCONV.EQ.0) THEN
            RE = RHO(I,J,K)*VEL*DEQ(I,J,K)/(VISC(I,J,K)*PRSTY(I,J,K))     !
           ELSE
            RE = RHO(I,J,K)*VEL*DEQ(I,J,K)/VISC(I,J,K)                    !
           ENDIF
           IF (PCONV.EQ.6) THEN    !FULL ASHRAF CORRELATION - constants found from his thesis
            HSF(I,J,K) = COND(I,JB,K)/COND0 + 0.82*(PR**(1.0/3.0))*
     C                   (RE**0.5)*
     C                   ((RHO0*VISC0/RHO(I,JB,K)/VISC(I,JB,K))**0.149)*
     C                   ((CP0/CP(I,JB,K))**1.424)
           ELSE
            HSF(I,J,K) = COND(I,J,K)*CT*(RE**MEXP)*(PR**NEXP)/DEQ(I,J,K)   !
           ENDIF 
         ENDIF
 10     CONTINUE
 20    CONTINUE
 30   CONTINUE
      RETURN
      END
