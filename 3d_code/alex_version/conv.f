*
*     file conv.f
************************************************************************
      SUBROUTINE CONV(HSF, U,V,W,PRSTY,PD,DEQ,
     C                RHO,CP,COND,VISC,
     C                PCONV,CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
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
      REAL*8 PRSTY(ID,JD,KD),PD,DEQ(ID,JD,KD),RHO,CP,COND,VISC
      REAL*8 PR,RE,CT,VEL,MEXP,NEXP
      INTEGER IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB,CVTYPE(ID,JD,KD,NNB+1)
      INTEGER PCONV,I,J,K
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
      ENDIF
*
      DO 30 K=KB,KE
       DO 20 J=JB,JE
        DO 10 I=IB,IE
         IF(CVTYPE(I,J,K,1).EQ.2) THEN
           VEL = (U(I,J,K)**2+V(I,J,K)**2+W(I,J,K)**2)**0.5
           PR = CP*VISC/COND
           IF(PCONV.EQ.0) THEN
            RE = RHO*VEL*DEQ(I,J,K)/(VISC*PRSTY(I,J,K))
           ELSE
            RE = RHO*VEL*DEQ(I,J,K)/VISC
           ENDIF
           HSF(I,J,K) = COND*CT*(RE**MEXP)*(PR**NEXP)/DEQ(I,J,K)
         ENDIF
 10     CONTINUE
 20    CONTINUE
 30   CONTINUE
      RETURN
      END
