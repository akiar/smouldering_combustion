*     file addisp.f
***********************************************************************
*
      SUBROUTINE ADDISP(KEFF, KEFF0,U,V,W,PRSTY,KPERM,PD,LD,
     C                  RHO,CP,COND,VISC,
     C                  CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
*
*     Routine to add dispersion conductivity.
*
*     KEFF(ID,JD,KD)  effective fluid conductivity array; output
*     KEFF0(ID,JD,KD)  stagnant effective fluid conductivity array;
*                      input
*     U,V,W  fluid velocity components at each volume; input
*     PRSTY(ID,JD,KD) porosity array; input
*     KPERM(ID,JD,KD) permeability array; input
*     PD pore diameter of porous medium; input
*     LD ligament diameter of porous medium; input
*     RHO  fluid density; input
*     CP  fluid specific heat capacity; input
*     COND  fluid thermal conductivity; input
*     VISC  fluid viscosity; input
*     CVTYPE(ID,JD,KD,NNB+1)  indicates type of each CV; input
*     IB,IE  first and last interior indices in i; input
*     JB,JE  first and last interior indices in j; input
*     KB,KE  first and last interior indices in k; input
*     ID,JD,KD  array dimensions; input
*     NNB  number of neighbors for interior nodes; input
*
***********************************************************************
*
      IMPLICIT NONE
      REAL*8 KEFF(ID,JD,KD),KEFF0(ID,JD,KD)
      REAL*8 U(ID,JD,KD),V(ID,JD,KD),W(ID,JD,KD)
      REAL*8 PRSTY(ID,JD,KD),KPERM(ID,JD,KD),PD,LD,RHO,CP,COND,VISC
      REAL*8 VEL,REK,PREKE,CD,KDISP
      INTEGER IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB,CVTYPE(ID,JD,KD,NNB+1)
      INTEGER I,J,K
*
*     Interior volumes
*
      CD = 0.06
      PREKE = VISC*CP
      DO 20 I=IB,IE      
       DO 10 J=JB,JE
        DO 5 K=KB,KE
         IF(CVTYPE(I,J,K,1).EQ.2) THEN
          VEL = (U(I,J,K)**2+V(I,J,K)**2+W(I,J,K)**2)**0.5
          REK = RHO*VEL*(KPERM(I,J,K)**0.5)/VISC
          KDISP = CD*REK*PREKE
          KEFF(I,J,K) = KEFF0(I,J,K) + KDISP
         ELSE
          KEFF(I,J,K) = KEFF0(I,J,K)
         ENDIF
    5   CONTINUE
   10  CONTINUE
   20 CONTINUE
*
*     Boundary nodes
*
      DO 30 I=IB,IE
       DO 25 K=KB,KE
        J = JB-1
        IF(CVTYPE(I,J,K,5).EQ.2) THEN
         VEL = (U(I,J,K)**2+V(I,J,K)**2+W(I,J,K)**2)**0.5
         REK = RHO*VEL*(KPERM(I,J,K)**0.5)/VISC
         KDISP = CD*REK*PREKE
         KEFF(I,J,K) = KEFF0(I,J,K) + KDISP
        ELSE
         KEFF(I,J,K) = KEFF0(I,J,K)
        ENDIF
*
        J = JE+1        
        IF(CVTYPE(I,J,K,4).EQ.2) THEN
         VEL = (U(I,J,K)**2+V(I,J,K)**2+W(I,J,K)**2)**0.5
         REK = RHO*VEL*(KPERM(I,J,K)**0.5)/VISC
         KDISP = CD*REK*PREKE
         KEFF(I,J,K) = KEFF0(I,J,K) + KDISP
        ELSE
         KEFF(I,J,K) = KEFF0(I,J,K)
        ENDIF
   25  CONTINUE
   30 CONTINUE
      DO 40 J=JB,JE
       DO 35 K=KB,KE
        I = IB-1
        IF(CVTYPE(I,J,K,3).EQ.2) THEN
         VEL = (U(I,J,K)**2+V(I,J,K)**2+W(I,J,K)**2)**0.5
         REK = RHO*VEL*(KPERM(I,J,K)**0.5)/VISC
         KDISP = CD*REK*PREKE
         KEFF(I,J,K) = KEFF0(I,J,K) + KDISP
        ELSE
         KEFF(I,J,K) = KEFF0(I,J,K)
        ENDIF
*
        I = IE+1        
        IF(CVTYPE(I,J,K,2).EQ.2) THEN
         VEL = (U(I,J,K)**2+V(I,J,K)**2+W(I,J,K)**2)**0.5
         REK = RHO*VEL*(KPERM(I,J,K)**0.5)/VISC
         KDISP = CD*REK*PREKE
         KEFF(I,J,K) = KEFF0(I,J,K) + KDISP
        ELSE
         KEFF(I,J,K) = KEFF0(I,J,K)
        ENDIF
   35  CONTINUE
   40 CONTINUE   
      DO 50 I=IB,IE
       DO 45 J=JB,JE
        K = KB-1
        IF(CVTYPE(I,J,K,7).EQ.2) THEN
         VEL = (U(I,J,K)**2+V(I,J,K)**2+W(I,J,K)**2)**0.5
         REK = RHO*VEL*(KPERM(I,J,K)**0.5)/VISC
         KDISP = CD*REK*PREKE
         KEFF(I,J,K) = KEFF0(I,J,K) + KDISP
        ELSE
         KEFF(I,J,K) = KEFF0(I,J,K)
        ENDIF
*
        K = KE+1        
        IF(CVTYPE(I,J,K,6).EQ.2) THEN
         VEL = (U(I,J,K)**2+V(I,J,K)**2+W(I,J,K)**2)**0.5
         REK = RHO*VEL*(KPERM(I,J,K)**0.5)/VISC
         KDISP = CD*REK*PREKE
         KEFF(I,J,K) = KEFF0(I,J,K) + KDISP
        ELSE
         KEFF(I,J,K) = KEFF0(I,J,K)
        ENDIF
   45  CONTINUE
   50 CONTINUE  
*
      RETURN
      END