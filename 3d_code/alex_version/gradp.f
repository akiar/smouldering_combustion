*
*     This file contains 2 subroutines: GRADP and PINT
*
************************************************************************
*
      SUBROUTINE GRADP(DPDX,DPDXP,DPDXM,DPDY,DPDYP,DPDYM,
     C                 DPDZ,DPDZP,DPDZM, P,PEAST,PNORTH,PTOP,
     C                 DIEP,DJNP,DKTP,DISE,DISN,DIST,
     C                 CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
*
*     Subroutine to calculate the pressure gradient at each nodal
*     location using the most current pressure field.
*
*     DPDX(ID,JD,KD) x-derivative of pressure at nodal location; output
*     DPDY(ID,JD,KD) y-derivative of pressure at nodal location; output
*     DPDZ(ID,JD,KD) z-derivative of pressure at nodal location; output
*     DPDXP(ID,JD,KD)  x-derivative from P to e; output
*     DPDYP(ID,JD,KD)  y-derivative from P to n; output
*     DPDZP(ID,JD,KD)  z-derivative from P to t; output
*     DPDXM(ID,JD,KD)  x-derivative from w to P; output
*     DPDYM(ID,JD,KD)  y-derivative from s to P; output
*     DPDZM(ID,JD,KD)  z-derivative from b to P; output
*     P(ID,JD,KD)  current pressure field; input
*     PEAST(ID,JD,KD)  current pressures at x-interface; input
*     PNORTH(ID,JD,KD)  current pressure at y-interface; input
*     PTOP(ID,JD,KD)  current pressure at z-interface; input
*     DIEP(ID)  distance from P to E nodes; input
*     DJNP(JD)  distance from P to N nodes; input
*     DKTP(KD)  distance from P to T nodes; input
*     DISE(ID)  distance from P to e; input
*     DISN(JD)  distance from P to n; input
*     DIST(KD)  distance from P to t; input
*     CVTYPE(ID,JD,KD,NNB+1)  control volume types; input
*     IB,IE     first and last interior indices in i; input
*     JB,JE     first and last interior indices in j; input
*     JB,JE     first and last interior indices in k; input
*     ID,JD,KD  array dimensions; input     
*     NNB       number of neighbors for interior nodes; input
*
***********************************************************************
*
      IMPLICIT NONE
      REAL*8 DPDX(ID,JD,KD),DPDXP(ID,JD,KD),DPDXM(ID,JD,KD)
      REAL*8 DPDY(ID,JD,KD),DPDYP(ID,JD,KD),DPDYM(ID,JD,KD)
      REAL*8 DPDZ(ID,JD,KD),DPDZP(ID,JD,KD),DPDZM(ID,JD,KD)
      REAL*8 P(ID,JD,KD),PEAST(ID,JD,KD),PNORTH(ID,JD,KD),PTOP(ID,JD,KD)
      REAL*8 DIEP(ID),DJNP(JD),DKTP(KD),DISE(ID),DISN(JD),DIST(KD)
      INTEGER IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB,CVTYPE(ID,JD,KD,NNB+1)
      INTEGER I,J,K
*
*--Derivatives of P at interior nodes
*
      DO 20 I=IB,IE
       DO 10 J=JB,JE
        DO 5 K=KB,KE
         IF(CVTYPE(I,J,K,1).NE.3) THEN 
*
*        Cases for dP/dx:
*
*         Clear fluid/porous media interface at east face 
*
          IF((CVTYPE(I,J,K,1).EQ.0 .AND. CVTYPE(I,J,K,3).EQ.2)
     C    .OR.(CVTYPE(I,J,K,1).EQ.2 .AND. CVTYPE(I,J,K,3).EQ.0)) THEN
           DPDXP(I,J,K)= (PEAST(I,J,K)-P(I,J,K))/DISE(I)
           DPDXM(I,J,K)= (P(I,J,K)-P(I-1,J,K))/DIEP(I-1)
           DPDX(I,J,K)= 0.5*(DPDXP(I,J,K)+DPDXM(I,J,K))
*
*         Clear fluid/porous media interface at west face 
*
          ELSEIF((CVTYPE(I,J,K,1).EQ.0 .AND. CVTYPE(I,J,K,2).EQ.2)
     C    .OR.(CVTYPE(I,J,K,1).EQ.2 .AND. CVTYPE(I,J,K,2).EQ.0)) THEN
           DPDXP(I,J,K)= (P(I+1,J,K)-P(I,J,K))/DIEP(I)
           DPDXM(I,J,K)= (P(I,J,K)-PEAST(I-1,J,K))/DISE(I)
           DPDX(I,J,K)= 0.5*(DPDXP(I,J,K)+DPDXM(I,J,K))          
*
*         Solid interface at east face 
*
          ELSEIF(CVTYPE(I,J,K,3).EQ.3) THEN
           DPDX(I,J,K)= (P(I,J,K)-P(I-1,J,K))/DIEP(I-1)
           DPDXP(I,J,K)= DPDX(I,J,K)
           DPDXM(I,J,K)= DPDX(I,J,K)
* 
*         Solid interface at west face 
*      
          ELSEIF(CVTYPE(I,J,K,2).EQ.3) THEN
           DPDX(I,J,K)= (P(I+1,J,K)-P(I,J,K))/DIEP(I)
           DPDXP(I,J,K)= DPDX(I,J,K)
           DPDXM(I,J,K)= DPDX(I,J,K)
*
*         No interface
*     
          ELSE
           DPDX(I,J,K)= 0.5*((P(I+1,J,K)-P(I,J,K))/DIEP(I)
     C                  +(P(I,J,K)-P(I-1,J,K))/DIEP(I-1))
           DPDXP(I,J,K)= DPDX(I,J,K)
           DPDXM(I,J,K)= DPDX(I,J,K)
          ENDIF
*
*        Cases for dP/dy:
*
*         Clear fluid/porous media interface at north face 
*
          IF((CVTYPE(I,J,K,1).EQ.0 .AND. CVTYPE(I,J,K,5).EQ.2)
     C    .OR.(CVTYPE(I,J,K,1).EQ.2 .AND. CVTYPE(I,J,K,5).EQ.0)) THEN
           DPDYP(I,J,K)= (PNORTH(I,J,K)-P(I,J,K))/DISN(J)
           DPDYM(I,J,K)= (P(I,J,K)-P(I,J-1,K))/DJNP(J-1)
           DPDY(I,J,K)= 0.5*(DPDYP(I,J,K)+DPDYM(I,J,K))
*
*         Clear fluid/porous media interface at south face 
*
          ELSEIF((CVTYPE(I,J,K,1).EQ.0 .AND. CVTYPE(I,J,K,4).EQ.2)
     C    .OR.(CVTYPE(I,J,K,1).EQ.2 .AND. CVTYPE(I,J,K,4).EQ.0)) THEN
           DPDYP(I,J,K)= (P(I,J+1,K)-P(I,J,K))/DJNP(J)
           DPDYM(I,J,K)= (P(I,J,K)-PNORTH(I,J-1,K))/DISN(J)
           DPDY(I,J,K)= 0.5*(DPDYP(I,J,K)+DPDYM(I,J,K))     
*
*         Solid interface at north face 
*
          ELSEIF(CVTYPE(I,J,K,5).EQ.3) THEN
           DPDY(I,J,K)= (P(I,J,K)-P(I,J-1,K))/DJNP(J-1)
           DPDYP(I,J,K)= DPDY(I,J,K)
           DPDYM(I,J,K)= DPDY(I,J,K)               
*
*         Solid interface at south face 
*     
          ELSEIF(CVTYPE(I,J,K,4).EQ.3) THEN
           DPDY(I,J,K)= (P(I,J+1,K)-P(I,J,K))/DJNP(J)
           DPDYP(I,J,K)= DPDY(I,J,K)
           DPDYM(I,J,K)= DPDY(I,J,K)     
*
*         No interface
*     
          ELSE
           DPDY(I,J,K)= 0.5*((P(I,J+1,K)-P(I,J,K))/DJNP(J)
     C                  +(P(I,J,K)-P(I,J-1,K))/DJNP(J-1))
           DPDYP(I,J,K)= DPDY(I,J,K)
           DPDYM(I,J,K)= DPDY(I,J,K)     
          ENDIF
*
*        Cases for dP/dz:
*
*         Clear fluid/porous media interface at top face 
*
          IF((CVTYPE(I,J,K,1).EQ.0 .AND. CVTYPE(I,J,K,7).EQ.2)
     C    .OR.(CVTYPE(I,J,K,1).EQ.2 .AND. CVTYPE(I,J,K,7).EQ.0)) THEN
           DPDZP(I,J,K)= (PTOP(I,J,K)-P(I,J,K))/DIST(K)
           DPDZM(I,J,K)= (P(I,J,K)-P(I,J,K-1))/DKTP(K-1)
           DPDZ(I,J,K)= 0.5*(DPDZP(I,J,K)+DPDZM(I,J,K))
*
*         Clear fluid/porous media interface at bottom face 
*
          ELSEIF((CVTYPE(I,J,K,1).EQ.0 .AND. CVTYPE(I,J,K,6).EQ.2)
     C    .OR.(CVTYPE(I,J,K,1).EQ.2 .AND. CVTYPE(I,J,K,6).EQ.0)) THEN
           DPDZP(I,J,K)= (P(I,J,K+1)-P(I,J,K))/DKTP(K)
           DPDZM(I,J,K)= (P(I,J,K)-PTOP(I,J,K-1))/DIST(K)
           DPDZ(I,J,K)= 0.5*(DPDZP(I,J,K)+DPDZM(I,J,K))          
*
*         Solid interface at top face 
*
          ELSEIF(CVTYPE(I,J,K,7).EQ.3) THEN
           DPDZ(I,J,K)= (P(I,J,K)-P(I,J,K-1))/DKTP(K-1)
           DPDZP(I,J,K)= DPDZ(I,J,K)
           DPDZM(I,J,K)= DPDZ(I,J,K)
* 
*         Solid interface at bottom face 
*      
          ELSEIF(CVTYPE(I,J,K,6).EQ.3) THEN
           DPDZ(I,J,K)= (P(I,J,K+1)-P(I,J,K))/DKTP(K)
           DPDZP(I,J,K)= DPDZ(I,J,K)
           DPDZM(I,J,K)= DPDZ(I,J,K)
*
*         No interface
*     
          ELSE
           DPDZ(I,J,K)= 0.5*((P(I,J,K+1)-P(I,J,K))/DKTP(K)
     C                  +(P(I,J,K)-P(I,J,K-1))/DKTP(K-1))
           DPDZP(I,J,K)= DPDZ(I,J,K)
           DPDZM(I,J,K)= DPDZ(I,J,K)
          ENDIF
         ENDIF
   5    CONTINUE
  10   CONTINUE
  20  CONTINUE
*
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE PINT(PEAST,PNORTH,PTOP, 
     C                P,DIEP,DJNP,DKTP,DISE,DISN,DIST,
     C                ME,MN,MT,U,V,W,PRSTY,DE,DN,DT,MUEFF,
     C                AREP,ARNP,ARTP,CONU,BLEND,ALFAE,ALFAN,ALFAT,
     C                CLVL,CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
*
*     Subroutine to calculate the pressure at porous/clear fluid
*     interfaces using the most current pressure and velocity fields.
*
*     Notes: CLVL controls the level of correction as follows:
*
*                CLVL=0: Interpolation 
*                CLVL=1: Extrapolation with correction for effects of
*                        contraction/expansion
*
*            For CLVL = 1, interface estimate comes from average
*            of second-order Taylor series extrapolation from porous 
*            region and second-order Taylor series extrapolation
*            with dynamic corrections from pure fluid region.
*
*     PEAST(ID,JD,KD)  current pressures at x-interface; output
*     PNORTH(ID,JD,KD)  current pressures at y-interface; output
*     PTOP(ID,JD,KD)  current pressures at z-interface; output
*
***********************************************************************
*
      IMPLICIT NONE
      REAL*8 PEAST(ID,JD,KD),PNORTH(ID,JD,KD),PTOP(ID,JD,KD)
      REAL*8 P(ID,JD,KD),DIEP(ID),DJNP(JD),DKTP(KD)
      REAL*8 DISE(ID),DISN(JD),DIST(KD)
      REAL*8 ME(ID,JD,KD),MN(ID,JD,KD),MT(ID,JD,KD)
      REAL*8 U(ID,JD,KD),V(ID,JD,KD),W(ID,JD,KD),PRSTY(ID,JD,KD)
      REAL*8 DE(ID,JD,KD),DN(ID,JD,KD),DT(ID,JD,KD)
      REAL*8 MUEFF(ID,JD,KD),AREP(JD,KD),ARNP(ID,KD),ARTP(ID,JD)
      REAL*8 BLEND,ALFAE(ID,JD,KD),ALFAN(ID,JD,KD),ALFAT(ID,JD,KD)
      REAL*8 CONU(ID,JD,KD,NNB+1),UBYEPS,VBYEPS,WBYEPS
      REAL*8 CSUBE,CSUBN,CSUBT
      INTEGER IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB,CVTYPE(ID,JD,KD,NNB+1)
      INTEGER CLVL,I,J,K
*
      DO 20 I=IB,IE
       DO 10 J=JB,JE
        DO 5 K=KB,KE
*
*        Estimate pressure at porous/clear fluid interfaces 
*        perpendicular to x-direction. 
*        
         IF(CVTYPE(I,J,K,1).EQ.0 .AND. CVTYPE(I,J,K,3).EQ.2) THEN 
          IF(CLVL.EQ.0) THEN
           PEAST(I,J,K)= (DISE(I+1)*P(I,J,K)
     C                   +DISE(I)*P(I+1,J,K))/DIEP(I)
          ELSE        
           CSUBE= DISE(I)*DE(I,J,K)/(MUEFF(I,J,K)*AREP(J,K))
           UBYEPS= BLEND*((CSUBE/PRSTY(I+1,J,K)-0.5*CONU(I,J,K,3)
     C        *(1.0-ALFAE(I,J,K)))*U(I+1,J,K)+((1.0-CSUBE)
     C        /PRSTY(I+1,J,K)-0.5*CONU(I,J,K,1)*(1.0
     C        +ALFAE(I,J,K)))*U(I,J,K))+0.5*CONU(I,J,K,1)*(1.0
     C        +ALFAE(I,J,K))*U(I,J,K)+0.5*CONU(I,J,K,3)*(1.0
     C        -ALFAE(I,J,K))*U(I+1,J,K)
           PEAST(I,J,K)= 0.5*( (P(I,J,K)+DISE(I)*(P(I,J,K)
     C      -P(I-1,J,K))/DIEP(I-1)+0.5*(((P(I,J,K)-P(I-1,J,K))
     C      /DIEP(I-1)-(P(I-1,J,K)-P(I-2,J,K))/DIEP(I-2))
     C      /DIEP(I-1))*(DISE(I))**2-ME(I,J,K)*UBYEPS*(1.0
     C      -PRSTY(I+1,J,K))/AREP(J,K)) + P(I+1,J,K)+DISE(I+1)
     C      *(P(I+1,J,K)-P(I+2,J,K))/DIEP(I+1)+0.5*(((P(I+1,J,K)
     C      -P(I+2,J,K))/DIEP(I+1)-(P(I+2,J,K)-P(I+3,J,K))
     C      /DIEP(I+2))/DIEP(I+1))*(DISE(I+1))**2 )
          ENDIF
         ELSEIF(CVTYPE(I,J,K,1).EQ.2 .AND. CVTYPE(I,J,K,3).EQ.0) THEN
          IF(CLVL.EQ.0) THEN
           PEAST(I,J,K)= (DISE(I+1)*P(I,J,K)
     C                   +DISE(I)*P(I+1,J,K))/DIEP(I)
          ELSE         
           CSUBE= DISE(I)*DE(I,J,K)/(MUEFF(I,J,K)*AREP(J,K))
           UBYEPS= BLEND*((CSUBE/PRSTY(I,J,K)-0.5*CONU(I,J,K,3)*(1.0
     C        -ALFAE(I,J,K)))*U(I+1,J,K)+((1.0-CSUBE)/PRSTY(I,J,K)
     C        -0.5*CONU(I,J,K,1)*(1.0+ALFAE(I,J,K)))*U(I,J,K))
     C        +0.5*CONU(I,J,K,1)*(1.0+ALFAE(I,J,K))*U(I,J,K)
     C        +0.5*CONU(I,J,K,3)*(1.0-ALFAE(I,J,K))*U(I+1,J,K)
           PEAST(I,J,K)= 0.5*( P(I+1,J,K)+DISE(I+1)*(P(I+1,J,K)
     C      -P(I+2,J,K))/DIEP(I+1)+0.5*(((P(I+1,J,K)-P(I+2,J,K))
     C      /DIEP(I+1)-(P(I+2,J,K)-P(I+3,J,K))/DIEP(I+2))
     C      /DIEP(I+1))*(DISE(I+1))**2-ME(I,J,K)*UBYEPS*(1.0
     C      -PRSTY(I,J,K))/AREP(J,K) + P(I,J,K)+DISE(I)*(P(I,J,K)
     C      -P(I-1,J,K))/DIEP(I-1)+0.5*(((P(I,J,K)-P(I-1,J,K))
     C      /DIEP(I-1)-(P(I-1,J,K)-P(I-2,J,K))/DIEP(I-2))
     C      /DIEP(I-1))*(DISE(I))**2 )    
          ENDIF         
         ENDIF
*
*        Estimate pressure at porous/clear fluid interfaces 
*        perpendicular to y-direction
*                
         IF(CVTYPE(I,J,K,1).EQ.0 .AND. CVTYPE(I,J,K,5).EQ.2) THEN 
          IF(CLVL.EQ.0) THEN
           PNORTH(I,J,K)= (DISN(J+1)*P(I,J,K)
     C                    +DISN(J)*P(I,J+1,K))/DJNP(J)
          ELSE
           CSUBN= DISN(J)*DN(I,J,K)/(MUEFF(I,J,K)*ARNP(I,K))
           VBYEPS= BLEND*((CSUBN/PRSTY(I,J+1,K)-0.5*CONU(I,J,K,5)
     C        *(1.0-ALFAN(I,J,K)))*V(I,J+1,K)+((1.0-CSUBN)
     C        /PRSTY(I,J+1,K)-0.5*CONU(I,J,K,1)*(1.0
     C        +ALFAN(I,J,K)))*V(I,J,K))+0.5*CONU(I,J,K,1)*(1.0
     C        +ALFAN(I,J,K))*V(I,J,K)+0.5*CONU(I,J,K,5)*(1.0
     C        -ALFAN(I,J,K))*V(I,J+1,K)
           PNORTH(I,J,K)= 0.5*( P(I,J,K)+DISN(J)*(P(I,J,K)
     C      -P(I,J-1,K))/DJNP(J-1)+0.5*(((P(I,J,K)-P(I,J-1,K))
     C      /DJNP(J-1)-(P(I,J-1,K)-P(I,J-2,K))/DJNP(J-2))
     C      /DJNP(J-1))*(DISN(J))**2-MN(I,J,K)*VBYEPS*(1.0
     C      -PRSTY(I,J+1,K))/ARNP(I,K) + P(I,J+1,K)+DISN(J+1)
     C      *(P(I,J+1,K)-P(I,J+2,K))/DJNP(J+1)+0.5*(((P(I,J+1,K)
     C      -P(I,J+2,K))/DJNP(J+1)-(P(I,J+2,K)-P(I,J+3,K))
     C      /DJNP(J+2))/DJNP(J+1))*(DISN(J+1))**2 )
          ENDIF
         ELSEIF(CVTYPE(I,J,K,1).EQ.2 .AND. CVTYPE(I,J,K,5).EQ.0) THEN
          IF(CLVL.EQ.0) THEN
           PNORTH(I,J,K)= (DISN(J+1)*P(I,J,K)
     C                    +DISN(J)*P(I,J+1,K))/DJNP(J)
          ELSE
           CSUBN= DISN(J)*DN(I,J,K)/(MUEFF(I,J,K)*ARNP(I,K))
           VBYEPS= BLEND*((CSUBN/PRSTY(I,J,K)-0.5*CONU(I,J,K,5)*(1.0
     C        -ALFAN(I,J,K)))*V(I,J+1,K)+((1.0-CSUBN)/PRSTY(I,J,K)
     C        -0.5*CONU(I,J,K,1)*(1.0+ALFAN(I,J,K)))*V(I,J,K))
     C        +0.5*CONU(I,J,K,1)*(1.0+ALFAN(I,J,K))*V(I,J,K)
     C        +0.5*CONU(I,J,K,5)*(1.0-ALFAN(I,J,K))*V(I,J+1,K)
           PNORTH(I,J,K)= 0.5*( P(I,J+1,K)+DISN(J+1)*(P(I,J+1,K)
     C      -P(I,J+2,K))/DJNP(J+1)+0.5*(((P(I,J+1,K)-P(I,J+2,K))
     C      /DJNP(J+1)-(P(I,J+2,K)-P(I,J+3,K))/DJNP(J+2))
     C      /DJNP(J+1))*(DISN(J+1))**2-MN(I,J,K)*VBYEPS*(1.0
     C      -PRSTY(I,J,K))/ARNP(I,K) + P(I,J,K)+DISN(J)*(P(I,J,K)
     C      -P(I,J-1,K))/DJNP(J-1)+0.5*(((P(I,J,K)-P(I,J-1,K))
     C      /DJNP(J-1)-(P(I,J-1,K)-P(I,J-2,K))/DJNP(J-2))
     C      /DJNP(J-1))*(DISN(J))**2 )
          ENDIF
         ENDIF 
*
*        Estimate pressure at porous/clear fluid interfaces 
*        perpendicular to z-direction
*        
         IF(CVTYPE(I,J,K,1).EQ.0 .AND. CVTYPE(I,J,K,7).EQ.2) THEN 
          IF(CLVL.EQ.0) THEN
           PTOP(I,J,K)= (DIST(K+1)*P(I,J,K)
     C                  +DIST(K)*P(I,J,K+1))/DKTP(K)
          ELSE
           CSUBT= DIST(K)*DT(I,J,K)/(MUEFF(I,J,K)*ARTP(I,J))
           WBYEPS= BLEND*((CSUBT/PRSTY(I,J,K+1)-0.5*CONU(I,J,K,7)
     C        *(1.0-ALFAT(I,J,K)))*W(I,J,K+1)+((1.0-CSUBT)
     C        /PRSTY(I,J,K+1)-0.5*CONU(I,J,K,1)*(1.0
     C        +ALFAT(I,J,K)))*W(I,J,K))+0.5*CONU(I,J,K,1)*(1.0
     C        +ALFAT(I,J,K))*W(I,J,K)+0.5*CONU(I,J,K,7)*(1.0
     C        -ALFAT(I,J,K))*W(I,J,K+1)
           PTOP(I,J,K)= 0.5*( (P(I,J,K)+DIST(K)*(P(I,J,K)
     C      -P(I,J,K-1))/DKTP(K-1)+0.5*(((P(I,J,K)-P(I,J,K-1))
     C      /DKTP(K-1)-(P(I,J,K-1)-P(I,J,K-2))/DKTP(K-2))
     C      /DKTP(K-1))*(DIST(K))**2-MT(I,J,K)*WBYEPS*(1.0
     C      -PRSTY(I,J,K+1))/ARTP(I,J)) + P(I,J,K+1)+DIST(K+1)
     C      *(P(I,J,K+1)-P(I,J,K+2))/DKTP(K+1)+0.5*(((P(I,J,K+1)
     C      -P(I,J,K+2))/DKTP(K+1)-(P(I,J,K+2)-P(I,J,K+3))
     C      /DKTP(K+2))/DKTP(K+1))*(DIST(K+1))**2 )
          ENDIF
         ELSEIF(CVTYPE(I,J,K,1).EQ.2 .AND. CVTYPE(I,J,K,7).EQ.0) THEN
          IF(CLVL.EQ.0) THEN
           PTOP(I,J,K)= (DIST(K+1)*P(I,J,K)
     C                  +DIST(K)*P(I,J,K+1))/DKTP(K)
          ELSE
           CSUBT= DIST(K)*DT(I,J,K)/(MUEFF(I,J,K)*ARTP(I,J))
           WBYEPS= BLEND*((CSUBT/PRSTY(I,J,K)-0.5*CONU(I,J,K,7)*(1.0
     C        -ALFAT(I,J,K)))*W(I,J,K+1)+((1.0-CSUBT)/PRSTY(I,J,K)
     C        -0.5*CONU(I,J,K,1)*(1.0+ALFAT(I,J,K)))*W(I,J,K))
     C        +0.5*CONU(I,J,K,1)*(1.0+ALFAT(I,J,K))*W(I,J,K)
     C        +0.5*CONU(I,J,K,7)*(1.0-ALFAT(I,J,K))*W(I,J,K+1)
           PTOP(I,J,K)= 0.5*( P(I,J,K+1)+DIST(K+1)*(P(I,J,K+1)
     C      -P(I,J,K+2))/DKTP(K+1)+0.5*(((P(I,J,K+1)-P(I,J,K+2))
     C      /DKTP(K+1)-(P(I,J,K+2)-P(I,J,K+3))/DKTP(K+2))
     C      /DKTP(K+1))*(DIST(K+1))**2-MT(I,J,K)*WBYEPS*(1.0
     C      -PRSTY(I,J,K))/ARTP(I,J) + P(I,J,K)+DIST(K)*(P(I,J,K)
     C      -P(I,J,K-1))/DKTP(K-1)+0.5*(((P(I,J,K)-P(I,J,K-1))
     C      /DKTP(K-1)-(P(I,J,K-1)-P(I,J,K-2))/DKTP(K-2))
     C      /DKTP(K-1))*(DIST(K))**2 )             
          ENDIF
         ENDIF         
   5    CONTINUE
  10   CONTINUE
  20  CONTINUE
*
      RETURN
      END

