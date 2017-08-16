*
*             file adcont.f
*********************************************************************
*
      SUBROUTINE ADCONT(AP,AW,AE,AS,AN,AB,AT,B,
     C                  ACUW,ACUE,ACVS,ACVN,ACWB,ACWT,
     C                  AEAST,ANORTH,ATOP,DHUE,DHVN,DHWT,
     C                  DPDXP,DPDXM,DPDYP,DPDYM,DPDZP,DPDZM,
     C                  P,PEAST,PNORTH,PTOP,RHO,
     C                  VOLP,VOLE,VOLN,VOLT,PRSTY,
     C                  DIEP,DISE,DJNP,DISN,DKTP,DIST,
     C                  CVTYPE,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB)
*
*     Subroutine to construct coefficients for all terms in the
*     conservation of mass equation for fluid and porous volumes.  
*
*     The coefficients are stored in the same matrix arrays as the
*     other equations.  
*
*     Note that for each face, three possibilities are considered:
*
*         1) The neighboring volume is a solid volume, or the
*            face is a boundary face.
*         2) The neighboring volume is of a different type (i.e. the
*            face is a pure fluid/porous interface).
*         3) The volumes are of the same type.
*
*********************************************************************
*
      IMPLICIT NONE
      REAL*8 AP(N,N,ID,JD,KD),AW(N,N,ID,JD,KD),AE(N,N,ID,JD,KD)
      REAL*8 AS(N,N,ID,JD,KD),AN(N,N,ID,JD,KD),AB(N,N,ID,JD,KD)
      REAL*8 AT(N,N,ID,JD,KD),B(N,ID,JD,KD)
      REAL*8 ACUW(ID,JD,KD),ACUE(ID,JD,KD)
      REAL*8 ACVS(ID,JD,KD),ACVN(ID,JD,KD)
      REAL*8 ACWB(ID,JD,KD),ACWT(ID,JD,KD)
      REAL*8 AEAST(ID,JD,KD),ANORTH(ID,JD,KD),ATOP(ID,JD,KD) 
      REAL*8 DHUE(ID,JD,KD),DHVN(ID,JD,KD),DHWT(ID,JD,KD)                 
      REAL*8 DPDXP(ID,JD,KD),DPDXM(ID,JD,KD)
      REAL*8 DPDYP(ID,JD,KD),DPDYM(ID,JD,KD)
      REAL*8 DPDZP(ID,JD,KD),DPDZM(ID,JD,KD)            
      REAL*8 P(ID,JD,KD),PEAST(ID,JD,KD),PNORTH(ID,JD,KD)
      REAL*8 PTOP(ID,JD,KD),RHO
      REAL*8 VOLP(ID,JD,KD),VOLE(ID,JD,KD),VOLN(ID,JD,KD)
      REAL*8 VOLT(ID,JD,KD),PRSTY(ID,JD,KD)
      REAL*8 DIEP(ID),DISE(ID),DJNP(JD),DISN(JD),DKTP(KD),DIST(KD)
      REAL*8 PCORR,PFACT,FACEPS
      INTEGER IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB
      INTEGER CVTYPE(ID,JD,KD,NNB+1),I,J,K
*  
      DO 3 I=IB,IE
      DO 2 J=JB,JE
      DO 1 K=KB,KE
*
*--Add components for west face of control-volume
*
       IF(CVTYPE(I,J,K,1).NE.3) THEN
*       
       IF(CVTYPE(I,J,K,2).EQ.1 .OR. CVTYPE(I,J,K,2).EQ.3) THEN
	AW(1,2,I,J,K)= AW(1,2,I,J,K)+ACUW(I,J,K)
       ELSEIF(CVTYPE(I,J,K,1).NE.CVTYPE(I,J,K,2)) THEN
        AP(1,2,I,J,K)= AP(1,2,I,J,K)
     C         -ACUW(I,J,K)*0.5*AP(2,2,I,J,K)/AEAST(I-1,J,K)
        AW(1,2,I,J,K)= AW(1,2,I,J,K)
     C         +ACUW(I,J,K)*0.5*AP(2,2,I-1,J,K)/AEAST(I-1,J,K)
        PFACT= (PRSTY(I-1,J,K)*VOLP(I-1,J,K)
     C    +PRSTY(I,J,K)*VOLP(I,J,K))/(2.0*VOLE(I-1,J,K)*DIEP(I-1))
        PCORR= (PRSTY(I-1,J,K)*VOLP(I-1,J,K)*(PEAST(I-1,J,K)
     C         -P(I-1,J,K))/DISE(I-1)+PRSTY(I,J,K)*VOLP(I,J,K)
     C         *(P(I,J,K)-PEAST(I-1,J,K))/DISE(I))
     C         /(2.0*VOLE(I-1,J,K))-PFACT*(P(I,J,K)-P(I-1,J,K))
        AP(1,1,I,J,K)= AP(1,1,I,J,K)
     C                 +ACUW(I,J,K)*DHUE(I-1,J,K)*PFACT
        AW(1,1,I,J,K)= AW(1,1,I,J,K)
     C                 +ACUW(I,J,K)*DHUE(I-1,J,K)*PFACT
        B(1,I,J,K)= B(1,I,J,K)+ACUW(I,J,K)*( DHUE(I-1,J,K)
     C       *((VOLP(I-1,J,K)*PRSTY(I-1,J,K)*DPDXP(I-1,J,K)
     C       +VOLP(I,J,K)*PRSTY(I,J,K)*DPDXM(I,J,K))
     C       /(2.0*VOLE(I-1,J,K)))-DHUE(I-1,J,K)*PCORR )
       ELSE
        FACEPS = (DISE(I-1)*PRSTY(I,J,K)
     C           +DISE(I)*PRSTY(I-1,J,K))/DIEP(I-1)
        AP(1,1,I,J,K)= AP(1,1,I,J,K)
     C         +ACUW(I,J,K)*FACEPS*DHUE(I-1,J,K)/DIEP(I-1)
        AW(1,1,I,J,K)= AW(1,1,I,J,K)
     C         +ACUW(I,J,K)*FACEPS*DHUE(I-1,J,K)/DIEP(I-1)
        AP(1,2,I,J,K)= AP(1,2,I,J,K)
     C                 -ACUW(I,J,K)*DISE(I-1)/DIEP(I-1)
        AW(1,2,I,J,K)= AW(1,2,I,J,K)
     C                 +ACUW(I,J,K)*DISE(I)/DIEP(I-1)
        B(1,I,J,K)= B(1,I,J,K)+ACUW(I,J,K)*( DHUE(I-1,J,K)
     C       *(DISE(I)*PRSTY(I-1,J,K)*DPDXP(I-1,J,K)
     C       +DISE(I-1)*PRSTY(I,J,K)*DPDXM(I,J,K))/DIEP(I-1) )
       ENDIF
*
*--Add components for east face of control-volume
*
       IF(CVTYPE(I,J,K,3).EQ.1 .OR. CVTYPE(I,J,K,3).EQ.3) THEN
	    AE(1,2,I,J,K)= AE(1,2,I,J,K)+ACUE(I,J,K)
       ELSEIF(CVTYPE(I,J,K,1).NE.CVTYPE(I,J,K,3)) THEN
        AP(1,2,I,J,K)= AP(1,2,I,J,K)
     C         -ACUE(I,J,K)*0.5*AP(2,2,I,J,K)/AEAST(I,J,K)
        AE(1,2,I,J,K)= AE(1,2,I,J,K)
     C         +ACUE(I,J,K)*0.5*AP(2,2,I+1,J,K)/AEAST(I,J,K)
        PFACT= (PRSTY(I,J,K)*VOLP(I,J,K)
     C    +PRSTY(I+1,J,K)*VOLP(I+1,J,K))/(2.0*VOLE(I,J,K)*DIEP(I))
        PCORR= (PRSTY(I,J,K)*VOLP(I,J,K)*(PEAST(I,J,K)
     C         -P(I,J,K))/DISE(I)+PRSTY(I+1,J,K)*VOLP(I+1,J,K)
     C         *(P(I+1,J,K)-PEAST(I,J,K))/DISE(I+1))
     C         /(2.0*VOLE(I,J,K))-PFACT*(P(I+1,J,K)-P(I,J,K))     
        AP(1,1,I,J,K)= AP(1,1,I,J,K)
     C                 -ACUE(I,J,K)*DHUE(I,J,K)*PFACT
        AE(1,1,I,J,K)= AE(1,1,I,J,K)
     C                 -ACUE(I,J,K)*DHUE(I,J,K)*PFACT
        B(1,I,J,K)= B(1,I,J,K)+ACUE(I,J,K)*( DHUE(I,J,K)
     C       *((VOLP(I,J,K)*PRSTY(I,J,K)*DPDXP(I,J,K)
     C       +VOLP(I+1,J,K)*PRSTY(I+1,J,K)*DPDXM(I+1,J,K))
     C       /(2.0*VOLE(I,J,K)))-DHUE(I,J,K)*PCORR )
       ELSE
        FACEPS = (DISE(I+1)*PRSTY(I,J,K)
     C           +DISE(I)*PRSTY(I+1,J,K))/DIEP(I)
	    AP(1,1,I,J,K)= AP(1,1,I,J,K)
     C         -ACUE(I,J,K)*FACEPS*DHUE(I,J,K)/DIEP(I)
	    AE(1,1,I,J,K)= AE(1,1,I,J,K)
     C         -ACUE(I,J,K)*FACEPS*DHUE(I,J,K)/DIEP(I)
        AP(1,2,I,J,K)= AP(1,2,I,J,K)
     C                 -ACUE(I,J,K)*DISE(I+1)/DIEP(I)
        AE(1,2,I,J,K)= AE(1,2,I,J,K)
     C                 +ACUE(I,J,K)*DISE(I)/DIEP(I)
        B(1,I,J,K)= B(1,I,J,K)+ACUE(I,J,K)*( DHUE(I,J,K)
     C       *(DISE(I+1)*PRSTY(I,J,K)*DPDXP(I,J,K)
     C       +DISE(I)*PRSTY(I+1,J,K)*DPDXM(I+1,J,K))/DIEP(I) )
       ENDIF
*
*--Add components for south face of control-volume
*
       IF(CVTYPE(I,J,K,4).EQ.1 .OR. CVTYPE(I,J,K,4).EQ.3) THEN
        AS(1,3,I,J,K)= AS(1,3,I,J,K)+ACVS(I,J,K)
       ELSEIF(CVTYPE(I,J,K,1).NE.CVTYPE(I,J,K,4)) THEN
        AP(1,3,I,J,K)= AP(1,3,I,J,K)
     C         -ACVS(I,J,K)*0.5*AP(3,3,I,J,K)/ANORTH(I,J-1,K)
        AS(1,3,I,J,K)= AS(1,3,I,J,K)
     C         +ACVS(I,J,K)*0.5*AP(3,3,I,J-1,K)/ANORTH(I,J-1,K)
        PFACT= (PRSTY(I,J-1,K)*VOLP(I,J-1,K)
     C    +PRSTY(I,J,K)*VOLP(I,J,K))/(2.0*VOLN(I,J-1,K)*DJNP(J-1))
        PCORR= (PRSTY(I,J-1,K)*VOLP(I,J-1,K)*(PNORTH(I,J-1,K)
     C         -P(I,J-1,K))/DISN(J-1)+PRSTY(I,J,K)*VOLP(I,J,K)
     C         *(P(I,J,K)-PNORTH(I,J-1,K))/DISN(J))
     C         /(2.0*VOLN(I,J-1,K))-PFACT*(P(I,J,K)-P(I,J-1,K))     
        AP(1,1,I,J,K)= AP(1,1,I,J,K)
     C                 +ACVS(I,J,K)*DHVN(I,J-1,K)*PFACT
        AS(1,1,I,J,K)= AS(1,1,I,J,K)
     C                 +ACVS(I,J,K)*DHVN(I,J-1,K)*PFACT
        B(1,I,J,K)= B(1,I,J,K)+ACVS(I,J,K)*( DHVN(I,J-1,K)
     C       *((VOLP(I,J-1,K)*PRSTY(I,J-1,K)*DPDYP(I,J-1,K)
     C       +VOLP(I,J,K)*PRSTY(I,J,K)*DPDYM(I,J,K))
     C       /(2.0*VOLN(I,J-1,K)))-DHVN(I,J-1,K)*PCORR )
       ELSE
        FACEPS = (DISN(J-1)*PRSTY(I,J,K)
     C           +DISN(J)*PRSTY(I,J-1,K))/DJNP(J-1)
	    AP(1,1,I,J,K)= AP(1,1,I,J,K)
     C         +ACVS(I,J,K)*FACEPS*DHVN(I,J-1,K)/DJNP(J-1)
	    AS(1,1,I,J,K)= AS(1,1,I,J,K)
     C         +ACVS(I,J,K)*FACEPS*DHVN(I,J-1,K)/DJNP(J-1)
        AP(1,3,I,J,K)= AP(1,3,I,J,K)
     C                 -ACVS(I,J,K)*DISN(J-1)/DJNP(J-1)
        AS(1,3,I,J,K)= AS(1,3,I,J,K)
     C                 +ACVS(I,J,K)*DISN(J)/DJNP(J-1)
        B(1,I,J,K)= B(1,I,J,K)+ACVS(I,J,K)*( DHVN(I,J-1,K)
     C       *(DISN(J)*PRSTY(I,J-1,K)*DPDYP(I,J-1,K)
     C       +DISN(J-1)*PRSTY(I,J,K)*DPDYM(I,J,K))/DJNP(J-1) )
       ENDIF
*
*--Add components for north face of control-volume
*
       IF(CVTYPE(I,J,K,5).EQ.1 .OR. CVTYPE(I,J,K,5).EQ.3) THEN
	    AN(1,3,I,J,K)= AN(1,3,I,J,K)+ACVN(I,J,K)
       ELSEIF(CVTYPE(I,J,K,1).NE.CVTYPE(I,J,K,5)) THEN 
        AP(1,3,I,J,K)= AP(1,3,I,J,K)
     C         -ACVN(I,J,K)*0.5*AP(3,3,I,J,K)/ANORTH(I,J,K)
        AN(1,3,I,J,K)= AN(1,3,I,J,K)
     C         +ACVN(I,J,K)*0.5*AP(3,3,I,J+1,K)/ANORTH(I,J,K)
        PFACT= (PRSTY(I,J,K)*VOLP(I,J,K)
     C    +PRSTY(I,J+1,K)*VOLP(I,J+1,K))/(2.0*VOLN(I,J,K)*DJNP(J))
        PCORR= (PRSTY(I,J,K)*VOLP(I,J,K)*(PNORTH(I,J,K)
     C         -P(I,J,K))/DISN(J)+PRSTY(I,J+1,K)*VOLP(I,J+1,K)
     C         *(P(I,J+1,K)-PNORTH(I,J,K))/DISN(J+1))
     C         /(2.0*VOLN(I,J,K))-PFACT*(P(I,J+1,K)-P(I,J,K))
        AP(1,1,I,J,K)= AP(1,1,I,J,K)
     C                 -ACVN(I,J,K)*DHVN(I,J,K)*PFACT
        AN(1,1,I,J,K)= AN(1,1,I,J,K)
     C                 -ACVN(I,J,K)*DHVN(I,J,K)*PFACT
        B(1,I,J,K)= B(1,I,J,K)+ACVN(I,J,K)*( DHVN(I,J,K)
     C       *((VOLP(I,J,K)*PRSTY(I,J,K)*DPDYP(I,J,K)
     C       +VOLP(I,J+1,K)*PRSTY(I,J+1,K)*DPDYM(I,J+1,K))
     C       /(2.0*VOLN(I,J,K)))-DHVN(I,J,K)*PCORR )
       ELSE
        FACEPS = (DISN(J+1)*PRSTY(I,J,K)
     C           +DISN(J)*PRSTY(I,J+1,K))/DJNP(J)
	    AP(1,1,I,J,K)= AP(1,1,I,J,K)
     C         -ACVN(I,J,K)*FACEPS*DHVN(I,J,K)/DJNP(J)
	    AN(1,1,I,J,K)= AN(1,1,I,J,K)
     C         -ACVN(I,J,K)*FACEPS*DHVN(I,J,K)/DJNP(J)
        AP(1,3,I,J,K)= AP(1,3,I,J,K)
     C                 -ACVN(I,J,K)*DISN(J+1)/DJNP(J)
        AN(1,3,I,J,K)= AN(1,3,I,J,K)
     C                 +ACVN(I,J,K)*DISN(J)/DJNP(J)
        B(1,I,J,K)= B(1,I,J,K)+ACVN(I,J,K)*( DHVN(I,J,K)
     C       *(DISN(J+1)*PRSTY(I,J,K)*DPDYP(I,J,K)
     C       +DISN(J)*PRSTY(I,J+1,K)*DPDYM(I,J+1,K))/DJNP(J) )       
       ENDIF
*
*--Add components for bottom face of control-volume
*
       IF(CVTYPE(I,J,K,6).EQ.1 .OR. CVTYPE(I,J,K,6).EQ.3) THEN
	AB(1,4,I,J,K)= AB(1,4,I,J,K)+ACWB(I,J,K)
       ELSEIF(CVTYPE(I,J,K,1).NE.CVTYPE(I,J,K,6)) THEN
        AP(1,4,I,J,K)= AP(1,4,I,J,K)
     C         -ACWB(I,J,K)*0.5*AP(4,4,I,J,K)/ATOP(I,J,K-1)
        AB(1,4,I,J,K)= AB(1,4,I,J,K)
     C         +ACWB(I,J,K)*0.5*AP(4,4,I,J,K-1)/ATOP(I,J,K-1)
        PFACT= (PRSTY(I,J,K-1)*VOLP(I,J,K-1)
     C    +PRSTY(I,J,K)*VOLP(I,J,K))/(2.0*VOLT(I,J,K-1)*DKTP(K-1))
        PCORR= (PRSTY(I,J,K-1)*VOLP(I,J,K-1)*(PTOP(I,J,K-1)
     C         -P(I,J,K-1))/DIST(K-1)+PRSTY(I,J,K)*VOLP(I,J,K)
     C         *(P(I,J,K)-PTOP(I,J,K-1))/DIST(K))
     C         /(2.0*VOLT(I,J,K-1))-PFACT*(P(I,J,K)-P(I,J-1,K))     
        AP(1,1,I,J,K)= AP(1,1,I,J,K)
     C                 +ACWB(I,J,K)*DHWT(I,J,K-1)*PFACT
        AB(1,1,I,J,K)= AB(1,1,I,J,K)
     C                 +ACWB(I,J,K)*DHWT(I,J,K-1)*PFACT
        B(1,I,J,K)= B(1,I,J,K)+ACWB(I,J,K)*( DHWT(I,J,K-1)
     C       *((VOLP(I,J,K-1)*PRSTY(I,J,K-1)*DPDZP(I,J,K-1)
     C       +VOLP(I,J,K)*PRSTY(I,J,K)*DPDZM(I,J,K))
     C       /(2.0*VOLT(I,J,K-1)))-DHWT(I,J,K-1)*PCORR )
       ELSE
        FACEPS = (DIST(K-1)*PRSTY(I,J,K)
     C           +DIST(K)*PRSTY(I,J,K-1))/DKTP(K-1)
	AP(1,1,I,J,K)= AP(1,1,I,J,K)
     C         +ACWB(I,J,K)*FACEPS*DHWT(I,J,K-1)/DKTP(K-1)
	AB(1,1,I,J,K)= AB(1,1,I,J,K)
     C         +ACWB(I,J,K)*FACEPS*DHWT(I,J,K-1)/DKTP(K-1)
        AP(1,4,I,J,K)= AP(1,4,I,J,K)
     C                 -ACWB(I,J,K)*DIST(K-1)/DKTP(K-1)
        AB(1,4,I,J,K)= AB(1,4,I,J,K)
     C                 +ACWB(I,J,K)*DIST(K)/DKTP(K-1)
        B(1,I,J,K)= B(1,I,J,K)+ACWB(I,J,K)*( DHWT(I,J,K-1)
     C       *(DIST(K)*PRSTY(I,J,K-1)*DPDZP(I,J,K-1)
     C       +DIST(K-1)*PRSTY(I,J,K)*DPDZM(I,J,K))/DKTP(K-1) )
       ENDIF
*
*--Add components for top face of control-volume
*
       IF(CVTYPE(I,J,K,7).EQ.1 .OR. CVTYPE(I,J,K,7).EQ.3) THEN
	AT(1,4,I,J,K)= AT(1,4,I,J,K)+ACWT(I,J,K)
       ELSEIF(CVTYPE(I,J,K,1).NE.CVTYPE(I,J,K,7)) THEN 
        AP(1,4,I,J,K)= AP(1,4,I,J,K)
     C         -ACWT(I,J,K)*0.5*AP(4,4,I,J,K)/ATOP(I,J,K)
        AT(1,4,I,J,K)= AT(1,4,I,J,K)
     C         +ACWT(I,J,K)*0.5*AP(4,4,I,J,K+1)/ATOP(I,J,K)
        PFACT= (PRSTY(I,J,K)*VOLP(I,J,K)
     C    +PRSTY(I,J,K+1)*VOLP(I,J,K+1))/(2.0*VOLT(I,J,K)*DKTP(K))
        PCORR= (PRSTY(I,J,K)*VOLP(I,J,K)*(PTOP(I,J,K)
     C         -P(I,J,K))/DIST(K)+PRSTY(I,J,K+1)*VOLP(I,J,K+1)
     C         *(P(I,J,K+1)-PTOP(I,J,K))/DIST(K+1))
     C         /(2.0*VOLT(I,J,K))-PFACT*(P(I,J,K+1)-P(I,J,K))
        AP(1,1,I,J,K)= AP(1,1,I,J,K)
     C                 -ACWT(I,J,K)*DHWT(I,J,K)*PFACT
        AT(1,1,I,J,K)= AT(1,1,I,J,K)
     C                 -ACWT(I,J,K)*DHWT(I,J,K)*PFACT
        B(1,I,J,K)= B(1,I,J,K)+ACWT(I,J,K)*( DHWT(I,J,K)
     C       *((VOLP(I,J,K)*PRSTY(I,J,K)*DPDZP(I,J,K)
     C       +VOLP(I,J,K+1)*PRSTY(I,J,K+1)*DPDZM(I,J,K+1))
     C       /(2.0*VOLT(I,J,K)))-DHWT(I,J,K)*PCORR )
       ELSE
        FACEPS = (DIST(K+1)*PRSTY(I,J,K)
     C           +DIST(K)*PRSTY(I,J,K+1))/DKTP(K)
	AP(1,1,I,J,K)= AP(1,1,I,J,K)
     C         -ACWT(I,J,K)*FACEPS*DHWT(I,J,K)/DKTP(K)
	AT(1,1,I,J,K)= AT(1,1,I,J,K)
     C         -ACWT(I,J,K)*FACEPS*DHWT(I,J,K)/DKTP(K)
        AP(1,4,I,J,K)= AP(1,4,I,J,K)
     C                 -ACWT(I,J,K)*DIST(K+1)/DKTP(K)
        AT(1,4,I,J,K)= AT(1,4,I,J,K)
     C                 +ACWT(I,J,K)*DIST(K)/DKTP(K)
        B(1,I,J,K)= B(1,I,J,K)+ACWT(I,J,K)*( DHWT(I,J,K)
     C       *(DIST(K+1)*PRSTY(I,J,K)*DPDZP(I,J,K)
     C       +DIST(K)*PRSTY(I,J,K+1)*DPDZM(I,J,K+1))/DKTP(K) )       
       ENDIF
*
       ENDIF
*
  1   CONTINUE
  2   CONTINUE
  3   CONTINUE
*
      RETURN
      END
