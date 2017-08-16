*     file cfppty.f
***********************************************************************
*
      SUBROUTINE CFPPTY(PRSTY,SLDTY,DEQ,CFORCH,KPERM,MUEFF,
     C                  KEFF0,KEFFX,KEFFY,KEFFZ,SPECSA, 
     C                  EPS,PD,DEPSDZ,VISC,VISCB,
     C                  PERM,FORCH,KS,COND,CONDS,ZP,GRDZ,PI,
     C                  CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
*
*     Routine to create property arrays for carbon foam with variable
*     porosity.
*
*     PRSTY(ID,JD,KD)  porosity array; output
*     SLDTY(ID,JD,KD)  solidity array (SLDTY = 1 - PRSTY); output
*     DEQ(ID,JD,KD)  equivalent solid particle diameter array; output
*     KPERM(ID,JD,KD)  permeability array; output
*     CFORCH(ID,JD,KD)  Forchheimer coefficient array; output
*     MUEFF(ID,JD,KD)  effective viscosity array; output
*     KEFF0(ID,JD,KD)  stagnant effective fluid conductivity array; 
*                      output
*     KEFFX(ID,JD,KD)  solid phase effective conductivity in x; output
*     KEFFY(ID,JD,KD)  solid phase effective conductivity in y; output
*     KEFFZ(ID,JD,KD)  solid phase effective conductivity in z; output
*     SPECSA(ID,JD,KD)  specific surface area array; output
*     EPS  mean foam porosity; input
*     PD  pore diameter of porous medium; input
*     DEPSDZ  porosity gradient in z-direction; input
*     VISC  fluid viscosity; input
*     VISCB  Brinkman effective viscosity; input
*     RHO  fluid density; input
*     CP  fluid specific heat capacity; input
*     PERM  coefficient for permeability correlation; input
*     FORCH  coefficient for Forchheimer correlation; input
*     KS  material conductivity of solid constituent of foam; input
*     COND  fluid thermal conductivity; input
*     CONDS  solid thermal condctivity in pure solid region; input
*     ZP(KD)  z-locations of P nodes; input
*     GRDZ  Total domain height in z-direction; input
*     CVTYPE(ID,JD,KD,NNB+1)  indicates type of each CV; input
*     IB,IE  first and last interior indices in i; input
*     JB,JE  first and last interior indices in j; input
*     KB,KE  first and last interior indices in k; input
*     ID,JD,KD  array dimensions; input
*     NNB  number of neighbors for interior nodes; input
*
*     Note: As presently constructed, routine is applicable for 
*           constant porosity, DEPSDZ = 0, or for variable porosity
*           in z direction with porous region height equal to full
*           height of domain.
*
***********************************************************************
*
      IMPLICIT NONE
      REAL*8 PRSTY(ID,JD,KD),SLDTY(ID,JD,KD),DEQ(ID,JD,KD)
      REAL*8 CFORCH(ID,JD,KD),KPERM(ID,JD,KD),MUEFF(ID,JD,KD)
      REAL*8 KEFF0(ID,JD,KD),SPECSA(ID,JD,KD)
      REAL*8 KEFFX(ID,JD,KD),KEFFY(ID,JD,KD),KEFFZ(ID,JD,KD)
      REAL*8 EPS,PD,DEPSDZ,VISC,VISCB
      REAL*8 PERM,FORCH,KS,COND,CONDS,ZP(KD),GRDZ,PI
      REAL*8 PDBYA,PDBYAO,DELTA,TOL,A,BBYA,BBYAO
      PARAMETER(TOL=1.0E-10)
      INTEGER IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB,CVTYPE(ID,JD,KD,NNB+1)
      INTEGER I,J,K
*
      DO 15 I=IB,IE
       DO 10 J=JB,JE
        DO 5 K=KB,KE
*
*        Fluid CVs
*         
         IF(CVTYPE(I,J,K,1).EQ.0) THEN
          PRSTY(I,J,K) = 1.0
          SLDTY(I,J,K) = 0.0
          MUEFF(I,J,K) = VISC
          KEFF0(I,J,K) = COND
          KEFFX(I,J,K) = COND          
          KEFFY(I,J,K) = COND          
          KEFFZ(I,J,K) = COND          
*
*        Porous CVs
*          
         ELSEIF(CVTYPE(I,J,K,1).EQ.2) THEN
*
*         Find porosity, assuming porosity is a linear function of z.
*         
          PRSTY(I,J,K) = DEPSDZ*(ZP(K)-0.5*GRDZ)+EPS
          SLDTY(I,J,K) = 1.0 - PRSTY(I,J,K)
*
*         Find effective viscosity.
*
          MUEFF(I,J,K) = VISCB/PRSTY(I,J,K)        
*
*         Find ratio of unit cube half-height to pore diameter using
*         Newton's method.
* 
          PDBYAO = 2.0
    2     CONTINUE
          PDBYA = PDBYAO-((PDBYAO**3)/24.0-3.0*(PDBYAO**2)/16.0
     C            +0.25+PRSTY(I,J,K)/PI)/((PDBYAO**2)/8.0
     C            -3.0*PDBYAO/8.0)
          DELTA = ABS(PDBYA-PDBYAO)
          IF(DELTA.GE.TOL) THEN
           PDBYAO = PDBYA
           GOTO 2          
          ENDIF
*
*         Find specific surface area.
*
          A = PD/PDBYA
          SPECSA(I,J,K) = 0.25*PI*(3.0*PDBYA-PDBYA**2)/A
*
*         Find equivalent diameter.
*
          DEQ(I,J,K) = 6.0*(1.0-PRSTY(I,J,K))/SPECSA(I,J,K)
*
*         Find Forchheimer coefficient and permeability.
*
          CFORCH(I,J,K) = FORCH/((PERM**0.5)*(PRSTY(I,J,K))**1.5)
          KPERM(I,J,K) = ((PRSTY(I,J,K))**3)*((DEQ(I,J,K))**2)
     C                   /(PERM*(1.0-PRSTY(I,J,K))**2)     
*
*         Find fluid phase stagnant effective conductivity.
*
          BBYAO = 0.1
    3     CONTINUE
          BBYA = BBYAO-(BBYAO**3-1.5*BBYAO**2+0.5*(1.0-PRSTY(I,J,K)))
     C           /(3.0*(BBYAO**2-BBYAO))
          DELTA = ABS(BBYA-BBYAO)
          IF(DELTA.GE.TOL) THEN
           BBYAO = BBYA
           GOTO 3
          ENDIF
*          
          KEFF0(I,J,K) = (1.0-2.0*BBYA+2.0*BBYA**2)*COND*((1.0/BBYA)
     C                   -1.0)**2/(((1.0/BBYA)-1.0)**2+1.0)
*
*         Find solid phase effective conductivity, assuming 
*         conductivity in x and y is half that in z.
*
          KEFFZ(I,J,K) = (1.0-2.0*BBYA+2.0*BBYA**2)*KS
     C                   /(((1.0/BBYA)-1.0)**2+1.0)
          KEFFX(I,J,K) = 0.5*KEFFZ(I,J,K)
          KEFFY(I,J,K) = 0.5*KEFFZ(I,J,K)          
*
*        Solid CVs
*          
         ELSEIF(CVTYPE(I,J,K,1).EQ.3) THEN
          PRSTY(I,J,K) = 0.0
          SLDTY(I,J,K) = 1.0          
          MUEFF(I,J,K) = 0.0
          KEFF0(I,J,K) = CONDS          
          KEFFX(I,J,K) = CONDS          
          KEFFY(I,J,K) = CONDS          
          KEFFZ(I,J,K) = CONDS                                        
         ENDIF
    5   CONTINUE
   10  CONTINUE
   15 CONTINUE
*
*
*     Boundary nodes
*
      DO 30 I=IB,IE
       DO 25 K=KB,KE
        J = JB-1
        PRSTY(I,J,K) = PRSTY(I,J+1,K)
        SLDTY(I,J,K) = SLDTY(I,J+1,K)        
        MUEFF(I,J,K) = MUEFF(I,J+1,K)        
        KEFF0(I,J,K) = KEFF0(I,J+1,K)
        KEFFX(I,J,K) = KEFFX(I,J+1,K)
        KEFFY(I,J,K) = KEFFY(I,J+1,K)
        KEFFZ(I,J,K) = KEFFZ(I,J+1,K)    
        KPERM(I,J,K) = KPERM(I,J+1,K)
*
        J = JE+1        
        PRSTY(I,J,K) = PRSTY(I,J-1,K)
        SLDTY(I,J,K) = SLDTY(I,J-1,K)                
        MUEFF(I,J,K) = MUEFF(I,J-1,K)                
        KEFF0(I,J,K) = KEFF0(I,J-1,K)
        KEFFX(I,J,K) = KEFFX(I,J-1,K)
        KEFFY(I,J,K) = KEFFY(I,J-1,K)
        KEFFZ(I,J,K) = KEFFZ(I,J-1,K)                                
        KPERM(I,J,K) = KPERM(I,J-1,K)
   25  CONTINUE
   30 CONTINUE
      DO 40 J=JB,JE
       DO 35 K=KB,KE
        I = IB-1
        PRSTY(I,J,K) = PRSTY(I+1,J,K)
        SLDTY(I,J,K) = SLDTY(I+1,J,K)                
        MUEFF(I,J,K) = MUEFF(I+1,J,K)                
        KEFF0(I,J,K) = KEFF0(I+1,J,K)
        KEFFX(I,J,K) = KEFFX(I+1,J,K)
        KEFFY(I,J,K) = KEFFY(I+1,J,K)
        KEFFZ(I,J,K) = KEFFZ(I+1,J,K)                                
        KPERM(I,J,K) = KPERM(I+1,J,K)        
*
        I = IE+1        
        PRSTY(I,J,K) = PRSTY(I-1,J,K)
        SLDTY(I,J,K) = SLDTY(I-1,J,K)                
        MUEFF(I,J,K) = MUEFF(I-1,J,K)        
        KEFF0(I,J,K) = KEFF0(I-1,J,K)
        KEFFX(I,J,K) = KEFFX(I-1,J,K)
        KEFFY(I,J,K) = KEFFY(I-1,J,K)
        KEFFZ(I,J,K) = KEFFZ(I-1,J,K)                                
        KPERM(I,J,K) = KPERM(I-1,J,K)
   35  CONTINUE
   40 CONTINUE   
      DO 50 I=IB,IE
       DO 45 J=JB,JE
        K = KB-1
        PRSTY(I,J,K) = PRSTY(I,J,K+1)
        SLDTY(I,J,K) = SLDTY(I,J,K+1)                
        MUEFF(I,J,K) = MUEFF(I,J,K+1)                
        KEFF0(I,J,K) = KEFF0(I,J,K+1)
        KEFFX(I,J,K) = KEFFX(I,J,K+1)
        KEFFY(I,J,K) = KEFFY(I,J,K+1)
        KEFFZ(I,J,K) = KEFFZ(I,J,K+1)                                
        KPERM(I,J,K) = KPERM(I,J,K+1)
*
        K = KE+1        
        PRSTY(I,J,K) = PRSTY(I,J,K-1)
        SLDTY(I,J,K) = SLDTY(I,J,K-1)                
        MUEFF(I,J,K) = MUEFF(I,J,K-1)        
        KEFF0(I,J,K) = KEFF0(I,J,K-1)
        KEFFX(I,J,K) = KEFFX(I,J,K-1)
        KEFFY(I,J,K) = KEFFY(I,J,K-1)
        KEFFZ(I,J,K) = KEFFZ(I,J,K-1)                                
        KPERM(I,J,K) = KPERM(I,J,K-1)
   45  CONTINUE
   50 CONTINUE    
*
      RETURN
      END