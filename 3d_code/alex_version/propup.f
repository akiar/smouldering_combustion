*
*     file propup.f
************************************************************************
*
*     Subroutine to update material properties in the porous region
      SUBROUTINE PROPUP(CPS,CPF,KS,KF,MUF,
     C                  ID,JD,KD,IB,IE,JB,JE,KB,KE,TS,TF) 
*
*     Update material properties using correlations developed 
*      by Zanoni et al.
*     
************************************************************************
*
*     Declare Variables
*
      INTEGER ID,JD,KD,IE,IB,JB,JE,KB,KE
      INTEGER I,J,K
      REAL CPS(ID,JD,KD),CPF(ID,JD,KD),KS(ID,JD,KD),KF(ID,JD,KD)
      REAL MUF(ID,JD,KD),TS(ID,JD,KD),TF(ID,JD,KD)
*
      DO 30 I=IB,IE
       DO 20 J=JB,JE
        DO 10 K=KB,KE
         CPS(I,J,K) = 2.49*TS(I,J,K) + 39.06
         CPF(I,J,K) = (-3.0E-5)*TF(I,J,K)**2 + 0.2261*TF(I,J,K) + 940.35
         KS(I,J,K)  = 0.000541*TS(I,J,K) + 0.1044
         KF(I,J,K)  = (-1E-8)*TF(I,J,K)**2 + (8E-5)*TF(I,J,K) + 4.3E-3
         MUF(I,J,K) = (-9E-12)*TF(I,J,K)**2 + (4E-8)*TF(I,J,K) + 6E-6
 10     CONTINUE
 20    CONTINUE
 30   CONTINUE
*
      RETURN
      END