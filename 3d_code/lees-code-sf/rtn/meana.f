*
*       file meana.f
*****************************************************************************
*
      SUBROUTINE MEANA(MEAN, 
     C           PHI,FL,POR,SOL,CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
*
*     Subroutine to find the magnitude of the average of PHI.
*
*****************************************************************************
*
      IMPLICIT NONE
      REAL*8 MEAN,PHI(ID,JD,KD),SUM
      INTEGER FL,POR,SOL,IB,IE,JB,JE,KB,KE,NPTS,I,J,K,ID,JD,KD,NNB
      INTEGER CVTYPE(ID,JD,KD,NNB+1)
*
      SUM= 0.0
      NPTS= 0
*
      DO 3 I=IB,IE
       DO 2 J=JB,JE
        DO 1 K=KB,KE
         IF(CVTYPE(I,J,K,1).EQ.0) THEN
          SUM= SUM+DBLE(FL)*ABS(PHI(I,J,K))
          NPTS= NPTS+FL
         ELSEIF(CVTYPE(I,J,K,1).EQ.2) THEN
          SUM= SUM+DBLE(POR)*ABS(PHI(I,J,K))
          NPTS= NPTS+POR
         ELSEIF(CVTYPE(I,J,K,1).EQ.3) THEN
          SUM= SUM+DBLE(SOL)*ABS(PHI(I,J,K))
          NPTS= NPTS+SOL
         ENDIF
  1     CONTINUE
  2    CONTINUE
  3   CONTINUE
      MEAN= SUM/DBLE(NPTS)
      IF(ABS(MEAN).GE.1.0E-20) THEN
       MEAN= MEAN
      ELSE
       MEAN= 1.0
      ENDIF
*
      RETURN
      END
