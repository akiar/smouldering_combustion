*
************************************************************************
       SUBROUTINE NULLM(A, ISTRT,IFIN,JSTRT,JFIN,KSTRT,KFIN,ID,JD,KD)
*
*      Subroutine to null a matrix.
*
*      A(ID,JD,KD) matrix
*      ISTRT,IFIN  beginning and ending indices in I direction
*      JSTRT,JFIN  beginning and ending indices in J direction
*      KSTRT,KFIN  beginning and ending indices in K direction
*      ID,JD,KD    array dimensions
*
************************************************************************
*
      IMPLICIT NONE
      REAL*8 A(ID,JD,KD)
      INTEGER ISTRT,IFIN,JSTRT,JFIN,KSTRT,KFIN,I,J,K,ID,JD,KD
*
      DO 3 I=ISTRT,IFIN
        DO 2 J=JSTRT,JFIN
          DO 1 K=KSTRT,KFIN
            A(I,J,K)= 0.0
  1       CONTINUE
  2     CONTINUE
  3   CONTINUE
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE NULLMM(A, ISTRT,IFIN,JSTRT,JFIN,KSTRT,KFIN,
     C                  N,ID,JD,KD)
*
*     Subroutine to null a block matrix.
*
*     A(N,ID,JD,KD)  matrix
*     ISTRT,IFIN  beginning and ending indices in I direction
*     JSTRT,JFIN  beginning and ending indices in J direction
*     KSTRT,KFIN  beginning and ending indices in K direction
*     ID,JD,KD    array dimensions
*
************************************************************************
*
      IMPLICIT NONE
      REAL*8 A(N,ID,JD,KD)
      INTEGER ISTRT,IFIN,JSTRT,JFIN,KSTRT,KFIN,I,J,K,L,N,ID,JD,KD
*
      DO 4 K=KSTRT,KFIN
       DO 3 I=ISTRT,IFIN
        DO 2 J=JSTRT,JFIN
         DO 1 L=1,N
           A(L,I,J,K)= 0.0
  1      CONTINUE
  2     CONTINUE
  3    CONTINUE
  4   CONTINUE
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE NULLV(A, ISTRT,IFIN,ID)
*
*     Subroutine to null a vector.
*
*     A(ID)     vector
*     ISRT,IFIN beginning and ending indices
*     ID        array dimension
************************************************************************
*
      IMPLICIT NONE
      REAL*8 A(ID)
      INTEGER ISTRT,IFIN,I,ID
*
      DO 1 I=ISTRT,IFIN
        A(I)= 0.0
  1   CONTINUE
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE NULLMN(A, ISTRT,IFIN,JSTRT,JFIN,KSTRT,KFIN,
     C                  N,ID,JD,KD)
*
*      Subroutine to null a block matrix.
*
*      A(N,N,ID,JD,KD)    matrix
*      ISTRT,IFIN  beginning and ending indices in I direction
*      JSTRT,JFIN  beginning and ending indices in J direction
*      KSTRT,KFIN  beginning and ending indices in K direction
*      ID,JD,KD    array dimensions
*
************************************************************************
*
      IMPLICIT NONE
      REAL*8 A(N,N,ID,JD,KD)
      INTEGER ISTRT,IFIN,JSTRT,JFIN,KSTRT,KFIN,N,ID,JD,KD
      INTEGER I,J,K,L,M
*
      DO 5 L=1,N
       DO 4 M=1,N
        DO 3 I=ISTRT,IFIN
         DO 2 J=JSTRT,JFIN
          DO 1 K=KSTRT,KFIN
            A(L,M,I,J,K)= 0.0
  1       CONTINUE
  2      CONTINUE
  3     CONTINUE
  4    CONTINUE
  5   CONTINUE
*
      RETURN
      END

