*      
*     file tecplt.f
*********************************************************************
*
      SUBROUTINE TECPLT(XP,YP,ZP,IB,IE,JB,JE,KB,KE,ID,JD,KD,
     C                  U,V,W,P,TS,TF,KNTOUT)
*     Subroutine to write solution to tecplot formatted file. 
*      
      IMPLICIT NONE
      INTEGER ID,JD,KD,IB,IE,JB,JE,KB,KE
      INTEGER I,J,K
      REAL*8 XP(ID),YP(JD),ZP(KD)
      REAL*8 U(ID,JD,KD),V(ID,JD,KD),W(ID,JD,KD)
      REAL*8 P(ID,JD,KD),TF(ID,JD,KD),TS(ID,JD,KD)
      CHARACTER*8 FILENUM
      INTEGER KNTOUT
*
      WRITE (FILENUM,'(I8)')KNTOUT
      OPEN(UNIT=20,FILE='tec_'//TRIM(FILENUM)//'.dat')
*
*  Write output in TECplot format for postprocessing
*
      WRITE(20,320)
      WRITE(20,322) 0,0,1,(JE-JB+3),(IE-IB+3)!,(KE-KB+3)
      K=10
      DO I=IB-1,IE+1
       DO J=JB-1,JE+1
        !DO K=KB-1,KE+1
         WRITE(20,325) XP(I),YP(J),P(I,J,K),U(I,J,K),V(I,J,K), !,ZP(K)
     C                 TS(I,J,K),TF(I,J,K)   !W(I,J,K),
        !END DO  
       END DO
      END DO
 320  FORMAT('VARIABLES="XP","YP","P","U","V","TS","TF"') !,"ZP","W"
 322  FORMAT('ZONE T="Zone ',I1,I1,I1,'"I=',I6,' J=',I6,' F=POINT') !' K=',I6,
 325  FORMAT(7(1X,1PE10.3))   !9
      CLOSE(20)
      RETURN
      END