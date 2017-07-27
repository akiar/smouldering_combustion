*
*     file htchar.f
*     
***********************************************************************
*
      SUBROUTINE HTCHAR(ID,JD,KD,IB,IE,JB,JE,KB,KE,GRDX,YP,
     C                  TF,TS,TIN,
     C                  CONDFE,CONDSE)
*
*     HEATFLUX    Heatflux through the first cells in x,y, or z.
*     NUSSELT     Nusselt number
*     TAVG        Average temperature over heatflux area
*     DTAVG       Difference between TIN and TAVG
*
***********************************************************************
*
*     Declare variables
*
      IMPLICIT NONE
      INTEGER ID,JD,KD,IB,IE,JB,JE,KB,KE,GRDX
      INTEGER YP(ID),I
      REAL*8 HFS(ID),HFF(ID),HF,NUSSELT      
      REAL*8 CONDFE,CONDSE                   
      REAL*8 TIN,TF(ID,JD,KD),TS(ID,JD,KD)   
      REAL*8 DY,DTDYF(ID),DTDYS(ID),WALLAV
*      
*     Calculate dy 
*
      DY = YP(11) - YP(JB-1) !SET Y LOCATION
*
*     Calculate fluid heat flux
*
      CALL TGRAD(DTDYF, TF,ID,JD,KD,IB,IE,JB,JE,KB,KE,DY)
      DO 1000 I=IB,IE
       HFF(I) = CONDFE*GRDX*DTDYF(I)
 1000 CONTINUE  
*
*     Calculate solid heat flux
*
      CALL TGRAD(DTDYS, TS,ID,JD,KD,IB,IE,JB,JE,KB,KE,DY)
      DO 1100 I=IB,IE
       HFS(I) = CONDSE*GRDX*DTDYS(I)
 1100 CONTINUE      
*
*     Calculate total heat flux
*
      DO 1200 I=IB,IE
       HF = HFS(I) + HFF(I) + HF
 1200 CONTINUE      
*
*     Calculate Nusselt Number
*      
      CALL WALLDIFF(WALLAV, TF,ID,JD,KD,IB,IE,JB,JE,KB,KE,TIN)
      NUSSELT = HF/CONDFE/WALLAV
*
*     Print results
*      
      
*
      RETURN
      END
*
***********************************************************************
* 
      SUBROUTINE TGRAD(DTDY, T,ID,JD,KD,IB,IE,JB,JE,KB,KE,DY)
*
*     T   Temperature solution 
*     Y   yp locations of the domain
*
***********************************************************************
*
*     Declare Variables
*
      INTEGER ID,JD,KD,I,J,K,IB,IE,JB,JE,KB,KE
      REAL*8 T(ID,JD,KD),DTDY(ID),DY
*
*     Calculate average tempearture
*
      CALL AVERAGET(TAVG, T,ID,JD,KD,IB,IE,JB,JE,KB,KE)
*
*     Calculate temperature gradient
*     
      DO 2000 I=IB,IE
       DTDY(I) = (T(I,JB-1,KB)-TAVG)/DY
 2000 CONTINUE
*
      RETURN 
      END
*
***********************************************************************
* 
      SUBROUTINE AVERAGET(TAVG, T,ID,JD,KD,IB,IE,JB,JE,KB,KE)
*
*     Subroutine to calculate the average temperature at a height y in the
*         domain
*
*     T   Temperature solution 
*
*********************************************************************** 
*
      INTEGER IB,IE,JB,JE,KB,KE,ID,JD,KD,I,J,K
      REAL*8 TAVG,SUMT,T(ID,JD,KD)
*
      J=11
      K=KB
      SUMT = 0.0 
*
      DO 3000 I=IB,IE
       SUMT = SUMT + T(I,J,K)
 3000 CONTINUE
      TAVG = SUMT/(IE-IB)
*
      RETURN
      END
*
***********************************************************************
* 
      SUBROUTINE WALLDIFF(WALLAV, T,ID,JD,KD,IB,IE,JB,JE,KB,KE,TIN)
*
*********************************************************************** 
*
      INTEGER ID,JD,KD,IB,IE,JB,JE,KB,KE
      REAL*8 T(ID,JD,KD),WALLAVG,WALLAV
*
      CALL AVERAGET(WALLAVG, T,ID,JD,KD,IB,IE,JB,JE,KB,KE)
      WALLAV = WALLAVG - TIN
*
      RETURN
      END