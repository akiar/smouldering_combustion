*
*     file interp.f
************************************************************************
*
      SUBROUTINE INTERP(PHI, PROPERTY,TF,ID,JD,KD,IB,IE,JB,JE,KB,KE)
*
*     Subroutine to interpolate variable PHI based on property table A.4
*         in Incropera & Hewitt
*     PROPERTY    1-Density 2-Specific Heat 3-Viscosity 4-Conductivity
*
************************************************************************
*
      IMPLICIT NONE
      INTEGER PROPERTY,ID,JD,KD,IE,IB,JB,JE,KB,KE,I,J,K
      REAL*8 PHI(ID,JD,KD),TF(ID,JD,KD),TMAX,TMIN,PMAX,PMIN
*
      DO 30 I=IB-1,IE+1
      DO 20 J=JB-1,JE+1
      DO 10 K=KB-1,KE+1
*
       CALL TRANGE(TMAX,TMIN,PMAX,PMIN, PROPERTY,TF(I,J,K))
       PHI(I,J,K) = (TF(I,J,K)-TMIN)/(TMAX-TMIN)*(PMAX-PMIN)+PMIN
 10   CONTINUE
 20   CONTINUE
 30   CONTINUE
*
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE TRANGE(TMAX,TMIN,PMAX,PMIN, PROPERTY,TTEST)
*
*     ALL PROPERTIES ARE FOR AIR AT AMOSPHERIC PRESSURE 
*         TEMPERATURE IN KELVIN
*
*     Determine range of temperature to interpolate over
*     PROPERTY    1-Density 2-Specific Heat 3-Viscosity 4-Conductivity
*     TMAX,TMIN   Temperature range to interpolate between
*     PMAX,PMIN   Property values to interpolate on 
*     TTEST       Known fluid temperature
*
************************************************************************
*
      IMPLICIT NONE
      REAL*8 TTEST,TMAX,TMIN,PMAX,PMIN
      INTEGER PROPERTY
*
      IF (TTEST <= 300) THEN
       TMAX  = 300.0
       TMIN  = 250.0
       IF (PROPERTY == 1) THEN
        PMIN = 1.3947
        PMAX = 1.1614
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.006E+03
        PMAX = 1.007E+03
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 159.6E-07
        PMAX = 184.6E-07
       ELSE 
        PMIN = 22.3E-03
        PMAX = 26.3E-03
       ENDIF 
*
      ELSEIF (TTEST <= 350) THEN
       TMAX  = 350.0
       TMIN  = 300.0
       IF (PROPERTY == 1) THEN
        PMIN = 1.1614
        PMAX = 0.9950
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.007E+03
        PMAX = 1.009E+03
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 184.6E-07
        PMAX = 208.2E-07
       ELSE 
        PMIN = 26.3E-03
        PMAX = 30.0E-03
       ENDIF
*
      ELSEIF (TTEST  <= 400) THEN
       TMAX  = 400.0
       TMIN  = 350.0 
       IF (PROPERTY == 1) THEN
        PMIN = 0.9950
        PMAX = 0.8711
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.009E+03
        PMAX = 1.014E+03
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 208.2E-07
        PMAX = 230.1E-07
       ELSE 
        PMIN = 30.0E-03
        PMAX = 33.8E-03
       ENDIF
*
      ELSEIF (TTEST  <= 450) THEN
       TMAX  = 450.0
       TMIN  = 400.0 
       IF (PROPERTY == 1) THEN
        PMIN = 0.8711
        PMAX = 0.7740
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.014E+03
        PMAX = 1.021E+03
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 230.1E-07
        PMAX = 250.7E-07
       ELSE 
        PMIN = 33.8E-03
        PMAX = 37.3E-03
       ENDIF
*
      ELSEIF (TTEST  <= 500) THEN
       TMAX  = 500.0
       TMIN  = 450.0 
       IF (PROPERTY == 1) THEN
        PMIN = 0.7740
        PMAX = 0.6964
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.021E+03
        PMAX = 1.030E+03
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 250.7E-07
        PMAX = 270.1E-07
       ELSE 
        PMIN = 37.3E-03
        PMAX = 40.7E-03
       ENDIF
*
      ELSEIF (TTEST  <= 550) THEN
       TMAX  = 550.0
       TMIN  = 500.0 
       IF (PROPERTY == 1) THEN
        PMIN = 0.6964
        PMAX = 0.6329
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.030E+03
        PMAX = 1.040E+03
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 270.1E-07
        PMAX = 288.4E-07
       ELSE 
        PMIN = 40.7E-03
        PMAX = 43.9E-03
       ENDIF
*
      ELSEIF (TTEST  <= 600) THEN
       TMAX  = 600.0
       TMIN  = 550.0 
       IF (PROPERTY == 1) THEN
        PMIN = 0.6329
        PMAX = 0.5804
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.040E+03
        PMAX = 1.051E+03
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 288.4E-07
        PMAX = 305.8E-07
       ELSE 
        PMIN = 43.9E-03
        PMAX = 46.9E-03
       ENDIF
*
      ELSEIF (TTEST  <= 650) THEN
       TMAX  = 650.0
       TMIN  = 600.0
       IF (PROPERTY == 1) THEN
        PMIN = 0.5804
        PMAX = 0.5356
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.051E+03
        PMAX = 1.063E+03
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 305.8E-07
        PMAX = 322.5E-07
       ELSE 
        PMIN = 46.9E-03
        PMAX = 49.7E-03
       ENDIF
*
      ELSEIF (TTEST  <= 700) THEN
       TMAX  = 700.0
       TMIN  = 650.0 
       IF (PROPERTY == 1) THEN
        PMIN = 0.5356
        PMAX = 0.4975
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.063E+03
        PMAX = 1.075E+03
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 322.5E-07
        PMAX = 338.8E-07
       ELSE 
        PMIN = 49.7E-03
        PMAX = 52.4E-03
       ENDIF
*
      ELSEIF (TTEST <= 750) THEN 
       TMAX  = 750.0
       TMIN  = 700.0
       IF (PROPERTY == 1) THEN
        PMIN = 0.4975
        PMAX = 0.4643
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.075E+03
        PMAX = 1.087E+03
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 338.8E-07
        PMAX = 354.6E-07
       ELSE 
        PMIN = 52.4E-03
        PMAX = 54.9E-03
       ENDIF 
      ELSEIF (TTEST <= 800) THEN
       TMAX  = 800.0
       TMIN  = 750.0
       IF (PROPERTY == 1) THEN
        PMIN = 0.4643
        PMAX = 0.4354
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.087E+03
        PMAX = 1.099E+03
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 354.6E-07
        PMAX = 369.8E-07
       ELSE 
        PMIN = 54.9E-03
        PMAX = 57.3E-03
       ENDIF
      ELSEIF (TTEST <= 850) THEN
       TMAX  = 850.0
       TMIN  = 800.0
       IF (PROPERTY == 1) THEN
        PMIN = 0.4354
        PMAX = 0.4097
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.099E+03
        PMAX = 1.110E+03
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 384.3E-07
        PMAX = 369.8E-07
       ELSE 
        PMIN = 59.6E-03
        PMAX = 57.3E-03
       ENDIF 
      ELSEIF (TTEST <= 900) THEN
       TMAX  = 900.0
       TMIN  = 850.0
       IF (PROPERTY == 1) THEN
        PMIN = 0.4097
        PMAX = 0.3868
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.110E+03
        PMAX = 1.121E+03
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 384.3E-07
        PMAX = 398.1E-07
       ELSE 
        PMIN = 59.6E-03
        PMAX = 62.0E-03
       ENDIF  
      ELSEIF (TTEST <= 950) THEN
       TMAX  = 950.0
       TMIN  = 900.0
       IF (PROPERTY == 1) THEN
        PMIN = 0.3868
        PMAX = 0.3666
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.121E+03
        PMAX = 1.131E+03
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 398.1E-07
        PMAX = 411.3E-07
       ELSE 
        PMIN = 62.0E-03
        PMAX = 64.3E-03
       ENDIF 
      ELSEIF (TTEST <= 1000) THEN
       TMAX  = 1000.0
       TMIN  = 950.0
       IF (PROPERTY == 1) THEN
        PMIN = 0.3666
        PMAX = 0.3482
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.131E+03
        PMAX = 1.141E+03
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 411.3E-07
        PMAX = 424.4E-07
       ELSE 
        PMIN = 64.3E-03
        PMAX = 66.7E-03
       ENDIF  
      ELSEIF (TTEST <= 1100) THEN
       TMAX  = 1100.0
       TMIN  = 1000.0
       IF (PROPERTY == 1) THEN
        PMIN = 0.3482
        PMAX = 0.3166
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.141E+03
        PMAX = 1.159E+03
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 424.4E-07
        PMAX = 449.0E-07
       ELSE 
        PMIN = 66.7E-03
        PMAX = 71.5E-03
       ENDIF  
      ELSEIF (TTEST <= 1200) THEN
       TMAX  = 1200.0
       TMIN  = 1100.0
       IF (PROPERTY == 1) THEN
        PMIN = 0.3166
        PMAX = 0.2902
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.159E+03
        PMAX = 1.175E+03
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 449.0E-07
        PMAX = 473.0E-07
       ELSE 
        PMIN = 71.5E-03
        PMAX = 76.3E-03
       ENDIF
      ELSEIF (TTEST <= 1300) THEN
       TMAX  = 1300.0
       TMIN  = 1200.0
       IF (PROPERTY == 1) THEN
        PMIN = 0.2902
        PMAX = 0.2679
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.175E+03
        PMAX = 1.189E+03
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 473.0E-07
        PMAX = 496.0E-07
       ELSE 
        PMIN = 76.3E-03
        PMAX = 82.0E-03
       ENDIF 
      ELSEIF (TTEST <= 1400) THEN
       TMAX  = 1400.0
       TMIN  = 1300.0
       IF (PROPERTY == 1) THEN
        PMIN = 0.2679
        PMAX = 0.2488
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.189E+03
        PMAX = 1.207E+03
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 496.0E-07
        PMAX = 530.0E-07
       ELSE 
        PMIN = 82.0E-03
        PMAX = 91.0E-03
       ENDIF 
      ELSEIF (TTEST <= 1500) THEN
       TMAX  = 1500.0
       TMIN  = 1400.0
       IF (PROPERTY == 1) THEN
        PMIN = 0.2488
        PMAX = 0.2322
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.207E+03
        PMAX = 1.230E+03
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 530.0E-07
        PMAX = 557.0E-07
       ELSE 
        PMIN = 91.0E-03
        PMAX = 100.0E-03
       ENDIF 
      ELSEIF (TTEST <= 1600) THEN
       TMAX  = 1600.0
       TMIN  = 1500.0
       IF (PROPERTY == 1) THEN
        PMIN = 0.2322
        PMAX = 0.2177
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.230E+03
        PMAX = 1.248E+03
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 557.0E-07
        PMAX = 584.0E-07
       ELSE 
        PMIN = 100.0E-03
        PMAX = 106.0E-03
       ENDIF
      ELSEIF (TTEST <= 1700) THEN
       TMAX  = 1700.0
       TMIN  = 1600.0
       IF (PROPERTY == 1) THEN
        PMIN = 0.2177
        PMAX = 0.2049
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.248E+03
        PMAX = 1.267E+03
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 584.0E-07
        PMAX = 611.0E-07
       ELSE 
        PMIN = 106.0E-03
        PMAX = 113.0E-03
       ENDIF
      ELSEIF (TTEST <= 1800) THEN
       TMAX  = 1800.0
       TMIN  = 1700.0
       IF (PROPERTY == 1) THEN
        PMIN = 0.2049
        PMAX = 0.1935
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.267E+03
        PMAX = 1.286E+03
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 611.0E-07
        PMAX = 637.0E-07
       ELSE 
        PMIN = 113.0E-03
        PMAX = 120.0E-03
       ENDIF 
      ELSE
       TMAX  = 1900.0
       TMIN  = 1800.0
       IF (PROPERTY == 1) THEN
        PMIN = 0.1935
        PMAX = 0.1833
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.286E+03
        PMAX = 1.307E+03
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 637.0E-07
        PMAX = 663.0E-07
       ELSE 
        PMIN = 120.0E-03
        PMAX = 128.0E-03
       ENDIF 
*
      ENDIF
*
      RETURN
      END