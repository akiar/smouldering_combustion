*
*     file interp.f
************************************************************************
*
      SUBROUTINE INTERP(PHI, PROPERTY,TF,ID,JD,KD,IE,IB,JB,JE,KB,KE)
*
*     Subroutine to interpolate variable PHI based on property table A.4
*         in Incropera & Hewitt
*     PROPERTY    1-Density 2-Specific Heat 3-Viscosity 4-Conductivity
*
************************************************************************
*
      INTEGER PROPERTY,ID,JD,KD,IE,IB,JB,JE,KB,KE,I,J,K
      REAL PHI(ID,JD,KD),TF(ID,JD,KD),TMAX,TMIN
      
*
      DO 30 I=IB,IE
      DO 20 J=JB,JE
      DO 10 K=KB,KE
*
       CALL TRANGE(TMAX,TMIN,PMAX,PMIN, PROPERTY,TF(I,J,K))
       PHI(I,J,K) = (TF(I,J,K)-TMIN)/(TMAX-TMIN)*(PMAX-PMIN)+PMIN
*
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
      REAL TTEST,TMAX,TMIN,PMAX,PMIN
      INTEGER PROPERTY
*
      IF (TTEST <= 300) THEN
       TMAX  = 300.0
       TMIN  = 250.0
       IF (PROPERTY == 1) THEN
        PMIN = 1.3947
        PMAX = 1.1614
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.006E+3
        PMAX = 1.007E+3
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 159.6E-7
        PMAX = 184.6E-7
       ELSE 
        PMIN = 22.3E-3
        PMAX = 26.3E-3
       ENDIF 
*
      ELSEIF (TTEST <= 350) THEN
       TMAX  = 350.0
       TMIN  = 300.0
       IF (PROPERTY == 1) THEN
        PMIN = 1.1614
        PMAX = 0.9950
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.007E+3
        PMAX = 1.009E+3
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 184.6E-7
        PMAX = 208.2E-7
       ELSE 
        PMIN = 26.3E-3
        PMAX = 30.0E-3
       ENDIF
*
      ELSEIF (TTEST  <= 400) THEN
       TMAX  = 400.0
       TMIN  = 350.0 
       IF (PROPERTY == 1) THEN
        PMIN = 0.9950
        PMAX = 0.8711
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.009E+3
        PMAX = 1.014E+3
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 208.2E-7
        PMAX = 230.1E-7
       ELSE 
        PMIN = 30.0E-3
        PMAX = 33.8E-3
       ENDIF
*
      ELSEIF (TTEST  <= 450) THEN
       TMAX  = 450.0
       TMIN  = 400.0 
       IF (PROPERTY == 1) THEN
        PMIN = 0.8711
        PMAX = 0.7740
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.014E+3
        PMAX = 1.021E+3
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 230.1E-7
        PMAX = 250.7E-7
       ELSE 
        PMIN = 33.8E-3
        PMAX = 37.3E-3
       ENDIF
*
      ELSEIF (TTEST  <= 500) THEN
       TMAX  = 500.0
       TMIN  = 450.0 
       IF (PROPERTY == 1) THEN
        PMIN = 0.7740
        PMAX = 0.6964
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.021E+3
        PMAX = 1.030E+3
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 250.7E-7
        PMAX = 270.1E-7
       ELSE 
        PMIN = 37.3E-3
        PMAX = 40.7E-3
       ENDIF
*
      ELSEIF (TTEST  <= 550) THEN
       TMAX  = 550.0
       TMIN  = 500.0 
       IF (PROPERTY == 1) THEN
        PMIN = 0.6964
        PMAX = 0.6329
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.030E+3
        PMAX = 1.040E+3
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 270.1E-7
        PMAX = 288.4E-7
       ELSE 
        PMIN = 40.7E-3
        PMAX = 43.9E-3
       ENDIF
*
      ELSEIF (TTEST  <= 600) THEN
       TMAX  = 600.0
       TMIN  = 550.0 
       IF (PROPERTY == 1) THEN
        PMIN = 0.6329
        PMAX = 0.5804
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.040E+3
        PMAX = 1.051E+3
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 288.4E-7
        PMAX = 305.8E-7
       ELSE 
        PMIN = 43.9E-3
        PMAX = 46.9E-3
       ENDIF
*
      ELSEIF (TTEST  <= 650) THEN
       TMAX  = 650.0
       TMIN  = 600.0
       IF (PROPERTY == 1) THEN
        PMIN = 0.5804
        PMAX = 0.5356
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.051E+3
        PMAX = 1.063E+3
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 305.8E-7
        PMAX = 322.5E-7
       ELSE 
        PMIN = 46.9E-3
        PMAX = 49.7E-3
       ENDIF
*
      ELSEIF (TTEST  <= 700) THEN
       TMAX  = 700.0
       TMIN  = 650.0 
       IF (PROPERTY == 1) THEN
        PMIN = 0.5356
        PMAX = 0.4975
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.063E+3
        PMAX = 1.075E+3
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 322.5E-7
        PMAX = 338.8E-7
       ELSE 
        PMIN = 49.7E-3
        PMAX = 52.4E-3
       ENDIF
*
      ELSE 
       TMAX  = 750.0
       TMIN  = 700.0
       IF (PROPERTY == 1) THEN
        PMIN = 0.4975
        PMAX = 0.4643
       ELSEIF (PROPERTY == 2) THEN
        PMIN = 1.075E+3
        PMAX = 1.087E+3
       ELSEIF (PROPERTY == 3) THEN
        PMIN = 338.8E-7
        PMAX = 354.6E-7
       ELSE 
        PMIN = 52.4E-3
        PMAX = 54.9E-3
       ENDIF
*
      ENDIF
*
      RETURN
      END