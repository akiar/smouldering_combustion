*     file acosh.f
***********************************************************************
*
      SUBROUTINE ACOSH(XFINAL, Y,X0,TOL)
*
*     Routine to find the inverse hyperbolic cosine of Y. XFINAL is
*     the final estimate to the value, X0 is the inital guess, and TOL
*     is the convergence criteria for the Newton's method iterations.
*
***********************************************************************
*
      IMPLICIT NONE
      REAL*8 XFINAL,Y,X0,TOL,XN,XNP1,DELTAX
*
      XN = X0
   10 CONTINUE
      XNP1 = XN+(Y-COSH(XN))/SINH(XN)
      DELTAX = ABS(XNP1-XN)
      IF(DELTAX.LT.TOL) THEN
       XFINAL = XNP1
      ELSE
       XN = XNP1
       GOTO 10
      ENDIF
*
      RETURN
      END