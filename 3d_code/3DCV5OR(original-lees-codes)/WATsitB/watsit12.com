C=================================================================
C
C WATSIT-1  v 2.0   by Simon S. Clift
C
C Derived from Mat1 v. 12 by P. A. Forsyth, J. R. Knightly and
C                            E. F. D'Azevedo
C
C "watsit12.com"
C
C This INCLUDE file lists the COMMON blocks, which must
C contain only informational variables (i.e. no parameters or
C control variables).
C
C=================================================================
C
C WATSI1 -- Definition of the common block controlling output.
C
C     wgunit     Unit number for output
C     wginfo     Level of infomation messages printed
C                <= -2 No printing
C                 = -1 Error messages only
C                 =  0 Warning and error messages (default)
C                 =  1 Convergence history
C                 =  2 Convergence history and fill
C     wgdebg     Debug information level
C                 =  0 No printing (default)
C                 =  1 Scalars printed
C                 =  2 Scalars and vectors printed
C
      INTEGER wgunit,wginfo,wgdebg
      COMMON /WATSI1/ wgunit,wginfo,wgdebg
C
C=================================================================
C
C WATSI2 -- Definition of the informational common block
C     
C These values are _not_ to be used as working values, but are
C only to be stored to by the solver routines.  These are only
C accessed by the user.
C
C     wgdiag    Smallest diagonal in the preconditioner
C               (Normalized by max. in row). WATSIT-1 only.
C     wgrow     Row where the small diagonal is (orig. ordering)
C               WATSIT-1 only.
C
C     wgrmsi    Initial l_2 norm of linear residual
C     wgrmsf    Final l_2 norm of linear residual
C     wgohm     Coefficient in last update of solution
C
C     wgiwku    Amount of IWork actually used.
C     wgrwku    Amount of RWork actually used.
C
C     wgfill    Number of blocks in the preconditioning matrix
C     wgasiz    Size of A() matrix required (=wgfill for WATSIT-1)
C               for preconditioner.
C
C     wgblak    Number of black points allocated
C     wgbleq    Total number of black point equations (=wgblak in
C               WATSIT-1)
C
C     wgnitr    Number of iterations required for the solve
C     wgrest    Number of restarts required
C     wgcnvc    Convergence criteria satisfied at convg.
C               May be summed if more than one is 
C               satisfied.
C                =   1  Residual reduction
C                =  10  Change < Abs. tolerance
C
      INTEGER wgrow,wgiwku,wgrwku,wgfill,wgasiz
      INTEGER wgblak,wgbleq,wgnitr,wgrest,wgcnvc
      DOUBLE PRECISION wgdiag,wgrmsi,wgrmsf,wgohm
C
      COMMON /WATSI2/ wgdiag,wgrmsi,wgrmsf,wgohm
      COMMON /WATSI2/ wgrow,wgiwku,wgrwku,wgfill,wgasiz
      COMMON /WATSI2/ wgblak,wgbleq,wgnitr,wgrest,wgcnvc
