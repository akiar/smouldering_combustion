C=================================================================
C
C WATSIT-1  v 2.0   by Simon S. Clift
C
C Derived from Mat1 v. 12 by P. A. Forsyth, J. R. Knightly and
C                            E. F. D'Azevedo
C
C "watsitmw.com"
C
C This file lists parameters indexing into the MWork array, the
C master work array.
C
C All parameters referenceing the mwork array must begin with "m".
C
C=================================================================
C
C The MWORK workspace is used to store INTEGER pointers and
C control values.
C
C NOTE:  A certain amount of space has been left unused in each
C        section of the MWORK array to allow for expansion.
C
C     mwlen       Length of MWork
C
      INTEGER mwlen
      PARAMETER(mwlen=50)
C
C=================================================================
C
C Parameters that can be read by the W1DATA routine.  These return
C internal storage allocation information
C 
C     MWork(mliwkt)    Amount of iwk total, i.e. pointer to end
C                       of iwk space.
C     MWork(mliwks)    Start of available iwk space
C     MWork(mliwke)    End of available iwk space
C
C     MWork(mlrwkt)    Amount of rwk total, i.e. pointer to end
C                       of rwk space.
C     MWork(mlrwks)    Start of available rwk space
C     MWork(mlrwke)    End of available rwk space
C
      INTEGER mliwkt
      PARAMETER ( mliwkt = 1 )
      INTEGER mliwks
      PARAMETER ( mliwks = 2 )
      INTEGER mliwke
      PARAMETER ( mliwke = 3 )
      INTEGER mlrwkt
      PARAMETER ( mlrwkt = 4 )
      INTEGER mlrwks
      PARAMETER ( mlrwks = 5 )
      INTEGER mlrwke
      PARAMETER ( mlrwke = 6 )
C
C=================================================================
C
C General parameters stored in MWORK
C
C     MWork(mneq)      Number of equations in full system 
C     MWork(mnelem)    Number of element in block system
C     MWork(mnja)      Length of ja() array
C     MWork(mndima)    Size of a() array
C
C     MWork(mnblak)       Number of black points
C     MWork(mnred)      + Number of red points
C
C     MWork(mnbunk)    Max. number of black unknowns
C
C     MWork(mstage)    Stage of processing of this matrix
C                       = 0 Nothing done.
C                       = 1 Basic initialization and check done
C                       = 2 Data structure analysis done
C                       = 3 Numeric factorization complete
C                       = 4 Solve complete and converged
C     MWork(mredbl)    Do red/black elimination on this thing
C                       = 0 No
C                       = 1 Yes
C
C     MWork(mnfill)    Number of fill elements, also length of 
C                      jaf() array.
C     MWork(mnafsz)    Amount of scalar fill, length of af() array.
C
C     MWork(mneqam)    Max. value of neqa()
C
      INTEGER mneq
      PARAMETER ( mneq = 10 )
      INTEGER mnelem
      PARAMETER ( mnelem = 11 )
      INTEGER mnja
      PARAMETER ( mnja = 12 )
      INTEGER mndima
      PARAMETER ( mndima = 13 )
      INTEGER mnblak
      PARAMETER ( mnblak = 14 )
      INTEGER mnred
      PARAMETER ( mnred = 15 )
      INTEGER mnbunk
      PARAMETER ( mnbunk = 16 )
      INTEGER mstage
      PARAMETER ( mstage = 17 )
      INTEGER mredbl
      PARAMETER ( mredbl = 18 )
      INTEGER mnfill
      PARAMETER ( mnfill = 19 )
      INTEGER mnafsz
      PARAMETER ( mnafsz = 20 )
      INTEGER mneqam
      PARAMETER ( mneqam = 21 )
C
C=================================================================
C
C Pointers set by allocation of the given IWork and RWork 
C workspaces.
C
C     IWork Allocation Pointers
C     -------------------------
C     MWork(mpatrb)    Ptr. from original matrix to red/black
C                       ordering: iatorb(neq)
C     MWork(mpbtoa)    Ptr. to array indexing black points 
C                       to original matrix: iblata(nblack)
C     MWork(mprtoa)    Ptr. to array indexing red points
C                       to original matrix: iredta(nred)
C
C     MWork(mpiaf)     Ptr. to iaf() ('ia' for preconditioner)
C     MWork(mpjaf)     Ptr. to jaf() ('ja' for preconditioner)
C     MWork(mpkaf)     Ptr. to block kaf() array for preconditioner
C
C     MWork(mpdiag)    Ptr. to diagf() ('diag' for preconditioner)
C
C     MWork(mpkbls)    Ptr. to black unknowns kvlst() type array
C
C     RWork Allocation Pointers
C     -------------------------
C     MWork(mpaf)      Ptr. to af() ('a' for preconditioner)
C
      INTEGER mpatrb
      PARAMETER ( mpatrb = 30 )
      INTEGER mpbtoa
      PARAMETER ( mpbtoa = 31 )
      INTEGER mprtoa
      PARAMETER ( mprtoa = 32 )
      INTEGER mpiaf
      PARAMETER ( mpiaf = 33 )
      INTEGER mpjaf
      PARAMETER ( mpjaf = 34 )
      INTEGER mpkaf
      PARAMETER ( mpkaf = 35 )
      INTEGER mpdiag
      PARAMETER ( mpdiag = 36 )
      INTEGER mpkbls
      PARAMETER ( mpkbls = 37 )
C
      INTEGER mpaf
      PARAMETER ( mpaf = 40 )
C
C------------------------------
C
C Sundry useful parameters
C
C Acceleration methods
C
      INTEGER zcg,zortho,zcgs,zgmres,zcgstb
      PARAMETER( zcg=0,zortho=1,zcgs=2,zcgstb=3,zgmres=4)
C
C------------------------------
C
C Parameter for ratio of REAL to DOUBLE PRECISION 
C
      REAL zrealr
      PARAMETER (zrealr=2.0)
