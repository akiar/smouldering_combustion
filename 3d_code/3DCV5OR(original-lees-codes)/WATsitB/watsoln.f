C
C file watsoln.f
C=================================================================
C
      SUBROUTINE WATSOLN(x, AP,AW,AE,AS,AN,AB,AT,BP,
     C                   IB,IE,JB,JE,KB,KE,ID,JD,KD,N)
C
C This routine constructs the global stiffness matrix and uses
C the WATSIT-B suite of routines to obtain the solution to:
C                           [A]{x}={B}
C where [A] is a sparse matrix with a block structure.
C
C No row rearrangement is used because the subroutine does funny
C things!  Thus, the matrix must be supplied, via GLOBMS, in the
C correct row order.
C=================================================================
C
C Active coefficients and indices for N simultaneous equations
C at each node
C
      IMPLICIT NONE
      REAL*8 AW(N,N,ID,JD,KD),AE(N,N,ID,JD,KD),AS(N,N,ID,JD,KD)
      REAL*8 AN(N,N,ID,JD,KD),AP(N,N,ID,JD,KD),BP(N,ID,JD,KD)
      REAL*8 AB(N,N,ID,JD,KD),AT(N,N,ID,JD,KD)
      INTEGER IB,IE,JB,JE,KB,KE,ID,JD,KD,N,i,prnt
C
C Set up basic parameters for array sizes
C
      INTEGER neqamx
      PARAMETER (neqamx=4)
      INTEGER nbprow
      PARAMETER (nbprow=7)
      INTEGER nelemx
      PARAMETER (nelemx=1000000)
      INTEGER njamax
      PARAMETER (njamax=nelemx*nbprow)
      INTEGER ndimax
      PARAMETER (ndimax=nelemx*nbprow*neqamx*neqamx)
C
      INTEGER liwork
      INTEGER lrwork
      PARAMETER (liwork=njamax*7,lrwork=njamax*12)
C
C------------------------------
C Set basic storage allocations
C
      INTEGER neqa(nelemx)
      INTEGER ia(nelemx+1)
      INTEGER ja(njamax)
      INTEGER ka(njamax+1)
      INTEGER kvlst(nelemx)
C
      INTEGER lorder(nelemx)
C
      DOUBLE PRECISION a(ndimax)
      DOUBLE PRECISION b(nelemx*neqamx)
      DOUBLE PRECISION x(nelemx*neqamx)
C
      INTEGER mwork(50)
      INTEGER iwork(liwork)
      DOUBLE PRECISION rwork(lrwork)
      INTEGER iwktmp(nelemx)
C
C These parameters are standard
C
      INTEGER neq,nelem,nja,ndima,ierror
      INTEGER redblk,level,blkupd
C
      INTEGER iaccel,north,itmax,ischek
      DOUBLE PRECISION rstall,rprec,rtoler(nelemx*neqamx)
C
      INTEGER rliwkt,rliwks,rliwke,rlrwkt,rlrwks,rlrwke
      INTEGER rnblak,rnbleq,rnjaf,rnaf,rpiaf,rpjaf,rpkaf
      INTEGER rpdiag,rpaf,rpkbls
C
C Misc. local variables
C
      CHARACTER*40 CODDIR
      CHARACTER*80 NAME
      DOUBLE PRECISION toler
C
C------------------------------
C
      INCLUDE 'watsit12.com'
C
C------------------------------
C External declarations
C
      EXTERNAL w1init
C
C------------------------------
C Read the solver parameters each time solver is used.  This
C way, parameters can be modified between iterations without
C halting the computational run.
C

      CODDIR='/home/hpstor2/mech/mech6/cfische/3DCV5OR'
      NAME=CODDIR//'/WATsitB/watsol.in'
      OPEN(unit=88,file=NAME)
      READ(88,8800) redblk,level,blkupd,iaccel,north,itmax,
     &              ischek,rstall,rprec,toler
      CLOSE(unit=88)
 8800 FORMAT(I5/I5/I5/I5/I5/I5/I5/1PD12.5/1PD12.5/1PD12.5)
C
C------------------------------
C Initialize the solver
C
      CALL w1init
C
C Announce startup (optional)
C
cc    WRITE (wgunit,9000)
cc    WRITE (wgunit,9010)
cc    WRITE (wgunit,9030)
C
C------------------------------
C Get the matrix
C
      CALL GLOBMS(a,b,ia,ja,ka,kvlst,neqa,
     C            AP,AW,AE,AS,AN,AB,AT,BP,
     C            IB,IE,JB,JE,KB,KE,ID,JD,KD,N,
     C            ndimax,nelemx*neqamx,nelemx,njamax)
C
C Calculate array lengths for processing
C
      nelem = (IE-IB+1)*(JE-JB+1)*(KE-KB+1)
      neq   = kvlst(nelem+1)-1
      nja   = ia(nelem+1)-1
      ndima = ka(nja+1)-1
C
C Set tolerance vector, and block update parameter
C
      DO 10 i = 1,neq
        rtoler(i) = toler
   10 CONTINUE
      blkupd = 0
C
C------------------------------
C NOTE: Do not use w2rowf; it does funny things!
C The w2rowf routine will fix the row order if the rows do not
C have the diagonal at the start, as they should.
C
cc    WRITE (wgunit,9020)
C
cc    CALL w2rowf(nelem,nja,ndima,neqa,ia,ja,ka,a,400,rwork(1),ierror)
C
cc    WRITE (wgunit,9090) 'w2rowf',ierror
C
C Check the structure
C
      CALL w2chek(nelem,nja,ndima,neqa,ia,ja,ka,lorder,iwktmp,ierror)
      IF(ierror.ne.0) THEN
        WRITE(*,*) 'Error in w2chek: ',ierror
	STOP
      ENDIF
C
C------------------------------
C Do the preconditioner.
C
      CALL w2anly(nelem,nja,ndima,neqa,ia,ja,ka,redblk,level,blkupd,
     1            lorder,liwork,lrwork,mwork,iwork,rwork,ierror)
      IF(ierror.ne.0) THEN
        WRITE(*,*) 'Error in w2anly: ',ierror
	STOP
      ENDIF
C
      CALL w2fctr(neqa,ia,ja,ka,a,kvlst,b,mwork,iwork,rwork,ierror)
      IF(ierror.ne.0) THEN
        WRITE(*,*) 'Error in w2fctr: ',ierror
	STOP
      ENDIF
C
C------------------------------
C Solve the problem
C
      CALL w2solv(neqa,ia,ja,ka,a,kvlst,b,iaccel,north,itmax,ischek,
     1            rstall,rprec,rtoler,x,mwork,iwork,rwork,ierror)
      IF(ierror.ne.0) THEN
        WRITE(*,*) 'Error in w2solv: ',ierror
	STOP
      ENDIF
C
C Print general information about the solve.  Setting prnt=1 will
C print a number of parameters to the screen which may be useful
C when trouble shooting a particular run; in general, you don't
C need to see it.  The run will stop if an error occurs.
C
      prnt=0
      IF(prnt.eq.1) THEN
      wgunit=6
      WRITE (wgunit,9030)
      WRITE (wgunit,9030)
      WRITE (wgunit,*) '****  WATSI2 COMMON block information ****'
      WRITE (wgunit,9030)
C
      WRITE (wgunit,9040) 'Smallest scaled matrix diagonal',wgdiag,
     1  wgrow
C
      WRITE (wgunit,9030)
C
      WRITE (wgunit,9080) 'Number of iterations in linear solve',wgnitr
      WRITE (wgunit,9080) 'Number of stall/restarts required',wgrest
C
      WRITE (wgunit,9030)
C
      IF (mod(wgcnvc,10).EQ.1) THEN
          WRITE (wgunit,9030)
          WRITE (wgunit,9060) 'Residual reduction convergence'//
     1      ' criterion satisfied'
      END IF
C
      IF (wgcnvc.GE.10) THEN
          WRITE (wgunit,9060) 'Update convergence criterion satisfied'
      END IF
C
      WRITE (wgunit,9030)
C
      WRITE (wgunit,9050) 'Initial l2 norm of linear residual',wgrmsi
      WRITE (wgunit,9050) 'Final l2 norm of linear residual',wgrmsf
      WRITE (wgunit,9050) 'Coefficient of last solution update',wgohm
C
      WRITE (wgunit,9030)
C
      WRITE (wgunit,9080) 'Number of equations',wgbleq
cc    WRITE (wgunit,9080) 'Number of black points',wgblak
C
cc    WRITE (wgunit,9030)
C
      WRITE (wgunit,9080) 'Amount of IWORK actually used',wgiwku
      WRITE (wgunit,9080) 'Amount of RWORK actually used',wgrwku
C
cc    WRITE (wgunit,9080) 'Number of fill elements in precond.',wgfill
cc    WRITE (wgunit,9080) 'Scalar fill in preconditioner',wgasiz
C
      WRITE (wgunit,9030)
C
C------------------------------
C
C Internal storage information
C
      CALL w1data(rliwkt,rliwks,rliwke,rlrwkt,rlrwks,rlrwke,rnblak,
     1            rnbleq,rnjaf,rnaf,rpiaf,rpjaf,rpkaf,rpdiag,rpaf,
     2            rpkbls,mwork)
C
      WRITE (wgunit,9070)
      WRITE (wgunit,9030)
C
      WRITE (wgunit,9080) 'Total IWORK space',rliwkt
      WRITE (wgunit,9080) 'Start of available IWORK space',rliwks
      WRITE (wgunit,9080) 'End of available IWORK space',rliwke
      WRITE (wgunit,9030)
C
      WRITE (wgunit,9080) 'Total RWORK space',rlrwkt
      WRITE (wgunit,9080) 'Start of available RWORK space',rlrwks
      WRITE (wgunit,9080) 'End of available RWORK space',rlrwke
      WRITE (wgunit,9030)
C
      WRITE (wgunit,9080) 'Number of black elements',rnblak
      WRITE (wgunit,9080) 'Number of black scalar equations',rnbleq
      WRITE (wgunit,9080) 'Size of preconditioner jaf array',rnjaf
      WRITE (wgunit,9080) 'Size of preconditioner af array',rnaf
      WRITE (wgunit,9030)
C
      WRITE (wgunit,9080) 'Start of preconditioner iaf array',rpiaf
      WRITE (wgunit,9080) 'Start of preconditioner jaf array',rpjaf
      WRITE (wgunit,9080) 'Start of preconditioner kaf array',rpkaf
      WRITE (wgunit,9080) 'Start of preconditioner diagf array',rpdiag
      WRITE (wgunit,9080) 'Start of preconditioner af array',rpaf
      WRITE (wgunit,9080) 'Start of preconditioner kvlst() array',
     1  rpkbls
      WRITE (wgunit,9030)
C
C------------------------------
C
      WRITE (wgunit,9110) ierror
      WRITE (wgunit,9030)
      WRITE (wgunit,9030)
      ENDIF
C
C------------------------------
C
 9000 FORMAT ('**** WATSIT-B v 2.0 Block Matrix Solver ',
     1       ' ****')
 9010 FORMAT ('**** by Simon S. Clift  (c) University of Waterloo 1993')
 9020 FORMAT ('****  Starting solution sequence ****')
 9030 FORMAT (1X)
 9040 FORMAT (4X,A39,T46,E14.8E2,' at row ',I8)
 9050 FORMAT (4X,A39,T46,E14.8E2)
 9060 FORMAT (4X,A54)
 9070 FORMAT ('****  Work array storage information ****')
 9080 FORMAT (4X,A39,T46,I14)
 9090 FORMAT (/17X,'Error after routine ',A6,T46,I14)
 9100 FORMAT ('**** Writing matrix solution to ',A40)
 9110 FORMAT ('**** Matrix test ends with error code ',I3)
 9120 FORMAT (1X,'Solution'/ (1P,8D10.2))
*
      RETURN
      END
