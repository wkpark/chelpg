C       CHELP-NET ATOMIC CHARGES FROM AB INITIO WAVE FUNCTIONS
C       Modified for Grid Operations by Curt Breneman, Yale University
C       Department of Chemistry, 3/88  (Currently of Rensselaer
C       Polytechnic Institute, Troy, NY 12180.)
C
C        CHELPG
C
C       (NET ATOMIC) CHARGES FIT TO ELECTROSTATIC POTENTIALS
C
C        Original CHELP code by:
C        M.M. FRANCL
C        L.E. CHIRLIAN
C
C        OCTOBER 1985
C        PRINCETON CHEMISTRY DEPARTMENT VAX 11/780
C        VMS 3.7
C
C        FEBRUARY 1988
C        MODIFIED TO USE GAUSSIAN86 CHECKPOINT FILES
C        Modified to use G88/90 checkpoint files 1/89
C        YALE UNIVERSITY DEPARTMENT OF CHEMISTRY
C        WIBERG GROUP VMS 4.5
C        CURT BRENEMAN
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*4 HANDLE1
      HANDLE1=0
C***    TRACE-7
C      ISTAT1=LIB$INIT_TIMER(HANDLE1)
C***
C
C        READ IN DATA FROM CHECKPOINT FILE
C
      CALL READIN
C
C
C        SELECT POINTS FOR FITTING, BEGIN WITH SHELL OF RADIUS 2A AND
C        INCREASING BY .5A SELECT POINTS FROM THE ROUGHLY RADIAL DISTRIBUTION
C        WHICH ARE NOT ENCLOSED BY THE VAN DER WAALS ENVELOPE OF THE MOLECULE
C        UNTIL A PREDETERMINED MAXIMUM NUMBER OF POINTS HAS BEEN REACHED
C
      CALL BALL
C
C        CALCULATE THE ELECTROSTATIC POTENTIAL USING FIRST ORDER HARTREE-FOCK
C        PERTURBATION THEORY
C
      CALL EP
C
C        USING METHOD OF LAGRANGE MULTIPLIERS, FIT BY LEAST SQUARES THE CHARGES
C        TO THE ELECTROSTATIC POTENTIAL, CONSTRAINING THE FIT TO REPRODUCE THE
C        TOTAL MOLECULAR CHARGE
C
      CALL FIT
C
C        PRINT OUT TABLE OF RESULTS
C
      CALL OUTPUT
C
C***    TRACE-7
C      ISTAT1=LIB$SHOW_TIMER(HANDLE1)
C***
      END
C
C
      SUBROUTINE BALL
C
C        ROUTINE TO SELECT POINTS FOR FITTING TO THE ELECTROSTATIC POTENTIAL.
C
C        POINTS WHICH LIE WITHIN THE VAN DER WAALS ENVELOPE OF THE MOLECULE
C        ARE EXCLUDED.
C
C        POINTS ARE INITIALLY SELECTED IN A CUBE AROUND THE MOLECULE WHICH
C        IS SCALED TO THE SIZE OF THE MOLECULE+RMAX. THIS IS PRESENTLY AN INPUT
C        PARAMETER.  POINTS ARE THEN EXCLUDED IF THEY FALL WITHIN THE INPUT
C        VDW RADIUS OF ANY OF THE ATOMS, OR, IF THEY FALL OUTSIDE
C        A DESIGNATED DISTANCE (RMAX) FROM ALL OF THE ATOMS.   THE REMAINING
C        POINTS ARE PACKED IN A SET OF THREE (X,Y,Z) VECTORS, AND SENT TO THE
C        LAGRANGE LEAST-SQUARES FITTING ROUTINE.  THE ORIGINAL CHELP INPUT
C        DECK IS AUGMENTED BY ADDING TWO FREE-FORMAT VARIABLES AT THE END.
C        THE TWO NEW INPUT VARIABLES ARE 'RMAX' AND 'DELR', WHERE RMAX
C        IS THE MAXIMUM DISTANCE A POINT CAN BE FROM ANY ATOM AND STILL
C        BE CONSIDERED IN THE FIT, AND DELR IS THE DISTANCE BETWEEN POINTS
C        IN THE GRID.  BOTH RMAX AND DELR ARE IN ANGSTROMS.
C
C        CURT BRENEMAN AND TERESA LEPAGE
C        YALE UNIVERSITY DEPARTMENT OF CHEMISTRY 3/88
C
C        ORIGINAL CODE BY:
C
C        L.E. CHIRLIAN
C        M.M. FRANCL
C        APRIL 1985
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (NPOINTS = 50000)
      COMMON /IO/ IN,IOUT
C+++
      COMMON /MOL/    NATOMS,ICHARG,MULTIP,NAE,NBE,NEL,NBASIS,
     $                IAN(401),ATMCHG(400),C(3,400)
C+++
      COMMON /IPO/ IPO(5)
      COMMON /SPHERE/ RADII(400),NTOTP
      COMMON /POINTS/ P(3,NPOINTS), MAXPNTS
C
      DATA ANG2AU /1.889726878D0/
C
C***    READ IN THE THE RMAX AND DELR VALUES IN ANGSTROMS.
C
        read(IN,*) RMAX, DELR
        write(IOUT,*) ' RMAX = ',RMAX,' (ANGS), DELR = ',DELR,' (ANGS).'
C***
C
C        CONVERT RADII TO AU
C
      DELR = DELR * ANG2AU
      RMAX = RMAX * ANG2AU
C
C        WHILE CONVERTING THE VDW RADII TO AU, FIND THE EXTREMA OF THE
C        MOLECULAR GEOMETRY.
C
      XMAX=-50.0D0
      XMIN=50.0D0
      YMAX=-50.0D0
      YMIN=50.0D0
      ZMAX=-50.0D0
      ZMIN=50.0D0
C
      WRITE(IOUT,*) ' THERE ARE ',NATOMS,' ATOMS TO CONSIDER.'
      DO 10 I=1,NATOMS
      RADII(I) = RADII(I) * ANG2AU
C
      IF (C(1,I) .GT. XMAX) XMAX = C(1,I)
      IF (C(1,I) .LT. XMIN) XMIN = C(1,I)
      IF (C(2,I) .GT. YMAX) YMAX = C(2,I)
      IF (C(2,I) .LT. YMIN) YMIN = C(2,I)
      IF (C(3,I) .GT. ZMAX) ZMAX = C(3,I)
      IF (C(3,I) .LT. ZMIN) ZMIN = C(3,I)
   10 CONTINUE
C
      WRITE(IOUT,*) ' XMAX = ',XMAX,' (AU), XMIN = ',XMIN,' (AU).'
      WRITE(IOUT,*) ' YMAX = ',YMAX,' (AU), YMIN = ',YMIN,' (AU).'
      WRITE(IOUT,*) ' ZMAX = ',ZMAX,' (AU), ZMIN = ',ZMIN,' (AU).'
C
C        DETERMINE THE MINIMUM CUBE DIMENSIONS REQUIRED TO CONTAIN
C        THE MOLECULE, INCLUDING THE MAXIMUM SELECTION RADIUS (RMAX)
C        ON BOTH SIDES.
C
      XRANGE = XMAX - XMIN + 2.0D0 * RMAX
      YRANGE = YMAX - YMIN + 2.0D0 * RMAX
      ZRANGE = ZMAX - ZMIN + 2.0D0 * RMAX
C
      NXPTS = INT(XRANGE/DELR)
      NYPTS = INT(YRANGE/DELR)
      NZPTS = INT(ZRANGE/DELR)
C
      WRITE(IOUT,*) ' NUMBER OF X POINTS REQUIRED = ',NXPTS
      WRITE(IOUT,*) ' NUMBER OF Y POINTS REQUIRED = ',NYPTS
      WRITE(IOUT,*) ' NUMBER OF Z POINTS REQUIRED = ',NZPTS
      MAXPOSS = NXPTS * NYPTS * NZPTS
      WRITE(IOUT,*) ' TOTAL NUMBER OF POINTS CONSIDERED = ',MAXPOSS
C
C
C        RESET POINT COUNTER FOR NUMBER OF SELECTED POINTS
C
      IPOINT = 0
C
C       LOOP OVER POSSIBLE POINTS
C
      DO 200 II = 1,NXPTS + 1
C
      P1 = XMIN - RMAX + DBLE(II-1) * DELR
C
      DO 200 JJ = 1,NYPTS + 1
C
      P2 = YMIN - RMAX + DBLE(JJ-1) * DELR
C
      DO 200 KK = 1,NZPTS + 1
C
      P3 = ZMIN - RMAX + DBLE(KK-1) * DELR
C
C
C        IS THIS POINT WITHIN A VAN DER WAALS SPHERE OR OUTSIDE THE
C        RMAX DISTANCE FROM ALL ATOMS?
C
      RADMIN=50.0D0
      DO 100 I=1,NATOMS
      VRAD = RADII(I)
      DIST = (P1 - C(1,I))**2 + (P2 - C(2,I))**2 + (P3 - C(3,I))**2
      DIST = DSQRT(DIST)
      IF (DIST .LT. VRAD) GOTO 210
      IF (DIST .LT. RADMIN) RADMIN = DIST
  100 CONTINUE
      IF (RADMIN .GT. RMAX) GOTO 210
C
C        STORE POINTS (IN ATOMIC UNITS)
C
      IPOINT = IPOINT + 1
      P(1,IPOINT) = P1
      P(2,IPOINT) = P2
      P(3,IPOINT) = P3
      IF (IPO(2) .EQ. 1)
     $ WRITE(IOUT,*) 'POINT ',IPOINT,' X,Y,Z ',P1,P2,P3
  210 CONTINUE
  200 CONTINUE
C
      MAXPNTS = IPOINT
      WRITE(IOUT,*) ' NUMBER OF POINTS SELECTED FOR FITTING : ',MAXPNTS
      RETURN
      END
C
      SUBROUTINE EP
C
C        ROUTINE TO CALCULATE THE ELECTROSTATIC POTENTIAL FROM FIRST ORDER
C        PERTURBATION THEORY
C
C        M.M. FRANCL    APRIL 1985
C        MODIFIED VERSION OF A MEPHISTO ROUTINE
C        RESTRICTED TO CLOSED SHELL MOLECULES
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NPOINTS = 50000)
      INTEGER*4 SHELLA,SHELLN,SHELLT,AOS,SHELLC,AON,HANDLE
C
      COMMON /IO/ IN,IOUT
      COMMON /IPO/ IPO(5)
C+++
      COMMON /MOL/    NATOMS,ICHARG,MULTIP,NAE,NBE,NEL,NBASIS,
     $                IAN(401),ATMCHG(400),C(3,400)
C
C===    Gaussian88 Modification for enlarged common /b/.
      Common/B/EXX(6000),C1(6000),C2(6000),C3(6000),X(2000),Y(2000),
     $Z(2000),JAN(2000),ShellA(2000),ShellN(2000),ShellT(2000),
     $ShellC(2000),AOS(2000),AON(2000),NShell,MaxTyp
C==== Old G86 Version of common /b/
c      COMMON/B/EXX(1200),C1(1200),C2(1200),C3(1200),
c     $         X(400),Y(400),Z(400),JAN(400),SHELLA(400),SHELLN(400),
c     $         SHELLT(400),SHELLC(400),AOS(400),AON(400),NSHELL,MAXTYP
C+++
C      COMMON /B/ EXX(240),C1(240),C2(240),C3(240),X(80),Y(80),Z(80),
C     $           JAN(80),SHELLA(80),SHELLN(80),SHELLT(80),SHELLC(80)
C     $          ,AOS(80),AON(80),NSHELL,MAXTYP
C      COMMON /MOL/ NATOMS,ICHARG,MULTIP,NAE,NBE,NEL,NBASIS,IAN(101),
C     $             ATMCHG(100),C(3,100)
      COMMON /POINTS/ P(3,NPOINTS),MAXPNTS
      COMMON /ELP/ ELECP(NPOINTS)
      COMMON /CHARGE/ COEF_ALPHA(100000),COEF_BETA(100000),IUHF
      COMMON /OUT/ NTITLE(20,3),CHKFIL,Q(400),I6TO5,RMS,PERCENT,
     1             NLIN,NEND(3)
C
      DIMENSION HPERT(100000),INDEX(1280)
C
      DATA IPTCHG/1.0/
      DATA ZERO/0.0/, TWO/2.0/, VNUCMAX/30.0/
C        DIVERT TO ROUTINE UEP IF WAVEFUNCTION IS UNRESTRICTED
C        HARTREE-FOCK WAVEFUNCTION
C
      IF (IUHF .EQ. 1) THEN
      CALL UEP
      RETURN
      END IF
C
      HANDLE = 0
C
C        SET UP THE INDEXING TABLE FOR HPERT
C
      DO 100 I=1,NBASIS
      INDEX(I) = (I-1)*I/2
  100 CONTINUE
C
C        BEGIN LOOP TO CALCULATE ELECTROSTATIC POTENTIAL
C
      NOCC = NEL / 2
      MVIR = NOCC + 1
C
C        START OF LOOP
C
      DO 200 NPNT=1,MAXPNTS
      X1 = P(1,NPNT)
      X2 = P(2,NPNT)
      X3 = P(3,NPNT)
C
C     CALCULATE THE ONE-ELECTRON INTEGRALS
C
      IF (IPO(5).EQ.1) THEN
      WRITE(IOUT,3010)
 3010 FORMAT(1X,'TIME FOR INTEGRALS')
C***
C      ISTAT = LIB$INIT_TIMER(HANDLE)
C***
      END IF
C
      CALL INTGRL (HPERT,X1,X2,X3,IPTCHG,I6TO5)
C
C***
C      IF (IPO(5).EQ.1) ISTAT = LIB$SHOW_TIMER(HANDLE)
C***
C
      IF (IPO(4).EQ.1) CALL LINOUT (HPERT,NBASIS,0,0)
C
      IF (IPO(5).EQ.1) THEN
      WRITE(IOUT,3000)
 3000 FORMAT(1X,'TIME FOR TRANSFORM')
C***
C      ISTAT = LIB$INIT_TIMER(HANDLE)
C***
      END IF
C
C     FORM THE HPERT MATRIX ELEMENTS
C
      E = ZERO
      ICOEFI = -NBASIS
C
C        SUM OVER OCCUPIED MOS
C
      DO 220 II=1,NOCC
      ICOEFI = ICOEFI + NBASIS
C
C        CALCULATE ELECTROSTATIC POTENTIAL
C
      DO 221 IP=1,NBASIS
      CPI = COEF_ALPHA(ICOEFI+IP)
      IPDEX = INDEX(IP)
C
      DO 222 IQ=1,IP
      E = E + CPI * COEF_ALPHA(ICOEFI+IQ) * HPERT(IPDEX+IQ)
  222 CONTINUE
      DO 223 IQ=IP+1,NBASIS
      E = E + CPI * COEF_ALPHA(ICOEFI+IQ) * HPERT(IP+INDEX(IQ))
  223 CONTINUE
C
  221 CONTINUE
  220 CONTINUE
C
C***
C      IF (IPO(5).EQ.1) ISTAT = LIB$SHOW_TIMER(HANDLE)
C***
C
C        CALCULATE NUCLEAR PART OF ELECTROSTATIC POTENTIAL
C
      VNUC = ZERO
      DO 300 IATOM=1,NATOMS
      DEL1 = C(1,IATOM) - X1
      DEL2 = C(2,IATOM) - X2
      DEL3 = C(3,IATOM) - X3
      RA = DSQRT(DEL1*DEL1 + DEL2*DEL2 + DEL3*DEL3)
      IF (RA.EQ.ZERO) THEN
      VNUC=VNUCMAX
      GOTO 310
      END IF
      VNUC = VNUC + IAN(IATOM) / RA
  300 CONTINUE
  310 CONTINUE
C
      ELECP(NPNT) = (E * TWO + VNUC * IPTCHG)
      IF (IPO(5) .EQ. 1) WRITE(IOUT,*) 'E(',NPNT,') = ',E
  200 CONTINUE
      RETURN
      END
      SUBROUTINE FIT
C
C        ROUTINE TO USE METHOD OF LAGRANGE MULTIPLIERS TO OBTAIN BEST
C        LEAST SQUARE FIT WITH CONSTRAINTS
C
C        M.M. FRANCL
C        APRIL 1985
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NPOINTS = 50000)
      INTEGER*4 WHICH1
      CHARACTER*40 CHKFIL
C
      COMMON /IO/ IN,IOUT
      COMMON /IPO/ IPO(5)
      COMMON /ELP/ E(NPOINTS)
      COMMON /POINTS/ P(3,NPOINTS),MAXPNTS
      COMMON /OUT/ NTITLE(20,3),CHKFIL,X(400),I6TO5,RMS,PERCENT,
     1             NLIN,NEND(3)
C+++
      COMMON /MOL/    NATOMS,ICHARG,MULTIP,NAE,NBE,NEL,NBASIS,
     $                IAN(401),ATMCHG(400),C(3,400)
C+++
C      COMMON /MOL/ NATOMS,ICHARG,MULTIP,NAE,NBE,NEL,NBASIS,IAN(101),
C     $             ATMCHG(100),C(3,100)
C
      DIMENSION A(400,400),Y(400),IS(2,400),IAD1(400),IAD2(400)
      DIMENSION D(400),WHICH1(3)
C
C        DEBYE = CONVERSION FROM DEBYES TO AU
C
      DATA ONE/1.0/, ZERO/0.0/, DEBYE/0.393427328/, MAXDIM/400/
      DATA AU2CAL/627.51/, HALF/0.5/, HUNDRED/100.0/,NCONSTR/1/
C
C        SET UP MATRIX OF LINEAR COEFFICIENTS, A
C
C        BEGIN LOOP OVER ROWS
C
      DO 100 K=1,NATOMS
C
C        BEGIN LOOP OVER COLUMNS
C
      DO 200 MU=1,NATOMS
C
      SUM = ZERO
      DO 400 I=1,MAXPNTS
      RIK = (P(1,I)-C(1,K))**2 + (P(2,I)-C(2,K))**2 + (P(3,I)-C(3,K))**2
      RIK = DSQRT(RIK)
      RIMU = (P(1,I)-C(1,MU))**2 + (P(2,I)-C(2,MU))**2 +
     $       (P(3,I)-C(3,MU))**2
      RIMU = DSQRT(RIMU)
      SUM = SUM + ONE / (RIK * RIMU)
  400 CONTINUE
C
      A(K,MU) = SUM
  200 CONTINUE
C
C        FILL OUT COLUMNS CORRESPONDING TO LAGRANGE MULTIPLIERS
C
      A(K,NATOMS+1) = HALF
C
C
  100 CONTINUE
C
C        FILL OUT THE ROWS CORRESPONDING TO CONSTRAINTS
C
      DO 500 MU=1,NATOMS
      A(NATOMS+1,MU) = ONE
C
  500 CONTINUE
C
C        FILL OUT THE BLOCK WHICH CONNECTS LAGRANGE MULTIPLIERS TO
C        CONSTRAINTS
C
      DO 600 K=NATOMS+1,NATOMS+NCONSTR
      DO 600 MU=NATOMS+1,NATOMS+NCONSTR
      A(K,MU) = ZERO
  600 CONTINUE
C
C****DEBUG*****
C
      IF (IPO(3) .EQ. 1) THEN
      WRITE(IOUT,*) 'A MATRIX'
      DO 699 K=1,NATOMS+NCONSTR
      WRITE(IOUT,1699) (A(K,MU),MU=1,NATOMS+NCONSTR)
 1699 FORMAT(1X,10F10.4)
  699 CONTINUE
      END IF
C***************
C
C        CONSTRUCT COLUMN VECTOR, Y
C
      DO 700 K=1,NATOMS
      SUM = ZERO
      DO 710 I=1,MAXPNTS
      RIK = (P(1,I)-C(1,K))**2 + (P(2,I)-C(2,K))**2 +
     $      (P(3,I)-C(3,K))**2
      RIK = DSQRT(RIK)
      SUM = SUM + E(I) / RIK
  710 CONTINUE
      Y(K) = SUM
      IF (IPO(3) .EQ. 1) WRITE(IOUT,*) K,Y(K)
  700 CONTINUE
C
C        CONSTRUCT THE PORTION OF Y CORRESPONDING TO LAGRANGE MULTIPLIERS
C
C
      Y(NATOMS+1) = DFLOAT(ICHARG)
C
C
      IF (IPO(3) .EQ. 1)
     $  WRITE(IOUT,*) 'COL VECTR Y', (Y(KK),KK=1,NATOMS+NCONSTR)
C
C        SOLVE MATRIX EQUATION AX = Y;
C        WHERE X = (Q1,Q2, ... QN,L1,L2, ... ,LN)
C
C        X = A(INV)Y
C
C        INVERT A
C
      CALL INV(A,NATOMS+NCONSTR,IS,IAD1,IAD2,D,MAXDIM)
C
C****DEBUG*****
C
      IF (IPO(3) .EQ. 1) THEN
      WRITE(IOUT,*) 'A INVERSE'
      DO 799 K=1,NATOMS+NCONSTR
      WRITE(IOUT,1699) (A(K,MU),MU=1,NATOMS+NCONSTR)
  799 CONTINUE
      END IF
C**************
C
C        PERFORM MATRIX MULTIPLICATION A(INV)Y
C
      CALL MULTAY(A,Y,X,NATOMS+NCONSTR,MAXDIM)
C
      IF (IPO(3) .EQ. 1) THEN
      WRITE(IOUT,*) 'CHARGES:  '
      DO 899 I=1,NATOMS
      WRITE(IOUT,*) IAN(I),X(I)
  899 CONTINUE
      END IF
C
C        COMPUTE RMS DEVIATION AND MEAN ABSOLUTE % DEVIATION
C
      RMS = ZERO
      PERCENT = ZERO
      DO 800 I=1,MAXPNTS
      EQ = ZERO
      DO 810 J=1,NATOMS
      DIST = (P(1,I)-C(1,J))**2 + (P(2,I)-C(2,J))**2 +
     $       (P(3,I)-C(3,J))**2
      DIST = DSQRT(DIST)
      EQ = EQ + X(J) / DIST
  810 CONTINUE
      RMS = RMS + (E(I) - EQ)**2
      PERCENT = PERCENT + DABS((E(I) - EQ) / E(I) * HUNDRED)
      IF (IPO(3) .EQ. 1) WRITE(IOUT,*) 'ACTUAL,CALC ',E(I),EQ
  800 CONTINUE
      IF (IPO(3) .EQ. 1) WRITE(IOUT,*) 'SUM OF SQUARES ',RMS
      RMS = DSQRT(RMS) * AU2CAL / MAXPNTS
      PERCENT = PERCENT / MAXPNTS
      IF (IPO(3) .EQ. 1) WRITE(IOUT,*) 'RMS, %',RMS,PERCENT
      RETURN
      END
      SUBROUTINE FMGEN(F,T,M)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/IO/IN,IOUT
C
      DIMENSION F(M)
      DIMENSION GA(35)
C
      EQUIVALENCE (APPROX,OLDSUM)
C
      DATA ZERO/0.0E0/, HALF/0.5E0/, ONE/1.0E0/, TWO/2.0E0/, TEN/10.0E0/
     $     ,PI/3.14159265358979E0/, F42/42.0E0/, F80/80.0E0/
C
 2001 FORMAT(42H1FAILURE IN FMGEN FOR SMALL T:  IX.GT.50, /
     $ 6H IX = ,I3,7H,  T = ,E20.14)
 2002 FORMAT(37H1FAILURE IN FMGEN FOR INTERMEDIATE T,/
     $ 6H  T = ,E20.14)
C
      TEXP=ZERO
      IF(T-F80)2,3,3
    2 TEXP=EXP(-T)
    3 CONTINUE
      IF(T-TEN)10,70,70
C***********************************************************************
C        0 .LT. T .LT. 10
C***********************************************************************
   10 TERM=HALF*GA(M)*TEXP
      TX=ONE
      IX=M+1
      SUM=TX/GA(IX)
      OLDSUM=SUM
   20 IX=IX+1
      TX=TX*T
      IF(IX - 35) 40,40,30
   30 WRITE(IOUT,2001)IX,T
      STOP 'FMGEN'
   40 SUM=SUM+TX/GA(IX)
      IF(TOL-ABS(OLDSUM/SUM-ONE))50,60,60
   50 OLDSUM=SUM
      GO TO 20
   60 F(M)=SUM*TERM
      GO TO 160
C
   70 IF(T-F42)80,150,150
C***********************************************************************
C        10 .LE. T .LT. 42
C***********************************************************************
   80 A=FLOAT(M-1)
      B=A+HALF
      A=A-HALF
      TX=ONE/T
      MM1=M-1
      APPROX=RPITWO*SQRT(TX)*(TX**MM1)
      IF(MM1)90,110,90
   90 DO 100 IX=1,MM1
      B=B-ONE
  100 APPROX=APPROX*B
  110 FIMULT=HALF*TEXP*TX
      SUM=ZERO
      IF(FIMULT)120,140,120
  120 FIPROP=FIMULT/APPROX
      TERM=ONE
      SUM =ONE
      NOTRMS=INT(T)+MM1
      DO 130 IX=2,NOTRMS
      TERM=TERM*A*TX
      SUM=SUM+TERM
      IF(ABS(TERM*FIPROP/SUM)-TOL)140,140,130
  130 A=A-ONE
      WRITE(IOUT,2002)T
      STOP 'FMGEN'
  140 F(M)=APPROX-FIMULT*SUM
      GO TO 160
C***********************************************************************
C        T .GE. 42
C***********************************************************************
  150 TX=FLOAT(M)-HALF
      F(M)=HALF*GA(M)/(T**TX)
C***********************************************************************
C        RECUR DOWNWARDS TO F(1)
C***********************************************************************
  160 TX=T+T
      SUM=FLOAT(M+M-3)
      MM1=M-1
      IF(MM1)170,190,170
  170 DO 180 IX=1,MM1
      F(M-IX)=(TX*F(M-IX+1)+TEXP)/SUM
  180 SUM=SUM-TWO
  190 RETURN
C
      ENTRY FMSET
C
      GA(1)=SQRT(PI)
      TOL=HALF
      DO 200 I=2,35
      GA(I)=GA(I-1)*TOL
  200 TOL=TOL+ONE
      TOL = 5.0E-09
      RPITWO=HALF*GA(1)
      RETURN
      END

      SUBROUTINE INTGRL (H,X1,X2,X3,ICHARG,I6TO5)
C
C     ROUTINE TO CALCULATE THE ELECTRON-CHARGE MATRIX ELEMENTS FOR THE
C     POLARIZATION POTENTIAL. CODE REVISED FROM THE ONE ELECTRON PACKAGE
C     AS IT EXISTED AUGUST, 1983.
C
C
C        REVISED BY M.M. FRANCL JANUARY 1984 FOR PRINCETON CHEMISTRY
C        DEPARTMENT VAX 11/780
C
C        REVISED TO BE COMPATIBLE WITH COMMON /B/ FROM GAUSSIAN 82
C       MAY 1984 M.M. FRANCL
C
C        REVISED TO USE ** BASIS SETS AND THOSE HAVING P ONLY SHELLS
C        JANUARY 1986 M.M. FRANCL
C
C        REVISED FOR GAUSSIAN 86 CHECKPOINT FILES FOR YALE UNIVERSITY
C       FEBRUARY 1988 CURT BRENEMAN
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*4 SHELLA,SHELLN,SHELLT,SHELLC,AOS,AON,SHLADF
C
C+++
      COMMON /MOL/    NATOMS,JCHARG,MULTIP,NAE,NBE,NEL,NBASIS,
     $                IAN(401),ATMCHG(400),C(3,400)
C
C===  Gaussian88 Modification.  New Common /b/ size.
      Common/B/EXX(6000),C1(6000),C2(6000),C3(2000),CF(2000),
     $SHLADF(4000),X(2000),Y(2000),
     $Z(2000),JAN(2000),ShellA(2000),ShellN(2000),ShellT(2000),
     $ShellC(2000),AOS(2000),AON(2000),NShell,MaxTyp
C
C===    Old G86 common /b/
c      COMMON/B/EXX(1200),C1(1200),C2(1200),C3(400),CF(400),SHLADF(800),
c     $         X(400),Y(400),Z(400),JAN(400),SHELLA(400),SHELLN(400),
c     $         SHELLT(400),SHELLC(400),AOS(400),AON(400),NSHELL,MAXTYP
c
C+++
C      COMMON /B/ EXX(240),C1(240),C2(240),C3(80),CF(80),SHLADF(160),
C     $           X(80),Y(80),Z(80),
C     $           JAN(80),SHELLA(80),SHELLN(80),SHELLT(80),SHELLC(80)
C     $          ,AOS(80),AON(80),NSHELL,MAXTYP
C      COMMON /MOL/ NATOMS,JCHARG,MULTIP,NAE,NBE,NEL,NBASIS,IAN(101),
C     $             ATMCHG(100),C(3,100)
      COMMON /IPO/ IPO(10)
      COMMON/IO/ IN,IOUT
C
      DIMENSION H(1)
      DIMENSION RENORM(10)
      DIMENSION OF(9),OX(9),TX(13),ABX(5),ABY(5),ABZ(5),ABSQ(5),
     *A(5),B(5),F(5),APB(5),CPX(5),CPY(5),CPZ(5),FM(5)
      DIMENSION EPN(100)
C
      COMMON/H100/
     $EP00,EP10,EP20,EP30,EP40,EP50,EP60,EP70,EP80,EP90,                 SPDSTV
     $EP01,EP11,EP21,EP31,EP41,EP51,EP61,EP71,EP81,EP91,                 SPDSTV
     $EP02,EP12,EP22,EP32,EP42,EP52,EP62,EP72,EP82,EP92,                 SPDSTV
     $EP03,EP13,EP23,EP33,EP43,EP53,EP63,EP73,EP83,EP93,                 SPDSTV
     $EP04,EP14,EP24,EP34,EP44,EP54,EP64,EP74,EP84,EP94,                 SPDSTV
     $EP05,EP15,EP25,EP35,EP45,EP55,EP65,EP75,EP85,EP95,                 SPDSTV
     $EP06,EP16,EP26,EP36,EP46,EP56,EP66,EP76,EP86,EP96,                 SPDSTV
     $EP07,EP17,EP27,EP37,EP47,EP57,EP67,EP77,EP87,EP97,                 SPDSTV
     $EP08,EP18,EP28,EP38,EP48,EP58,EP68,EP78,EP88,EP98,                 SPDSTV
     $EP09,EP19,EP29,EP39,EP49,EP59,EP69,EP79,EP89,EP99                  SPDSTV
C
      DIMENSION EEP(100)                                                 SPDSTV
      DIMENSION MAX(6)
C                                                                        SPDSTV
C     LOCAL VARIABLES.                                                   SPDSTV
C
      DIMENSION AG(6),CSA(6),CPA(6),CDA(6),                              SPDSTV
     $          BG(6),CSB(6),CPB(6),CDB(6),                              SPDSTV
     $          DPP(9)                                                   SPDSTV
      EQUIVALENCE(OF0,OF(1)),(OF1,OF(2)),(OF2,OF(3)),                    SPDSTV
     $           (OF3,OF(4)),(OF4,OF(5)),(OF5,OF(6)),                    SPDSTV
     $           (OF6,OF(7)),(OF7,OF(8)),(OF8,OF(9))                     SPDSTV
      EQUIVALENCE(OX0,OX(1)),(OX1,OX(2)),(OX2,OX(3)),                    SPDSTV
     $           (OX3,OX(4)),(OX4,OX(5)),(OX5,OX(6)),                    SPDSTV
     $           (OX6,OX(7)),(OX7,OX(8)),(OX8,OX(9))                     SPDSTV
      EQUIVALENCE(A1,A(2)),(A2,A(3)),(A3,A(4)),(A4,A(5))                 SPDSTV
      EQUIVALENCE(B1,B(2)),(B2,B(3)),(B3,B(4)),(B4,B(5))                 SPDSTV
      EQUIVALENCE(T01,T0),(T02,T1),(T03,T2),                             SPDSTV
     $           (T04,T3),(T05,T4),(T06,T5),                             SPDSTV
     $           (T07,T6),(T08,T7),(T09,T8)                              SPDSTV
      EQUIVALENCE(T10,TX(10)),(T11,TX(11)),(T12,TX(12)),(T13,TX(13))     SPDSTV
      EQUIVALENCE(T0,TX(1)),(T1,TX(2)),(T2,TX(3)),                       SPDSTV
     $           (T3,TX(4)),(T4,TX(5)),(T5,TX(6)),                       SPDSTV
     $           (T6,TX(7)),(T7,TX(8)),(T8,TX(9))                        SPDSTV
      EQUIVALENCE(C001,T01),(C050,T02),(C054,T09),                       SPDSTV
     $           (C067,T13),(C068,T08),(C074,T03)                        SPDSTV
      EQUIVALENCE(ABX1,ABX(2)),(ABX2,ABX(3)),                            SPDSTV
     $           (ABX3,ABX(4)),(ABX4,ABX(5))                             SPDSTV
      EQUIVALENCE(AB004,ABX1),(AB006,ABX2),(AB023,ABX3),(AB029,ABX4)     SPDSTV
      EQUIVALENCE(ABY1,ABY(2)),(ABY2,ABY(3)),                            SPDSTV
     $           (ABY3,ABY(4)),(ABY4,ABY(5))                             SPDSTV
      EQUIVALENCE(AB007,ABY1),(AB010,ABY2),(AB032,ABY3),(AB035,ABY4)     SPDSTV
      EQUIVALENCE(ABZ1,ABZ(2)),(ABZ2,ABZ(3)),                            SPDSTV
     $           (ABZ3,ABZ(4)),(ABZ4,ABZ(5))                             SPDSTV
      EQUIVALENCE(AB002,ABZ1),(AB003,ABZ2),(AB011,ABZ3),(AB017,ABZ4)     SPDSTV
      EQUIVALENCE(ABSQ1,ABSQ(2)),(ABSQ2,ABSQ(3)),                        SPDSTV
     $           (ABSQ3,ABSQ(4)),(ABSQ4,ABSQ(5))                         SPDSTV
      EQUIVALENCE(APB1,APB(2)),(APB2,APB(3)),                            SPDSTV
     $           (APB3,APB(4)),(APB4,APB(5))                             SPDSTV
      EQUIVALENCE(CPX1,CPX(2)),(CPX2,CPX(3)),                            SPDSTV
     $           (CPX3,CPX(4)),(CPX4,CPX(5))                             SPDSTV
      EQUIVALENCE(CPY1,CPY(2)),(CPY2,CPY(3)),                            SPDSTV
     $           (CPY3,CPY(4)),(CPY4,CPY(5))                             SPDSTV
      EQUIVALENCE(CPZ1,CPZ(2)),(CPZ2,CPZ(3)),                            SPDSTV
     $           (CPZ3,CPZ(4)),(CPZ4,CPZ(5))                             SPDSTV
      EQUIVALENCE(F1,F(2)),(F2,F(3)),(F3,F(4)),(F4,F(5))                 SPDSTV
      EQUIVALENCE(FM0,FM(1)),(FM1,FM(2)),(FM2,FM(3)),(FM3,FM(4)),        SPDSTV
     $           (FM4,FM(5))                                             SPDSTV
      EQUIVALENCE (D001,FM0)                                             SPDSTV
      EQUIVALENCE(EP00,EEP(1))                                           SPDSTV
C
      DATA MAX/1,4,9,1,4,10/                                             SPDSTV
      DATA TX/1.0E0,0.5E0,0.25E0,0.125E0,0.375E0,0.625E-01,0.1875E0,
     $ 0.75E0,1.5E0,2.25E0,1.125E0,0.0E0,3.0E0/
      DATA ZERO/0.0/,HALF/0.5/,ONE/1.0/,ONEPT5/1.5/,TWO/2.0/,THREE/3.0/,
     *ROOT3/1.732050808/,PI/3.14159265358979/
      DATA ANTOAU /1.889726878D0/
C
 2010 FORMAT(/1X,'ELECTRON-CHARGE MATRIX ELEMENTS'/)
C
C        CALL ROUTINE TO MODIFY COMMON /B/ IF P ONLY SHELLS ARE PRESENT
C
      CALL STAR (NBASIS,SHELLT,SHELLC,AOS,NSHELL,NOSTAR)
C
C*********************************************************************** SPDSTV
C        INITIALIZE THIS SEGMENT.                                        SPDSTV
C*********************************************************************** SPDSTV
C                                                                        SPDSTV
C     ****************************************************************** SPDSTV
C     COMPUTE SIZE OF S T AND V ARRAYS
C     ****************************************************************** SPDSTV
      NTT=(NBASIS*(NBASIS+1))/2
      I5OR6=3
CC    IF(IGO(4) .NE. 0) I5OR6 = 0
C     ****************************************************************** SPDSTV
C     INITIALIZE RENORM  USED TO NORMALIZE D FUNCTIONS
C     ****************************************************************** SPDSTV
      DO 100 I=1,10                                                      SPDSTV
  100 RENORM(I)=ONE                                                      SPDSTV
      RENORM(5)=ROOT3                                                    SPDSTV
      RENORM(8)=ROOT3                                                    SPDSTV
      RENORM(9)=ROOT3                                                    SPDSTV
C     ****************************************************************** SPDSTV
C     CLEAR H ARRAY
C     ****************************************************************** SPDSTV
      DO 50 I=1,NTT                                                      SPDSTV
   50 H(I)=ZERO
C     ****************************************************************** SPDSTV
C     *  INITIALIZE THE VARIABLES USED BY ROUTINE FMGEN.               * SPDSTV
C     ****************************************************************** SPDSTV
      CALL FMSET
      DO 95 I=1,5
   95 FM(I)=ZERO
      ABX(1)=ONE                                                         SPDSTV
      ABY(1)=ONE                                                         SPDSTV
      ABZ(1)=ONE                                                         SPDSTV
      A(1)=ONE                                                           SPDSTV
      B(1)=ONE                                                           SPDSTV
      F(1)=ONE                                                           SPDSTV
      CPX(1)=ONE                                                         SPDSTV
      CPY(1)=ONE                                                         SPDSTV
      CPZ(1)=ONE                                                         SPDSTV
      APB(1)=ONE                                                         SPDSTV
      ABSQ(1)=ONE                                                        SPDSTV
C*********************************************************************** SPDSTV
C        LOOP OVER SHELLS ISHELL AND JSHELL.                             SPDSTV
C*********************************************************************** SPDSTV
      DO 1000 ISHELL=1,NSHELL                                            SPDSTV
      DO 1000 JSHELL=1,ISHELL                                            SPDSTV
      SYMFAC = ONE
C     ****************************************************************** SPDSTV
C     ZERO LOCATIONS                                                     SPDSTV
C     ****************************************************************** SPDSTV
   80 CONTINUE
      DO 9447 JI=1,100                                                   SPDSTV
      EPN(JI)=ZERO                                                       SPDSTV
 9447 CONTINUE                                                           SPDSTV
      IF(SHELLT(ISHELL)-SHELLT(JSHELL))120,120,110                       SPDSTV
  110 INEW=JSHELL                                                        SPDSTV
      JNEW=ISHELL                                                        SPDSTV
      LA=SHELLT(JSHELL)                                                  SPDSTV
      LB=SHELLT(ISHELL)                                                  SPDSTV
      GO TO 200                                                          SPDSTV
  120 INEW=ISHELL                                                        SPDSTV
      JNEW=JSHELL                                                        SPDSTV
      LA=SHELLT(ISHELL)                                                  SPDSTV
      LB=SHELLT(JSHELL)                                                  SPDSTV
  200 CONTINUE                                                           SPDSTV
      LAP1=LA+1                                                          SPDSTV
      LBP1=LB+1                                                          SPDSTV
      LAMAX=MAX(LAP1+I5OR6)                                              SPDSTV
      LBMAX=MAX(LBP1+I5OR6)                                              SPDSTV
      ITYPE=3*LB+LA                                                      SPDSTV
      M=LA+LB+1                                                          SPDSTV
      NGA=SHELLN(INEW)                                                   SPDSTV
      NGB=SHELLN(JNEW)                                                   SPDSTV
      AX=X(INEW)                                                         SPDSTV
      BX=X(JNEW)                                                         SPDSTV
      AY=Y(INEW)                                                         SPDSTV
      BY=Y(JNEW)                                                         SPDSTV
      AZ=Z(INEW)                                                         SPDSTV
      BZ=Z(JNEW)                                                         SPDSTV
      ISHA=SHELLA(INEW)                                                  SPDSTV
      ISHB=SHELLA(JNEW)                                                  SPDSTV
      ISHAD = SHLADF(INEW)
      ISHBD = SHLADF(JNEW)
      IAOS=AOS(INEW)                                                     SPDSTV
      JAOS=AOS(JNEW)                                                     SPDSTV
C     ****************************************************************** SPDSTV
C     OBTAIN INFORMATION ABOUT SHELLS INEW AND JNEW                      SPDSTV
C     ****************************************************************** SPDSTV
      DO 101 I=1,NGA                                                     SPDSTV
      N=ISHA+I-1                                                         SPDSTV
      ND = ISHAD + I -1
      IF (MAXTYP .LE. 1) ND=1
      AG(I)=EXX(N)                                                       SPDSTV
      CSA(I)=C1(N)                                                       SPDSTV
      CPA(I)=C2(N)                                                       SPDSTV
  101 CDA(I)=C3(ND)                                                       SPDSTV

      DO 102 I=1,NGB                                                     SPDSTV
      N=ISHB+I-1                                                         SPDSTV
      ND = ISHBD + I -1
      BG(I)=EXX(N)                                                       SPDSTV
      CSB(I)=C1(N)                                                       SPDSTV
      CPB(I)=C2(N)                                                       SPDSTV
  102 CDB(I)=C3(ND)                                                       SPDSTV

      ABX(2)=BX-AX                                                       SPDSTV
      ABY(2)=BY-AY                                                       SPDSTV
      ABZ(2)=BZ-AZ                                                       SPDSTV
      RABSQ=ABX(2)*ABX(2)+ABY(2)*ABY(2)+ABZ(2)*ABZ(2)                    SPDSTV
      ABSQ(2)=RABSQ                                                      SPDSTV
      DO 103 I=3,5                                                       SPDSTV
      ABX(I)=ABX(I-1)*ABX(2)                                             SPDSTV
      ABY(I)=ABY(I-1)*ABY(2)                                             SPDSTV
      ABZ(I)=ABZ(I-1)*ABZ(2)                                             SPDSTV
  103 ABSQ(I)=ABSQ(I-1)*ABSQ(2)                                          SPDSTV
      AB001=ONE                                                          SPDSTV
      AB005=ABX1*ABZ1                                                    SPDSTV
      AB008=ABY1*ABZ1                                                    SPDSTV
      AB009=ABX1*ABY1                                                    SPDSTV
      AB012=ABX1*ABZ2                                                    SPDSTV
      AB013=ABX2*ABZ1                                                    SPDSTV
      AB014=ABY1*ABZ2                                                    SPDSTV
      AB015=ABX1*ABY1*ABZ1                                               SPDSTV
      AB016=ABY2*ABZ1                                                    SPDSTV
      AB018=ABX1*ABZ3                                                    SPDSTV
      AB019=ABX2*ABZ2                                                    SPDSTV
      AB020=ABY1*ABZ3                                                    SPDSTV
      AB021=ABX1*ABY1*ABZ2                                               SPDSTV
      AB022=ABY2*ABZ2                                                    SPDSTV
      AB024=ABX2*ABY1                                                    SPDSTV
      AB025=ABX1*ABY2                                                    SPDSTV
      AB026=ABX3*ABZ1                                                    SPDSTV
      AB027=ABX2*ABY1*ABZ1                                               SPDSTV
      AB028=ABX1*ABY2*ABZ1                                               SPDSTV
      AB030=ABX3*ABY1                                                    SPDSTV
      AB031=ABX2*ABY2                                                    SPDSTV
      AB033=ABY3*ABZ1                                                    SPDSTV
      AB034=ABX1*ABY3                                                    SPDSTV
C*********************************************************************** SPDSTV
C        LOOP OVER GAUSSIANS  (CONTRACTION LOOP).                        SPDSTV
C*********************************************************************** SPDSTV
      DO 105 IGAUSS=1,NGA                                                SPDSTV
      AA=AG(IGAUSS)                                                      SPDSTV
      DO 105 JGAUSS=1,NGB                                                SPDSTV
      BB=BG(JGAUSS)                                                      SPDSTV
      AAPBB=AA+BB                                                        SPDSTV
      APBB=ONE/AAPBB                                                     SPDSTV
      F2=TWO*AA*BB*APBB                                                  SPDSTV
      PX=(AA*AX+BB*BX)*APBB                                              SPDSTV
      PY=(AA*AY+BB*BY)*APBB                                              SPDSTV
      PZ=(AA*AZ+BB*BZ)*APBB                                              SPDSTV
      A(2)=ONE/AA                                                        SPDSTV
      B(2)=ONE/BB                                                        SPDSTV
      F(2)=F2                                                            SPDSTV
      APB(2)=APBB                                                        SPDSTV
      YX=PI*APBB                                                         SPDSTV
      EXX1=HALF*F2*RABSQ                                                 SPDSTV
      IF(EXX1-80.0E0)4172,4173,4173
 4173 EXX1=ZERO
      GO TO 4714                                                         SPDSTV
 4172 EXX1=EXP(-EXX1)
 4714 CONTINUE                                                           SPDSTV
      OV=(YX**ONEPT5)*EXX1                                               SPDSTV
      OVEK=THREE*AA*BB*APBB                                              SPDSTV
      EK=F2*AA*BB*APBB*OV                                                SPDSTV
      EP=TWO*YX*EXX1                                                     SPDSTV
      DO 119 I=3,5                                                       SPDSTV
      A(I)=A(I-1)*A(2)                                                   SPDSTV
      B(I)=B(I-1)*B(2)                                                   SPDSTV
      APB(I)=APB(I-1)*APB(2)                                             SPDSTV
  119 F(I)=F(I-1)*F(2)                                                   SPDSTV
      DPP(1)=CSA(IGAUSS)*CSB(JGAUSS)                                     SPDSTV
      DPP(2)=CPA(IGAUSS)*CSB(JGAUSS)                                     SPDSTV
      DPP(3)=CDA(IGAUSS)*CSB(JGAUSS)                                     SPDSTV
      DPP(4)=CSA(IGAUSS)*CPB(JGAUSS)                                     SPDSTV
      DPP(5)=CPA(IGAUSS)*CPB(JGAUSS)                                     SPDSTV
      DPP(6)=CDA(IGAUSS)*CPB(JGAUSS)                                     SPDSTV
      DPP(7)=CSA(IGAUSS)*CDB(JGAUSS)                                     SPDSTV
      DPP(8)=CPA(IGAUSS)*CDB(JGAUSS)                                     SPDSTV
      DPP(9)=CDA(IGAUSS)*CDB(JGAUSS)                                     SPDSTV
      DO 2132 I=1,9                                                      SPDSTV
      OF(I)=DPP(I)*OV                                                    SPDSTV
 2132 OX(I)=DPP(I)*EK                                                    SPDSTV
      DO 2139 I=1,100
 2139 EEP(I)=ZERO
      C002=T02*A1*F1                                                     SPDSTV
      C006=T02*B1*F1                                                     SPDSTV
      C007=T03*A1*B1*F2                                                  SPDSTV
      C008=T03*A1*B1*F1                                                  SPDSTV
      C027=T01*A1                                                        SPDSTV
      C031=T01*A1*B1*F1                                                  SPDSTV
      C032=T02*A1*B1                                                     SPDSTV
      C051=T02*A1*B1*F2                                                  SPDSTV
      C012=T02*B1                                                        SPDSTV
      C013=T03*B2*F2                                                     SPDSTV
      C014=T03*B2*F1                                                     SPDSTV
      C036=T01*B2*F1                                                     SPDSTV
      C037=T02*B2                                                        SPDSTV
      C056=T01*B1*F1                                                     SPDSTV
      C030=T01*B1                                                        SPDSTV
      C018=T04*A1*B2*F2                                                  SPDSTV
      IF(ITYPE-7)3060,3040,3041                                          SPDSTV
 3041 CONTINUE                                                           SPDSTV
      C003=T02*A1                                                        SPDSTV
      C004=T03*A2*F2                                                     SPDSTV
      C005=T03*A2*F1                                                     SPDSTV
      C009=T04*A2*B1*F3                                                  SPDSTV
      C010=T05*A2*B1*F2                                                  SPDSTV
      C011=T04*A2*B1*F2                                                  SPDSTV
      C017=T03*A1*B1                                                     SPDSTV
      C019=T04*A1*B2*F1                                                  SPDSTV
      C020=T04*A2*B1*F1                                                  SPDSTV
      C021=T06*A2*B2*F4                                                  SPDSTV
      C022=T05*A2*B2*F3                                                  SPDSTV
      C023=T07*A2*B2*F2                                                  SPDSTV
      C024=T07*A2*B2*F3                                                  SPDSTV
      C025=T06*A2*B2*F3                                                  SPDSTV
      C026=T06*A2*B2*F2                                                  SPDSTV
      C028=T01*A2*F1                                                     SPDSTV
      C029=T02*A2                                                        SPDSTV
      C033=T08*A2*B1*F2                                                  SPDSTV
      C034=T09*A2*B1*F1                                                  SPDSTV
      C035=T02*A2*B1*F1                                                  SPDSTV
      C040=T02*A1*B2*F1                                                  SPDSTV
      C041=T03*A1*B2                                                     SPDSTV
      C042=T03*A2*B1                                                     SPDSTV
      C043=T02*A2*B2*F3                                                  SPDSTV
      C044=T10*A2*B2*F2                                                  SPDSTV
      C045=T08*A2*B2*F1                                                  SPDSTV
      C046=T11*A2*B2*F2                                                  SPDSTV
      C047=T05*A2*B2*F2                                                  SPDSTV
      C048=T03*A2*B2*F1                                                  SPDSTV
      C049=T01*A1*F1                                                     SPDSTV
      C057=T12*A1*B1*F1                                                  SPDSTV
      C058=T03*A1                                                        SPDSTV
      C059=T03*B1                                                        SPDSTV
      C060=T03*A1*B2*F3                                                  SPDSTV
      C061=T04*B2*F2                                                     SPDSTV
      C062=T04*B2*F1                                                     SPDSTV
      C063=T03*A2*B1*F3                                                  SPDSTV
      C064=T01*A1*B1*F2                                                  SPDSTV
      C065=T09*B1*F1                                                     SPDSTV
      C066=T09*A1*F1                                                     SPDSTV
      C069=T04*A2*F2                                                     SPDSTV
      C070=T04*A2*F1                                                     SPDSTV
      C071=T03*A2*B1*F2                                                  SPDSTV
      C072=T08*A1*F1                                                     SPDSTV
      C073=T03*A1*B2*F2                                                  SPDSTV
      C075=T08*B1*F1                                                     SPDSTV
      C076=T04*A1*B1*F2                                                  SPDSTV
 3040 CONTINUE                                                           SPDSTV
      C015=T04*A1*B2*F3                                                  SPDSTV
      C016=T05*A1*B2*F2                                                  SPDSTV
      C038=T08*A1*B2*F2                                                  SPDSTV
      C039=T09*A1*B2*F1                                                  SPDSTV
      C040=T02*A1*B2*F1                                                  SPDSTV
      C052=T02*A1*B1*F1                                                  SPDSTV
      C053=T03*B1*F1                                                     SPDSTV
      C055=T03*A1*F1                                                     SPDSTV
 3060 CONTINUE
      CX=X1
      CY=X2
      CZ=X3
      CPX(2)=PX-CX                                                       SPDSTV
      CPY(2)=PY-CY                                                       SPDSTV
      CPZ(2)=PZ-CZ                                                       SPDSTV
      CP2=CPX(2)*CPX(2)+CPY(2)*CPY(2)+CPZ(2)*CPZ(2)                      SPDSTV
      CALL FMGEN(FM,AAPBB*CP2,M)                                         SPDSTV
      DO 108 I=3,5                                                       SPDSTV
      CPX(I)=CPX(I-1)*CPX(2)                                             SPDSTV
      CPY(I)=CPY(I-1)*CPY(2)                                             SPDSTV
  108 CPZ(I)=CPZ(I-1)*CPZ(2)                                             SPDSTV
      EPAN=EP*FLOAT(-ICHARG)
      DO 2136 I=1,9                                                      SPDSTV
 2136 OF(I)=DPP(I)*EPAN                                                  SPDSTV
      D002=CPZ1*FM1                                                      SPDSTV
      D003=CPZ2*FM2                                                      SPDSTV
      D004=APB1*FM1                                                      SPDSTV
      D005=CPX1*FM1                                                      SPDSTV
      D006=CPX1*CPZ1*FM2                                                 SPDSTV
      D007=CPX2*FM2                                                      SPDSTV
      D008=CPY1*FM1                                                      SPDSTV
      D009=CPY1*CPZ1*FM2                                                 SPDSTV
      D010=CPX1*CPY1*FM2                                                 SPDSTV
      D011=CPY2*FM2                                                      SPDSTV
      D012=CPZ3*FM3                                                      SPDSTV
      D013=APB1*CPZ1*FM2                                                 SPDSTV
      D014=CPX1*CPZ2*FM3                                                 SPDSTV
      D015=APB1*CPX1*FM2                                                 SPDSTV
      D016=CPX2*CPZ1*FM3                                                 SPDSTV
      D017=CPY1*CPZ2*FM3                                                 SPDSTV
      D018=APB1*CPY1*FM2                                                 SPDSTV
      D019=CPX1*CPY1*CPZ1*FM3                                            SPDSTV
      D020=CPY2*CPZ1*FM3                                                 SPDSTV
      D034=CPX3*FM3                                                      SPDSTV
      D035=CPX2*CPY1*FM3                                                 SPDSTV
      D036=CPX1*CPY2*FM3                                                 SPDSTV
      D043=CPY3*FM3                                                      SPDSTV
C     ****************************************************************** SPDSTV
C     *                                SS                              * SPDSTV
C     ****************************************************************** SPDSTV
      EP00=OF0*(+C001*AB001*D001)                                        SPDSTV
      IF(ITYPE)3230,3262,3230                                            SPDSTV
C     ****************************************************************** SPDSTV
C     *                                SP                              * SPDSTV
C     ****************************************************************** SPDSTV
 3230 CONTINUE                                                           SPDSTV
      EP01=OF3*(-C006*AB002*D001-C001*AB001*D002)                        SPDSTV
      EP03=OF3*(-C006*AB004*D001-C001*AB001*D005)                        SPDSTV
      EP06=OF3*(-C006*AB007*D001-C001*AB001*D008)                        SPDSTV
      IF(ITYPE-7)3240,3242,3241                                          SPDSTV
 3240 IF(ITYPE-4)3262,3261,3260                                          SPDSTV
C     ****************************************************************** SPDSTV
C     *                                DD                              * SPDSTV
C     ****************************************************************** SPDSTV
 3241 CONTINUE                                                           SPDSTV
      D021=CPZ4*FM4                                                      SPDSTV
      D022=APB1*CPZ2*FM3                                                 SPDSTV
      D023=APB2*FM2                                                      SPDSTV
      D024=CPX1*CPZ3*FM4                                                 SPDSTV
      D025=APB1*CPX1*CPZ1*FM3                                            SPDSTV
      D026=CPX2*CPZ2*FM4                                                 SPDSTV
      D027=APB1*CPX2*FM3                                                 SPDSTV
      D028=CPY1*CPZ3*FM4                                                 SPDSTV
      D029=APB1*CPY1*CPZ1*FM3                                            SPDSTV
      D030=CPX1*CPY1*CPZ2*FM4                                            SPDSTV
      D031=APB1*CPX1*CPY1*FM3                                            SPDSTV
      D032=CPY2*CPZ2*FM4                                                 SPDSTV
      D033=APB1*CPY2*FM3                                                 SPDSTV
      D037=CPX3*CPZ1*FM4                                                 SPDSTV
      D038=CPX2*CPY1*CPZ1*FM4                                            SPDSTV
      D039=CPX1*CPY2*CPZ1*FM4                                            SPDSTV
      D040=CPX4*FM4                                                      SPDSTV
      D041=CPX3*CPY1*FM4                                                 SPDSTV
      D042=CPX2*CPY2*FM4                                                 SPDSTV
      D044=CPY3*CPZ1*FM4                                                 SPDSTV
      D045=CPX1*CPY3*FM4                                                 SPDSTV
      D046=CPY4*FM4                                                      SPDSTV
      EP20=OF2*(+C003*AB001*D001+C004*AB003*D001-C005*AB001*D001-C049*AB SPDSTV
     $002*D002+C001*AB001*D003-C050*AB001*D004)                          SPDSTV
      EP40=OF2*(+C004*AB005*D001-C002*AB004*D002-C002*AB002*D005+C001*AB SPDSTV
     $001*D006)                                                          SPDSTV
      EP50=OF2*(+C003*AB001*D001+C004*AB006*D001-C005*AB001*D001-C049*AB SPDSTV
     $004*D005+C001*AB001*D007-C050*AB001*D004)                          SPDSTV
      EP70=OF2*(+C004*AB008*D001-C002*AB007*D002-C002*AB002*D008+C001*AB SPDSTV
     $001*D009)                                                          SPDSTV
      EP80=OF2*(+C004*AB009*D001-C002*AB004*D008-C002*AB007*D005+C001*AB SPDSTV
     $001*D010)                                                          SPDSTV
      EP90=OF2*(+C003*AB001*D001+C004*AB010*D001-C005*AB001*D001-C049*AB SPDSTV
     $007*D008+C001*AB001*D011-C050*AB001*D004)                          SPDSTV
      EP21=OF5*(-C008*AB002*D001-C003*AB001*D002-C009*AB011*D001+C010*AB SPDSTV
     $002*D001+C051*AB003*D002-C052*AB001*D002-C006*AB002*D003+C053*AB00 SPDSTV
     $2*D004-C004*AB003*D002+C005*AB001*D002+C049*AB002*D003-C002*AB002* SPDSTV
     $D004-C001*AB001*D012+C054*AB001*D013)                              SPDSTV
      EP41=OF5*(-C009*AB012*D001+C011*AB004*D001+C007*AB005*D002+C007*AB SPDSTV
     $003*D005-C008*AB001*D005-C006*AB002*D006-C004*AB005*D002+C002*AB00 SPDSTV
     $4*D003-C055*AB004*D004+C002*AB002*D006-C001*AB001*D014+C050*AB001* SPDSTV
     $D015)                                                              SPDSTV
      EP51=OF5*(-C008*AB002*D001-C003*AB001*D002-C009*AB013*D001+C011*AB SPDSTV
     $002*D001+C051*AB005*D005-C006*AB002*D007+C053*AB002*D004-C004*AB00 SPDSTV
     $6*D002+C005*AB001*D002+C049*AB004*D006-C001*AB001*D016+C050*AB001* SPDSTV
     $D013)                                                              SPDSTV
      EP71=OF5*(-C009*AB014*D001+C011*AB007*D001+C007*AB008*D002+C007*AB SPDSTV
     $003*D008-C008*AB001*D008-C006*AB002*D009-C004*AB008*D002+C002*AB00 SPDSTV
     $7*D003-C055*AB007*D004+C002*AB002*D009-C001*AB001*D017+C050*AB001* SPDSTV
     $D018)                                                              SPDSTV
      EP81=OF5*(-C009*AB015*D001+C007*AB005*D008+C007*AB008*D005-C006*AB SPDSTV
     $002*D010-C004*AB009*D002+C002*AB004*D009+C002*AB007*D006-C001*AB00 SPDSTV
     $1*D019)                                                            SPDSTV
      EP91=OF5*(-C008*AB002*D001-C003*AB001*D002-C009*AB016*D001+C011*AB SPDSTV
     $002*D001+C051*AB008*D008-C006*AB002*D011+C053*AB002*D004-C004*AB01 SPDSTV
     $0*D002+C005*AB001*D002+C049*AB007*D009-C001*AB001*D020+C050*AB001* SPDSTV
     $D013)                                                              SPDSTV
      EP22=OF8*(+C017*AB001*D001+C018*AB003*D001-C019*AB001*D001+C057*AB SPDSTV
     $002*D002+C003*AB001*D003-C058*AB001*D004+C011*AB003*D001-C020*AB00 SPDSTV
     $1*D001+C012*AB001*D003-C059*AB001*D004+C021*AB017*D001-C022*AB003* SPDSTV
     $D001-C060*AB011*D002+C023*AB001*D001+C038*AB002*D002+C013*AB003*D0 SPDSTV
     $03-C061*AB003*D004-C014*AB001*D003+C062*AB001*D004+C063*AB011*D002 SPDSTV
     $-C033*AB002*D002-C064*AB003*D003+C051*AB003*D004+C031*AB001*D003-C SPDSTV
     $052*AB001*D004+C056*AB002*D012-C065*AB002*D013+C004*AB003*D003-C00 SPDSTV
     $5*AB001*D003-C049*AB002*D012+C066*AB002*D013+C001*AB001*D021-C067* SPDSTV
     $AB001*D022+C068*AB001*D023-C069*AB003*D004+C070*AB001*D004)        SPDSTV
      EP42=OF8*(+C011*AB005*D001-C008*AB004*D002-C008*AB002*D005+C012*AB SPDSTV
     $001*D006+C021*AB018*D001-C024*AB005*D001-C015*AB012*D002-C015*AB01 SPDSTV
     $1*D005+C016*AB002*D005+C013*AB003*D006+C018*AB004*D002-C014*AB001* SPDSTV
     $D006+C063*AB012*D002-C071*AB004*D002-C051*AB005*D003+C007*AB005*D0 SPDSTV
     $04-C051*AB003*D006+C052*AB001*D006+C056*AB002*D014-C006*AB002*D015 SPDSTV
     $+C004*AB005*D003-C002*AB004*D012+C072*AB004*D013-C002*AB002*D014+C SPDSTV
     $001*AB001*D024-C054*AB001*D025-C069*AB005*D004+C055*AB002*D015)    SPDSTV
      EP52=OF8*(+C017*AB001*D001+C018*AB003*D001-C019*AB001*D001+C052*AB SPDSTV
     $002*D002+C003*AB001*D003-C058*AB001*D004+C011*AB006*D001-C020*AB00 SPDSTV
     $1*D001-C052*AB004*D005+C012*AB001*D007-C059*AB001*D004+C021*AB019* SPDSTV
     $D001-C025*AB003*D001-C060*AB012*D005+C013*AB003*D007-C061*AB003*D0 SPDSTV
     $04-C025*AB006*D001+C026*AB001*D001+C073*AB004*D005-C014*AB001*D007 SPDSTV
     $+C062*AB001*D004+C063*AB013*D002-C071*AB002*D002-C064*AB005*D006+C SPDSTV
     $056*AB002*D016-C006*AB002*D013+C004*AB006*D003-C005*AB001*D003-C04 SPDSTV
     $9*AB004*D014+C001*AB001*D026-C050*AB001*D022-C069*AB006*D004+C070* SPDSTV
     $AB001*D004+C002*AB004*D015-C050*AB001*D027+C074*AB001*D023)        SPDSTV
      EP72=OF8*(+C011*AB008*D001-C008*AB007*D002-C008*AB002*D008+C012*AB SPDSTV
     $001*D009+C021*AB020*D001-C024*AB008*D001-C015*AB014*D002-C015*AB01 SPDSTV
     $1*D008+C016*AB002*D008+C013*AB003*D009+C018*AB007*D002-C014*AB001* SPDSTV
     $D009+C063*AB014*D002-C071*AB007*D002-C051*AB008*D003+C007*AB008*D0 SPDSTV
     $04-C051*AB003*D009+C052*AB001*D009+C056*AB002*D017-C006*AB002*D018 SPDSTV
     $+C004*AB008*D003-C002*AB007*D012+C072*AB007*D013-C002*AB002*D017+C SPDSTV
     $001*AB001*D028-C054*AB001*D029-C069*AB008*D004+C055*AB002*D018)    SPDSTV
      EP82=OF8*(+C011*AB009*D001-C008*AB004*D008-C008*AB007*D005+C012*AB SPDSTV
     $001*D010+C021*AB021*D001-C015*AB012*D008-C015*AB014*D005+C013*AB00 SPDSTV
     $3*D010-C025*AB009*D001+C018*AB004*D008+C018*AB007*D005-C014*AB001* SPDSTV
     $D010+C063*AB015*D002-C051*AB005*D009-C051*AB008*D006+C056*AB002*D0 SPDSTV
     $19+C004*AB009*D003-C002*AB004*D017-C002*AB007*D014+C001*AB001*D030 SPDSTV
     $-C069*AB009*D004+C055*AB004*D018+C055*AB007*D015-C050*AB001*D031)  SPDSTV
      EP92=OF8*(+C017*AB001*D001+C018*AB003*D001-C019*AB001*D001+C052*AB SPDSTV
     $002*D002+C003*AB001*D003-C058*AB001*D004+C011*AB010*D001-C020*AB00 SPDSTV
     $1*D001-C052*AB007*D008+C012*AB001*D011-C059*AB001*D004+C021*AB022* SPDSTV
     $D001-C025*AB003*D001-C060*AB014*D008+C013*AB003*D011-C061*AB003*D0 SPDSTV
     $04-C025*AB010*D001+C026*AB001*D001+C073*AB007*D008-C014*AB001*D011 SPDSTV
     $+C062*AB001*D004+C063*AB016*D002-C071*AB002*D002-C064*AB008*D009+C SPDSTV
     $056*AB002*D020-C006*AB002*D013+C004*AB010*D003-C005*AB001*D003-C04 SPDSTV
     $9*AB007*D017+C001*AB001*D032-C050*AB001*D022-C069*AB010*D004+C070* SPDSTV
     $AB001*D004+C002*AB007*D018-C050*AB001*D033+C074*AB001*D023)        SPDSTV
      EP23=OF5*(-C008*AB004*D001-C003*AB001*D005-C009*AB012*D001+C011*AB SPDSTV
     $004*D001+C051*AB005*D002-C006*AB004*D003+C053*AB004*D004-C004*AB00 SPDSTV
     $3*D005+C005*AB001*D005+C049*AB002*D006-C001*AB001*D014+C050*AB001* SPDSTV
     $D015)                                                              SPDSTV
      EP43=OF5*(-C009*AB013*D001+C007*AB006*D002+C011*AB002*D001-C008*AB SPDSTV
     $001*D002+C007*AB005*D005-C006*AB004*D006-C004*AB005*D005+C002*AB00 SPDSTV
     $4*D006+C002*AB002*D007-C001*AB001*D016-C055*AB002*D004+C050*AB001* SPDSTV
     $D013)                                                              SPDSTV
      EP53=OF5*(-C008*AB004*D001-C003*AB001*D005-C009*AB023*D001+C010*AB SPDSTV
     $004*D001+C051*AB006*D005-C052*AB001*D005-C006*AB004*D007+C053*AB00 SPDSTV
     $4*D004-C004*AB006*D005+C005*AB001*D005+C049*AB004*D007-C002*AB004* SPDSTV
     $D004-C001*AB001*D034+C054*AB001*D015)                              SPDSTV
      EP73=OF5*(-C009*AB015*D001+C007*AB009*D002+C007*AB005*D008-C006*AB SPDSTV
     $004*D009-C004*AB008*D005+C002*AB007*D006+C002*AB002*D010-C001*AB00 SPDSTV
     $1*D019)                                                            SPDSTV
      EP83=OF5*(-C009*AB024*D001+C007*AB006*D008+C011*AB007*D001-C008*AB SPDSTV
     $001*D008+C007*AB009*D005-C006*AB004*D010-C004*AB009*D005+C002*AB00 SPDSTV
     $4*D010+C002*AB007*D007-C001*AB001*D035-C055*AB007*D004+C050*AB001* SPDSTV
     $D018)                                                              SPDSTV
      EP93=OF5*(-C008*AB004*D001-C003*AB001*D005-C009*AB025*D001+C011*AB SPDSTV
     $004*D001+C051*AB009*D008-C006*AB004*D011+C053*AB004*D004-C004*AB01 SPDSTV
     $0*D005+C005*AB001*D005+C049*AB007*D010-C001*AB001*D036+C050*AB001* SPDSTV
     $D015)                                                              SPDSTV
      EP24=OF8*(+C018*AB005*D001+C008*AB004*D002+C008*AB002*D005+C003*AB SPDSTV
     $001*D006+C021*AB018*D001-C024*AB005*D001-C060*AB012*D002+C073*AB00 SPDSTV
     $4*D002+C013*AB005*D003-C061*AB005*D004+C009*AB012*D002-C011*AB004* SPDSTV
     $D002-C051*AB005*D003+C007*AB005*D004+C006*AB004*D012-C075*AB004*D0 SPDSTV
     $13+C009*AB011*D005-C010*AB002*D005-C051*AB003*D006+C052*AB001*D006 SPDSTV
     $+C006*AB002*D014-C053*AB002*D015+C004*AB003*D006-C005*AB001*D006-C SPDSTV
     $049*AB002*D014+C002*AB002*D015+C001*AB001*D024-C054*AB001*D025)    SPDSTV
      EP44=OF8*(+C021*AB019*D001-C025*AB006*D001-C015*AB013*D002-C025*AB SPDSTV
     $003*D001+C026*AB001*D001+C018*AB002*D002-C015*AB012*D005+C018*AB00 SPDSTV
     $4*D005+C013*AB005*D006+C009*AB013*D002-C007*AB006*D003+C076*AB006* SPDSTV
     $D004-C011*AB002*D002+C008*AB001*D003-C008*AB001*D004-C051*AB005*D0 SPDSTV
     $06+C006*AB004*D014-C053*AB004*D015+C009*AB012*D005-C011*AB004*D005 SPDSTV
     $-C007*AB003*D007+C008*AB001*D007+C006*AB002*D016+C076*AB003*D004-C SPDSTV
     $053*AB002*D013+C004*AB005*D006-C002*AB004*D014+C055*AB004*D015-C00 SPDSTV
     $2*AB002*D016+C001*AB001*D026-C050*AB001*D027+C055*AB002*D013-C050* SPDSTV
     $AB001*D022+C074*AB001*D023)                                        SPDSTV
      EP54=OF8*(+C018*AB005*D001+C008*AB004*D002+C008*AB002*D005+C003*AB SPDSTV
     $001*D006+C021*AB026*D001-C024*AB005*D001-C060*AB013*D005+C073*AB00 SPDSTV
     $2*D005+C013*AB005*D007-C061*AB005*D004+C009*AB023*D002-C010*AB004* SPDSTV
     $D002-C051*AB006*D006+C052*AB001*D006+C006*AB004*D016-C053*AB004*D0 SPDSTV
     $13+C009*AB013*D005-C011*AB002*D005-C051*AB005*D007+C007*AB005*D004 SPDSTV
     $+C006*AB002*D034-C075*AB002*D015+C004*AB006*D006-C005*AB001*D006-C SPDSTV
     $049*AB004*D016+C002*AB004*D013+C001*AB001*D037-C054*AB001*D025)    SPDSTV
      EP74=OF8*(+C021*AB021*D001-C025*AB009*D001-C015*AB015*D002-C015*AB SPDSTV
     $012*D008+C018*AB004*D008+C013*AB005*D009+C009*AB015*D002-C007*AB00 SPDSTV
     $9*D003+C076*AB009*D004-C007*AB005*D009+C006*AB004*D017-C053*AB004* SPDSTV
     $D018+C009*AB014*D005-C011*AB007*D005-C007*AB008*D006-C007*AB003*D0 SPDSTV
     $10+C008*AB001*D010+C006*AB002*D019+C004*AB008*D006-C002*AB007*D014 SPDSTV
     $+C055*AB007*D015-C002*AB002*D019+C001*AB001*D030-C050*AB001*D031)  SPDSTV
      EP84=OF8*(+C021*AB027*D001-C015*AB013*D008-C025*AB008*D001+C018*AB SPDSTV
     $002*D008-C015*AB015*D005+C013*AB005*D010+C009*AB024*D002-C007*AB00 SPDSTV
     $6*D009-C011*AB007*D002+C008*AB001*D009-C007*AB009*D006+C006*AB004* SPDSTV
     $D019+C009*AB015*D005-C007*AB005*D010-C007*AB008*D007+C006*AB002*D0 SPDSTV
     $35+C076*AB008*D004-C053*AB002*D018+C004*AB009*D006-C002*AB004*D019 SPDSTV
     $-C002*AB007*D016+C001*AB001*D038+C055*AB007*D013-C050*AB001*D029)  SPDSTV
      EP94=OF8*(+C018*AB005*D001+C008*AB004*D002+C008*AB002*D005+C003*AB SPDSTV
     $001*D006+C021*AB028*D001-C025*AB005*D001-C060*AB015*D008+C013*AB00 SPDSTV
     $5*D011-C061*AB005*D004+C009*AB025*D002-C011*AB004*D002-C051*AB009* SPDSTV
     $D009+C006*AB004*D020-C053*AB004*D013+C009*AB016*D005-C011*AB002*D0 SPDSTV
     $05-C051*AB008*D010+C006*AB002*D036-C053*AB002*D015+C004*AB010*D006 SPDSTV
     $-C005*AB001*D006-C049*AB007*D019+C001*AB001*D039-C050*AB001*D025)  SPDSTV
      EP25=OF8*(+C017*AB001*D001+C018*AB006*D001-C019*AB001*D001+C052*AB SPDSTV
     $004*D005+C003*AB001*D007-C058*AB001*D004+C011*AB003*D001-C020*AB00 SPDSTV
     $1*D001-C052*AB002*D002+C012*AB001*D003-C059*AB001*D004+C021*AB019* SPDSTV
     $D001-C025*AB006*D001-C060*AB013*D002+C013*AB006*D003-C061*AB006*D0 SPDSTV
     $04-C025*AB003*D001+C026*AB001*D001+C073*AB002*D002-C014*AB001*D003 SPDSTV
     $+C062*AB001*D004+C063*AB012*D005-C071*AB004*D005-C064*AB005*D006+C SPDSTV
     $056*AB004*D014-C006*AB004*D015+C004*AB003*D007-C005*AB001*D007-C04 SPDSTV
     $9*AB002*D016+C001*AB001*D026-C050*AB001*D027-C069*AB003*D004+C070* SPDSTV
     $AB001*D004+C002*AB002*D013-C050*AB001*D022+C074*AB001*D023)        SPDSTV
      EP45=OF8*(+C011*AB005*D001-C008*AB004*D002-C008*AB002*D005+C012*AB SPDSTV
     $001*D006+C021*AB026*D001-C015*AB023*D002-C024*AB005*D001+C016*AB00 SPDSTV
     $4*D002-C015*AB013*D005+C013*AB006*D006+C018*AB002*D005-C014*AB001* SPDSTV
     $D006+C063*AB013*D005-C051*AB006*D006-C071*AB002*D005+C052*AB001*D0 SPDSTV
     $06-C051*AB005*D007+C056*AB004*D016+C007*AB005*D004-C006*AB004*D013 SPDSTV
     $+C004*AB005*D007-C002*AB004*D016-C002*AB002*D034+C001*AB001*D037+C SPDSTV
     $072*AB002*D015-C054*AB001*D025-C069*AB005*D004+C055*AB004*D013)    SPDSTV
      EP55=OF8*(+C017*AB001*D001+C018*AB006*D001-C019*AB001*D001+C057*AB SPDSTV
     $004*D005+C003*AB001*D007-C058*AB001*D004+C011*AB006*D001-C020*AB00 SPDSTV
     $1*D001+C012*AB001*D007-C059*AB001*D004+C021*AB029*D001-C022*AB006* SPDSTV
     $D001-C060*AB023*D005+C023*AB001*D001+C038*AB004*D005+C013*AB006*D0 SPDSTV
     $07-C061*AB006*D004-C014*AB001*D007+C062*AB001*D004+C063*AB023*D005 SPDSTV
     $-C033*AB004*D005-C064*AB006*D007+C051*AB006*D004+C031*AB001*D007-C SPDSTV
     $052*AB001*D004+C056*AB004*D034-C065*AB004*D015+C004*AB006*D007-C00 SPDSTV
     $5*AB001*D007-C049*AB004*D034+C066*AB004*D015+C001*AB001*D040-C067* SPDSTV
     $AB001*D027+C068*AB001*D023-C069*AB006*D004+C070*AB001*D004)        SPDSTV
      EP75=OF8*(+C011*AB008*D001-C008*AB007*D002-C008*AB002*D008+C012*AB SPDSTV
     $001*D009+C021*AB027*D001-C015*AB024*D002-C015*AB013*D008+C013*AB00 SPDSTV
     $6*D009-C025*AB008*D001+C018*AB007*D002+C018*AB002*D008-C014*AB001* SPDSTV
     $D009+C063*AB015*D005-C051*AB009*D006-C051*AB005*D010+C056*AB004*D0 SPDSTV
     $19+C004*AB008*D007-C002*AB007*D016-C002*AB002*D035+C001*AB001*D038 SPDSTV
     $-C069*AB008*D004+C055*AB007*D013+C055*AB002*D018-C050*AB001*D029)  SPDSTV
      EP85=OF8*(+C011*AB009*D001-C008*AB004*D008-C008*AB007*D005+C012*AB SPDSTV
     $001*D010+C021*AB030*D001-C015*AB023*D008-C024*AB009*D001+C016*AB00 SPDSTV
     $4*D008-C015*AB024*D005+C013*AB006*D010+C018*AB007*D005-C014*AB001* SPDSTV
     $D010+C063*AB024*D005-C051*AB006*D010-C071*AB007*D005+C052*AB001*D0 SPDSTV
     $10-C051*AB009*D007+C056*AB004*D035+C007*AB009*D004-C006*AB004*D018 SPDSTV
     $+C004*AB009*D007-C002*AB004*D035-C002*AB007*D034+C001*AB001*D041+C SPDSTV
     $072*AB007*D015-C054*AB001*D031-C069*AB009*D004+C055*AB004*D018)    SPDSTV
      EP95=OF8*(+C017*AB001*D001+C018*AB006*D001-C019*AB001*D001+C052*AB SPDSTV
     $004*D005+C003*AB001*D007-C058*AB001*D004+C011*AB010*D001-C020*AB00 SPDSTV
     $1*D001-C052*AB007*D008+C012*AB001*D011-C059*AB001*D004+C021*AB031* SPDSTV
     $D001-C025*AB006*D001-C060*AB024*D008+C013*AB006*D011-C061*AB006*D0 SPDSTV
     $04-C025*AB010*D001+C026*AB001*D001+C073*AB007*D008-C014*AB001*D011 SPDSTV
     $+C062*AB001*D004+C063*AB025*D005-C071*AB004*D005-C064*AB009*D010+C SPDSTV
     $056*AB004*D036-C006*AB004*D015+C004*AB010*D007-C005*AB001*D007-C04 SPDSTV
     $9*AB007*D035+C001*AB001*D042-C050*AB001*D027-C069*AB010*D004+C070* SPDSTV
     $AB001*D004+C002*AB007*D018-C050*AB001*D033+C074*AB001*D023)        SPDSTV
      EP26=OF5*(-C008*AB007*D001-C003*AB001*D008-C009*AB014*D001+C011*AB SPDSTV
     $007*D001+C051*AB008*D002-C006*AB007*D003+C053*AB007*D004-C004*AB00 SPDSTV
     $3*D008+C005*AB001*D008+C049*AB002*D009-C001*AB001*D017+C050*AB001* SPDSTV
     $D018)                                                              SPDSTV
      EP46=OF5*(-C009*AB015*D001+C007*AB009*D002+C007*AB008*D005-C006*AB SPDSTV
     $007*D006-C004*AB005*D008+C002*AB004*D009+C002*AB002*D010-C001*AB00 SPDSTV
     $1*D019)                                                            SPDSTV
      EP56=OF5*(-C008*AB007*D001-C003*AB001*D008-C009*AB024*D001+C011*AB SPDSTV
     $007*D001+C051*AB009*D005-C006*AB007*D007+C053*AB007*D004-C004*AB00 SPDSTV
     $6*D008+C005*AB001*D008+C049*AB004*D010-C001*AB001*D035+C050*AB001* SPDSTV
     $D018)                                                              SPDSTV
      EP76=OF5*(-C009*AB016*D001+C007*AB010*D002+C011*AB002*D001-C008*AB SPDSTV
     $001*D002+C007*AB008*D008-C006*AB007*D009-C004*AB008*D008+C002*AB00 SPDSTV
     $7*D009+C002*AB002*D011-C001*AB001*D020-C055*AB002*D004+C050*AB001* SPDSTV
     $D013)                                                              SPDSTV
      EP86=OF5*(-C009*AB025*D001+C011*AB004*D001+C007*AB009*D008+C007*AB SPDSTV
     $010*D005-C008*AB001*D005-C006*AB007*D010-C004*AB009*D008+C002*AB00 SPDSTV
     $4*D011-C055*AB004*D004+C002*AB007*D010-C001*AB001*D036+C050*AB001* SPDSTV
     $D015)                                                              SPDSTV
      EP96=OF5*(-C008*AB007*D001-C003*AB001*D008-C009*AB032*D001+C010*AB SPDSTV
     $007*D001+C051*AB010*D008-C052*AB001*D008-C006*AB007*D011+C053*AB00 SPDSTV
     $7*D004-C004*AB010*D008+C005*AB001*D008+C049*AB007*D011-C002*AB007* SPDSTV
     $D004-C001*AB001*D043+C054*AB001*D018)                              SPDSTV
      EP27=OF8*(+C018*AB008*D001+C008*AB007*D002+C008*AB002*D008+C003*AB SPDSTV
     $001*D009+C021*AB020*D001-C024*AB008*D001-C060*AB014*D002+C073*AB00 SPDSTV
     $7*D002+C013*AB008*D003-C061*AB008*D004+C009*AB014*D002-C011*AB007* SPDSTV
     $D002-C051*AB008*D003+C007*AB008*D004+C006*AB007*D012-C075*AB007*D0 SPDSTV
     $13+C009*AB011*D008-C010*AB002*D008-C051*AB003*D009+C052*AB001*D009 SPDSTV
     $+C006*AB002*D017-C053*AB002*D018+C004*AB003*D009-C005*AB001*D009-C SPDSTV
     $049*AB002*D017+C002*AB002*D018+C001*AB001*D028-C054*AB001*D029)    SPDSTV
      EP47=OF8*(+C021*AB021*D001-C025*AB009*D001-C015*AB015*D002-C015*AB SPDSTV
     $014*D005+C018*AB007*D005+C013*AB008*D006+C009*AB015*D002-C007*AB00 SPDSTV
     $9*D003+C076*AB009*D004-C007*AB008*D006+C006*AB007*D014-C053*AB007* SPDSTV
     $D015+C009*AB012*D008-C011*AB004*D008-C007*AB005*D009-C007*AB003*D0 SPDSTV
     $10+C008*AB001*D010+C006*AB002*D019+C004*AB005*D009-C002*AB004*D017 SPDSTV
     $+C055*AB004*D018-C002*AB002*D019+C001*AB001*D030-C050*AB001*D031)  SPDSTV
      EP57=OF8*(+C018*AB008*D001+C008*AB007*D002+C008*AB002*D008+C003*AB SPDSTV
     $001*D009+C021*AB027*D001-C025*AB008*D001-C060*AB015*D005+C013*AB00 SPDSTV
     $8*D007-C061*AB008*D004+C009*AB024*D002-C011*AB007*D002-C051*AB009* SPDSTV
     $D006+C006*AB007*D016-C053*AB007*D013+C009*AB013*D008-C011*AB002*D0 SPDSTV
     $08-C051*AB005*D010+C006*AB002*D035-C053*AB002*D018+C004*AB006*D009 SPDSTV
     $-C005*AB001*D009-C049*AB004*D019+C001*AB001*D038-C050*AB001*D029)  SPDSTV
      EP77=OF8*(+C021*AB022*D001-C025*AB010*D001-C015*AB016*D002-C025*AB SPDSTV
     $003*D001+C026*AB001*D001+C018*AB002*D002-C015*AB014*D008+C018*AB00 SPDSTV
     $7*D008+C013*AB008*D009+C009*AB016*D002-C007*AB010*D003+C076*AB010* SPDSTV
     $D004-C011*AB002*D002+C008*AB001*D003-C008*AB001*D004-C051*AB008*D0 SPDSTV
     $09+C006*AB007*D017-C053*AB007*D018+C009*AB014*D008-C011*AB007*D008 SPDSTV
     $-C007*AB003*D011+C008*AB001*D011+C006*AB002*D020+C076*AB003*D004-C SPDSTV
     $053*AB002*D013+C004*AB008*D009-C002*AB007*D017+C055*AB007*D018-C00 SPDSTV
     $2*AB002*D020+C001*AB001*D032-C050*AB001*D033+C055*AB002*D013-C050* SPDSTV
     $AB001*D022+C074*AB001*D023)                                        SPDSTV
      EP87=OF8*(+C021*AB028*D001-C025*AB005*D001-C015*AB015*D008-C015*AB SPDSTV
     $016*D005+C018*AB002*D005+C013*AB008*D010+C009*AB025*D002-C011*AB00 SPDSTV
     $4*D002-C007*AB009*D009-C007*AB010*D006+C008*AB001*D006+C006*AB007* SPDSTV
     $D019+C009*AB015*D008-C007*AB005*D011+C076*AB005*D004-C007*AB008*D0 SPDSTV
     $10+C006*AB002*D036-C053*AB002*D015+C004*AB009*D009-C002*AB004*D020 SPDSTV
     $+C055*AB004*D013-C002*AB007*D019+C001*AB001*D039-C050*AB001*D025)  SPDSTV
      EP97=OF8*(+C018*AB008*D001+C008*AB007*D002+C008*AB002*D008+C003*AB SPDSTV
     $001*D009+C021*AB033*D001-C024*AB008*D001-C060*AB016*D008+C073*AB00 SPDSTV
     $2*D008+C013*AB008*D011-C061*AB008*D004+C009*AB032*D002-C010*AB007* SPDSTV
     $D002-C051*AB010*D009+C052*AB001*D009+C006*AB007*D020-C053*AB007*D0 SPDSTV
     $13+C009*AB016*D008-C011*AB002*D008-C051*AB008*D011+C007*AB008*D004 SPDSTV
     $+C006*AB002*D043-C075*AB002*D018+C004*AB010*D009-C005*AB001*D009-C SPDSTV
     $049*AB007*D020+C002*AB007*D013+C001*AB001*D044-C054*AB001*D029)    SPDSTV
      EP28=OF8*(+C018*AB009*D001+C008*AB004*D008+C008*AB007*D005+C003*AB SPDSTV
     $001*D010+C021*AB021*D001-C025*AB009*D001-C060*AB015*D002+C013*AB00 SPDSTV
     $9*D003-C061*AB009*D004+C009*AB012*D008-C011*AB004*D008-C051*AB005* SPDSTV
     $D009+C006*AB004*D017-C053*AB004*D018+C009*AB014*D005-C011*AB007*D0 SPDSTV
     $05-C051*AB008*D006+C006*AB007*D014-C053*AB007*D015+C004*AB003*D010 SPDSTV
     $-C005*AB001*D010-C049*AB002*D019+C001*AB001*D030-C050*AB001*D031)  SPDSTV
      EP48=OF8*(+C021*AB027*D001-C015*AB024*D002-C025*AB008*D001+C018*AB SPDSTV
     $007*D002-C015*AB015*D005+C013*AB009*D006+C009*AB013*D008-C007*AB00 SPDSTV
     $6*D009-C011*AB002*D008+C008*AB001*D009-C007*AB005*D010+C006*AB004* SPDSTV
     $D019+C009*AB015*D005-C007*AB009*D006-C007*AB008*D007+C006*AB007*D0 SPDSTV
     $16+C076*AB008*D004-C053*AB007*D013+C004*AB005*D010-C002*AB004*D019 SPDSTV
     $-C002*AB002*D035+C001*AB001*D038+C055*AB002*D018-C050*AB001*D029)  SPDSTV
      EP58=OF8*(+C018*AB009*D001+C008*AB004*D008+C008*AB007*D005+C003*AB SPDSTV
     $001*D010+C021*AB030*D001-C024*AB009*D001-C060*AB024*D005+C073*AB00 SPDSTV
     $7*D005+C013*AB009*D007-C061*AB009*D004+C009*AB023*D008-C010*AB004* SPDSTV
     $D008-C051*AB006*D010+C052*AB001*D010+C006*AB004*D035-C053*AB004*D0 SPDSTV
     $18+C009*AB024*D005-C011*AB007*D005-C051*AB009*D007+C007*AB009*D004 SPDSTV
     $+C006*AB007*D034-C075*AB007*D015+C004*AB006*D010-C005*AB001*D010-C SPDSTV
     $049*AB004*D035+C002*AB004*D018+C001*AB001*D041-C054*AB001*D031)    SPDSTV
      EP78=OF8*(+C021*AB028*D001-C015*AB025*D002-C025*AB005*D001+C018*AB SPDSTV
     $004*D002-C015*AB015*D008+C013*AB009*D009+C009*AB015*D008-C007*AB00 SPDSTV
     $9*D009-C007*AB005*D011+C006*AB004*D020+C076*AB005*D004-C053*AB004* SPDSTV
     $D013+C009*AB016*D005-C007*AB010*D006-C011*AB002*D005+C008*AB001*D0 SPDSTV
     $06-C007*AB008*D010+C006*AB007*D019+C004*AB008*D010-C002*AB007*D019 SPDSTV
     $-C002*AB002*D036+C001*AB001*D039+C055*AB002*D015-C050*AB001*D025)  SPDSTV
      EP88=OF8*(+C021*AB031*D001-C025*AB006*D001-C015*AB024*D008-C025*AB SPDSTV
     $010*D001+C026*AB001*D001+C018*AB007*D008-C015*AB025*D005+C018*AB00 SPDSTV
     $4*D005+C013*AB009*D010+C009*AB024*D008-C007*AB006*D011+C076*AB006* SPDSTV
     $D004-C011*AB007*D008+C008*AB001*D011-C008*AB001*D004-C051*AB009*D0 SPDSTV
     $10+C006*AB004*D036-C053*AB004*D015+C009*AB025*D005-C011*AB004*D005 SPDSTV
     $-C007*AB010*D007+C008*AB001*D007+C006*AB007*D035+C076*AB010*D004-C SPDSTV
     $053*AB007*D018+C004*AB009*D010-C002*AB004*D036+C055*AB004*D015-C00 SPDSTV
     $2*AB007*D035+C001*AB001*D042-C050*AB001*D027+C055*AB007*D018-C050* SPDSTV
     $AB001*D033+C074*AB001*D023)                                        SPDSTV
      EP98=OF8*(+C018*AB009*D001+C008*AB004*D008+C008*AB007*D005+C003*AB SPDSTV
     $001*D010+C021*AB034*D001-C024*AB009*D001-C060*AB025*D008+C073*AB00 SPDSTV
     $4*D008+C013*AB009*D011-C061*AB009*D004+C009*AB025*D008-C011*AB004* SPDSTV
     $D008-C051*AB009*D011+C007*AB009*D004+C006*AB004*D043-C075*AB004*D0 SPDSTV
     $18+C009*AB032*D005-C010*AB007*D005-C051*AB010*D010+C052*AB001*D010 SPDSTV
     $+C006*AB007*D036-C053*AB007*D015+C004*AB010*D010-C005*AB001*D010-C SPDSTV
     $049*AB007*D036+C002*AB007*D015+C001*AB001*D045-C054*AB001*D031)    SPDSTV
      EP29=OF8*(+C017*AB001*D001+C018*AB010*D001-C019*AB001*D001+C052*AB SPDSTV
     $007*D008+C003*AB001*D011-C058*AB001*D004+C011*AB003*D001-C020*AB00 SPDSTV
     $1*D001-C052*AB002*D002+C012*AB001*D003-C059*AB001*D004+C021*AB022* SPDSTV
     $D001-C025*AB010*D001-C060*AB016*D002+C013*AB010*D003-C061*AB010*D0 SPDSTV
     $04-C025*AB003*D001+C026*AB001*D001+C073*AB002*D002-C014*AB001*D003 SPDSTV
     $+C062*AB001*D004+C063*AB014*D008-C071*AB007*D008-C064*AB008*D009+C SPDSTV
     $056*AB007*D017-C006*AB007*D018+C004*AB003*D011-C005*AB001*D011-C04 SPDSTV
     $9*AB002*D020+C001*AB001*D032-C050*AB001*D033-C069*AB003*D004+C070* SPDSTV
     $AB001*D004+C002*AB002*D013-C050*AB001*D022+C074*AB001*D023)        SPDSTV
      EP49=OF8*(+C011*AB005*D001-C008*AB004*D002-C008*AB002*D005+C012*AB SPDSTV
     $001*D006+C021*AB028*D001-C015*AB025*D002-C015*AB016*D005+C013*AB01 SPDSTV
     $0*D006-C025*AB005*D001+C018*AB004*D002+C018*AB002*D005-C014*AB001* SPDSTV
     $D006+C063*AB015*D008-C051*AB009*D009-C051*AB008*D010+C056*AB007*D0 SPDSTV
     $19+C004*AB005*D011-C002*AB004*D020-C002*AB002*D036+C001*AB001*D039 SPDSTV
     $-C069*AB005*D004+C055*AB004*D013+C055*AB002*D015-C050*AB001*D025)  SPDSTV
      EP59=OF8*(+C017*AB001*D001+C018*AB010*D001-C019*AB001*D001+C052*AB SPDSTV
     $007*D008+C003*AB001*D011-C058*AB001*D004+C011*AB006*D001-C020*AB00 SPDSTV
     $1*D001-C052*AB004*D005+C012*AB001*D007-C059*AB001*D004+C021*AB031* SPDSTV
     $D001-C025*AB010*D001-C060*AB025*D005+C013*AB010*D007-C061*AB010*D0 SPDSTV
     $04-C025*AB006*D001+C026*AB001*D001+C073*AB004*D005-C014*AB001*D007 SPDSTV
     $+C062*AB001*D004+C063*AB024*D008-C071*AB007*D008-C064*AB009*D010+C SPDSTV
     $056*AB007*D035-C006*AB007*D018+C004*AB006*D011-C005*AB001*D011-C04 SPDSTV
     $9*AB004*D036+C001*AB001*D042-C050*AB001*D033-C069*AB006*D004+C070* SPDSTV
     $AB001*D004+C002*AB004*D015-C050*AB001*D027+C074*AB001*D023)        SPDSTV
      EP79=OF8*(+C011*AB008*D001-C008*AB007*D002-C008*AB002*D008+C012*AB SPDSTV
     $001*D009+C021*AB033*D001-C015*AB032*D002-C024*AB008*D001+C016*AB00 SPDSTV
     $7*D002-C015*AB016*D008+C013*AB010*D009+C018*AB002*D008-C014*AB001* SPDSTV
     $D009+C063*AB016*D008-C051*AB010*D009-C071*AB002*D008+C052*AB001*D0 SPDSTV
     $09-C051*AB008*D011+C056*AB007*D020+C007*AB008*D004-C006*AB007*D013 SPDSTV
     $+C004*AB008*D011-C002*AB007*D020-C002*AB002*D043+C001*AB001*D044+C SPDSTV
     $072*AB002*D018-C054*AB001*D029-C069*AB008*D004+C055*AB007*D013)    SPDSTV
      EP89=OF8*(+C011*AB009*D001-C008*AB004*D008-C008*AB007*D005+C012*AB SPDSTV
     $001*D010+C021*AB034*D001-C024*AB009*D001-C015*AB025*D008-C015*AB03 SPDSTV
     $2*D005+C016*AB007*D005+C013*AB010*D010+C018*AB004*D008-C014*AB001* SPDSTV
     $D010+C063*AB025*D008-C071*AB004*D008-C051*AB009*D011+C007*AB009*D0 SPDSTV
     $04-C051*AB010*D010+C052*AB001*D010+C056*AB007*D036-C006*AB007*D015 SPDSTV
     $+C004*AB009*D011-C002*AB004*D043+C072*AB004*D018-C002*AB007*D036+C SPDSTV
     $001*AB001*D045-C054*AB001*D031-C069*AB009*D004+C055*AB007*D015)    SPDSTV
      EP99=OF8*(+C017*AB001*D001+C018*AB010*D001-C019*AB001*D001+C057*AB SPDSTV
     $007*D008+C003*AB001*D011-C058*AB001*D004+C011*AB010*D001-C020*AB00 SPDSTV
     $1*D001+C012*AB001*D011-C059*AB001*D004+C021*AB035*D001-C022*AB010* SPDSTV
     $D001-C060*AB032*D008+C023*AB001*D001+C038*AB007*D008+C013*AB010*D0 SPDSTV
     $11-C061*AB010*D004-C014*AB001*D011+C062*AB001*D004+C063*AB032*D008 SPDSTV
     $-C033*AB007*D008-C064*AB010*D011+C051*AB010*D004+C031*AB001*D011-C SPDSTV
     $052*AB001*D004+C056*AB007*D043-C065*AB007*D018+C004*AB010*D011-C00 SPDSTV
     $5*AB001*D011-C049*AB007*D043+C066*AB007*D018+C001*AB001*D046-C067* SPDSTV
     $AB001*D033+C068*AB001*D023-C069*AB010*D004+C070*AB001*D004)        SPDSTV
C     ****************************************************************** SPDSTV
C     *                                PD                              * SPDSTV
C     ****************************************************************** SPDSTV
 3242 CONTINUE                                                           SPDSTV
      EP12=OF7*(+C008*AB002*D001-C012*AB001*D002+C015*AB011*D001-C016*AB SPDSTV
     $002*D001-C013*AB003*D002+C014*AB001*D002+C051*AB003*D002-C052*AB00 SPDSTV
     $1*D002-C056*AB002*D003+C006*AB002*D004+C002*AB002*D003-C001*AB001* SPDSTV
     $D012+C054*AB001*D013-C055*AB002*D004)                              SPDSTV
      EP32=OF7*(+C008*AB004*D001-C012*AB001*D005+C015*AB012*D001-C013*AB SPDSTV
     $003*D005-C018*AB004*D001+C014*AB001*D005+C051*AB005*D002-C056*AB00 SPDSTV
     $2*D006+C002*AB004*D003-C001*AB001*D014-C055*AB004*D004+C050*AB001* SPDSTV
     $D015)                                                              SPDSTV
      EP62=OF7*(+C008*AB007*D001-C012*AB001*D008+C015*AB014*D001-C013*AB SPDSTV
     $003*D008-C018*AB007*D001+C014*AB001*D008+C051*AB008*D002-C056*AB00 SPDSTV
     $2*D009+C002*AB007*D003-C001*AB001*D017-C055*AB007*D004+C050*AB001* SPDSTV
     $D018)                                                              SPDSTV
      EP14=OF7*(+C015*AB012*D001-C018*AB004*D001-C013*AB005*D002+C007*AB SPDSTV
     $005*D002-C006*AB004*D003+C053*AB004*D004+C007*AB003*D005-C008*AB00 SPDSTV
     $1*D005-C006*AB002*D006+C002*AB002*D006-C001*AB001*D014+C050*AB001* SPDSTV
     $D015)                                                              SPDSTV
      EP34=OF7*(+C015*AB013*D001-C018*AB002*D001-C013*AB005*D005+C007*AB SPDSTV
     $006*D002-C008*AB001*D002-C006*AB004*D006+C007*AB005*D005-C006*AB00 SPDSTV
     $2*D007+C053*AB002*D004+C002*AB004*D006-C001*AB001*D016+C050*AB001* SPDSTV
     $D013)                                                              SPDSTV
      EP64=OF7*(+C015*AB015*D001-C013*AB005*D008+C007*AB009*D002-C006*AB SPDSTV
     $004*D009+C007*AB008*D005-C006*AB002*D010+C002*AB007*D006-C001*AB00 SPDSTV
     $1*D019)                                                            SPDSTV
      EP15=OF7*(+C008*AB002*D001-C012*AB001*D002+C015*AB013*D001-C013*AB SPDSTV
     $006*D002-C018*AB002*D001+C014*AB001*D002+C051*AB005*D005-C056*AB00 SPDSTV
     $4*D006+C002*AB002*D007-C001*AB001*D016-C055*AB002*D004+C050*AB001* SPDSTV
     $D013)                                                              SPDSTV
      EP35=OF7*(+C008*AB004*D001-C012*AB001*D005+C015*AB023*D001-C016*AB SPDSTV
     $004*D001-C013*AB006*D005+C014*AB001*D005+C051*AB006*D005-C052*AB00 SPDSTV
     $1*D005-C056*AB004*D007+C006*AB004*D004+C002*AB004*D007-C001*AB001* SPDSTV
     $D034+C054*AB001*D015-C055*AB004*D004)                              SPDSTV
      EP65=OF7*(+C008*AB007*D001-C012*AB001*D008+C015*AB024*D001-C013*AB SPDSTV
     $006*D008-C018*AB007*D001+C014*AB001*D008+C051*AB009*D005-C056*AB00 SPDSTV
     $4*D010+C002*AB007*D007-C001*AB001*D035-C055*AB007*D004+C050*AB001* SPDSTV
     $D018)                                                              SPDSTV
      EP17=OF7*(+C015*AB014*D001-C018*AB007*D001-C013*AB008*D002+C007*AB SPDSTV
     $008*D002-C006*AB007*D003+C053*AB007*D004+C007*AB003*D008-C008*AB00 SPDSTV
     $1*D008-C006*AB002*D009+C002*AB002*D009-C001*AB001*D017+C050*AB001* SPDSTV
     $D018)                                                              SPDSTV
      EP37=OF7*(+C015*AB015*D001-C013*AB008*D005+C007*AB009*D002-C006*AB SPDSTV
     $007*D006+C007*AB005*D008-C006*AB002*D010+C002*AB004*D009-C001*AB00 SPDSTV
     $1*D019)                                                            SPDSTV
      EP67=OF7*(+C015*AB016*D001-C018*AB002*D001-C013*AB008*D008+C007*AB SPDSTV
     $010*D002-C008*AB001*D002-C006*AB007*D009+C007*AB008*D008-C006*AB00 SPDSTV
     $2*D011+C053*AB002*D004+C002*AB007*D009-C001*AB001*D020+C050*AB001* SPDSTV
     $D013)                                                              SPDSTV
      EP18=OF7*(+C015*AB015*D001-C013*AB009*D002+C007*AB005*D008-C006*AB SPDSTV
     $004*D009+C007*AB008*D005-C006*AB007*D006+C002*AB002*D010-C001*AB00 SPDSTV
     $1*D019)                                                            SPDSTV
      EP38=OF7*(+C015*AB024*D001-C018*AB007*D001-C013*AB009*D005+C007*AB SPDSTV
     $006*D008-C008*AB001*D008-C006*AB004*D010+C007*AB009*D005-C006*AB00 SPDSTV
     $7*D007+C053*AB007*D004+C002*AB004*D010-C001*AB001*D035+C050*AB001* SPDSTV
     $D018)                                                              SPDSTV
      EP68=OF7*(+C015*AB025*D001-C018*AB004*D001-C013*AB009*D008+C007*AB SPDSTV
     $009*D008-C006*AB004*D011+C053*AB004*D004+C007*AB010*D005-C008*AB00 SPDSTV
     $1*D005-C006*AB007*D010+C002*AB007*D010-C001*AB001*D036+C050*AB001* SPDSTV
     $D015)                                                              SPDSTV
      EP19=OF7*(+C008*AB002*D001-C012*AB001*D002+C015*AB016*D001-C013*AB SPDSTV
     $010*D002-C018*AB002*D001+C014*AB001*D002+C051*AB008*D008-C056*AB00 SPDSTV
     $7*D009+C002*AB002*D011-C001*AB001*D020-C055*AB002*D004+C050*AB001* SPDSTV
     $D013)                                                              SPDSTV
      EP39=OF7*(+C008*AB004*D001-C012*AB001*D005+C015*AB025*D001-C013*AB SPDSTV
     $010*D005-C018*AB004*D001+C014*AB001*D005+C051*AB009*D008-C056*AB00 SPDSTV
     $7*D010+C002*AB004*D011-C001*AB001*D036-C055*AB004*D004+C050*AB001* SPDSTV
     $D015)                                                              SPDSTV
      EP69=OF7*(+C008*AB007*D001-C012*AB001*D008+C015*AB032*D001-C016*AB SPDSTV
     $007*D001-C013*AB010*D008+C014*AB001*D008+C051*AB010*D008-C052*AB00 SPDSTV
     $1*D008-C056*AB007*D011+C006*AB007*D004+C002*AB007*D011-C001*AB001* SPDSTV
     $D043+C054*AB001*D018-C055*AB007*D004)                              SPDSTV
C     ****************************************************************** SPDSTV
C     *                                SD                              * SPDSTV
C     ****************************************************************** SPDSTV
 3260 CONTINUE                                                           SPDSTV
      EP02=OF6*(+C012*AB001*D001+C013*AB003*D001-C014*AB001*D001+C056*AB SPDSTV
     $002*D002+C001*AB001*D003-C050*AB001*D004)                          SPDSTV
      EP04=OF6*(+C013*AB005*D001+C006*AB004*D002+C006*AB002*D005+C001*AB SPDSTV
     $001*D006)                                                          SPDSTV
      EP05=OF6*(+C012*AB001*D001+C013*AB006*D001-C014*AB001*D001+C056*AB SPDSTV
     $004*D005+C001*AB001*D007-C050*AB001*D004)                          SPDSTV
      EP07=OF6*(+C013*AB008*D001+C006*AB007*D002+C006*AB002*D008+C001*AB SPDSTV
     $001*D009)                                                          SPDSTV
      EP08=OF6*(+C013*AB009*D001+C006*AB004*D008+C006*AB007*D005+C001*AB SPDSTV
     $001*D010)                                                          SPDSTV
      EP09=OF6*(+C012*AB001*D001+C013*AB010*D001-C014*AB001*D001+C056*AB SPDSTV
     $007*D008+C001*AB001*D011-C050*AB001*D004)                          SPDSTV
      IF(ITYPE-6)3261,3262,3261                                          SPDSTV
C     ****************************************************************** SPDSTV
C     *                                PP                              * SPDSTV
C     ****************************************************************** SPDSTV
 3261 CONTINUE                                                           SPDSTV
      EP10=OF1*(+C002*AB002*D001-C001*AB001*D002)                        SPDSTV
      EP30=OF1*(+C002*AB004*D001-C001*AB001*D005)                        SPDSTV
      EP60=OF1*(+C002*AB007*D001-C001*AB001*D008)                        SPDSTV
      EP11=OF4*(-C007*AB003*D001+C008*AB001*D001+C006*AB002*D002-C002*AB SPDSTV
     $002*D002+C001*AB001*D003-C050*AB001*D004)                          SPDSTV
      EP31=OF4*(-C007*AB005*D001+C006*AB002*D005-C002*AB004*D002+C001*AB SPDSTV
     $001*D006)                                                          SPDSTV
      EP61=OF4*(-C007*AB008*D001+C006*AB002*D008-C002*AB007*D002+C001*AB SPDSTV
     $001*D009)                                                          SPDSTV
      EP13=OF4*(-C007*AB005*D001+C006*AB004*D002-C002*AB002*D005+C001*AB SPDSTV
     $001*D006)                                                          SPDSTV
      EP33=OF4*(-C007*AB006*D001+C008*AB001*D001+C006*AB004*D005-C002*AB SPDSTV
     $004*D005+C001*AB001*D007-C050*AB001*D004)                          SPDSTV
      EP63=OF4*(-C007*AB009*D001+C006*AB004*D008-C002*AB007*D005+C001*AB SPDSTV
     $001*D010)                                                          SPDSTV
      EP16=OF4*(-C007*AB008*D001+C006*AB007*D002-C002*AB002*D008+C001*AB SPDSTV
     $001*D009)                                                          SPDSTV
      EP36=OF4*(-C007*AB009*D001+C006*AB007*D005-C002*AB004*D008+C001*AB SPDSTV
     $001*D010)                                                          SPDSTV
      EP66=OF4*(-C007*AB010*D001+C008*AB001*D001+C006*AB007*D008-C002*AB SPDSTV
     $007*D008+C001*AB001*D011-C050*AB001*D004)                          SPDSTV
 3262 CONTINUE
      DO 2137 I=1,100                                                    SPDSTV
 2137 EPN(I)=EPN(I)+EEP(I)                                               SPDSTV
  105 CONTINUE                                                           SPDSTV
C     ****************************************************************** SPDSTV
C     END OF LOOP OVER GAUSSIANS                                         SPDSTV
C     STORE IN ARRAYS                                                    SPDSTV
C     ****************************************************************** SPDSTV
      INTC=0                                                             SPDSTV
      DO 500 J=1,10                                                      SPDSTV
      R3B=RENORM(J)                                                      SPDSTV
      DO 500 I=1,10                                                      SPDSTV
      R3A=R3B*RENORM(I)                                                  SPDSTV
      INTC=INTC+1                                                        SPDSTV
  500 EPN(INTC) = ( EPN(INTC) )*R3A*SYMFAC
      CALL REDUC1(EPN,LAMAX,LBMAX,I6TO5)
      CALL MATFIL(H,EPN,AOS,SHELLT,INEW,JNEW,LAMAX,LBMAX,LA,LB)
 1000 CONTINUE                                                           SPDSTV
C
C        REFORMAT COMMON /B/ AND THE H ARRAY IF THIS BASIS CONTAINS
C        P ONLY SHELLS
C
      IF (IPO(4) .EQ. 0) GOTO 1285
      WRITE(IOUT,*) 'DEBUG OF UNSTAR'
      CALL LINOUT (H,NBASIS,0,0)
 1285 CONTINUE
C
      CALL UNSTAR (NBASIS,SHELLT,SHELLC,AOS,NSHELL,H,NOSTAR)
C
      IF(IPO(4).EQ.0) GOTO 1500
      WRITE(IOUT,2010)
      CALL LINOUT(H,NBASIS,0,0)
 1500 CONTINUE
      RETURN
      END                                                                SPDSTV
      SUBROUTINE INV(A,N,IS,IAD1,IAD2,D,MDM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ******************************************************************
C     INVERSION OF SQUARE MATRIX A BY MEANS OF THE GAUSS-JORDAN
C     ALGORITHM
C
C     APRIL 72/RS9B
C     ******************************************************************
      DIMENSION A(MDM,MDM),IS(2,MDM),IAD1(MDM),IAD2(MDM),D(MDM)
C
      COMMON/IO/IN,IOUT,IPUNCH
C
      DATA ZERO/0.0D0/, ONE/1.0D0/, SMALL/1.0D-20/
C
 2000 FORMAT(' WARNING FROM INV: MATRIX IS SINGULAR')
C     ******************************************************************
      DO 1 L=1,N
      IS(1,L)=0
  1   IS(2,L)=0
      DO 9 IMA=1,N
      B= ZERO
      DO 2 L=1,N
      DO 2 M=1,N
      IF(IS(1,L).EQ.1.OR.IS(2,M).EQ.1) GOTO 2
      E=DABS(A(L,M))
      IF(E.LT.B) GOTO 8
      I=L
      K=M
    8 B=DMAX1(B,E)
  2   CONTINUE
      IS(1,I)=1
      IS(2,K)=1
      IAD1(K)=I
      IAD2(I)=K
      B=A(I,K)
C.....PIVOT
      IF(DABS(B).LT. SMALL) GOTO 20
      A(I,K)=ONE/B
      DO 6 L=1,N
      IF(L.EQ.K) GOTO 6
C.....KELLERZEILE
      A(I,L)=-A(I,L)/B
  6   CONTINUE
      DO 5 L=1,N
      DO 5 M=1,N
      IF(L.EQ.I.OR.M.EQ.K) GOTO 5
C.....RECHTECK-REGEL
      A(L,M)=A(L,M)+A(L,K)*A(I,M)
  5   CONTINUE
      DO 11 L=1,N
      IF(L.EQ.I) GOTO 11
C.....PIVOT-SPALTE
      A(L,K)=A(L,K)/B
  11  CONTINUE
  9   CONTINUE
C.....PERMUTATION DER ZEILEN, UM DIE NATUERLICHE ORDNUNG WIEDER HERZUSTE
      DO 15 L=1,N
      DO 13 J=1,N
      K=IAD1(J)
  13  D(J)=A(K,L)
      DO 14 J=1,N
  14  A(J,L)=D(J)
  15  CONTINUE
C.....PERMUTATION DER SPALTEN
      DO 16 L=1,N
      DO 17 J=1,N
      K=IAD2(J)
  17  D(J)=A(L,K)
      DO 18 J=1,N
  18  A(L,J)=D(J)
  16  CONTINUE
      RETURN
C
C     ERROR EXIT: MATRIX IS SINGULAR
  20  WRITE(IOUT,2000)
      STOP 'INV IN POLAR'
      END
      SUBROUTINE LINOUT(X,N,KEY,IZERO)
C
C     GENERAL LINEAR MATRIX OUTPUT ROUTINE
C
C     KEY=0  MATRIX SYMMETRIC
C     KEY=1  MATRIX SQUARE ASYMMETRIC
C
C     IZERO=0  ZERO MATRIX ELEMENTS LESS THAN CUTOFF
C     IZERO=1  DO NOT ZERO MATRIX ELEMENTS
C
C     CUTOFF=1.0E-06
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/IO/IN,IOUT
C
      DIMENSION S(9),X(1)
C
      DATA CUTOFF/1.0E-06/
      DATA ZERO/0.0E0/
C
      IA(I)=(I*(I-1))/2
C
C
      ILOWER=1
  100 IUPPER=MIN0(ILOWER+8,N)
      IRANGE=MIN0(IUPPER-ILOWER+1,9)
      WRITE (IOUT,9000) (J,J=ILOWER,IUPPER)
      WRITE (IOUT,9010)
      DO 160 I=1,N
      K=1
      DO 150 J=ILOWER,IUPPER
      IF(KEY)110,120,110
  110 IJ=N*(J-1)+I
      GO TO 140
  120 IJ=IA(I)+J
      IF(I-J)130,140,140
  130 IJ=IA(J)+I
  140 S(K)=X(IJ)
      IF(IZERO.EQ.0.AND.ABS(S(K)).LE.CUTOFF) S(K)=ZERO
  150 K=K+1
  160 WRITE (IOUT,9020) I,(S(J),J=1,IRANGE)
      WRITE (IOUT,9010)
      ILOWER=ILOWER+9
      IF(N-IUPPER)170,170,100
  170 RETURN
 9000 FORMAT(12X,8(I3,11X),I3)
 9010 FORMAT(/)
 9020 FORMAT(1X,I3,2X,9E14.6)
      END
      SUBROUTINE MATFIL(A,AA,AOS,SHELLT,INEW,JNEW,LAMAX,LBMAX,LA,LB)
C
C     GAUSSIAN 77/UCI
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER AOS(1), SHELLT(1)
C
      DIMENSION A(1),AA(1)
C
      LIND(I)=(I*(I-1))/2
C
      ISTART=AOS(INEW)                                                   MATFIL
      JSTART=AOS(JNEW)                                                   MATFIL
      IAL = 0
      IAU = 5
      IBL = 0
      IBU = 5
      IMA = 0
      IMB = 0
      IF(SHELLT(INEW) .EQ. 2) IMA = 1
      IF(SHELLT(JNEW) .EQ. 2) IMB = 1
C                                                                        MATFIL
  120 INTC=0                                                             MATFIL
      DO 170 J=1,LBMAX                                                   MATFIL
      DO 170 I=1,LAMAX                                                   MATFIL
      INTC=INTC+1                                                        MATFIL
      IF( LA.GT.1 .AND. I.GT.IAL .AND. I.LT.IAU )GO TO 170               MATFIL
      IF( LA.EQ.1 .AND. I.EQ.IMA                )GO TO 170               MATFIL
      IF( LB.GT.1 .AND. J.GT.IBL .AND. J.LT.IBU )GO TO 170               MATFIL
      IF( LB.EQ.1 .AND. J.EQ.IMB                )GO TO 170               MATFIL
      IND=ISTART+I-1                                                     MATFIL
      JND=JSTART+J-1                                                     MATFIL
      IF(IND-JND)130,140,150                                             MATFIL
  130 IJ=LIND(JND)+IND                                                   MATFIL
      GO TO 160                                                          MATFIL
  140 IJ=LIND(IND+1)                                                     MATFIL
      GO TO 160                                                          MATFIL
  150 IJ=LIND(IND)+JND                                                   MATFIL
  160 A(IJ)=AA(INTC)                                                     MATFIL
  170 CONTINUE                                                           MATFIL
      RETURN                                                             MATFIL
      END                                                                MATFIL
      SUBROUTINE MULTAY(A,Y,X,N,MAXDIM)
C
C        MATRIX MULTIPLICATION ROUTINE
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(MAXDIM,MAXDIM),Y(MAXDIM),X(MAXDIM)
C
      DATA ZERO/0.0/
C
      DO 200 IROW=1,N
      SUM = ZERO
      DO 100 JCOL=1,N
      SUM = SUM + A(IROW,JCOL) * Y(JCOL)
  100 CONTINUE
      X(IROW) = SUM
  200 CONTINUE
      RETURN
      END
C
      SUBROUTINE OUTPUT
C
C
C        L.E. CHIRLIAN
C        APRIL 1985
C
C        A SUBROUTINE TO OUTPUT THE CHARGES AND OTHER PERTINANT
C        INFORMATION FROM THE CHELP PROGRAM
C
C        Slightly Modified for CHELPG operations by Curt Breneman
C        Yale University Department of Chemistry, 2/88
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NPOINTS = 50000)
      INTEGER*4 SHELLA,SHELLN,SHELLT,SHELLC,AOS,AON,SHLADF
      CHARACTER*40 CHKFIL
C
      COMMON /IO/ IN,IOUT
      COMMON /IPO/ IPO(5)
C+++
      COMMON /MOL/    NATOMS,ICHARG,MULTIP,NAE,NBE,NE,NBASIS,
     $                IAN(401),ATMCHG(400),C(3,400)
C+++
C      COMMON /MOL/ NATOMS,ICHARG,MULTIP,NAE,NBE,NE,NBASIS,IAN(101),
C     1             ATMCHG(100),C(3,100)
      COMMON /OUT/ NTITLE(20,3),CHKFIL,Q(400),ND,RMS,PD,
     1             NLIN,NEND(3)
      COMMON /POINTS/ P(3,NPOINTS), NP
      DATA DEB/0.393427328/
c
c
c       write(6,*) 'Debug --', nlin,nend(1),nend(2),nend(3),nwords
C
C        CALCULATE THE DIPOLE MOMENT FROM THE FITTED CHARGES
C
C
      DIPX=0.
      DIPY=0.
      DIPZ=0.
      DO 99 I=1,NATOMS
      DIPX=DIPX+(Q(I)*C(1,I))
      DIPY=DIPY+(Q(I)*C(2,I))
      DIPZ=DIPZ+(Q(I)*C(3,I))
99    CONTINUE
      DIPX=DIPX/DEB
      DIPY=DIPY/DEB
      DIPZ=DIPZ/DEB
C
C     CALCULATE TOTAL DIPOLE MOMENT
C
      DIPTOT = DSQRT(DIPX**2+DIPY**2+DIPZ**2)
C
C        CREATE OUTPUT
C
      WRITE (IOUT,100)
100   FORMAT (/,/,17X,'CHARGES FROM ELECTROSTATIC POTENTIAL GRID')
      WRITE (IOUT,110)
      WRITE (IOUT,111)
110   FORMAT (/,36X,'CHELPGrid',/)
111   FORMAT (/,15X,'Grid Modification.')
      DO 24 I=1,NLIN
      WRITE (6,1200)(NTITLE(J,I),J=1,NEND(I))
1200  FORMAT(2X,19A4)
24    CONTINUE
C
C        WRITE CHECKPOINT FILE NAME
C
      WRITE (IOUT,150)CHKFIL
150   FORMAT(/2X,'CHECKPOINT FILE:  ',A40)
C
C        PRINT DATE
C
C***    Take out for Trace-7
C      CALL FOR$JDATE(IMONTH,IDATE,IYEAR)
C***
      WRITE (IOUT,170)IMONTH,IDATE,IYEAR
170   FORMAT (/2X,I2,'-'I2,'-'I2)
C
C**************************************************************************
C                WRITE GEOMETRY
C**************************************************************************
C
      WRITE (IOUT,180)
180   FORMAT (/2X,36X,'MOLECULAR GEOMETRY')
      WRITE (IOUT,190)
190   FORMAT (/,/,17X,'ATOMIC NUMBER',8X,'X',12X,'Y',12X,'Z')
      DO 30 I=1,NATOMS
      WRITE (IOUT,200)IAN(I),C(1,I),C(2,I),C(3,I)
200   FORMAT  (/,23X,I2,8X,F10.7,3X,F10.7,3X,F10.7)
30    CONTINUE
      WRITE (IOUT,210)ICHARG
210   FORMAT (/2X,'THE TOTAL CHARGE IS CONSTRAINED TO:  ',I3)
      WRITE  (IOUT,240)
240   FORMAT (/,36X,'NET CHARGES')
      WRITE (IOUT,250)
250   FORMAT (/,28X,'ATOMIC NUMBER',5X,'CHARGE')
      WRITE (IOUT,260)(IAN(I),Q(I),I=1,NATOMS)
      WRITE (6,101) DIPTOT
101   FORMAT (/,2X,'THE DIPOLE MOMENT OF THESE CHARGES IS:  ',F8.5)
260   FORMAT (/,34X,I2,10X,F8.4)
      WRITE (IOUT,270)NP
270   FORMAT(/,2X,'FIT TO ELECTROSTATIC POTENTIAL AT ',I6,' POINTS')
      WRITE (IOUT,280)RMS
280   FORMAT (/,2X,'ROOT MEAN SQUARE DEVIATION IS ',F6.4,' KCAL/MOLE')
      RETURN
      END

      SUBROUTINE READIN
C
C       WRITTEN BY M.M. FRANCL FOR THE
C        PRINCETON CHEMISTRY DEPARTMENT VAX 11/780 UNDER VMS 3.4.
C        MODIFIED BY L.E. CHIRLIAN UNDER VMS 3.7.
C
C        MODIFIED FOR GAUSSIAN 86 BY CURT BRENEMAN (VMS 4.5)
C        Modified for G88/90 by Curt Breneman, 2/89
C        YALE UNIVERSITY DEPARTMENT OF CHEMISTRY.
C
C        THIS VERSION IS COMPATIBLE WITH GAUSSIAN 82 FROM CARNEGIE-
C        MELLON UNIVERSITY AND IS DESIGNED FOR THE INPUT OF MO AND
C        BASIS INFORMATION FROM CHECKPOINT FILES.  THIS VERSION IS
C        TO BE USED FOR THE DETERMINATION OF ATOMIC CHARGES FROM
C        ELECTROSTATIC POTENTIALS DETERMINED BY FIRST ORDER HARTREE
C        FOCK PERURBATION THEORY.
C
C        OLD LIMITATIONS:   NO MORE THAN 256 BASIS FUNCTIONS
C                      NO MORE THAN 80 SHELLS
C
C        NEW LIMITATIONS:   NO MORE THAN 1280 BASIS FUNCTIONS
C                      NO MORE THAN 400 SHELLS
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*4 SHELLA,SHELLN,SHELLT,SHELLC,AOS,AON,SHLADF,FILNUM
      CHARACTER*40 CHKFIL
      PARAMETER (NUMPTS = 20)
      DIMENSION LINE(20)
C
      COMMON /IO/ IN,IOUT
c+++
c       Change for G86 : New Commons /MOL/ and /B/
c
      COMMON /MOL/    NATOMS,ICHARG,MULTIP,NAE,NBE,NE,NBASIS,
     $                IAN(401),ATMCHG(400),C(3,400)
C
C===  Gaussian88 Modification.  New Common /b/ size.
      Common/B/EXX(6000),C1(6000),C2(6000),C3(2000),CF(2000),
     $SHLADF(4000),X(2000),Y(2000),
     $Z(2000),JAN(2000),ShellA(2000),ShellN(2000),ShellT(2000),
     $ShellC(2000),AOS(2000),AON(2000),NShell,MaxTyp
C===    Old G86 Common /b/
c      COMMON/B/EXX(1200),C1(1200),C2(1200),C3(400),CF(400),SHLADF(800),
c     $         X(400),Y(400),Z(400),JAN(400),SHELLA(400),SHELLN(400),
c     $         SHELLT(400),SHELLC(400),AOS(400),AON(400),NSHELL,MAXTYP
c
c+++
C%%%
c       Original CHELP common /B/
c
c      COMMON /B/ EXX(240),C1(240),C2(240),C3(80),CF(80),SHLADF(160),
c     $           X(80),Y(80),Z(80),
c     $           JAN(80),SHELLA(80),SHELLN(80),SHELLT(80),SHELLC(80),
c     $           AOS(80),AON(80),NSHELL,MAXTYP
C      COMMON /MOL/ NATOMS,ICHARG,MULTIP,NAE,NBE,NE,NBASIS,IAN(101),
C     $             ATMCHG(100),C(3,100)
C%%%
      COMMON /IPO/ IPO(5)
      COMMON /OUT/ NTITLE(20,3),CHKFIL,Q(400),ND,RMS,PD,
     1             NLIN,NEND(3)
      COMMON /CHARGE/ COEF_ALPHA(100000),COEF_BETA(100000),IUHF
      COMMON /SPHERE/ VDWR(400), NPTS
C        VECT(3,NUMPTS)
C
      DATA MAXBAS /1280/
      DATA MAXPTS/ 50000/
      DATA NPO/ 5/
      DATA IUNIT/ 3/, IREAD/ 2/, IFIND/11/, IBLNK/4H    /
C      DATA VECT/0.00000000,     0.00000000,     1.00000000,
C     $  0.23807485,    -0.08801514,    -0.96725059,
C     $ -0.46055608,     0.17026540,     0.87114740,
C     $  0.65287142,    -0.24136347,    -0.71798508,
C     $ -0.80242446,     0.29665251,     0.51779559,
C     $  0.00000000,     0.00000000,     1.00000000,
C     $ -0.20401674,    -0.15100817,    -0.96725059,
C     $  0.39467063,     0.29212549,     0.87114740,
C     $ -0.55947405,    -0.41410893,    -0.71798508,
C     $  0.68763258,     0.50896872,     0.51779559,
C     $ -0.94978689,    -0.28553486,    -0.12796369,
C     $  0.88757697,     0.26683266,     0.37550960,
C     $ -0.76723180,    -0.23065324,    -0.59846007,
C     $  0.59663385,     0.17936630,     0.78221211,
C     $ -0.38695708,    -0.11633108,    -0.91473018,
C     $  0.28119199,     0.95108168,    -0.12796369,
C     $ -0.26277424,    -0.88878695,     0.37550960,
C     $  0.22714509,     0.76827772,    -0.59846007,
C     $ -0.17663821,    -0.59744720,     0.78221211,
C     $  0.11456173,     0.38748460,    -0.91473018/
C
C       Old "spherical" unit vectors (14 of them)
C
c        0.5773502691896258,0.5773502691896258,
c     $  0.5773502691896258,
c     1 -0.5773502691896258,-0.5773502691896258,0.5773502691896258,
c     2  0.5773502691896258,-0.5773502691896258,-0.5773502691896258,
c     3 -0.5773502691896258,0.5773502691896258,-0.5773502691896258,
c     4  0.0000000000000000E+00,-1.000000000000000,
c     $  0.0000000000000000E+00,
c     5  0.0000000000000000E+00,0.0000000000000000E+00,
c     $  -1.000000000000000,
c     6 -1.000000000000000,0.0000000000000000E+00,
c     $  0.0000000000000000E+00,
c     7 -0.5773502691896258,-0.5773502691896258,-0.5773502691896258,
c     8  1.000000000000000,0.0000000000000000E+00,
c     $  0.0000000000000000E+00,
c     9  0.0000000000000000E+00,1.000000000000000,
c     $  0.0000000000000000E+00,
c     $  0.5773502691896258,0.5773502691896258,-0.5773502691896258,
c     $  0.0000000000000000E+00,0.0000000000000000E+00,
c     $  1.000000000000000,
c     $ -0.5773502691896258,0.5773502691896258,0.5773502691896258,
c     $  0.5773502691896258,-0.5773502691896258,0.5773502691896258,

      IN = 5
      IOUT = 6
C
C        CHECKPOINT FILE NAME
C
      READ(IN,1000) CHKFIL
 1000 FORMAT(A40)
C
C     INPUT INFORMATION
C
 1010 FORMAT(20A4)
      DO 1921 I=1,3
      NLIN= I
      READ (5,1010) LINE
      IF (LINE(1) .EQ. IBLNK) THEN
              NLIN=NLIN-1
              GOTO 192
      END IF
      DO 1922 J=1,20
      L=20-J
      IF (LINE(L) .NE. IBLNK) THEN
              NEND(I) = L
              DO 1923 K=1,NEND(I)
              NTITLE(K,I) = LINE(K)
1923    CONTINUE
      GOTO 1921
      END IF
1922  CONTINUE
1921  CONTINUE
      NLIN=NLIN+1
192   CONTINUE
C
C
C         READ IN PRINT OPTIONS
C
3000  READ(IN,*) (IPO(I),I=1,NPO)
C
C        READ IN # OF D FUNCTIONS
C        NOTE: IF THE BASIS SET USES 5 D FUNCTION, ND MUST BE
C        SET EQUAL TO 1 TO ACCOMODATE THE INTEGRAL PACKAGE.
C        IF THE BASIS SET USES  6 D FUNCTIONS, ND IS SET EQUAL TO
C        0.
C
      READ(IN,*) ND
      IF (ND .NE. 5 .AND. ND .NE. 6) THEN
         STOP '# OF D FUNCTIONS MUST BE 5 OR 6'
      END IF
      IF (ND .EQ. 5) THEN
                      ND = 1
                      GOTO 15
      END IF
      ND = 0
15    CONTINUE
C
C
C        INITIATE FILEIO
C
C***   Different Fopen statement for Trace-7
       CALL FOPEN (IUNIT,5,CHKFIL,IALLOC,junk)
c       CALL FOPEN (IUNIT,'old',CHKFIL(1:linend(chkfil))//char(0))
C
      IWWRIT = IPO(1)
C***********************************************************************
C
C        READ IN COMMON /MOL/
C
c      NWORDS = 1804
c      IFILENO = 30997
c      CALL FILEIO (IREAD,IFILENO,NWORDS,NATOMS,IALLOC)
      IRwMol=997
      MaxAtm=400
      LenMol = 4*MaxAtm + InToWP(8+MaxAtm)
      Call FileIO(2,-FilNum(IRwMol,IUnit),LenMol,NAtoms,0)
C
C***********************************************************************
C
C        READ IN SPHERE DATA (VAN DER WAALS RADII, # OF POINTS TO
C        FIT
C
C
      READ (IN,*)(VDWR(I),I=1,NATOMS)
      READ (IN,*)NPTS
      IF (NPTS .GT. MAXPTS) THEN
         STOP 'MAXIMUM NUMBER OF POINTS MUST BE LESS THAN 50000'
      END IF
C
C        SET MAXIMUM VAN DER WAALS RADII TO 4
C
      VMAX=4.
      DO 20 I=1,NATOMS
      IF (VDWR(I) .GT. VMAX) THEN
            WRITE (IOUT, 2500) I
2500                  FORMAT (3X, 'THE VAN DER WAALS RADII OF ATOM', I3,
     1            'IS OUT OF RANGE')
                  STOP
      END IF
20    CONTINUE
C
      READ (IN,*) VFACT
      DO 21 I=1,NATOMS
      VDWR(I)=VDWR(I)*VFACT
21    CONTINUE
C***********************************************************************
C
C        READ IN BASIS SET INFORMATION (COMMON /B/)
C
      IFILENO = 30506
      CALL FILEIO (IFIND,IFILENO,NWORDS,EXX,0)
      CALL FILEIO (IREAD,-IFILENO,NWORDS,EXX,0)
C***********************************************************************
      IF(IWWRIT .NE. 1) GOTO 170
      WRITE(IOUT,8000)(C(1,I),C(2,I),C(3,I),I=1,NATOMS)
 8000 FORMAT(/1X,'COORDINATES'/(1X,3F12.6))
      WRITE(IOUT,8020) NATOMS,ICHARG,MULTIP,NAE,NBE,NE,NBASIS
     $ ,NSHELL,MAXTYP
 8020 FORMAT(/1X,'NATOMS   = ',I3
     $/1X,'ICHARG = ',I3
     $/1X,'MULTIP = ',I3
     $/1X,'NAE    = ',I3
     $/1X,'NBE    = ',I3
     $/1X,'NE     = ',I3
     $/1X,'NBASIS = ',I3
     $/1X,'NSHELL = ',I3
     $/1X,'MAXTYP = ',I3)
      WRITE(IOUT,8030) (IAN(I),I=1,NATOMS)
 8030 FORMAT(/1X,'IAN'/(1X,20I3))
      WRITE(IOUT,8050) (JAN(I),SHELLT(I),SHELLA(I),I=1,NSHELL)
 8050 FORMAT(/1X,'CENTER TYPE SHELLA'/(1X,3I7))
      WRITE(IOUT,8055) SHELLA(NSHELL+1)
 8055 FORMAT(1X,14X,I7)
      WRITE(IOUT,8060) (EXX(I),C1(I),C2(I),I=1,NSHELL)
 8060 FORMAT(/1X,12X,'EXPON',8X,'EXPCOF(S)',8X,'EXPCOF(P)',
     $/(1X,3E17.9))
      WRITE(IOUT,8070) (C3(I),CF(I),I=1,NSHELL)
 8070 FORMAT(/1X,12X,'EXPCOF(D)',8X,'EXPCOF(F)',
     $/(1X,2E17.9))
      WRITE (IOUT,3575)
3575  FORMAT(/,15X,'ATOM #',5X,'V.D.W. RADII (MULTIPLIED BY FACTOR)')
      DO 3500 I=1,NATOMS
      WRITE (IOUT,4000)I,VDWR(I)
4000  FORMAT (15X,I5,20X,F5.2)
3500  CONTINUE
c      WRITE (IOUT,4500)NPTS
4500  FORMAT(/2X,'NUMBER OF POINTS TO FIT',I6)
  170 CONTINUE
C***********************************************************************
C
C        READ IN ALPHA MO COEFFICIENTS
C
      IFILENO = 30524
      CALL FILEIO (IFIND,-IFILENO,NWORDS,COEF_ALPHA,0)
      CALL FILEIO (IREAD,-IFILENO,NWORDS,COEF_ALPHA,0)
C***********************************************************************
C
C        READ IN THE BETA MO COEFFICIENTS
C
      IFILENO = 30526
      CALL FILEIO (IFIND,IFILENO,NWORDS,COEF_BETA,0)
      IF (NWORDS.EQ.0) THEN
      IUHF = 0
      GOTO 300
      END IF
      IUHF = 1
      CALL FILEIO (IREAD,-IFILENO,NWORDS,COEF_BETA,0)
C***********************************************************************
  300 CONTINUE
C***********************************************************************
      RETURN
      END
      SUBROUTINE REDUC1(X,LAMAX,LBMAX,I6TO5)
C
C        MODIFIED FOR POLARIZATION POTENTIAL CALCULATIONS
C        M.M. FRANCL  FEBRUARY 1984
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION X(100),S(100),IND5(9),IND6(10)
C
      DATA PT5/0.5/
      DATA R3OV2/0.8660254040/
      DATA IND5/1,
     $          4,7,2,
     $          3,6,9,5,8/
      DATA IND6/1,
     $          4,7,2,
     $          6,10,3,9,5,8/
C
C     ******************************************************************
C     ROUTINE REORDERS FROM ARRANGEMENT: S,Z,ZZ,X,XZ,XX,Y,YZ,XY,YY
C     TO                                 S,X,Y,Z,XX,YY,ZZ,XY,XZ,YZ
C     OR FROM   S,Z,X,Y TO S,X,Y,Z
C     THIS ENSURES LABELING COMPATIBILITY BETWEEN THE SP AND D PACKAGES
C     AT THE SAME TIME THE INTEGRALS ARE MOVED TO THE FIRST 1,4,10,16,
C     40 OR 100 LOCATIONS, DEPENDING ON THE SHELL QUANTUM NUMBERS
C     ******************************************************************
      NWORD=LAMAX*LBMAX
      IF(NWORD-1)5,180,5
    5 INTC=0
      IF(I6TO5 .EQ. 1) GOTO 40
   10 DO 20 I=1,LBMAX
      ISB=10*(IND6(I)-1)
      DO 20 J=1,LAMAX
      ISA=ISB+IND6(J)
      INTC=INTC+1
   20 S(INTC)=X(ISA)
      GO TO 160
C     ******************************************************************
C     ROUTINE TO REDUCE SIX D FUNCTIONS TO FIVE
C     ALSO REORDERS FROM : S,Z,ZZ,X,XZ,XX,Y,YZ,YX,YY
C                     TO   S,X,Y,Z,ZZ,XX-YY,XY,XZ,YZ
C     OR FROM S,Z,X,Y TO S,X,Y,Z
C     FOR COMPATIBILITY WITH SP PACKAGE
C     ******************************************************************
   40 DO 150 I=1,LBMAX
      ISB=10*(IND5(I)-1)
C     IFB=0 FOR S,X,Y,Z,XY,XZ,YZ, IFB=1 FOR ZZ-RR, IFB=2 FOR XX-YY
      IFB = 0
      IF(I .EQ. 5) IFB = 1
      IF(I .EQ. 6) IFB = 2
   80 DO 150 J=1,LAMAX
      ISA=ISB+IND5(J)
      IFA = 0
      IF(J .EQ. 5) IFA = 1
      IF(J .EQ. 6) IFA = 2
  120 IHOP = 3*IFB + IFA + 1
      GOTO(130,122,123,124,125,126,127,128,129),IHOP
C
C     ******************************************************************
C     *                           (F,O,ZA2)                            *
C     ******************************************************************
  122 XX=ZZ1(X,ISA,3,7)
      GO TO 140
C
C     ******************************************************************
C     *                           (F,O,XA2-YA2)                        *
C     ******************************************************************
  123 XX=XY1(X,ISA,4)
      GO TO 140
C
C     ******************************************************************
C     *                           (ZB2,O,F)                            *
C     ******************************************************************
  124 XX=ZZ1(X,ISA,30,70)
      GO TO 140
C
C     ******************************************************************
C     *                           (ZB2,O,ZA2)                          *
C     ******************************************************************
  125 XX=ZZ1(X,ISA,30,70)-PT5*(ZZ1(X,ISA+3,30,70)+ZZ1(X,ISA+7,30,70))
      GO TO 140
C
C     ******************************************************************
C     *                           (ZB2,O,XA2-YA2)                      *
C     ******************************************************************
  126 XX=R3OV2*(ZZ1(X,ISA,30,70)-ZZ1(X,ISA+4,30,70))
      GO TO 140
C
C     ******************************************************************
C     *                           (XB2-YB2,O,F)                        *
C     ******************************************************************
  127 XX=XY1(X,ISA,40)
      GO TO 140
C
C     ******************************************************************
C     *                           (XB2-YB2,O,ZA2)                      *
C     ******************************************************************
  128 XX=R3OV2*(ZZ1(X,ISA,3,7)-ZZ1(X,ISA+40,3,7))
      GO TO 140
C
C     ******************************************************************
C     *                           (XB2-YB2,O,XA2-YA2)                  *
C     ******************************************************************
  129 XX=R3OV2*(XY1(X,ISA,4)-XY1(X,ISA+40,4))
      GO TO 140
C
C     ******************************************************************
C     *                           (F,O,F)                              *
C     ******************************************************************
  130 XX=X(ISA)
  140 INTC=INTC+1
  150 S(INTC)=XX
  160 DO 170 I=1,NWORD
  170 X(I)=S(I)
  180 RETURN
      END
      SUBROUTINE STAR(NBASIS,SHELLT,SHELLC,AOS,NSHELL,NOSTAR)
C
C        ROUTINE TO MODIFY COMMON /B/ TO THE EXPECTED FORMAT FOR INTGRL
C        FOR BASIS SETS HAVING P ONLY SHELLS, SUCH AS THE 6-31G** BASIS
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*4 SHELLC,SHELLT,AOS
C
      DIMENSION SHELLC(2000),SHELLT(2000),AOS(2000)
C
C        LOOP OVER SHELLS
C
      DO 100 I=1,NSHELL
      IF (SHELLT(I).EQ.1 .AND. SHELLC(I).EQ.1) THEN
      NBASIS = NBASIS + 1
      NOSTAR = 1
      DO 200 J=I,NSHELL
      AOS(J) = AOS(J) + 1
  200 CONTINUE
      END IF
  100 CONTINUE
      RETURN
      END
      SUBROUTINE UEP
C
C        ROUTINE TO CALCULATE THE ELECTROSTATIC POTENTIAL FROM FIRST ORDER
C        PERTURBATION THEORY
C
C        M.M. FRANCL    JULY 1985
C        MODIFIED VERSION OF A MEPHISTO ROUTINE
C        RESTRICTED TO UNRESTRICTED HARTREE-FOCK WAVEFUNCTIONS
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NPOINTS = 50000)
      INTEGER*4 SHELLA,SHELLN,SHELLT,AOS,SHELLC,AON,HANDLE
C
      COMMON /IO/ IN,IOUT
      COMMON /IPO/ IPO(5)
C+++
      COMMON /MOL/    NATOMS,ICHARG,MULTIP,NAE,NBE,NEL,NBASIS,
     $                IAN(401),ATMCHG(400),C(3,400)
C
C===    Gaussian88 Modification for enlarged common /b/.
      Common/BB/EXX(6000),C1(6000),C2(6000),C3(6000),X(2000),Y(2000),
     $Z(2000),JAN(2000),ShellA(2000),ShellN(2000),ShellT(2000),
     $ShellC(2000),AOS(2000),AON(2000),NShell,MaxTyp
C
C===    Old G86 Common /b/
c      COMMON/B/EXX(1200),C1(1200),C2(1200),C3(1200),
c     $         X(400),Y(400),Z(400),JAN(400),SHELLA(400),SHELLN(400),
c     $         SHELLT(400),SHELLC(400),AOS(400),AON(400),NSHELL,MAXTYP
C+++
C      COMMON /B/ EXX(240),C1(240),C2(240),C3(240),X(80),Y(80),Z(80),
C     $           JAN(80),SHELLA(80),SHELLN(80),SHELLT(80),SHELLC(80)
C     $          ,AOS(80),AON(80),NSHELL,MAXTYP
C      COMMON /MOL/ NATOMS,ICHARG,MULTIP,NAE,NBE,NEL,NBASIS,IAN(101),
C     $             ATMCHG(100),C(3,100)
      COMMON /POINTS/ P(3,NPOINTS),MAXPNTS
      COMMON /ELP/ ELECP(NPOINTS)
      COMMON /CHARGE/ COEF_ALPHA(100000),COEF_BETA(100000),IUHF
      COMMON /OUT/ NTITLE(20,3),CHKFIL,Q(400),I6TO5,RMS,PERCENT,
     1             NLIN,NEND(3)
C
      DIMENSION HPERT(100000),INDEX(1280)
C
      DATA IPTCHG/1.0/
      DATA ZERO/0.0/, TWO/2.0/, VNUCMAX/30.0/
C
C        INTIALIZE TIMING
C
      HANDLE = 0
C
C        SET UP THE INDEXING TABLE FOR HPERT
C
      DO 100 I=1,NBASIS
      INDEX(I) = (I-1)*I/2
  100 CONTINUE
C
C        BEGIN LOOP TO CALCULATE ELECTROSTATIC POTENTIAL
C
      DO 200 NPNT=1,MAXPNTS
      X1 = P(1,NPNT)
      X2 = P(2,NPNT)
      X3 = P(3,NPNT)
C
C     CALCULATE THE ONE-ELECTRON INTEGRALS
C
      IF (IPO(5).EQ.1) THEN
      WRITE(IOUT,3010)
 3010 FORMAT(1X,'TIME FOR INTEGRALS')
C***
C      ISTAT = LIB$INIT_TIMER(HANDLE)
C***
      END IF
C
      CALL INTGRL (HPERT,X1,X2,X3,IPTCHG,I6TO5)
C
C***
C      IF (IPO(5).EQ.1) ISTAT = LIB$SHOW_TIMER(HANDLE)
C***
C
      IF (IPO(4).EQ.1) CALL LINOUT (HPERT,NBASIS,0,0)
C
      IF (IPO(5).EQ.1) THEN
      WRITE(IOUT,3000)
 3000 FORMAT(1X,'TIME FOR TRANSFORM')
C***
C      ISTAT = LIB$INIT_TIMER(HANDLE)
C***
      END IF
C
C     FORM THE HPERT MATRIX ELEMENTS
C
C        ALPHA CODE
C
      E = ZERO
      ICOEFI = -NBASIS
C
C        SUM OVER OCCUPIED ALPHA MOS
C
      DO 220 II=1,NAE
      ICOEFI = ICOEFI + NBASIS
C
C        CALCULATE ELECTROSTATIC POTENTIAL
C
      DO 221 IP=1,NBASIS
      CPI = COEF_ALPHA(ICOEFI+IP)
      IPDEX = INDEX(IP)
C
      DO 222 IQ=1,IP
      E = E + CPI * COEF_ALPHA(ICOEFI+IQ) * HPERT(IPDEX+IQ)
  222 CONTINUE
      DO 223 IQ=IP+1,NBASIS
      E = E + CPI * COEF_ALPHA(ICOEFI+IQ) * HPERT(IP+INDEX(IQ))
  223 CONTINUE
C
  221 CONTINUE
  220 CONTINUE
C
C        BETA CODE
C
      ICOEFI = -NBASIS
C
C        SUM OVER OCCUPIED BETA MOS
C
      DO 420 II=1,NBE
      ICOEFI = ICOEFI + NBASIS
C
C        CALCULATE ELECTROSTATIC POTENTIAL
C
      DO 421 IP=1,NBASIS
      CPI = COEF_BETA(ICOEFI+IP)
      IPDEX = INDEX(IP)
C
      DO 422 IQ=1,IP
      E = E + CPI * COEF_BETA(ICOEFI+IQ) * HPERT(IPDEX+IQ)
  422 CONTINUE
      DO 423 IQ=IP+1,NBASIS
      E = E + CPI * COEF_BETA(ICOEFI+IQ) * HPERT(IP+INDEX(IQ))
      E = E + CPI * COEF_* HPERT(IP+INDEX(IQ))
  423 CONTINUE
C
  421 CONTINUE
  420 CONTINUE
C
C***
C      IF (IPO(5) .EQ. 1) ISTAT = LIB$SHOW_TIMER(HANDLE)
C***
C
C        CALCULATE NUCLEAR PART OF ELECTROSTATIC POTENTIAL
C
      VNUC = ZERO
      DO 300 IATOM=1,NATOMS
      DEL1 = C(1,IATOM) - X1
      DEL2 = C(2,IATOM) - X2
      DEL3 = C(3,IATOM) - X3
      RA = DSQRT(DEL1*DEL1 + DEL2*DEL2 + DEL3*DEL3)
      IF (RA.EQ.ZERO) THEN
      VNUC=VNUCMAX
      GOTO 310
      END IF
      VNUC = VNUC + IAN(IATOM) / RA
  300 CONTINUE
  310 CONTINUE
C
      ELECP(NPNT) = (E + VNUC * IPTCHG)
      IF (IPO(5) .EQ. 1) WRITE(IOUT,*) 'E(',NPNT,') = ',E
  200 CONTINUE
      RETURN
      END
      SUBROUTINE UNSTAR(NBASIS,SHELLT,SHELLC,AOS,NSHELL,H,NOSTAR)
C
C        ROUTINE TO REFORMAT COMMON/B/ AND THE H ARRAY FOR BASIS
C        SETS HAVING P ONLY SHELLS
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*4 SHELLC,SHELLT,AOS
C
      DIMENSION SHELLC(2000),SHELLT(2000),AOS(2000),H(1)
C
      IF (NOSTAR.EQ.0) RETURN
C
C        LOOP OVER SHELLS
C
      DO 100 I=1,NSHELL
      IF (SHELLT(I).EQ.1 .AND. SHELLC(I).EQ.1) THEN
C
C        REMOVE EXTRA ROWS AND COLUMNS
C
C        LOOP OVER ROWS
C
      IBASIS = AOS(I)
      ITEM = (IBASIS-1) * IBASIS / 2
      DO 200 J=IBASIS+1,NBASIS
C
C        LOOP OVER COLUMNS
C
      NEWITEM = (J-1) * J /2
      DO 250 K=1,IBASIS-1
      ITEM = ITEM + 1
      NEWITEM = NEWITEM + 1
      H(ITEM) = H(NEWITEM)
  250 CONTINUE
C
C        SKIP THE VALUE IN THE IBASIS TH COLUMN
C
      NEWITEM = NEWITEM + 1
C
      DO 260 K=IBASIS+1,J
      ITEM = ITEM + 1
      NEWITEM = NEWITEM + 1
      H(ITEM) = H(NEWITEM)
  260 CONTINUE
  200 CONTINUE
      NBASIS = NBASIS - 1
C
C        RESTRUCTURE AOS TO ACCOUNT FOR THE S SHELL REMOVED
C
      DO 300 IAOS = I,NSHELL
      AOS(IAOS) = AOS(IAOS) - 1
  300 CONTINUE
C
      END IF
  100 CONTINUE
      RETURN
      END
      FUNCTION XY1(X,I,IY)
C
      DIMENSION X(100)
C
      DATA HALFR3/0.8660254040/
C
      XY1=HALFR3*(X(I)-X(I+IY))
      RETURN
      END
      FUNCTION ZZ1(X,I,IX,IY)
C
      DIMENSION X(100)
C
      DATA HALF/0.5/
C
      ZZ1=X(I)-HALF*(X(I+IX)+X(I+IY))
      RETURN
      END

