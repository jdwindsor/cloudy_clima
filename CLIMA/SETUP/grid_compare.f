      SUBROUTINE GRID(P0,FAC,ZCON,Z,ZF,DZ)

C   THIS SUBROUTINE ESTABLISHES A LOG PRESSURE GRID WITH A SPACING OF
C   DZ0 AT THE BOTTOM AND FAC*DZ0 AT THE TOP.
C       PF AND ZF ARE DEFINED AT THE FLUX GRID POINTS, P AND Z AT THE
C   TEMPERATURE GRIP POINTS.  BOTH GRIDS HAVE A POINT AT THE GROUND.

      INCLUDE 'CLIMA/INCLUDE/header.inc'
      PARAMETER(NT=76)
      PARAMETER(NSOL=38, NGS=5)

      DIMENSION Z(ND),ZF(ND),DZ(ND)
      DIMENSION PML(ND)
      COMMON/SOLARBLK/AMU0,SRFALB,OMG0A(NSOL,ND-1),
     &  ASYA(NSOL,ND-1),TAUAER(NSOL),SIGERT(NSOL),FMA(NSOL),PF(ND),
     &  ALAMBDA(NSOL),CGAS(ND,NGS),FUPSOL(ND),FDNSOL(ND),
     &  NGAS(2,NSOL),WGHT(4,2,NSOL),NPR(2,NSOL),SOLINT(NSOL),
     &  TAULAM(ND-1),ASY(ND-1),OMG0(ND-1),FMT(ND-1),QEXT(NSOL,ND-1)
      COMMON/PRESSURE/P(ND),PLOG(ND)
      COMMON/EBLOK/PG,TG,PG0,IMW,RSURF,OMEGA,POCEAN,IMOIST,
     2  BETA1,BETA2,FVDRY,PDRY
      COMMON/FBLOK/TTAB(NT),PVAP(NT),DPVAP(NT),SVC(NT),DSV(NT),DSC(NT)
     2  ,RHOV(NT),DRHOV(NT),BETAM(70,75),TCP(75),PCP(70),DPDTL(70,75),
     3  DRDTL(70,75)
C

c-jdh Pressure data for comparision with Mlawer's run
      DATA PML / 0.690,1.420,2.610,4.410,6.950,10.370,14.810,20.400,
     2           27.260,35.510,45.290,56.730,69.970,85.180,102.050,
     3           122.040,143.840,167.950,194.360,222.940,253.710,
     3           286.600,321.500,358.280,396.810,436.950,478.540,
     4           521.460,565.540,610.600,656.430,702.730,749.120,
     5           795.090,839.950,882.800,922.460,957.440,985.880,
     6           1005.430,1013.250 /

      DN = ND
      DIV = FAC*DN - 1. - (FAC-1.)*(0.5*DN*(DN+1.) - 1.)/(DN-1.)
C
C   RECOMPUTE PG
      TGR = TG
      CALL SATRAT(TGR,PSAT)
      POCEAN = 270.
C   (THIS ASSUMES A FULL TERRESTRIAL OCEAN OF WATER.)
      PH2O = AMIN1(PSAT,POCEAN)
      PG = PG0 + PH2O
c     PG = PG0
      IMOIST = 1
      PDRY = 0.
      IF (PSAT.LT.POCEAN) GO TO 4
C
C   RUNAWAY GREENHOUSE (I.E. UNSATURATED LOWER ATMOSPHERE)
      IMOIST = 0
      FPH2O = POCEAN/PG
      TC = TG - 273.15
      N = TC/10. + 2
      IF (TC.GT.600.) N = (TC - 600.)/100. + 62
      PV = POCEAN
      M = PV/5. + 1
      M1 = M - 1
      IF (M.EQ.1) M1 = M

      FP = 0.
      IF (M.GT.1) FP = (PV - PCP(M))/(PCP(M) - PCP(M1))
      FT = (TC - TCP(N))/(TCP(N+1) - TCP(N))
      BETA2 = (1.+FP-FT)*BETAM(M,N) - FP*BETAM(M1,N) + FT*BETAM(M,N+1)
      FVDRY = 1./(1. + (1.-FPH2O)/BETA2/FPH2O)

   4  P0L = log(P0)
      PGL = log(PG)
      DZ0 = (PGL - P0L)/DIV

      PF(1) = P0
      ZF(1) = P0L + ZCON
      DO 1 J=2,ND
      DZ(J) = DZ0 + (ND-J)*(FAC-1.)*DZ0/FLOAT(ND-1)
      ZF(J) = ZF(J-1) + DZ(J)
   1  PF(J) = EXP(ZF(J) - ZCON)

C   FIND P AND Z
      ND1 = ND - 1
      DO 2 J=1,ND1
      Z(J) = 0.5*(ZF(J) + ZF(J+1))
      P(J) = EXP(Z(J) - ZCON)
   2  PLOG(J) = log10(P(J))
      Z(ND) = ZF(ND)
      P(ND) = PF(ND)
      PLOG(ND) = log10(P(ND))      

C 	CODE MODIFIED: modification for the purpose of exact duplication
C	of Mlawer's RRTM MLS T-P-Alt grid. 

C    4      Z(1) = log(0.067) + ZCON

c      DATA PML/ 0.067, 0.0775, 0.091, 0.107, 0.125, 0.149, 0.179, 0.216,
c     5		0.261, 0.334, 0.431, 0.589, 1.339, 1.782, 2.292, 2.964,
c     6		3.853, 5.035, 6.613, 8.786, 11.8, 15.0693, 18.791,
c     7		23.156, 28.12, 33.733, 40.535, 48.69, 57.694, 68.429,
c     8		81.2, 96.49, 114.564, 136.511, 162.91, 196.44, 218.668,
c     9		243., 269.015, 297.469, 328.507, 361.862, 398.085, 437.556,
c     1		480.526, 532.986, 589.84, 651.552, 718.70, 792.287,
c     2		891.46, 1013./


C	PF(ND) = P(ND)
C	do j = 1, (ND-1)
C	  PF(j) = (P(j) + P(j+1))/2.
C	end do

      P(ND) = PML(ND)/1000.
      PF(ND) = PML(ND)/1000.
      DO J=1,ND-1
        PF(J) = PML(J)/1000.
        P(J) = SQRT(PML(J)/1000.) * SQRT(PML(J+1)/1000.)
        Z(J) = log(P(J))
        ZF(J) = log(PF(J))
      END DO



C   RECOMPUTE DZ AS THE DISTANCE BETWEEN TEMPERATURE GRID POINTS
      DZ(1) = log(P(1)) - log(PF(1))
      DO J=2,ND
        DZ(J) = log(P(J)) - log(P(J-1))
      end do
      RETURN
      END
