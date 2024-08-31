cccccccccccccccccccccccccccc
      subroutine init_optics( do_optics, read_mie, ngas, gas_name, 
     $  nwave, wave, nrad, radius, dr, qscat, qext, cos_qscat,do_IR_wav )

c   --------------------------------------------------------------------
c
c   Setup up a particle size grid and calculate single-particle scattering
c   and absorption efficiencies and other parameters to be used by
c   calc_optics()
c
c   input scalars:
c
c     do_optics    .false. means do nothing
c     read_mie     .true. means read Mie coefficients from file
c                  'gas_name'.mieff
c     ngas         number of condensing gases
c
c   input vectors:
c
c     gas_name     names of condensing gases
c
c
c ( intermediate vectors defined by this routine:
c
c     rup          upper bounds of radius bins (cm)
c     rmin         minimum radius of grid (cm)
c )
c
c   output vectors:
c
c     wave         wavelength bin centers (cm)
c     radius       radius bin centers (cm)
c     dr           widths of radius bins (cm)
c     qscat        scattering efficiency
c     qext         scattering efficiency
c     cos_qscat    qscat * average <cos (scattering angle)> 
c
c
c   A. Ackerman Feb-2000
c
c   --------------------------------------------------------------------

C      implicit none


c   Include common data shared with eddysed()

      include 'globals.h'
      include 'prog_params'
      PARAMETER (NSOL=38,NGS=8,NF=55, IK=8)
      COMMON /SPECTI/ BWNI(NSPC1IT),WNOI(NSPECIT),DWNI(NSPECIT),
     &                WLNI(NSPECIT)
      !COMMON/WAVE/AV(NF),LAM(NF),W(NF) ! Added WAVE common, need SOLWAVE damnit. JDW 2021
c   Declare local storage

      logical okay, read_mie, bad_mie, do_optics,do_IR_wav
      integer ngas, igas, nwave, iwave, nrad(MAXNGAS), irad
      integer mwave, mrad, iskip, isub, istatus
      integer n_thetd, nskip, ns, idum
      integer i
      character*(*) gas_name(ngas)
      character*(80) this_gas, filename
      character*1 dumstr
      double precision wave(MAXNWAVE)
      double precision radius(MAXNRAD,MAXNGAS)
      double precision dr(MAXNRAD,MAXNGAS)
      double precision rup(MAXNRAD,MAXNGAS)
      double precision qscat(MAXNWAVE,MAXNRAD,MAXNGAS)
      double precision qext(MAXNWAVE,MAXNRAD,MAXNGAS)
      double precision cos_qscat(MAXNWAVE,MAXNRAD,MAXNGAS)
      double precision rmin(MAXNGAS)
      double precision wvno, qs, qe, c_qs, corerad, corereal, coreimag
      double precision vrat, pw, f1, f2, r_in, wave_in
      double precision dr5, rr, nn, kk, thetd
      double precision qs_pass, qe_pass, c_qs_pass
      double precision bwni, wnoi, dwni, wlni

C   --------------------------------------------------------------------
C   Read in Specti files from egp JDW!

C      open(unit=001,file='init_optics_bwni.txt')
C      open(unit=002,file='init_optics_wnoi.txt')
C      open(unit=003,file='init_optics_dwni.txt')
C      open(unit=004,file='init_optics_wave.txt')



C      do i=1,NSPC1IT
C        read(001,*) bwni!bwni(i)
C        print*,'bwni',bwni(i)
C      enddo
C      read(002,*) wnoi
C      read(003,*) dwni
C      read(004,*) wlni


c   --------------------------------------------------------------------

c   Report input controls

      !write(*,*) ''
      !write(*,*) 'Begin init_optics'



c   If do_optics is false, zero nrad and nwave and return
c   (geometric scatterers only)

c      if( .not. do_optics )then
c        nrad = 0
c        nwave = 0
c        return
c      endif
C      print*, 'init_optics(): do_optics, read_mie = ',
C     & do_optics,read_mie
c   Define radius grids 
      do igas = 1, ngas
            nrad(igas) = 40
            rmin(igas) = 1e-5
            vrat = 2.2
 
        if( nrad(igas) .gt. MAXNRAD )then
          print*, 'init_optics(): nrad > MAXNRAD'
          stop 1
        endif

c      do igas = 1, ngas

c        rmin(igas) = 1e-7
       


c       if( gas_name(igas) .eq. 'CH4' )then
c       elseif( gas_name(igas) .eq. 'NH3' )then
c       elseif( gas_name(igas) .eq. 'H2O' )then
c       elseif( gas_name(igas) .eq. 'Fe' )then
c       elseif( gas_name(igas) .eq. 'KCl' )then
c       elseif( gas_name(igas) .eq. 'MgSiO3' )then
c       elseif( gas_name(igas) .eq. 'Al2O3' )then
c       endif

        pw = 1. / 3.
        f1 = ( 2*vrat / ( 1 + vrat) )**pw
        f2 = ( 2 / ( 1 + vrat ) )**pw * (vrat**pw-1)

        do irad = 1, nrad(igas)
          radius(irad,igas) = rmin(igas) * vrat**(float(irad-1)/3.)
          rup(irad,igas) = f1*radius(irad,igas)
          dr(irad,igas) = f2*radius(irad,igas)
        enddo

c hack: Banfield et al
c       radius(1,igas) = 0.2e-4
c       dr(1,igas) = f2*radius(1,igas)
c       rup(1,igas) = f1*radius(1,igas)

c hack: Brooke et al
c       radius(1,igas) = 1e-4
c       dr(1,igas) = f2*radius(1,igas)
c       rup(1,igas) = f1*radius(1,igas)
c       radius(2,igas) = 10e-4
c       dr(2,igas) = f2*radius(2,igas)
c       rup(2,igas) = f1*radius(2,igas)

      enddo

c   --------------------------------------------------------------------
c   Define optical properties

      !nwave = 196 !Why is this hardcoded? jdw
c hack
c     nwave = 1

      if( nwave .gt. MAXNWAVE )then
        print*, 'init_optics(): nwave > MAXNWAVE'
        stop 1
      endif

c      print*, 'init_optics(): nwave, nrad = ', nwave, nrad(igas)

c   --------------------------------------------------------------------
c   Read extinction and scattering coefficients for each condensing vapor

      if( read_mie )then

        do igas = 1, ngas
 
          this_gas = gas_name(igas)
          call dblank( this_gas, ns )
          IF(do_IR_wav)then
          filename = 'cloud_optics/' // this_gas(1:ns) // 'IR.mieff' !Make use of already calculated .mie arrays JDW
          print*,'filename',filename
          else 
          filename = 'cloud_optics/' // this_gas(1:ns) // 'SOL.mieff'
            print*,'filename',filename
          endif
          open( LUNIO, file=filename, form='formatted', status='old' ) !LUNIO is just like unit. jdw
 

c   Check that input file is consistent with radius and wavelength grids
 
          read( LUNIO, * ) mwave, mrad
          okay = (mwave .eq. nwave) .and. (mrad .eq. nrad(igas)) 
   
          if ( .not. okay )then
            print*,'init_optics(): input grid bad ',
     $        'for gas(', igas, ') = ', this_gas(1:ns), ':'
            print*,' mwave, mrad = ', mwave, mrad
            stop 1
          endif
 

          do irad = 1, nrad(igas)


c   Read and check input radii

            read( LUNIO, * ) r_in
C            print*,'is this reading in?'
C            print *,'radius=',r_in
            okay = abs( 1 - radius(irad,igas)/r_in ) .lt. 1d-6
            if ( .not. okay )then
              print*,'init_optics(): input radius grid bad ',
     $          'for gas(', igas, ') = ', this_gas(1:ns), ':'
              print*,' irad, radius, r_in = ', 
     $          irad, radius(irad,igas), r_in
              stop 1
            endif


c   Read and check wavelength and scattering efficiencies etc.

            do iwave = 1, nwave 

              read( LUNIO ,* ) wave_in, qscat(iwave,irad,igas),
     $          qext(iwave,irad,igas), cos_qscat(iwave,irad,igas)
c              print *, 'opt properties:', wave_in, qscat(iwave, irad, igas), qext(iwave, irad, igas)

        write(*,871)iwave,irad,igas,wlni(iwave),qext(iwave,irad,igas)
        
871   format (3i5,2f15.5)

              if(( igas .eq. 1 ).and.( irad .eq. 1)) then
                wave(iwave) = wave_in
                print*,'wave_in',wave_in
			    okay = abs( 1 - wlni(iwave)/(wave_in*1d4) ) .lt. 5d-2
                if ( .not. okay )then
                  print*,'init_optics(): input mieff wavelength grid ',
     $              'not same as gas opacity: '
                  print*,' iwave, wave, wave_in = ', 
     $              iwave, wlni(iwave), wave_in*1d4, abs( 1 - wlni(iwave)/(wave_in*1d4) )
                  stop 1
                endif

              else

                okay = abs( 1 - wave(iwave)/wave_in ) .lt. 1d-7  !jjf, was 1d-7
                if ( .not. okay )then
                  print*,'init_optics(): input wavelength grid bad ',
     $              'for gas(', igas, ') = ', this_gas(1:ns), ':'
                  print*,' iwave, wave, wave_in = ', 
     $              iwave, wave(iwave), wave_in, abs( 1 - wave(iwave)/wave_in )
c                  stop 1
                endif

              endif

            enddo
c            STOP
         enddo

          close( LUNIO )

        enddo

      else 

c   --------------------------------------------------------------------
c   Calculate single-scattering efficiencies etc from refractive indices
c   for each condensing vapor

C   This needs to be updated for the Clima wavelegnth grid JDW 2021
C   Specifically the H2O.refrind file needs to be updated




c   Mie parameters:
c   thetd is angle between incident and scattered radiation !Seems weird to let this
C   be a constant in the code. 
c   n_thetd is number of thetd values to consider

        thetd = 0.0
        n_thetd = 1

        do igas = 1, ngas

          this_gas = gas_name(igas)
          call dblank( this_gas, ns )
          IF(do_IR_wav)then
          filename = 'cloud_optics/' // this_gas(1:ns) // 'IR.refrind'
          else
          filename = 'cloud_optics/' // this_gas(1:ns) // 'SOL.refrind'
          endif
          open( LUNIO, file=filename, form='formatted', status='old' )


c   Skip filler and header lines

          nskip = 0
          do iskip = 1, nskip
            read( LUNIO, * ) dumstr
          enddo


          do iwave = 1, nwave

            read( LUNIO, * ) idum,wave_in, nn, kk

            if( igas .eq. 1 )then

              wave(iwave) =wave_in*1d-4!wave_in*1d-4 ! CAN I CHANGE WAVE_IN JDW
              ! Before the wave(iwave) was computed by using the 
c              wvno = 2*PI / wave(iwave)
              !print*,'wave',wave_in
            else

c   Consistency check

              okay = abs( 1 - wave(iwave)/(wave_in*1d-4) ) .lt. 1d-7
              if( .not. okay )then
                print*,'init_optics(): input wavelength grid bad ',
     $            'for gas(', igas, ') = ', gas_name(igas), ':'
                print*,' iwave, wave, wave_in = ', 
     $            iwave, wave(iwave), wave_in*1d-4
                stop 1
              endif

            endif
c
c*TDR       bug -- Mark says move this from above to here
c
            wvno = 2*PIE / wave(iwave)
c
   
            do irad = 1, nrad(igas)

c   --------------------------------------------------------------------
c   Subdivide radius grid into 6 bins (to avg out oscillations) and
c   call Mie code

              if( irad .eq. 1 )then
                dr5 = ( rup(1,igas) - radius(1,igas) ) / 5.
                rr  = radius(1,igas)
              else
                dr5 = ( rup(irad,igas) - rup(irad-1,igas) ) / 5.
                rr  = rup(irad-1,igas)
              endif

              qext(iwave,irad,igas) = 0.
              qscat(iwave,irad,igas) = 0.
              cos_qscat(iwave,irad,igas) = 0.

              corerad = 0.
              corereal = 1.
              coreimag = 0.

c    Only want one warning (not 6) message per radius bin
              bad_mie = .false.

              do isub = 1, 6

                call mie_calc( rr, nn, kk, thetd, n_thetd, 
     $            qe_pass, qs_pass, c_qs_pass, 
     $            corerad, corereal, coreimag, wvno,
     $            istatus )

                if( istatus .eq. 0 )then
                  qe = qe_pass
                  qs = qs_pass
                  c_qs = c_qs_pass
                else
                  if( .not. bad_mie )then
                    bad_mie = .true.
                    print*,'init_optics(): no Mie solution for '
     $                // 'irad, r(um), iwave, wave(um), n, k, gas = '
                  endif
                  print*,
     $              irad,rr*1e4,iwave,wave(iwave)*1e4,nn,kk,
     $              ' ',gas_name(igas)
                endif

                qext(iwave,irad,igas) = qext(iwave,irad,igas) + qe
                qscat(iwave,irad,igas) = qscat(iwave,irad,igas) + qs
                cos_qscat(iwave,irad,igas) = 
     $            cos_qscat(iwave,irad,igas) + c_qs

                rr = rr + dr5

              enddo
              
              qext(iwave,irad,igas) = qext(iwave,irad,igas) / 6.
              qscat(iwave,irad,igas) = qscat(iwave,irad,igas) / 6.
              cos_qscat(iwave,irad,igas) = 
     $          cos_qscat(iwave,irad,igas) / 6.

            enddo
          enddo
        enddo

        close(LUNIO)

c   --------------------------------------------------------------------
c   Write extinction and scattering coefficients
 
        do igas = 1, ngas

          this_gas = gas_name(igas)
          call dblank( this_gas, ns )
          IF(do_IR_wav)then
          filename =  'cloud_optics/' // this_gas(1:ns) // 'IR.mieff' !ISSUE HERE WITH MIE COEFFICIENTS NEEDED TO BE KEPT SEPARATE JDW 2021
          else
          filename =  'cloud_optics/' // this_gas(1:ns) // 'SOL.mieff' 
          endif
          open( LUNIO, file=filename, form='formatted', 
     $      status='replace' )
 
          write( LUNIO, * ) nwave, nrad(igas)
 
          do irad = 1, nrad(igas)
            write( LUNIO, * ) radius(irad,igas)
            do iwave = 1, nwave
              write( LUNIO ,8888 ) wave(iwave), qscat(iwave,irad,igas),
     $          qext(iwave,irad,igas), cos_qscat(iwave,irad,igas)
            enddo
          enddo
8888      format (4F15.11) 
          close( LUNIO )

        enddo
      endif
c
      !write(*,*) 'End init_optics'
      !write(*,*) ''

c
      return
      end

      subroutine dblank(s, ns)
c
c
c  @(#) dblank.f  McKie  Jun-1988
c  This routine finds index of last nonblank char in a string.
c
c  Argument list input:
c    s = string to be examined
c
c  Argument list output:
c    ns = # chars in s, up to & including last non-blank
c
c
c  Declare subprogram arg(s)
c
      character*(*) s
      integer ns
c
c
c  Find last non-blank char in string (beyond 1st char), return to caller
c
      ns = len(s)
 2100 if( ( s(ns:ns) .eq. ' ' ) .and. ( ns .gt. 1 ) )then
       ns = ns - 1
       goto 2100
      endif
      return
      end











            SUBROUTINE mie_calc( RO, RFR, RFI, THETD, JX, QEXT, QSCAT, CTBRQS,
     1                   R, RE2, TMAG2, WVNO, istatus )
C
C **********************************************************************
C    THIS SUBROUTINE COMPUTES MIE SCATTERING BY A STRATIFIED SPHERE,
C    I.E. A PARTICLE CONSISTING OF A SPHERICAL CORE SURROUNDED BY A
C    SPHERICAL SHELL.  THE BASIC CODE USED WAS THAT DESCRIBED IN THE
C    REPORT: " SUBROUTINES FOR COMPUTING THE PARAMETERS OF THE
C    ELECTROMAGNETIC RADIATION SCATTERED BY A SPHERE " J.V. DAVE,
C    I B M SCIENTIFIC CENTER, PALO ALTO , CALIFORNIA.
C    REPORT NO. 320 - 3236 .. MAY 1968 .
C
C    THE MODIFICATIONS FOR STRATIFIED SPHERES ARE DESCRIBED IN
C        TOON AND ACKERMAN, APPL. OPTICS, IN PRESS, 1981
C
C    THE PARAMETERS IN THE CALLING STATEMENT ARE DEFINED AS FOLLOWS :
C      RO IS THE OUTER (SHELL) RADIUS;
C      R  IS THE CORE RADIUS;
C      RFR, RFI  ARE THE REAL AND IMAGINARY PARTS OF THE SHELL INDEX
C          OF REFRACTION IN THE FORM (RFR - I* RFI);
C      RE2, TMAG2  ARE THE INDEX PARTS FOR THE CORE;
C          ( WE ASSUME SPACE HAS UNIT INDEX. )
C      THETD(J): ANGLE IN DEGREES BETWEEN THE DIRECTIONS OF THE INCIDENT
C          AND THE SCATTERED RADIATION.  THETD(J) IS< OR= 90.0
C          IF THETD(J) SHOULD HAPPEN TO BE GREATER THAN 90.0, ENTER WITH
C          SUPPLEMENTARY VALUE, SEE COMMENTS BELOW ON ELTRMX;
C      JX: TOTAL NUMBER OF THETD FOR WHICH THE COMPUTATIONS ARE
C          REQUIRED.  JX SHOULD NOT EXCEED IT UNLESS THE DIMENSIONS
C          STATEMENTS ARE APPROPRIATEDLY MODIFIED;
C
C      THE DEFINITIONS FOR THE FOLLOWING SYMBOLS CAN BE FOUND IN"LIGHT
C          SCATTERING BY SMALL PARTICLES,H.C.VAN DE HULST, JOHN WILEY '
C          SONS, INC., NEW YORK, 1957" .
C      QEXT: EFFICIENCY FACTOR FOR EXTINCTION,VAN DE HULST,P.14 ' 127.
C      QSCAT: EFFICIENCY FACTOR FOR SCATTERING,V.D. HULST,P.14 ' 127.
C      CTBRQS: AVERAGE(COSINE THETA) * QSCAT,VAN DE HULST,P.128
C      ELTRMX(I,J,K): ELEMENTS OF THE TRANSFORMATION MATRIX F,V.D.HULST
C          ,P.34,45 ' 125. I=1: ELEMENT M SUB 2..I=2: ELEMENT M SUB 1..
C          I = 3: ELEMENT S SUB 21.. I = 4: ELEMENT D SUB 21..
C      ELTRMX(I,J,1) REPRESENTS THE ITH ELEMENT OF THE MATRIX FOR
C          THE ANGLE THETD(J).. ELTRMX(I,J,2) REPRESENTS THE ITH ELEMENT
C          OF THE MATRIX FOR THE ANGLE 180.0 - THETD(J) ..
C      QBS IS THE BACK SCATTER CROSS SECTION.
C
C      IT: IS THE DIMENSION OF THETD, ELTRMX, CSTHT, PI, TAU, SI2THT,
C          IT MUST CORRESPOND EXACTLY TO THE SECOND DIMENSION OF ELTRMX.  
C      nacap IS THE DIMENSION OF ACAP 
C          IN THE ORIGINAL PROGRAM THE DIMENSION OF ACAP WAS 7000.
C          FOR CONSERVING SPACE THIS SHOULD BE NOT MUCH HIGHER THAN
C          THE VALUE, N=1.1*(NREAL**2 + NIMAG**2)**.5 * X + 1
C      WVNO: 2*PI / WAVELENGTH
C
C    ALSO THE SUBROUTINE COMPUTES THE CAPITAL A FUNCTION BY MAKING USE O
C    DOWNWARD RECURRENCE RELATIONSHIP.
C
C      TA(1): REAL PART OF WFN(1).  TA(2): IMAGINARY PART OF WFN(1).
C      TA(3): REAL PART OF WFN(2).  TA(4): IMAGINARY PART OF WFN(2).
C      TB(1): REAL PART OF FNA.     TB(2): IMAGINARY PART OF FNA.
C      TC(1): REAL PART OF FNB.     TC(2): IMAGINARY PART OF FNB.
C      TD(1): REAL PART OF FNAP.    TD(2): IMAGINARY PART OF FNAP.
C      TE(1): REAL PART OF FNBP.    TE(2): IMAGINARY PART OF FNBP.
C      FNAP, FNBP  ARE THE PRECEDING VALUES OF FNA, FNB RESPECTIVELY.
C **********************************************************************
c
c
c  Define implicit type for local floating point variables
c
c  This is for 32-bit single precision 
c
      implicit double precision ( a-h, o-z )
c
c  This is for 64-bit single precision
c
c     implicit real ( a-h, o-z )
c
c
c  Define implicit type for integer variables
c
      implicit integer ( i-n )
c
c
c  Define tolerance for series convergence in Mie calculations
c
c     parameter( EPSILON_MIE = 1d-14 )
      parameter( EPSILON_MIE = 1d-7 )
c
c
c  Explicitly declare floating point arguments as single-precision
c   (not done for calling routines in double precision)
c
c     real RO, RFR, RFI, THETD, QEXT, QSCAT, CTBRQS, R, RE2, TMAG2, WVNO
c
c
c   Define various dimensions of local arrays
c
      parameter( nacap = 1000000,IT=1)
      integer IT
c 
c
c   Declare local variables
c
c     COMPLEX        FNAP,   FNBP,   ACAP(nacap), W(3,nacap),
      double COMPLEX FNAP,   FNBP,   ACAP(nacap), W(3,nacap),
     2        FNA,    FNB,    RF,       RRF,
     3        RRFX,   WM1,    FN1,      FN2,
     4        TC1,    TC2,    WFN(2),   Z(4),
     5        K1,     K2,     K3,
     6        RC,     U(8),   DH1,
     7        DH2,    DH4,    P24H24,   P24H21,
     8        PSTORE, HSTORE, DUMMY,    DUMSQ
C
      COMMON / WARRAY / W
C
      DIMENSION     T(5),      TA(4),     TB(2),        TC(2),
     2              TD(2),     TE(2),     PI( 3,IT ),   TAU( 3,IT ),
     3              CSTHT(IT), THETD(IT), SI2THT(IT),   ELTRMX( 4,IT,2 )
C     

      ! DIMENSION     T(5),      TA(4),     TB(2),        TC(2),
      !2              TD(2),     TE(2),     PI( 3,IT ),   TAU( 3,IT ),
      !3              CSTHT(IT),SI2THT(IT),   ELTRMX( 4,IT,2 )


      !double precision THETD(IT)
      EQUIVALENCE   (FNA,TB(1)),(FNB,TC(1)),(FNAP,TD(1)),(FNBP,TE(1))
C
C   IF THE CORE IS SMALL SCATTERING IS COMPUTED FOR THE SHELL ONLY
C
      IFLAG = 1
      IF ( R/RO .LT. 1e-6 )   IFLAG = 2
      IF ( JX .LE. IT )   GO TO 20
         WRITE( *,7 )
         WRITE( *,6 )
         STOP 30
   20 RF =  CMPLX( RFR,  -RFI )
      RC =  CMPLX( RE2,-TMAG2 )
      X  =  RO * WVNO
      K1 =  RC * WVNO
      K2 =  RF * WVNO
      K3 =  CMPLX( WVNO, 0.0 )
      Z(1) =  K2 * RO
      Z(2) =  K3 * RO
      Z(3) =  K1 * R
      Z(4) =  K2 * R
      X1   =  REAL( Z(1) )
      X4   =  REAL( Z(4) )
c for AbSoft
c     Y1   =  dimag( Z(1) )
c     Y4   =  dimag( Z(4) )
      Y1   =  imag( Z(1) )
      Y4   =  imag( Z(4) )
      RRF  =  1.0 / RF
      RX   =  1.0 / X
      RRFX =  RRF * RX
      T(1) =  ( X**2 ) * ( RFR**2 + RFI**2 )
      T(1) =  SQRT( T(1) )
      NMX1 =  1.10 * T(1)
C
      IF ( NMX1 .LE. nacap-1 )   GO TO 21
         istatus = -1
c Ack hack April 2000 in response to Marley's RFR < 1
         return
c        WRITE(*,8)
c        STOP 32
   21 NMX2 = T(1)
      IF ( NMX1 .GT.  150 )   GO TO 22
         NMX1 = 150
         NMX2 = 135
C
   22 ACAP( NMX1+1 )  =  ( 0.0,0.0 )
      IF ( IFLAG .EQ. 2 )   GO TO 26
         DO 29   N = 1,3
   29    W( N,NMX1+1 )  =  ( 0.0,0.0 )
   26 CONTINUE
      DO 23   N = 1,NMX1
         NN = NMX1 - N + 1
         ACAP(NN) = (NN+1) * RRFX - 1.0 / ( (NN+1) * RRFX + ACAP(NN+1) )
         IF ( IFLAG .EQ. 2 )   GO TO 23
            DO 31   M = 1,3
   31       W( M,NN ) = (NN+1) / Z(M+1)  -
     1                   1.0 / (  (NN+1) / Z(M+1)  +  W( M,NN+1 )  )
   23 CONTINUE
C
      DO 30   J = 1,JX
      IF ( THETD(J) .LT. 0.0 )  THETD(J) =  ABS( THETD(J) )
      IF ( THETD(J) .GT. 0.0 )  GO TO 24
      CSTHT(J)  = 1.0
      SI2THT(J) = 0.0
      GO TO 30
   24 IF ( THETD(J) .GE. 90.0 )  GO TO 25
      T(1)      =  ( 3.14159265359 * THETD(J) ) / 180.0
      CSTHT(J)  =  COS( T(1) )
      SI2THT(J) =  1.0 - CSTHT(J)**2
      GO TO 30
   25 IF ( THETD(J) .GT. 90.0 )  GO TO 28
      CSTHT(J)  =  0.0
      SI2THT(J) =  1.0
      GO TO 30
   28 WRITE( *,5 )  THETD(J)
      WRITE( *,6 )
      STOP 34
   30 CONTINUE
C
      DO 35  J = 1,JX
      PI(1,J)  =  0.0
      PI(2,J)  =  1.0
      TAU(1,J) =  0.0
      TAU(2,J) =  CSTHT(J)
   35 CONTINUE
C
C INITIALIZATION OF HOMOGENEOUS SPHERE
C
      T(1)   =  COS(X)
      T(2)   =  SIN(X)
      WM1    =  CMPLX( T(1),-T(2) )
      WFN(1) =  CMPLX( T(2), T(1) )
      TA(1)  =  T(2)
      TA(2)  =  T(1)
      WFN(2) =  RX * WFN(1) - WM1
      TA(3)  =  REAL(WFN(2))
c dimag for AbSoft
c     TA(4)  =  dimag(WFN(2))
      TA(4)  =  imag(WFN(2))
C
      IF ( IFLAG .EQ. 2 )   GO TO 560
      N = 1
C
C INITIALIZATION PROCEDURE FOR STRATIFIED SPHERE BEGINS HERE
C
      SINX1   =  SIN( X1 )
      SINX4   =  SIN( X4 )
      COSX1   =  COS( X1 )
      COSX4   =  COS( X4 )
      EY1     =  EXP( Y1 )
      E2Y1    =  EY1 * EY1
      EY4     =  EXP( Y4 )
      EY1MY4  =  EXP( Y1 - Y4 )
      EY1PY4  =  EY1 * EY4
      EY1MY4  =  EXP( Y1 - Y4 )
      AA  =  SINX4 * ( EY1PY4 + EY1MY4 )
      BB  =  COSX4 * ( EY1PY4 - EY1MY4 )
      CC  =  SINX1 * ( E2Y1 + 1.0 )
      DD  =  COSX1 * ( E2Y1 - 1.0 )
      DENOM   =  1.0  +  E2Y1 * ( 4.0 * SINX1 * SINX1 - 2.0 + E2Y1 )
      REALP   =  ( AA * CC  +  BB * DD ) / DENOM
      AMAGP   =  ( BB * CC  -  AA * DD ) / DENOM
      DUMMY   =  CMPLX( REALP, AMAGP )
      AA  =  SINX4 * SINX4 - 0.5
      BB  =  COSX4 * SINX4
      P24H24  =  0.5 + CMPLX( AA,BB ) * EY4 * EY4
      AA  =  SINX1 * SINX4  -  COSX1 * COSX4
      BB  =  SINX1 * COSX4  +  COSX1 * SINX4
      CC  =  SINX1 * SINX4  +  COSX1 * COSX4
      DD  = -SINX1 * COSX4  +  COSX1 * SINX4
      P24H21  =  0.5 * CMPLX( AA,BB ) * EY1 * EY4  +
     2           0.5 * CMPLX( CC,DD ) * EY1MY4
      DH4  =  Z(4) / ( 1.0 + ( 0.0,1.0 ) * Z(4) )  -  1.0 / Z(4)
      DH1  =  Z(1) / ( 1.0 + ( 0.0,1.0 ) * Z(1) )  -  1.0 / Z(1)
      DH2  =  Z(2) / ( 1.0 + ( 0.0,1.0 ) * Z(2) )  -  1.0 / Z(2)
      PSTORE  =  ( DH4 + N / Z(4) )  *  ( W(3,N) + N / Z(4) )
      P24H24  =  P24H24 / PSTORE
      HSTORE  =  ( DH1 + N / Z(1) )  *  ( W(3,N) + N / Z(4) )
      P24H21  =  P24H21 / HSTORE
      PSTORE  =  ( ACAP(N) + N / Z(1) )  /  ( W(3,N) + N / Z(4) )
      DUMMY   =  DUMMY * PSTORE
      DUMSQ   =  DUMMY * DUMMY
C
C NOTE:  THE DEFINITIONS OF U(I) IN THIS PROGRAM ARE NOT THE SAME AS
C        THE USUBI DEFINED IN THE ARTICLE BY TOON AND ACKERMAN.  THE
C        CORRESPONDING TERMS ARE:
C          USUB1 = U(1)                       USUB2 = U(5)
C          USUB3 = U(7)                       USUB4 = DUMSQ
C          USUB5 = U(2)                       USUB6 = U(3)
C          USUB7 = U(6)                       USUB8 = U(4)
C          RATIO OF SPHERICAL BESSEL FTN TO SPHERICAL HENKAL FTN = U(8)
C
      U(1) =  K3 * ACAP(N)  -  K2 * W(1,N)
      U(2) =  K3 * ACAP(N)  -  K2 * DH2
      U(3) =  K2 * ACAP(N)  -  K3 * W(1,N)
      U(4) =  K2 * ACAP(N)  -  K3 * DH2
      U(5) =  K1 *  W(3,N)  -  K2 * W(2,N)
      U(6) =  K2 *  W(3,N)  -  K1 * W(2,N)
      U(7) =  ( 0.0,-1.0 )  *  ( DUMMY * P24H21 - P24H24 )
      U(8) =  TA(3) / WFN(2)
C
      FNA  =  U(8) * ( U(1)*U(5)*U(7)  +  K1*U(1)  -  DUMSQ*K3*U(5) ) /
     2               ( U(2)*U(5)*U(7)  +  K1*U(2)  -  DUMSQ*K3*U(5) )
      FNB  =  U(8) * ( U(3)*U(6)*U(7)  +  K2*U(3)  -  DUMSQ*K2*U(6) ) /
     2               ( U(4)*U(6)*U(7)  +  K2*U(4)  -  DUMSQ*K2*U(6) )
      GO TO 561
  560 TC1  =  ACAP(1) * RRF  +  RX
      TC2  =  ACAP(1) * RF   +  RX
      FNA  =  ( TC1 * TA(3)  -  TA(1) ) / ( TC1 * WFN(2)  -  WFN(1) )
      FNB  =  ( TC2 * TA(3)  -  TA(1) ) / ( TC2 * WFN(2)  -  WFN(1) )
C
  561 CONTINUE
      FNAP = FNA
      FNBP = FNB
      T(1) = 1.50
C
C    FROM HERE TO THE STATMENT NUMBER 90, ELTRMX(I,J,K) HAS
C    FOLLOWING MEANING:
C    ELTRMX(1,J,K): REAL PART OF THE FIRST COMPLEX AMPLITUDE.
C    ELTRMX(2,J,K): IMAGINARY PART OF THE FIRST COMPLEX AMPLITUDE.
C    ELTRMX(3,J,K): REAL PART OF THE SECOND COMPLEX AMPLITUDE.
C    ELTRMX(4,J,K): IMAGINARY PART OF THE SECOND COMPLEX AMPLITUDE.
C    K = 1 : FOR THETD(J) AND K = 2 : FOR 180.0 - THETD(J)
C    DEFINITION OF THE COMPLEX AMPLITUDE: VAN DE HULST,P.125.
C
      TB(1) = T(1) * TB(1)
      TB(2) = T(1) * TB(2)
      TC(1) = T(1) * TC(1)
      TC(2) = T(1) * TC(2)
      DO 60 J = 1,JX
          ELTRMX(1,J,1) = TB(1) * PI(2,J) + TC(1) * TAU(2,J)
          ELTRMX(2,J,1) = TB(2) * PI(2,J) + TC(2) * TAU(2,J)
          ELTRMX(3,J,1) = TC(1) * PI(2,J) + TB(1) * TAU(2,J)
          ELTRMX(4,J,1) = TC(2) * PI(2,J) + TB(2) * TAU(2,J)
          ELTRMX(1,J,2) = TB(1) * PI(2,J) - TC(1) * TAU(2,J)
          ELTRMX(2,J,2) = TB(2) * PI(2,J) - TC(2) * TAU(2,J)
          ELTRMX(3,J,2) = TC(1) * PI(2,J) - TB(1) * TAU(2,J)
          ELTRMX(4,J,2) = TC(2) * PI(2,J) - TB(2) * TAU(2,J)
   60 CONTINUE
C
      QEXT   = 2.0 * ( TB(1) + TC(1))
      QSCAT  = ( TB(1)**2 + TB(2)**2 + TC(1)**2 + TC(2)**2 ) / 0.75
      CTBRQS = 0.0
      QBSR   = -2.0*(TC(1) - TB(1))
      QBSI   = -2.0*(TC(2) - TB(2))
      RMM    = -1.0
      N = 2
   65 T(1) = 2*N - 1
      T(2) =   N - 1
      T(3) = 2*N + 1
      DO 70  J = 1,JX
          PI(3,J)  = ( T(1) * PI(2,J) * CSTHT(J) - N * PI(1,J) ) / T(2)
          TAU(3,J) = CSTHT(J) * ( PI(3,J) - PI(1,J) )  -
     1                          T(1) * SI2THT(J) * PI(2,J)  +  TAU(1,J)
   70 CONTINUE
C
C HERE SET UP HOMOGENEOUS SPHERE
C
      WM1    =  WFN(1)
      WFN(1) =  WFN(2)
      TA(1)  =  REAL(WFN(1))
c dimag for AbSoft
c     TA(2)  =  dimag(WFN(1))
c     TA(4)  =  dimag(WFN(2))
      TA(2)  =  imag(WFN(1))
      TA(4)  =  imag(WFN(2))
      WFN(2) =  T(1) * RX * WFN(1)  -  WM1
      TA(3)  =  REAL(WFN(2))
C
      IF ( IFLAG .EQ. 2 )   GO TO 1000
C
C HERE SET UP STRATIFIED SPHERE
C
      DH2  =  - N / Z(2)  +  1.0 / ( N / Z(2) - DH2 )
      DH4  =  - N / Z(4)  +  1.0 / ( N / Z(4) - DH4 )
      DH1  =  - N / Z(1)  +  1.0 / ( N / Z(1) - DH1 )
      PSTORE  =  ( DH4 + N / Z(4) )  *  ( W(3,N) + N / Z(4) )
      P24H24  =  P24H24 / PSTORE
      HSTORE  =  ( DH1 + N / Z(1) )  *  ( W(3,N) + N / Z(4) )
      P24H21  =  P24H21 / HSTORE
      PSTORE  =  ( ACAP(N) + N / Z(1) )  /  ( W(3,N) + N / Z(4) )
      DUMMY   =  DUMMY * PSTORE
      DUMSQ   =  DUMMY * DUMMY
C
      U(1) =  K3 * ACAP(N)  -  K2 * W(1,N)
      U(2) =  K3 * ACAP(N)  -  K2 * DH2
      U(3) =  K2 * ACAP(N)  -  K3 * W(1,N)
      U(4) =  K2 * ACAP(N)  -  K3 * DH2
      U(5) =  K1 *  W(3,N)  -  K2 * W(2,N)
      U(6) =  K2 *  W(3,N)  -  K1 * W(2,N)
      U(7) =  ( 0.0,-1.0 )  *  ( DUMMY * P24H21 - P24H24 )
      U(8) =  TA(3) / WFN(2)
C
      FNA  =  U(8) * ( U(1)*U(5)*U(7)  +  K1*U(1)  -  DUMSQ*K3*U(5) ) /
     2               ( U(2)*U(5)*U(7)  +  K1*U(2)  -  DUMSQ*K3*U(5) )
      FNB  =  U(8) * ( U(3)*U(6)*U(7)  +  K2*U(3)  -  DUMSQ*K2*U(6) ) /
     2               ( U(4)*U(6)*U(7)  +  K2*U(4)  -  DUMSQ*K2*U(6) )
C
 1000 CONTINUE
      TC1  =  ACAP(N) * RRF  +  N * RX
      TC2  =  ACAP(N) * RF   +  N * RX
      FN1  =  ( TC1 * TA(3)  -  TA(1) ) /  ( TC1 * WFN(2) - WFN(1) )
      FN2  =  ( TC2 * TA(3)  -  TA(1) ) /  ( TC2 * WFN(2) - WFN(1) )
      M    =  WVNO * R
      IF ( N .LT. M )   GO TO 1002
      IF ( IFLAG .EQ. 2 )   GO TO 1001
      IF ( abs(  ( FN1-FNA ) / FN1  )  .LT.  EPSILON_MIE   .AND.
     1     abs(  ( FN2-FNB ) / FN2  )  .LT . EPSILON_MIE  ) IFLAG = 2
      IF ( IFLAG .EQ. 1 )   GO TO 1002
 1001 FNA  =  FN1
      FNB  =  FN2
C
 1002 CONTINUE
      T(5)  =  N
      T(4)  =  T(1) / ( T(5) * T(2) )
      T(2)  =  (  T(2) * ( T(5) + 1.0 )  ) / T(5)
C
      CTBRQS  =  CTBRQS  +  T(2) * ( TD(1) * TB(1)  +  TD(2) * TB(2)
     1                   +           TE(1) * TC(1)  +  TE(2) * TC(2) )
     2                   +  T(4) * ( TD(1) * TE(1)  +  TD(2) * TE(2) )
      QEXT    =   QEXT  +  T(3) * ( TB(1) + TC(1) )
      T(4)    =  TB(1)**2 + TB(2)**2 + TC(1)**2 + TC(2)**2
      QSCAT   =  QSCAT  +  T(3) * T(4)
      RMM     =  -RMM
      QBSR    =  QBSR + T(3)*RMM*(TC(1) - TB(1))
      QBSI    =  QBSI + T(3)*RMM*(TC(2) - TB(2))
C
      T(2)    =  N * (N+1)
      T(1)    =  T(3) / T(2)
      K = (N/2)*2
      DO 80 J = 1,JX
       ELTRMX(1,J,1) = ELTRMX(1,J,1)+T(1)*(TB(1)*PI(3,J)+TC(1)*TAU(3,J))
       ELTRMX(2,J,1) = ELTRMX(2,J,1)+T(1)*(TB(2)*PI(3,J)+TC(2)*TAU(3,J))
       ELTRMX(3,J,1) = ELTRMX(3,J,1)+T(1)*(TC(1)*PI(3,J)+TB(1)*TAU(3,J))
       ELTRMX(4,J,1) = ELTRMX(4,J,1)+T(1)*(TC(2)*PI(3,J)+TB(2)*TAU(3,J))
      IF ( K .EQ. N )  THEN
       ELTRMX(1,J,2) =ELTRMX(1,J,2)+T(1)*(-TB(1)*PI(3,J)+TC(1)*TAU(3,J))
       ELTRMX(2,J,2) =ELTRMX(2,J,2)+T(1)*(-TB(2)*PI(3,J)+TC(2)*TAU(3,J))
       ELTRMX(3,J,2) =ELTRMX(3,J,2)+T(1)*(-TC(1)*PI(3,J)+TB(1)*TAU(3,J))
       ELTRMX(4,J,2) =ELTRMX(4,J,2)+T(1)*(-TC(2)*PI(3,J)+TB(2)*TAU(3,J))
      ELSE
       ELTRMX(1,J,2) = ELTRMX(1,J,2)+T(1)*(TB(1)*PI(3,J)-TC(1)*TAU(3,J))
       ELTRMX(2,J,2) = ELTRMX(2,J,2)+T(1)*(TB(2)*PI(3,J)-TC(2)*TAU(3,J))
       ELTRMX(3,J,2) = ELTRMX(3,J,2)+T(1)*(TC(1)*PI(3,J)-TB(1)*TAU(3,J))
       ELTRMX(4,J,2) = ELTRMX(4,J,2)+T(1)*(TC(2)*PI(3,J)-TB(2)*TAU(3,J))
      END IF
   80 CONTINUE
C
      IF ( T(4) .LT. EPSILON_MIE )   GO TO 100
      N = N + 1
      DO 90 J = 1,JX
         PI(1,J)   =   PI(2,J)
         PI(2,J)   =   PI(3,J)
         TAU(1,J)  =  TAU(2,J)
         TAU(2,J)  =  TAU(3,J)
   90 CONTINUE
      FNAP  =  FNA
      FNBP  =  FNB
      IF ( N .LE. NMX2 )   GO TO 65
         istatus = -1
c Ack hack April 2000 in response to Marley's RFR < 1
         return
c        WRITE( *,8 )
c        STOP 36
  100 DO 120 J = 1,JX
      DO 120 K = 1,2
         DO  115  I= 1,4
         T(I)  =  ELTRMX(I,J,K)
  115    CONTINUE
         ELTRMX(2,J,K)  =      T(1)**2  +  T(2)**2
         ELTRMX(1,J,K)  =      T(3)**2  +  T(4)**2
         ELTRMX(3,J,K)  =  T(1) * T(3)  +  T(2) * T(4)
         ELTRMX(4,J,K)  =  T(2) * T(3)  -  T(4) * T(1)
  120 CONTINUE
      T(1)    =    2.0 * RX**2
      QEXT    =   QEXT * T(1)
      QSCAT   =  QSCAT * T(1)
      CTBRQS  =  2.0 * CTBRQS * T(1)
C
C QBS IS THE BACK SCATTER CROSS SECTION
C
c     PIG   = ACOS(-1.0)
c     RXP4  = RX*RX/(4.0*PIG)
c     QBS   = RXP4*(QBSR**2 + QBSI**2)
C
    5 FORMAT( 10X,' THE VALUE OF THE SCATTERING ANGLE IS GREATER THAN
     1 90.0 DEGREES. IT IS ', E15.4 )
    6 FORMAT( // 10X, 'PLEASE READ COMMENTS.' // )
    7 FORMAT( // 10X, 'THE VALUE OF THE ARGUMENT JX IS GREATER THAN IT')
    8 FORMAT( // 10X, 'THE UPPER LIMIT FOR ACAP IS NOT ENOUGH. SUGGEST
     1 GET DETAILED OUTPUT AND MODIFY SUBROUTINE' // )
C
      istatus = 0
      RETURN
      END

      subroutine calc_optics( do_subcloud, nz, ngas, nwave, nrad,
     $  z, T, gas_name, radius, dr, qscat, qext, cos_qscat, 
     $  ndz, sig, rg, opd, opd_gas, w0, g0, opd_tot,do_IR_wav )

c   --------------------------------------------------------------------
c
c   Calculate spectrally-resolved profiles of optical depth, single-scattering
c   albedo, and asymmetry parameter.
c
c   input scalars:
c
c     do_subcloud  logical flag for Marley subcloud kludge 
c     nz           number of layers
c     ngas         number of condensing gases
c     nwave        number of wavelength bins
c     nrad         number of radius bins
c
c   input vectors:
c
c     z            layer altitude (cm)
c     T            layer temperature (K)
c     gas_name     names of condensing gases
c     radius       radius bin centers (cm)
c     dr           width of radius bins (cm)
c     qscat        scattering efficiency
c     qext         scattering efficiency
c     cos_qscat    qscat-weighted <cos (scattering angle)>
c     ndz          number column density of condensate (cm^-3) !JDW these units are wrong and should be (#/cm^2)
c     sig          geometric standard deviation of lognormal size distribution
c     rg           geometric mean radius of lognormal size distribution
c
c   output scalars:
c
c     opd_tot      total optical depth for geometric conservative scatterers
c
c   output vectors:
c
c     opd          extinction optical depth due to all condensates in layer
c     opd_gas      cumulative (from top) opd by condensing vapor as
c                  geometric conservative scatterers
c     w0           single scattering albedo
c     g0           asymmetry parameter = Q_scat wtd avg of <cos theta>
c
c
c   A. Ackerman Apr-2000
c
c   --------------------------------------------------------------------

      implicit none


c   Include common data shared with eddysed()

      include 'globals.h'
      include 'prog_params'


c   Declare local storage

      logical do_subcloud,do_IR_wav
      integer nz, ngas, iz, igas, nwave, nrad(MAXNGAS), iwave, irad
      integer itop, ibot, incr, ibot_subcloud, ibot_gas
      integer ibot_cloud(MAXNGAS)
      COMMON /SPECTI/ BWNI(NSPC1IT),WNOI(NSPECIT),DWNI(NSPECIT)
     &           ,WLNI(NSPECIT)
      common/smart/QextsIR,QExtsSOL,QscasIR,QscasSOL,g0sIR,g0sSOL
     &,dtauextsIR,dtauextsSOL,Qexts,Qscas,g0s,dtauexts !JDW+TDR

      character*(*) gas_name(MAXNGAS)
      double precision bwni,wnoI,DWNI,WLNI
      double precision z(nlevel),T(nlevel)
      double precision radius(MAXNRAD,MAXNGAS)
      double precision dr(MAXNRAD,MAXNGAS)
      double precision qscat(MAXNWAVE,MAXNRAD,MAXNGAS)
      double precision qext(MAXNWAVE,MAXNRAD,MAXNGAS)
      double precision cos_qscat(MAXNWAVE,MAXNRAD,MAXNGAS)
      double precision scat_gas(MAXNZ,MAXNWAVE,MAXNGAS)
      double precision ext_gas(MAXNZ,MAXNWAVE,MAXNGAS)
      double precision cqs_gas(MAXNZ,MAXNWAVE,MAXNGAS)
      double precision ndz(MAXNZ,MAXNGAS)
      double precision sig(MAXNZ,MAXNGAS)
      double precision rg(MAXNZ,MAXNGAS)
      double precision opd(MAXNZ,MAXNWAVE)
      double precision w0(MAXNZ,MAXNWAVE)
      double precision g0(MAXNZ,MAXNWAVE)
      double precision opd_gas(MAXNZ,MAXNGAS)
      double precision opd_layer(MAXNZ,MAXNGAS)
      double precision opd_tot, rr, r2, norm
      double precision opd_scat, opd_ext, cos_qs
      double precision rsig, pir2ndz, arg1, arg2
      double precision, dimension(MAXNZ,MAXNWAVE,MAXNGAS) :: QextsIR
     & ,QscasIR
     & ,g0sIR, dtauextsIR,QextsSOL,QscasSOL,g0sSOL,dtauextsSOL,Qexts,
     & Qscas,g0s,dtauexts !JDW
       double precision norms(MAXNWAVE) !JDW

      IF(do_IR_wav)then
c   Determine indices of top and bottom layers

      if( z(2) .gt. z(1) )then
        itop = nz
        ibot = 1
        incr = -1
      else
        itop = 1
        ibot = nz
        incr = 1
      endif


c   Initialize indices of bottoms of cloud layers for subcloud kludge

      if( do_subcloud )then
        do igas = 1, ngas
          ibot_cloud(igas) = ibot
        enddo
      endif

c   --------------------------------------------------------------------
c   Loop over layers and condensing vapors

      do iz = itop, ibot, incr
        do igas = 1, ngas

          opd_layer(iz,igas) = 0.


c   Zero spectral sums

          do iwave = 1, nwave            !Could be problem here with nwave JDW
            scat_gas(iz,iwave,igas) = 0.
            ext_gas(iz,iwave,igas) = 0.
            cqs_gas(iz,iwave,igas) = 0.


            Qexts(iz,iwave,igas) = 0. !JDW
            Qscas(iz,iwave,igas) = 0. !JDW
            g0s(iz,iwave,igas) = 0. !JDW

C            QextsSOL(iz,iwave,igas) = 0. !JDW
C            QscasSOL(iz,iwave,igas) = 0. !JDW
C            g0sSOL(iz,iwave,igas) = 0. !JDW
          enddo

  
c   Optical depth for conservative geometric scatterers 

          if( ndz(iz,igas) .gt. 0. )then

            r2 = rg(iz,igas)**2 * exp( 2*log( sig(iz,igas) )**2 )
            opd_layer(iz,igas) = 2.*PIE*r2*ndz(iz,igas)

c  Calculate normalization factor (forces lognormal sum = 1.0)

            rsig = sig(iz,igas)
            norm = 0.

            do irad = 1,nrad(igas)
              rr = radius(irad,igas)
              arg1 = dr(irad,igas) / ( sqrt(2.*PIE)*rr*log(rsig) )
              arg2 = -log( rr/rg(iz,igas) )**2 / ( 2*log(rsig)**2 )
              norm = norm + arg1*exp( arg2 )
            enddo

            norm = ndz(iz,igas) / norm



C
C*JDW+TDR for computing number density-weighted mean Qexts and Qscas

C When writing to the qextsIR, it is overwritten by the first 38 spectral values with zeros. 
C I must be able to initalize them atthe same time.
C
            do iwave =1,nwave
              norms(iwave) = 0.
            enddo

c   Loop over wavelength and radius

            do irad = 1, nrad(igas)

              rr = radius(irad,igas)
              arg1 = dr(irad,igas) / ( sqrt(2.*PIE)*log(rsig) )
              arg2 = -log( rr/rg(iz,igas) )**2 / ( 2*log(rsig)**2 )
              pir2ndz = norm*PIE*rr*arg1*exp( arg2 )

              do iwave = 1, nwave
                scat_gas(iz,iwave,igas) = scat_gas(iz,iwave,igas) + 
     $            qscat(iwave,irad,igas)*pir2ndz
                ext_gas(iz,iwave,igas) = ext_gas(iz,iwave,igas) + 
     $            qext(iwave,irad,igas)*pir2ndz
                cqs_gas(iz,iwave,igas) = cqs_gas(iz,iwave,igas) + 
     $            cos_qscat(iwave,irad,igas)*pir2ndz
C                IF(do_IR_wav)then
                Qexts(iz,iwave,igas) = Qexts(iz,iwave,igas) +
     &           qext(iwave,irad,igas)*pir2ndz !JDW+TDR
                Qscas(iz,iwave,igas) = Qscas(iz,iwave,igas) +
     &           qscat(iwave,irad,igas)*pir2ndz !JDW+TDR
C                else
C                QextsSOL(iz,iwave,igas) = QextsSOL(iz,iwave,igas) +
C     &           qext(iwave,irad,igas)*pir2ndz !JDW+TDR
C                QscasSOL(iz,iwave,igas) = QscasSOL(iz,iwave,igas) +
C     &           qscat(iwave,irad,igas)*pir2ndz !JDW+TDR
C                endif
                  norms(iwave)=norms(iwave) + pir2ndz !JDW + TDR
              enddo

            enddo   ! irad:1,nrad
C
C JDW + TDR, compute number density-weighted mean Qexts and Qscats
C

            do iwave=1,nwave
C              IF(do_IR_wav)then
              Qexts(iz,iwave,igas)=Qexts(iz,iwave,igas)/norms(iwave)
              Qscas(iz,iwave,igas)=Qscas(iz,iwave,igas)/norms(iwave)
C              else
C             QextsSOL(iz,iwave,igas)=QextsSOL(iz,iwave,igas)/norms(iwave)
C              QscasSOL(iz,iwave,igas)=QscasSOL(iz,iwave,igas)/norms(iwave)
C              endif
            enddo 




c   index of bottom of cloud layer for subcloud kludge

            if( do_subcloud )then
              ibot_cloud(igas) = iz
            endif

          endif     ! ndz > 0
        enddo       ! igas:1,ngas
      enddo         ! iz:itop,ibot,incr

c   --------------------------------------------------------------------
c   subcloud kludge to soften discontinuity at cloud base
c   (10% in first layer below, 5% in second layer)

      if( do_subcloud )then
        do igas = 1, ngas
          ibot_gas = ibot_cloud(igas)
          if( ibot_gas .ne. ibot )then


c    Choose lower index to be within grid

            if( ibot_gas+incr .eq. ibot )then
              ibot_subcloud = ibot_gas + incr
            else
              ibot_subcloud = ibot_gas + 2*incr
            endif


            do iz = ibot_gas+incr, ibot_subcloud, incr

              if( iz .eq. ibot_gas+incr )then
                norm = 0.10
              else
                norm = 0.05
              endif




C
C JDW + TDR for computing number density-weighted mean Qexts and Qscas
C
              do iwave=1,nwave
                norms(iwave)=0.
              enddo 

              opd_layer(iz,igas) = opd_layer(iz,igas) + 
     $          opd_layer(ibot_gas,igas)*norm

              do iwave = 1, nwave
                scat_gas(iz,iwave,igas) = scat_gas(iz,iwave,igas) +  !Likely there is a bug in scat_gas+ due to initialization with waves diffent JDW 2021
     $            scat_gas(ibot_gas,iwave,igas)*norm
                ext_gas(iz,iwave,igas) = ext_gas(iz,iwave,igas) + 
     $            ext_gas(ibot_gas,iwave,igas)*norm
                cqs_gas(iz,iwave,igas) = cqs_gas(iz,iwave,igas) + 
     $            cqs_gas(ibot_gas,iwave,igas)*norm
C                IF(do_IR_wav)then
                Qexts(iz,iwave,igas)=Qexts(iz,iwave,igas) + 
     &           qext(iwave,irad,igas)*norm !JDW+TDR
                Qscas(iz,iwave,igas)=Qscas(iz,iwave,igas) + 
     &           qscat(iwave,irad,igas)*norm !JDW+TDR
                norms(iwave)=norms(iwave) + norm !JDW + TDR
C                else  

C                QextsSOL(iz,iwave,igas)=QextsSOL(iz,iwave,igas) + 
C     &           qext(iwave,irad,igas)*norm !JDW+TDR
C                QscasSOL(iz,iwave,igas)=QscasSOL(iz,iwave,igas) + 
C     &           qscat(iwave,irad,igas)*norm !JDW+TDR
                !norms(iwave)=norms(iwave) + norm !JDW + TDR
C                endif


              enddo

C
C JDW + TDR compute number density-weighted mean Qexts and Qscas
C
              do iwave=1,nwave
C                IF(do_IR_wav)then
                Qexts(iz,iwave,igas)=Qexts(iz,iwave,igas)/norms(iwave)
                Qscas(iz,iwave,igas)=Qscas(iz,iwave,igas)/norms(iwave)
C                else 

C                QextsSOL(iz,iwave,igas)=QextsSOL(iz,iwave,igas)/norms(iwave)
C                QscasSOL(iz,iwave,igas)=QscasSOL(iz,iwave,igas)/norms(iwave)
C                endif
              enddo

            enddo
          endif
        enddo 
      endif 

c   --------------------------------------------------------------------
c   Sum over gases and compute spectral optical depth profile etc

      do iz = itop, ibot, incr
        do iwave = 1, nwave

          opd_scat = 0.
          opd_ext = 0.
          cos_qs = 0.

          do igas = 1, ngas
            opd_scat = opd_scat + scat_gas(iz,iwave,igas)
            opd_ext = opd_ext + ext_gas(iz,iwave,igas)
            cos_qs = cos_qs + cqs_gas(iz,iwave,igas)
            if( scat_gas(iz,iwave,igas) .gt. 0. )then !JDW+TDR
C              IF(do_IR_wav)then
              g0s(iz,iwave,igas) = cqs_gas(iz,iwave,igas)
     &        /scat_gas(iz,iwave,igas) !JDW+TDR
C              else
C              g0sSOL(iz,iwave,igas) = cqs_gas(iz,iwave,igas)
C     &        /scat_gas(iz,iwave,igas) !JDW+TDR
C              endif
            else !JDW+TDR
C              IF(do_IR_wav)then
              g0s(iz,iwave,igas) = 0. !JDW +TDR
C              else
C              g0sSOL(iz,iwave,igas) = 0. !JDW +TDR
C             endif
              
            endif !JDW+TDR
C            IF(do_IR_wav)then
            dtauexts(iz,iwave,igas) = ext_gas(iz,iwave,igas) !JDW+TDR
C            else
C            dtauextsSOL(iz,iwave,igas) = ext_gas(iz,iwave,igas) !JDW+TDR
C            endif
c			if (igas.eq.3) print *,iz,iwave,wlni(iwave),ext_gas(iz,iwave,igas)
          enddo

          if( opd_scat .gt. 0. )then
            
            opd(iz,iwave) = opd_ext
            w0(iz,iwave) = opd_scat / opd_ext
            g0(iz,iwave) = cos_qs / opd_scat
          else 
            opd(iz,iwave) = 0.
            w0(iz,iwave) = 0.
            g0(iz,iwave) = 0.
          endif
          !print*,opd(iz,1),iz
        enddo
      enddo


c   cumulative optical depths for conservative geometric scatterers
C   Add Hack to control opd values. i.e. if greater than 100, they are set to 100. JDW Hacked 2021
      opd_tot = 0.
      do igas = 1, ngas
        opd_gas(itop,igas) = opd_layer(itop,igas)
        do iz = itop+incr, ibot, incr
          opd_gas(iz,igas) = opd_gas(iz-incr,igas) + opd_layer(iz,igas)
        enddo
        opd_tot = opd_tot + opd_gas(ibot,igas)
      enddo
c		stop
      !print*,opd





C      IF(do_IR_wav)then 
      dtauextsIR=dtauexts
      g0sIR=g0s !JDW Hacked
      QextsIR=Qexts
      QscasIR=Qscas
C      dtauextsSOL=dtauexts
C      g0sSOL=g0s
C      QextsSOL=Qexts
C      QscasSOL=Qscas



      else 




c   Determine indices of top and bottom layers

      if( z(2) .gt. z(1) )then
        itop = nz
        ibot = 1
        incr = -1
      else
        itop = 1
        ibot = nz
        incr = 1
      endif


c   Initialize indices of bottoms of cloud layers for subcloud kludge

      if( do_subcloud )then
        do igas = 1, ngas
          ibot_cloud(igas) = ibot
        enddo
      endif

c   --------------------------------------------------------------------
c   Loop over layers and condensing vapors

      do iz = itop, ibot, incr
        do igas = 1, ngas

          opd_layer(iz,igas) = 0.


c   Zero spectral sums

          do iwave = 1, nwave            !Could be problem here with nwave JDW
            scat_gas(iz,iwave,igas) = 0.
            ext_gas(iz,iwave,igas) = 0.
            cqs_gas(iz,iwave,igas) = 0.


            Qexts(iz,iwave,igas) = 0. !JDW
            Qscas(iz,iwave,igas) = 0. !JDW
            g0s(iz,iwave,igas) = 0. !JDW

C            QextsSOL(iz,iwave,igas) = 0. !JDW
C            QscasSOL(iz,iwave,igas) = 0. !JDW
C            g0sSOL(iz,iwave,igas) = 0. !JDW
          enddo

  
c   Optical depth for conservative geometric scatterers 

          if( ndz(iz,igas) .gt. 0. )then

            r2 = rg(iz,igas)**2 * exp( 2*log( sig(iz,igas) )**2 )
            opd_layer(iz,igas) = 2.*PIE*r2*ndz(iz,igas)

c  Calculate normalization factor (forces lognormal sum = 1.0)

            rsig = sig(iz,igas)
            norm = 0.

            do irad = 1,nrad(igas)
              rr = radius(irad,igas)
              arg1 = dr(irad,igas) / ( sqrt(2.*PIE)*rr*log(rsig) )
              arg2 = -log( rr/rg(iz,igas) )**2 / ( 2*log(rsig)**2 )
              norm = norm + arg1*exp( arg2 )
            enddo

            norm = ndz(iz,igas) / norm



C
C*JDW+TDR for computing number density-weighted mean Qexts and Qscas

C When writing to the qextsIR, it is overwritten by the first 38 spectral values with zeros. 
C I must be able to initalize them atthe same time.
C
            do iwave =1,nwave
              norms(iwave) = 0.
            enddo

c   Loop over wavelength and radius

            do irad = 1, nrad(igas)

              rr = radius(irad,igas)
              arg1 = dr(irad,igas) / ( sqrt(2.*PIE)*log(rsig) )
              arg2 = -log( rr/rg(iz,igas) )**2 / ( 2*log(rsig)**2 )
              pir2ndz = norm*PIE*rr*arg1*exp( arg2 )

              do iwave = 1, nwave
                scat_gas(iz,iwave,igas) = scat_gas(iz,iwave,igas) + 
     $            qscat(iwave,irad,igas)*pir2ndz
                ext_gas(iz,iwave,igas) = ext_gas(iz,iwave,igas) + 
     $            qext(iwave,irad,igas)*pir2ndz
                cqs_gas(iz,iwave,igas) = cqs_gas(iz,iwave,igas) + 
     $            cos_qscat(iwave,irad,igas)*pir2ndz
C                IF(do_IR_wav)then
                Qexts(iz,iwave,igas) = Qexts(iz,iwave,igas) +
     &           qext(iwave,irad,igas)*pir2ndz !JDW+TDR
                Qscas(iz,iwave,igas) = Qscas(iz,iwave,igas) +
     &           qscat(iwave,irad,igas)*pir2ndz !JDW+TDR
C                else
C                QextsSOL(iz,iwave,igas) = QextsSOL(iz,iwave,igas) +
C     &           qext(iwave,irad,igas)*pir2ndz !JDW+TDR
C                QscasSOL(iz,iwave,igas) = QscasSOL(iz,iwave,igas) +
C     &           qscat(iwave,irad,igas)*pir2ndz !JDW+TDR
C                endif
                  norms(iwave)=norms(iwave) + pir2ndz !JDW + TDR
              enddo

            enddo   ! irad:1,nrad
C
C JDW + TDR, compute number density-weighted mean Qexts and Qscats
C

            do iwave=1,nwave
C              IF(do_IR_wav)then
              Qexts(iz,iwave,igas)=Qexts(iz,iwave,igas)/norms(iwave)
              Qscas(iz,iwave,igas)=Qscas(iz,iwave,igas)/norms(iwave)
C              else
C             QextsSOL(iz,iwave,igas)=QextsSOL(iz,iwave,igas)/norms(iwave)
C              QscasSOL(iz,iwave,igas)=QscasSOL(iz,iwave,igas)/norms(iwave)
C              endif
            enddo 




c   index of bottom of cloud layer for subcloud kludge

            if( do_subcloud )then
              ibot_cloud(igas) = iz
            endif

          endif     ! ndz > 0
        enddo       ! igas:1,ngas
      enddo         ! iz:itop,ibot,incr

c   --------------------------------------------------------------------
c   subcloud kludge to soften discontinuity at cloud base
c   (10% in first layer below, 5% in second layer)

      if( do_subcloud )then
        do igas = 1, ngas
          ibot_gas = ibot_cloud(igas)
          if( ibot_gas .ne. ibot )then


c    Choose lower index to be within grid

            if( ibot_gas+incr .eq. ibot )then
              ibot_subcloud = ibot_gas + incr
            else
              ibot_subcloud = ibot_gas + 2*incr
            endif


            do iz = ibot_gas+incr, ibot_subcloud, incr

              if( iz .eq. ibot_gas+incr )then
                norm = 0.10
              else
                norm = 0.05
              endif




C
C JDW + TDR for computing number density-weighted mean Qexts and Qscas
C
              do iwave=1,nwave
                norms(iwave)=0.
              enddo 

              opd_layer(iz,igas) = opd_layer(iz,igas) + 
     $          opd_layer(ibot_gas,igas)*norm

              do iwave = 1, nwave
                scat_gas(iz,iwave,igas) = scat_gas(iz,iwave,igas) +  !Likely there is a bug in scat_gas+ due to initialization with waves diffent JDW 2021
     $            scat_gas(ibot_gas,iwave,igas)*norm
                ext_gas(iz,iwave,igas) = ext_gas(iz,iwave,igas) + 
     $            ext_gas(ibot_gas,iwave,igas)*norm
                cqs_gas(iz,iwave,igas) = cqs_gas(iz,iwave,igas) + 
     $            cqs_gas(ibot_gas,iwave,igas)*norm
C                IF(do_IR_wav)then
                Qexts(iz,iwave,igas)=Qexts(iz,iwave,igas) + 
     &           qext(iwave,irad,igas)*norm !JDW+TDR
                Qscas(iz,iwave,igas)=Qscas(iz,iwave,igas) + 
     &           qscat(iwave,irad,igas)*norm !JDW+TDR
                norms(iwave)=norms(iwave) + norm !JDW + TDR
C                else  

C                QextsSOL(iz,iwave,igas)=QextsSOL(iz,iwave,igas) + 
C     &           qext(iwave,irad,igas)*norm !JDW+TDR
C                QscasSOL(iz,iwave,igas)=QscasSOL(iz,iwave,igas) + 
C     &           qscat(iwave,irad,igas)*norm !JDW+TDR
                !norms(iwave)=norms(iwave) + norm !JDW + TDR
C                endif


              enddo

C
C JDW + TDR compute number density-weighted mean Qexts and Qscas
C
              do iwave=1,nwave
C                IF(do_IR_wav)then
                Qexts(iz,iwave,igas)=Qexts(iz,iwave,igas)/norms(iwave)
                Qscas(iz,iwave,igas)=Qscas(iz,iwave,igas)/norms(iwave)
C                else 

C                QextsSOL(iz,iwave,igas)=QextsSOL(iz,iwave,igas)/norms(iwave)
C                QscasSOL(iz,iwave,igas)=QscasSOL(iz,iwave,igas)/norms(iwave)
C                endif
              enddo

            enddo
          endif
        enddo 
      endif 

c   --------------------------------------------------------------------
c   Sum over gases and compute spectral optical depth profile etc

      do iz = itop, ibot, incr
        do iwave = 1, nwave

          opd_scat = 0.
          opd_ext = 0.
          cos_qs = 0.

          do igas = 1, ngas
            opd_scat = opd_scat + scat_gas(iz,iwave,igas)
            opd_ext = opd_ext + ext_gas(iz,iwave,igas)
            cos_qs = cos_qs + cqs_gas(iz,iwave,igas)
            if( scat_gas(iz,iwave,igas) .gt. 0. )then !JDW+TDR
C              IF(do_IR_wav)then
              g0s(iz,iwave,igas) = cqs_gas(iz,iwave,igas)
     &        /scat_gas(iz,iwave,igas) !JDW+TDR
C              else
C              g0sSOL(iz,iwave,igas) = cqs_gas(iz,iwave,igas)
C     &        /scat_gas(iz,iwave,igas) !JDW+TDR
C              endif
            else !JDW+TDR
C              IF(do_IR_wav)then
              g0s(iz,iwave,igas) = 0. !JDW +TDR
C              else
C              g0sSOL(iz,iwave,igas) = 0. !JDW +TDR
C             endif
              
            endif !JDW+TDR
C            IF(do_IR_wav)then
            dtauexts(iz,iwave,igas) = ext_gas(iz,iwave,igas) !JDW+TDR
C            else
C            dtauextsSOL(iz,iwave,igas) = ext_gas(iz,iwave,igas) !JDW+TDR
C            endif
c			if (igas.eq.3) print *,iz,iwave,wlni(iwave),ext_gas(iz,iwave,igas)
          enddo

          if( opd_scat .gt. 0. )then
            
            opd(iz,iwave) = opd_ext
            w0(iz,iwave) = opd_scat / opd_ext
            g0(iz,iwave) = cos_qs / opd_scat
          else 
            opd(iz,iwave) = 0.
            w0(iz,iwave) = 0.
            g0(iz,iwave) = 0.
          endif
          !print*,opd(iz,1),iz
        enddo
      enddo


c   cumulative optical depths for conservative geometric scatterers

      opd_tot = 0.
      do igas = 1, ngas
        opd_gas(itop,igas) = opd_layer(itop,igas)
        do iz = itop+incr, ibot, incr
          opd_gas(iz,igas) = opd_gas(iz-incr,igas) + opd_layer(iz,igas)
        enddo
        opd_tot = opd_tot + opd_gas(ibot,igas)
      enddo
c		stop
      !print*,opd





C      IF(do_IR_wav)then 
C      dtauextsIR=dtauexts
C      g0sIR=g0s
C      QextsIR=Qexts
C      QscasIR=Qscas
C      dtauextsSOL=dtauexts
C      g0sSOL=g0s
C      QextsSOL=Qexts
C      QscasSOL=Qscas
      dtauextsSOL=dtauexts
      g0sSOL=g0s !JDW Hacked 2021
      QextsSOL=Qexts
      QscasSOL=Qscas

      endif 
      return
      end








