      FUNCTION RELHUM(P)
         
C
C   THIS FUNCTION CALCULATES THE RELATIVE HUMIDITY (RELHUM) AT A
C   GIVEN PRESSURE P.  IT IS CURRENTLY SET UP TO YIELD EITHER A
C   STANDARD MANABE/WETHERALD RH PROFILE OR A FULLY-SATURATED
C   ATMOSPHERE, DEPENDING UPON THE VALUE OF IMW.
C       
      INCLUDE 'CLIMA/INCLUDE/header.inc'
      INCLUDE 'globals.h'
      COMMON/EBLOK/PG,TG,PG0,IMW,RSURF,OMEGA,POCEAN,IMOIST,
     2  BETA1,BETA2,FVDRY,PDRY
      common/relhum_eddy/relhum_eddy_july_new !JDW
      common/cinputs/doEddy,doCloud,fcloud,kz_min,Crainf,Csig
     &,supsat,nsub_max,cld_hum1,clr_hum1,new_relhum(ND) !JDW Changed for separate atmospheric columns
      !COMMON/STEPS/NST

      COMMON/STEPS/NST
      common/jrelhum/ndz_jeddy(MAXNGAS),p_ndz_jeddy(MAXNGAS)
     &, qt_jeddy(MAXNGAS)
     &,p_relhum,eddy_relhum

      Q = P/PG
      Q2 = AMAX1(Q-0.02,1.E-10)

      RELHUMMW = Q2/0.98

      OMGA = 1  ! For normal Manabe-Wetherald
!--------------------------------------------Using Modified Cess et al. 1976 expression to calculation FH2O for all layers
       IF (IMW.eq.1)THEN
       CALL SATRAT(TG,PSAT)

       FSAT = PSAT/PG
       FP = 1.66e-02


      IF(FSAT.lt.FP)THEN
! If the saturation mixing ratio is less than saturation mixing ratio for present atmosphere (288 K)
! use normal Manabe-Wetherald.
       OMGA = 1.

      ELSEIF((FSAT.ge.FP).and.(FSAT.le.0.1))THEN ! If saturation mixing ratio is between
!present saturation mixing ratio and 0.1 use Kasting '86 sigma parameterization

       OMGA=1.-(FSAT-FP)/(0.1-FP)

      ELSEIF(FSAT.gt.0.1)THEN ! If the saturation mixing ratio is gt 0.1 then set RH = 80%
! throughout entire troposphere
        OMGA = 0.

      ENDIF
!S      print*, OMGA, '= this is OMEGA!'

       ENDIF   ! Only use this logic when Manabe-Wetherald is used. c-rr 10/23/2012
!---------------------------------------------------------------------------------------


      RELHUM = RSURF*(RELHUMMW**OMGA)

!      print 55,RELHUM, OMGA, P,PG,
!     &   FSAT

55    format(1x,'RELHUM=',1pe10.3,2x,'OMGA=',e10.3,2x,'P=',e10.3,2x,
     &  'PG=',e10.3, 3x, 'FSAT=', e10.3)


      IF (IMW.EQ.0) RELHUM = 1.0
      IF (IMW.EQ.3) RELHUM = 0.5
      IF (IMW.EQ.4) RELHUM = 0.1
C-KK	added for low-O2 environments, to prevent water from
C-KK	zeroing itself out. 8% rel hum is present-day atmosphere
C-KK	at approximately 15 km (cold trap level)
      !IF (RELHUM .LT. 0.08) RELHUM = 0.08

      RETURN
      END

! C     Specified relhum value only for IMW=2, otherwise follows conventions below
!       Q = P/PG 
!       Q2 = AMAX1(Q-0.02,1.E-10)
!       RELHUM = RSURF * Q2/0.98
!       IF (NST .gt. 1) then 
!             hum1=RSURF
!             cld_hum1=RSURF
!             !print*,'changed relhum to eddy_dependent',
!       !& relhum_eddy_july_new,hum1
!       else 
!             hum1=RSURF
!             cld_hum1=RSURF
!       endif

!       !print*,'hum1',hum1,nst
!       !print*,'eddy_relhum',relhum_eddy_july_new
!       RELHUM=hum1*Q2/0.98 !JDW uopdated for separate columns

!       OMGA = 1  ! For normal Manabe-Wetherald
! !--------------------------------------------Using Modified Cess et al. 1976 expression to calculation FH2O for all layers
!        IF (IMW.eq.1)THEN
!        CALL SATRAT(TG,PSAT)

!        FSAT = PSAT/PG
!        FP = 1.66e-02


!       IF(FSAT.lt.FP)THEN
! ! If the saturation mixing ratio is less than saturation mixing ratio for present atmosphere (288 K)
! ! use normal Manabe-Wetherald.
!        OMGA = 1.

!       ELSEIF((FSAT.ge.FP).and.(FSAT.le.0.1))THEN ! If saturation mixing ratio is between
! !present saturation mixing ratio and 0.1 use Kasting '86 sigma parameterization

!        OMGA=1.-(FSAT-FP)/(0.1-FP)

!       ELSEIF(FSAT.gt.0.1)THEN ! If the saturation mixing ratio is gt 0.1 then set RH = 80%
! ! throughout entire troposphere
!         OMGA = 0.

!       ENDIF
! !S      print*, OMGA, '= this is OMEGA!'

!        ENDIF   ! Only use this logic when Manabe-Wetherald is used. c-rr 10/23/2012
! !---------------------------------------------------------------------------------------




!       IF (IMW.EQ.0)then
!              RELHUM = 1.
!              cld_hum1=1.
!       endif
!       IF (IMW.EQ.3)then 
!             RELHUM = 0.5
!             cld_hum1=0.5
!       endif

!       IF (IMW.EQ.4) then 
!             RELHUM = 0.1
!             cld_hum1 = 0.1
!       endif
! C-KK  added for low-O2 environments, to prevent water from
! C-KK  zeroing itself out. 8% rel hum is present-day atmosphere
! C-KK  at approximately 15 km (cold trap level)
!       !IF (RELHUM .LT. 0.08) RELHUM = 0.08
!       IF (IMW.EQ.5) then 
!             RELHUM = 1.e-10
!             cld_hum1 = 1.e-10
!       endif
!       END
