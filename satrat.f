      SUBROUTINE SATRAT(T,PSAT)
      
C
C   THIS SUBROUTINE FINDS THE SATURATION VAPOR PRESSURE (PSAT) AT
C   A GIVEN TEMPERATURE T.
C
C     below 273.15 K this returns the Clausius Clapeyron relation for pvap
C     above 273.15 K this uses data from Haar et al. 1984      JDW 2022
C
C
C
      PARAMETER(NT=76)
      COMMON/FBLOK/TTAB(NT),PVAP(NT),DPVAP(NT),SVC(NT),DSV(NT),DSC(NT)
     2  ,RHOV(NT),DRHOV(NT),BETAM(70,75),TCP(75),PCP(70),DPDTL(70,75),
     3  DRDTL(70,75)
      COMMON/SBLOK/P0P,T0P,R,SUBL

      RV = R/18. !very specific to water 
      print*,'subl',SUBL
      print*,'R',R
      print*,'t0p',T0P
      print*,'P0P',p0p
      print*,'i','TTAB(i)','PVAP(i)'
      ! do i=1,NT
      !      print*,i,TTAB(i),PVAP(i)
      ! enddo

      !do i=1,NT
      !      print*,PVAP(i)
      !enddo
C
      IF(T.GT.646.96) GO TO 3
      IF(T.GT.273.16) GO TO 2 !273.16
      HL = SUBL
      PSAT = P0P * EXP(-HL/RV * (1./T - 1./T0P))
      ! What equation is this?
      RETURN
C
   2  TC = T - 273.15
      N = TC/5. + 1
      !print*,'N',N
      FR = (TC - TTAB(N))/(TTAB(N+1) - TTAB(N))
      PSAT = FR*PVAP(N+1) + (1.-FR)*PVAP(N)
      PSAT=PVAP(N)
      !print*,'FR',FR
      RETURN
C
   3  PSAT = 1.E38 !Someone kludging in the critical point.
      RETURN
      END
