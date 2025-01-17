      subroutine GASCON(T,PF,FO2,FH2,FI,FI_cloudy,FNC,FNC_cloudy,CGAS,CGAS_cloudy,NST)
      INCLUDE 'CLIMA/INCLUDE/header.inc'
      PARAMETER(NS=3, NS1=NS+2,NGS2=8) ! NGS2 to differentiate from NGS 5/30/2012
      !gna - changed ngas2 to 8 (ethane added) and NS1 to NS+2
c      NST = NST ! EWS - note that 'NST' isn't used in the subroutine
c      T = T   ! EWS - 'T' isn't used either
C
C  NGS is the number of gases in the solar code. The order of gases in
C  this subroutine is: AIR, CH4,O2, O3, CO2, H2O, and H2. 

      REAL CON
      DIMENSION T(ND),PF(ND),FI(NS1,ND),F(NGS2,ND),FF(NGS2,ND),
     2  CGAS(ND,NGS2),FNC(ND)  ! Added noncondensible mixing fraction array c-rr 5/14/2011 
      !DIMENSION T(ND),PF(ND),F(NGS2,ND),FF(NGS2,ND),
      !2  CGAS(ND,NGS2),FNC(ND)  ! Added noncondensible mixing fraction array c-rr 5/14/2011 
      Dimension FI_cloudy(NS1,ND),F_cloudy(NGS2,ND),FF_cloudy(NGS2,ND), !JDW
     & CGAS_cloudy(ND,NGS2),FNC_cloudy(ND)

      !Dimension F_cloudy(NGS2,ND),FF_cloudy(NGS2,ND), !JDW
      ! & CGAS_cloudy(ND,NGS2),FNC_cloudy(ND)
       COMMON/CONSS/C,BK,G,GNEW(ND),PI,SM,DM,DM2
      !common/cfi/cFI(NS1,ND)
       !COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4,FC2H6,FNO2, FI(NS1,ND),FH22
      !& ,FI_cloudy(NS1,ND)

C  
      NST = NST ! EWS - note that these are dummy arguments here
      T = T   
      ND1 = ND - 1 !gna - number of layers = ND I think
      CON = 2.687E19 ! to convert CGAS into molecules/cm^2.
      ! CON is taken by multiplying Loshcmidt's numnber 
C


      ! FNC is the volume amount of dry air/volume amount of total air
      ! With the exception of FCO2 and FH2O, all mixing ratios inputted here are volume amounts with respect to dry air (noncondensible species).
      ! So to get the mixing ratios with respect to total air---- (volume species/volume dry air)*(volume dry air/volume of total air) = (volume species/volume total air) c-rr 6/5/2012


     
c These are the two condensibles, water and CO2(gases 1 and 2 in FI). c-rr 5/14/2011
!      DO 1 I=1,2
!      DO 1 J=1,ND
!   1  F(I+1,J) = FI(I,J)
      
      
C

      DO 5 J = 1,ND
!      F(4,J) = FI(4,J)*FNC(J)   ! Ethane mixing ratio obtained by multiplying FI(4,J) by noncondensible mixing ratio c-rr 5/14/2011.. 8/31/2011 removed this
!   5  F(5,J) = FI(3,J)*FNC(J)   ! Methane mixing ratio obtained by multiplyling FI(3,J) by noncondensible mixing ratio c-rr 5/14/2011 8/31/2011 removed this
       F(1,J) = 1.  ! AIR 
       F(2,J) = FI(3,J)*FNC(J) ! methane
!       F(3,J) = FO2*(1. - FI(1,J)) !oxygen
       F(3,J) = FO2*FNC(J) ! oxygen
       F(4,J) = FI(4,J)*FNC(J)  ! ozone
       F(5,J) = FI(2,J) ! CO2
       !F(6,J) = cFI(1,J) ! H2O ! This is where I should incorperate the EddySed mixing ratio? JDW
       F(6,J) = FI(1,J)
!       print *, 'FH2O=', F(6,J),J
!       pause
       F(7,J) = FH2*FNC(J) ! hydrogen	
       F(8,J) = FI(5,J)*FNC(J) ! ethane
  5    CONTINUE
  


C JDW Added below for the Cloudy Column with EddySed.
      DO 100 J=1,ND
      F_cloudy(1,J) = 1.
      F_cloudy(2,J) = FI_cloudy(3,J)*FNC_cloudy(J) ! methane
      F_cloudy(3,J) = FO2*FNC_cloudy(J) ! oxygen
      F_cloudy(4,J) = FI_cloudy(4,J)*FNC_cloudy(J)  ! ozone
      F_cloudy(5,J) = FI_cloudy(2,J) ! CO2

      F_cloudy(6,J) = FI_cloudy(1,J)

      F_cloudy(7,J) = FH2*FNC_cloudy(J) !Hydrogen
      F_cloudy(8,J) = FI_cloudy(5,J)*FNC_cloudy(J) !ethane
  100 continue


      

      DO 3 I=1,NGS2
      FF(I,1) = F(I,1)
      FF(I,ND) = F(I,ND)

      FF_cloudy(I,1) = F_cloudy(I,1) !JDW added for cloudy column. 
      FF_cloudy(I,ND) = F_cloudy(I,ND)



      DO 3 J=2,ND1
      FF(I,J) = SQRT(F(I,J)*F(I,J-1))  ! FF is sqrt of F for all layers except for layer 1 and layer ND

      FF_cloudy(I,J) = SQRT(F_cloudy(I,J)*F_cloudy(I,J-1))

  3   CONTINUE
C

      DO 4 I=1,NGS2   ! 8 species
      
      DO 4 J=1,ND1   ! 100 layers
        BKMG = BK*273.15/(SM*GNEW(J)) ! scale height. BK, boltzman's constant. SM is the mass of a hydrogen atom (1.67e-24 g)
      FACT = CON*BKMG
      CGAS(J,I) = FACT * (PF(J+1) - PF(J))*1.013*(FF(I,J) + 4.*F(I,J)  ! Pressure is in bars converted to atm cuz Loschmidt # units are in molecules/(cm^3*atm). 
     2  + FF(I,J+1))/(6.*DM) 

      CGAS_cloudy(J,I)=FACT * (PF(J+1) - PF(J))*1.013*(FF_cloudy(I,J)           !Changed in 2021 JDW
     & + 4.*F_cloudy(I,J)  ! Pressure is in bars converted to atm cuz Loschmidt # units are in molecules/(cm^3*atm). 
     &  + FF_cloudy(I,J+1))/(6.*DM) 
      !do K=1,size(CGAS)
      !print*,'cgas', CGAS(1,k),CGAS_cloudy(1,k)
      !enddo
!     Use Simpson's rule for more accuracy on mixing ratios    
       
  4    CONTINUE


!          print *, CGAS(1,6), PF(51), PF(50)
!          pause

          !ASUMCGAS = 0.

           
!       print *, CGAS(ND,1)
       

!      DO J = 1,ND
 !            print 200, CGAS(J,6), CGAS(J,5)
 !     enddo
 !200  format(1p2e12.5)
!      stop

c 100  FORMAT(1X,1P10E12.5)
 
      RETURN
      END
