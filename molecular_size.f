      subroutine molecular_size(print_size)
C ==============================================================
C
C   Based on tables from Rosner 2000, Chapter 3. 
C
C ==============================================================
      INCLUDE 'CLIMA/INCLUDE/header.inc'
      INCLUDE 'globals.h'
      INCLUDE 'prog_params' !included from EddySed JDW
      PARAMETER(NS=3, NS1=NS+2, NS4=NS+5) !gna: changed NS1 from NS+1 to NS+2 to add ethane
      double precision cumul_molec_size
      logical print_size
      double precision mix(10,nd)
      COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4,FC2H6,FNO2, FI(NS1,ND),FH22
     & ,FI_cloudy(NS1,ND)
      COMMON/ROSNER/cumul_molec_viscosity(MAXNZ),cumul_molec_size(MAXNZ)

C         READ(114,*) FAR                  !Argon
C         READ(114,*) FCH4                 !Methane
C         READ(114,*) FC2H6                !Ethane        
C         READ(114,*) FCO2                 !Carbon dioxide
C         READ(114,*) FN2                  !Nitrogen - added Nitrogen mixing ratio c-rr 6/5/2012        
C         READ(114,*) FO2                  !Oxygen        
C         READ(114,*) FH22                 ! c-rr 5/29/2012 added H2 mixing ratio
C         READ(114,*) FNO2                 !Nitrogen dioxide
c C   FI = SPECIES MIXING RATIOS   1 = water, 2 = co2, 3 = ch4, 4 = o3, 5 = ethane
      do zz=1,nd
            mix(1,zz)=FI(1,zz) !water 1 !Numbered according to HITRAN
            !print*,'FI(1,zz)',FI(1,zz)
            mix(2,zz)=FI(2,zz) !CO2 2
            mix(3,zz)=FI(4,zz) !O3  3   !Need O3 molecular size data
            mix(4,zz)=FI(3,zz) !CH4 4
            mix(5,zz)=FO2  !O2      5
            mix(6,zz)=FNO2!Nitrogen Dioxide 6   !Need NO2 molecular size data
            mix(7,zz)=FN2 - FI(2,zz)-FI(3,zz)-FI(4,zz)
     & - FO2 - FNO2 -FI(5,zz)-FH22 - FAR  !Nitrogen 7
            mix(8,zz)=FI(5,zz) !C2H6 8
            mix(9,zz)=FH22 !H2      9
            mix(10,zz)=FAR  !Argon    10




      enddo

C Load in the data. sigma, units are in angstrom
      data Ar_s,He_s,CH4_s,CO_s,CO2_s,C2H6_s,H2O_s,NH3_s,N2_s,O2_s,H2_s
     &/3.542,2.551,3.758,3.690,3.941,4.443,2.641,2.900,3.798,3.467,2.827/


      do i=1,nd !Skipped O3   !Skipped NO2 
      cumul_molec_size(i)=(mix(1,i)*H2O_s) + (mix(2,i)*CO2_s) + 
     & (mix(4,i)*CH4_s) + (mix(5,i)*O2_s) + 
     & ((mix(7,i)+FI(4,i)+FNO2)*N2_s)!Backfill for N2
     & + (mix(8,i)*C2H6_s) + (mix(9,i)*H2_s) + (mix(10,i)*Ar_s)
      enddo

      !if(print_size) print*,'Within sibroutine',cumul_molec_size




      return
      end