      subroutine calleddy(eJCOLD,laststep,clim_ALT)
C        implicit none
C =====================================================================
C     This subroutine calls Eddysed from A&M 2001 with 
C     input vectors from CLIMA. The outputs are cloud 
C     optical properties & water vapor volume mixing ratios.
C =====================================================================  

      INCLUDE 'CLIMA/INCLUDE/header.inc'
      include 'globals.h'
      PARAMETER (NSOL=38,NGS=8,NF=55, IK=8)  !c-rr Adding IK=8, number of sums for CO2 and H2O absorption coefficients 8/27/2012
      PARAMETER(NS=3, NS1=NS+2, NS4=NS+5) ! Adding parameter statement needed for FI(NS1,ND) 5/23/2011
      REAL kmatrix_solh2o, kmatrix_solco2, weights, KAPPALAYERSOL_CO2,
     &   KAPPALAYERSOL_H2O,relhum_eddy_july_new ! EWS - BETA not used



      double precision eopd_gas(NF,MAXNGAS)!
      double precision eopd_gasIR(NF,MAXNGAS)
      double precision eopd_gasSOL(NSOL,MAXNGAS)


      double precision clim_ALT(MAXNZ)

      

      double precision interp_eopd(MAXNZ,NF)
      double precision interp_ew0(MAXNZ,NF)
      double precision interp_eg0(MAXNZ,NF)

      double precision interp_eopdIR(MAXNZ,NF)
      double precision interp_ew0IR(MAXNZ,NF)
      double precision interp_eg0IR(MAXNZ,NF)

      double precision interp_eopdSOL(MAXNZ,NSOL)
      double precision interp_ew0SOL(MAXNZ,NSOL)
      double precision interp_eg0SOL(MAXNZ,NSOL)


      double precision Re_GAS,rho_atmos


      double precision tran_int_eopd(NF,MAXNZ)
      double precision tran_int_ew0(NF,MAXNZ)
      double precision tran_int_eg0(NF,MAXNZ)

      double precision tran_int_eopdIR(NF,MAXNZ)
      double precision tran_int_ew0IR(NF,MAXNZ)
      double precision tran_int_eg0IR(NF,MAXNZ)


      double precision tran_int_eopdSOL(NSOL,MAXNZ)
      double precision tran_int_ew0SOL(NSOL,MAXNZ)
      double precision tran_int_eg0SOL(NSOL,MAXNZ)

      double precision eopd_tot
      double precision eopd_totSOL
      double precision eopd_totIR

      double precision eopd(MAXNZ,NF)!
      double precision ew0(MAXNZ,NF)!
      double precision eg0(MAXNZ,NF)!

      double precision eopdIR(MAXNZ,NF)
      double precision ew0IR(MAXNZ,NF)
      double precision eg0IR(MAXNZ,NF)

      double precision eopdSOL(MAXNZ,NSOL)
      double precision ew0SOL(MAXNZ,NSOL)
      double precision eg0SOL(MAXNZ,NSOL)



      double precision eddy_psat


      logical laststep
      logical do_IR_wav !JDW 2021
      logical do_subcloud
      integer i,j


C      double precision T(ND)
      double precision pe(MAXNZ) !p eddysed
      double precision te(MAXNZ)
C      double precision z !no given dimension in eddysed
      double precision ze(MAXNZ)
      double precision z_top(MAXNZ+1)
      double precision p_top(MAXNZ+1)
      double precision t_top(MAXNZ+1)
      double precision chf(MAXNZ)
      double precision cp_eddy(MAXNZ)
      double precision gas_mmr(MAXNGAS)
      !double precision gas_mmr_cold(MAXNGAS)

      !double precision gas_mmr(MAXNGAS,MAXNZ)

      double precision gas_mw(MAXNGAS)
      double precision sig(MAXNZ,MAXNGAS)
      double precision sigIR(MAXNZ,MAXNGAS)
      double precision sigSOL(MAXNZ,MAXNGAS)


      double precision rainf(MAXNZ,MAXNGAS)
      double precision cloudf(MAXNZ)
      double precision rho_p(MAXNGAS)
      double precision grav, teff, mw_atmos, kz_min, cloudf_min
      integer nsub_max
      double precision supsat

      logical do_virtual 
      integer nz,engas

      character*10 gas_name(MAXNGAS)
      character*10 sgas_name(MAXNGAS)
      integer nsteps
      double precision new_relhum


      double precision kz(MAXNZ)
      double precision qt(MAXNZ,MAXNGAS)
      double precision qc(MAXNZ,MAXNGAS)
      double precision ndz(MAXNZ,MAXNGAS)

      double precision ndzIR(MAXNZ,MAXNGAS)
      double precision ndzSOL(MAXNZ,MAXNGAS)


      double precision ndz1(MAXNZ)
      double precision rg(MAXNZ,MAXNGAS)


      double precision rgIR(MAXNZ,MAXNGAS)
      double precision rgSOL(MAXNZ,MAXNGAS)


      double precision reff(MAXNZ,MAXNGAS)

      double precision reffIR(MAXNZ,MAXNGAS)
      double precision reffSOL(MAXNZ,MAXNGAS)


      double precision eFSATUR
      double precision ePSAT
      CHARACTER*5 :: STARR   !Changed to make STARR hold up to 5 characters
      CHARACTER*11 :: AA
      CHARACTER :: DIRINOUT*8,DIRDATA*10

      logical doEddy



      integer sngas !Added for SMART_output
      integer snwave !Added for SMART_output
      integer snrad(MAXNGAS) !Added for SMART_output
      double precision sradius(MAXNRAD,MAXNGAS) !Added for SMART_output
      double precision sdr(MAXNRAD,MAXNGAS) !Added for SMART_output
      double precision sqscat(NF,MAXNRAD,MAXNGAS) !Added for SMART_output
      double precision scos_qscat(NF,MAXNRAD,MAXNGAS) !Added for SMART_output
      double precision sqext(NF,MAXNRAD,MAXNGAS) !Added for SMART_output


C     Highres optical properties vector build
      PARAMETER (NSOL_highres=1000,NF_highres=1000)
      integer snwave_highres !Added for SMART_output
      integer snrad_highres(MAXNGAS) !Added for SMART_output
      double precision sradius_highres(MAXNRAD,MAXNGAS) !Added for SMART_output
      double precision sdr_highres(MAXNRAD,MAXNGAS) !Added for SMART_output
      double precision sqscat_highres(NF_highres,MAXNRAD,MAXNGAS) !Added for SMART_output
      double precision scos_qscat_highres(NF_highres,MAXNRAD,MAXNGAS) !Added for SMART_output
      double precision sqext_highres(NF_highres,MAXNRAD,MAXNGAS) !Added for SMART_output
      
      double precision eopdIR_highres(MAXNZ,NF_highres)
      double precision ew0IR_highres(MAXNZ,NF_highres)
      double precision eg0IR_highres(MAXNZ,NF_highres)

      double precision eopdSOL_highres(MAXNZ,NSOL_highres)
      double precision ew0SOL_highres(MAXNZ,NSOL_highres)
      double precision eg0SOL_highres(MAXNZ,NSOL_highres)

      double precision waveIR_highres(NF_highres),
     $ waveSOL_highres(NSOL_highres)

      integer sngas_highres !Added for SMART_output
      logical do_highres

      double precision eopd_gasIR_highres(NF_highres,MAXNGAS)
      double precision eopd_gasSOL_highres(NSOL_highres,MAXNGAS)

      COMMON/WAVE_highres/AV_highres(NF_highres),LAM_highres(NF_highres)
     &,W_highres(NF_highres)

      common/smart_optics_highres/sngas_highres,sgas_name_highres,
     & snwave_highres,swave_highres,snrad_highres,
     & sradius_highres,sdr_highres   !added for the clima outputs.
     & ,sqscat_highres,sqext_highres,scos_qscat_highres


      common/smart_highres/QextsIR_highres,QExtsSOL_highres,
     & QscasIR_highres,QscasSOL_highres,g0sIR_highres,g0sSOL_highres
     &,dtauextsIR_highres,dtauextsSOL_highres,Qexts_highres,
     & Qscas_highres,g0s_highres,dtauexts_highres !JDW+TDR

      common/optical_properties_smart_highres/
     & eqscatIR_highres(NF_highres,MAXNRAD,MAXNGAS),
     & ecos_qscatIR_highres(NF_highres,MAXNRAD,MAXNGAS),
     & eqextIR_highres(NF_highres,MAXNRAD,MAXNGAS),
     & eqscatSOL_highres(NSOL_highres,MAXNRAD,MAXNGAS),
     & ecos_qscatSOL_highres(NSOL_highres,MAXNRAD,MAXNGAS),
     & eqextSOL_highres(NSOL_highres,MAXNRAD,MAXNGAS)

      common /optics_highres/engas_highres,gas_name_highres,
     & nwave_highres,wave_highres,
     & nrad_highres,eradius_highres,edr_highres,eqscat_highres,
     & eqext_highres,ecos_qscat_highres

      common/eddyblok_highres/eddyopdIR_highres(NF_highres,MAXNZ), 
     & eddyw0IR_highres(NF_highres,MAXNZ),
     & eddyopdSOL_highres(NSOL_highres,MAXNZ),
     & eddyw0SOL_highres(NSOL_highres,MAXNZ),
     & eddyg0SOL_highres(NSOL_highres,MAXNZ),
     & eddyg0IR_highres(NF_highres,MAXNZ),
     & eddyqt_highres(MAXNZ,MAXNGAS),eddyqc_highres(MAXNZ,MAXNGAS)


      double precision eddyopd1(NF,MAXNZ),eddyopd2(NF,MAXNZ)
     &,eddyopd3(NF,MAXNZ),eddyg01(NF,MAXNZ),eddyg02(NF,MAXNZ),
     & eddyg03(NF,MAXNZ),eddyw01(NF,MAXNZ),eddyw02(NF,MAXNZ),
     & eddyw03(NF,MAXNZ)


      double precision eddyqt1(MAXNZ,MAXNGAS),eddyqt2(MAXNZ,MAXNGAS),
     & eddyqt3(MAXNZ,MAXNGAS),eddyqc1(MAXNZ,MAXNGAS),
     & eddyqc2(MAXNZ,MAXNGAS),eddyqc3(MAXNZ,MAXNGAS)

      integer JCOLD1,JCOLD2,JCOLD3


      double precision bwni, wnoi, dwni, wlni
      integer nwave,eJCOLD

C      double precision wave(MAXNWAVE)
      double precision waveIR(NF),waveSOL(NSOL)



      integer nrad(MAXNGAS),foundloc
      integer CO2_cloud_flag,H2O_cloud_flag,CH4_cloud_flag
      double precision eradius(MAXNRAD,MAXNGAS)
      double precision eradiusIR(MAXNRAD,MAXNGAS)
      double precision eradiusSOL(MAXNRAD,MAXNGAS)


      double precision edr(MAXNRAD,MAXNGAS)
      double precision edrIR(MAXNRAD,MAXNGAS)
      double precision edrSOL(MAXNRAD,MAXNGAS)
      double precision s_eff

      double precision epvap_co2
     
      !double precision eqscat(NF,MAXNRAD,MAXNGAS) !NF used to be MAXNWAVE



      !double precision ecos_qscat(NF,MAXNRAD,MAXNGAS) !




      !double precision eqext(NF,MAXNRAD,MAXNGAS)!!
C

      COMMON/fh2o_ed/FH2O_e(ND)
      COMMON/PRESSURE/P(ND),PLOG(ND)

      COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4,FC2H6,FNO2, FI(NS1,ND),FH22
     & ,FI_cloudy(NS1,ND) !FI(1,101)
      COMMON/IRBLK/FUPIR(ND),FDNIR(ND),SRFALBIR,OMG0AIR(NF,ND-1),
     & ASYAIR(NF,ND-1),IO3,QEXTIR(NF,ND-1),FUPIR_cloudy(ND),
     & FDNIR_cloudy(ND),FUPIR_clear(ND),FDNIR_clear(ND) !JDW

      common/tclima/tclim(ND),alt_convec(ND)
      COMMON/CONSS/C,BK,G,GNEW(ND),PI,SM,DM,DM2

      COMMON/EBLOK/PG,TG,PG0,IMW,RSURF,OMEGA,POCEAN,IMOIST,
     2  BETA1,BETA2,FVDRY,PDRY

      common/cinputs/doEddy,doCloud,fcloud,kz_min,Crainf,Csig
     &,supsat,nsub_max,cld_hum1,clr_hum1,new_relhum(ND)


      common/smart_optics/sngas,sgas_name,snwave,swave,snrad,sradius,sdr   !added for the clima outputs.
     & ,sqscat,sqext,scos_qscat
      common/smart/QextsIR,QExtsSOL,QscasIR,QscasSOL,g0sIR,g0sSOL
     &,dtauextsIR,dtauextsSOL,Qexts,Qscas,g0s,dtauexts !JDW+TDR

      double precision ndz_jeddy,p_ndz_jeddy,qt_jeddy,
     &p_sat_relhum,eddy_relhum,
     & eddy_nst

      common/jrelhum/ndz_jeddy(MAXNGAS),p_ndz_jeddy(MAXNGAS)
     &, qt_jeddy(MAXNGAS)
     &,p_relhum,eddy_relhum
      common/fc_min/fc_minf

      COMMON/WAVE/AV(NF),LAM(NF),W(NF)
      common/molec_weight/mw_atmos
      !DATA C/3.E10/

      !common/colblok/eJCOLD,upatm_mix

      common/coldtrap/gas_mmr_cold(MAXNGAS)

        COMMON/SOLARBLK/AMU0,SRFALB,OMG0A(NSOL,ND-1),
     &  ASYA(NSOL,ND-1),TAUAER(NSOL),SIGERT(NSOL),FMA(NSOL),PF(ND),
     &  ALAMBDA(NSOL),CGAS(ND,NGS),FUPSOL(ND),FDNSOL(ND),
     &  NGAS(2,NSOL),WGHT(4,2,NSOL),NPR(2,NSOL),SOLINT(NSOL),
     &  TAULAM(ND-1),ASY(ND-1),OMG0(ND-1),FMT(ND-1),QEXT(NSOL,ND-1),
     &  fdnsol_cloudy(nd),fupsol_cloudy(nd),
     &  fdnsol_clear(nd),fupsol_clear(nd),ASY_cloudy(ND-1),
     &  OMG0_cloudy(ND-1),TAULAM_cloudy(ND-1),
     &  ASY_clear(ND-1),OMG0_clear(ND-1),TAULAM_clear(ND-1),
     &  ASYEDDY(ND-1),OMG0EDDY(ND-1),TAUEDDY(ND-1),CGAS_cloudy(ND,NGS)
     &  ,relhum_eddy_july(nd) !JDW 


      common/relhum_eddy/relhum_eddy_july_new,foundloc



      common/eddyblok/eddyopdIR(NF,MAXNZ), eddyw0IR(NF,MAXNZ),
     & eddyopdSOL(NSOL,MAXNZ),eddyw0SOL(NSOL,MAXNZ),
     & eddyg0SOL(NSOL,MAXNZ),
     & eddyg0IR(NF,MAXNZ),eddyqt(MAXNZ,MAXNGAS),eddyqc(MAXNZ,MAXNGAS),
     & JCOLD





      common/eddy_cp/cp_eddy
      COMMON/STEPS/NST


      include 'prog_params' !from EddySed
      common /optics/engas,gas_name,nwave,wave,
     & nrad,eradius,edr,eqscat,eqext,ecos_qscat


      COMMON /SPECTI/ BWNI(NSPC1IT),WNOI(NSPECIT),DWNI(NSPECIT),
     &                WLNI(NSPECIT)



      

      common/optical_properties_smart/eqscatIR(NF,MAXNRAD,MAXNGAS),
     & ecos_qscatIR(NF,MAXNRAD,MAXNGAS),eqextIR(NF,MAXNRAD,MAXNGAS),
     & eqscatSOL(NSOL,MAXNRAD,MAXNGAS),
     & ecos_qscatSOL(NSOL,MAXNRAD,MAXNGAS),
     & eqextSOL(NSOL,MAXNRAD,MAXNGAS)
      COMMON/DIR/DIRINOUT,DIRDATA

      COMMON/CPHEAT/CPO2(ND),CPCO2(ND), CPN2(ND), CPH2O(ND), 
     &  CPN(ND), CPNT(ND), CPH2(ND)  ! Added CPH2 5/31/2012 c-rr JDW 2021
      COMMON/ROSNER/cumul_molec_viscosity(MAXNZ),cumul_molec_size(MAXNZ)

      common/inverse/INVERSE


C ===================================================================
C
C     Set up text files to record useful output vectors JDW
C
C ====================================================================
 
      open(unit=1995,FILE= DIRINOUT//'/cloud_atmosphere.tab')
      open(unit=1996,FILE=DIRINOUT//'/cloud_optics_ew0IR.tab')
      open(unit=1997,FILE=DIRINOUT//'/cloud_optics_eg0IR.tab')
      open(unit=1998,FILE=DIRINOUT//'/cloud_optics_eopdIR.tab')

      open(unit=19961,FILE=DIRINOUT//'/cloud_optics_ew0SOL.tab')
      open(unit=19971,FILE=DIRINOUT//'/cloud_optics_eg0SOL.tab')
      open(unit=19981,FILE=DIRINOUT//'/cloud_optics_eopdSOL.tab')

      open(unit=2996,FILE=DIRINOUT//'/cloud_optics_ew0IR_highres.tab')
      open(unit=2997,FILE=DIRINOUT//'/cloud_optics_eg0IR_highres.tab')
      open(unit=2998,FILE=DIRINOUT//'/cloud_optics_eopdIR_highres.tab')

      open(unit=29961,FILE=DIRINOUT//'/cloud_optics_ew0SOL_highres.tab')
      open(unit=29971,FILE=DIRINOUT//'/cloud_optics_eg0SOL_highres.tab')
      open(unit=29981,FILE=DIRINOUT//'/cloud_optics_eopdSOL_highres.tab')


      open(unit=19951,FILE=DIRINOUT//'/cloud_radii.tab')
      open(unit=19953,FILE=DIRINOUT//'/cloud_radii_CO2.tab')

      open(unit=19952,FILE=DIRINOUT//'/cloud_kzz.tab')


      open(unit=19991,FILE=DIRINOUT//'/CO2_cloud_flag.tab')
      open(unit=19992,FILE=DIRINOUT//'/H2O_cloud_flag.tab')
      open(unit=19993,FILE=DIRINOUT//'/CH4_cloud_flag.tab')

C ================================================================
C     Zero all Eddysed Vectors
C ================================================================
      if(nst.eq.1)then !This changes the performance of the model. JDW 2022
      eddyopdSOL=0.0 
      eddyw0SOL=0.0
      eddyg0SOL=0.0
      eddyopdIR=0.0 
      eddyw0IR=0.0
      eddyg0IR=0.0
      kz=0
      qt=0 
      qc=0
      ndz=0
      rg=0
      reff=0 
      cloudf=0
      ndz_eddy=0.0
      nsub_max=  1!MAXNZ         !Maximum number of subducting layers
      endif



C Multiply the CHF by the S-Eff updated from within the model. Just for the shortwave fluxes.
C ==============================================================================================
C       Load stuff into Eddysed inputs. 
C ==============================================================================================

      cloudf_min=0.0                    !Not sure if this should be Hardcoded. JDW

      do_virtual  = .False.                   ! I expect for terrestrial planets this should not change. 
      gas_name(1) = 'H2O'
      
      
      gas_mw(1) = 18.01528          ! This is Hardcoded for Water 

      rho_p(1)=0.9970000000000000            ! Values for Water density taken from EddySed This is for water ice, we want liquid water. JDW 2022


        grav=G

        !kz=1e5
C 
C     Build the Convective Heat Flux vector from Upwelling and Downwelling fluxes -JDW
C


      do i=1,MAXNZ
            IF(INVERSE.eq.1)then
                  IF(NST.eq.1)then
                        s_eff=1.0
                  else
                  s_eff = ABS((FDNIR(1)-FUPIR(1))/(FDNSOL(1)-FUPSOL(1)))
                  endif
                  
                  !print*,'S_eff updated =',s_eff
                  chf(i)=(ABS((FUPSOL(i+1)-FDNSOL(i+1))*s_eff)-ABS((FUPIR(i+1)-FDNIR(i+1))))!ABS((ABS((FUPSOL(i+1)-FDNSOL(i+1)))-ABS((FUPIR(i+1)-FDNIR(i+1)))))!(ABS((FUPSOL(i+1)-FDNSOL(i+1)))-ABS((FUPIR(i+1)-FDNIR(i+1))))
            else
                  s_eff=1.0
                  chf(i)=(ABS((FUPSOL(i+1)-FDNSOL(i+1)))-ABS((FUPIR(i+1)-FDNIR(i+1))))!ABS((ABS((FUPSOL(i+1)-FDNSOL(i+1)))-ABS((FUPIR(i+1)-FDNIR(i+1)))))!(ABS((FUPSOL(i+1)-FDNSOL(i+1)))-ABS((FUPIR(i+1)-FDNIR(i+1))))
            endif
            !CHF should be zero above JCOLD
            if(i.lt.10)chf(i)=0.0 !Make sure to test this further JDW 2022
            if(chf(i).lt.0.0)chf(i)=0.0
            if(chf(i).lt.fc_minf*ABS((FUPSOL(1)-FDNSOL(1))))then
                   chf(i)=fc_minf*ABS((FUPSOL(1)-FDNSOL(1)))
            endif
            if(nst.eq.1)chf(i)=1000
            !print*,'chf',chf(i),i
            
      enddo

      
      do i=1,MAXNZ
            cp_eddy(i)=CPNT(i+1)
      enddo
C
C     Build the Pressure and Temperature Layer Tops,Bottoms, and Flux points -JDW
C

      do i=1,nd

        p_top(i)=PF(i)*1e6 !Convert bar to dyne 1000000
        t_top(i)=tclim(i)
        if(tclim(i).lt. 0.0) then
            print*,tclim
            stop '"tclim" is negative' !tclim is negative. 
        endif
        if(t_top(i).lt. 0.0) stop '"t_top" is negative, stop calleddy'
      enddo 

      nz=ND-1 !(nlayer-1)
      engas=1
      gas_mmr(1)=FI(1,MAXNZ+1)*(18.01528/mw_atmos)
      gas_mmr_cold(1)=FI(1,eJCOLD-1)*(18.01528/mw_atmos)

C ================================================================================================
C   Compute stuff like is done in the Ackerman Eddysed subroutine. 
C ================================================================================================

      r_atmos = 8.3143e7 / mw_atmos
      z_top(nz+1) = 0.

      do iz=nz,1,-1
        itop=iz
        ibot=iz+1
        dlnp=log( p_top(ibot)/p_top(itop))
        pe(iz) = 0.5*( p_top(itop) + p_top(ibot))

C
        dtdlnp = (t_top(itop) -t_top(ibot) )/dlnp
        te(iz) = t_top(ibot) + log( p_top(ibot)/pe(iz))*dtdlnp
        scale_h=r_atmos * te(iz)/grav
        dz_pmid = scale_h * log( p_top(ibot)/pe(iz) )
        dz_layer = scale_h * dlnp
        ze(iz) = z_top(ibot) + dz_pmid
        z_top(iz) = z_top(ibot) + dz_layer
C        DTDP(iz)=(DLOG( TC(iz)) - DLOG( TC(iz+1)))
C     &            /(DLOG(PF(iz)) - DLOG(PF(iz+1)))
      enddo      

      sig=Csig
      rainf=Crainf
      if( .not. doEddy)then
      print*, '   ALT(km)      ','           Temperature (k)    ','     Mixing Ratio(g/g)'  
      do i=1,MAXNZ
            print*,z_top(i)/1e5,t_top(i),FI(1,i)!,!ndz(i,1)
    
      enddo
      endif


      do i=1,MAXNZ
            z_top(i) = clim_ALT(i)*1e2*1e3
      enddo

      do i=nz,1,-1
            itop=iz
            ibot=iz+1
            
            ze(iz) = 0.5*( z_top(itop) + z_top(ibot))
      enddo

c   vapor pressure of condensate (dyne/cm^2)
      !pvap = pvap_gas( gas_name, t_layer , p_layer, metallicity,mw_atmos)
      call satrat(t_top(MAXNZ),eddy_psat)
      IF(EDDY_PSAT.GT.PG) eddy_psat = POCEAN
      Re_GAS=8.3143e7
c   atmospheric density (g/cm^3)
      rho_atmos = p_top(MAXNZ) / ( Re_GAS/mw_atmos * t_top(MAXNZ) )
c   mass mixing ratio of saturated vapor (g/g)
      gas_mmr(1)=eddy_psat*1e6/pe(MAXNZ)*RELHUM(P(MAXNZ+1))*(18.01528/mw_atmos)



      call satco2(te(MAXNZ),epvap_co2)
      gas_mmr(2)=FCO2*(44.01/mw_atmos)
      ! print*,'pvap_co2',epvap_co2,gas_mmr(2)




      if(doEddy)then 
      call molecular_size(.true.)
      !print*,cumul_molec_size
      call molecular_viscosity(.true.)



      !print*,'within eddysed',cumul_molec_size
      call eddysed( 
     $  grav, kz_min, cloudf_min, nsub_max, 
     $  supsat, mw_atmos, do_virtual,
     $  nz, ze, z_top, pe, p_top, te, t_top, chf, 
     $  engas, gas_name, gas_mmr, gas_mw, rho_p, sig, rainf,
     $  kz,cp_eddy, qt, qc, ndz, rg, reff, cloudf )
      !print*,'grav after eddy',grav
      eddyqt=qt
      eddyqc=qc

      ndz_jeddy(1)=ndz(eJCOLD,1)
      qt_jeddy(1)=qt(eJCOLD,1)


      !move fi_cloudy saved timesteps to here, before optical properties are calulated
      !eddyopd3=0.0
      !eddyw03=0.0
      !eddyg03=0.0
      !do j=1,NF
      !      do i=1,MAXNZ
      !             eddyopd3(j,i)=eddyopd2(j,i)
      !             eddyw03(j,i)=eddyw02(j,i)
      !             eddyg03(i,j)=eddyg02(j,i)
      !             if(NST.eq.1)eddyg03(i,J)=eddyg01(i,j)
      !             if(NST.eq.1)eddyw03(i,J)=eddyw01(i,j)
      !             if(NST.eq.1)eddyopd3(i,J)=eddyopd1(i,j)
      !             if(NST.eq.2)eddyg03(i,J)=eddyg02(i,j)
      !             if(NST.eq.2)eddyw03(i,J)=eddyw02(i,j)
      !             if(NST.eq.2)eddyopd3(i,J)=eddyopd2(i,j)
      !       enddo
      ! enddo

      eddyqt3=0.0
      eddyqc3=0.0
      JCOLD3=0

      JCOLD1=JCOLD
      JCOLD2=JCOLD
      JCOLD3=JCOLD

      do i=1,MAXNZ
            do j=1,MAXNGAS
                  eddyqt3(i,j)=eddyqt2(i,j)
                  eddyqc3(i,j)=eddyqc2(i,j)
                  if(NST.eq.1) eddyqt3(i,j)=eddyqt(i,j)
                  if(NST.eq.1) eddyqc3(i,j)=eddyqc(i,j)
                  if(NST.eq.2) eddyqt3(i,j)=eddyqt1(i,j)
                  if(NST.eq.2) eddyqc3(i,j)=eddyqc1(i,j)
            enddo
      enddo
      JCOLD3=JCOLD2

      ! eddyopd2=0.0
      ! eddyw02=0.0
      ! eddyg02=0.0
      ! do j=1,NF
      !       do i=1,MAXNZ
      !             eddyopd2(j,i)=eddyopd1(j,i)
      !             eddyw02(j,i)=eddyw01(j,i)
      !             eddyg02(i,j)=eddyg01(j,i)
      !             if(NST.eq.1)eddyg02(i,J)=eddyg0(i,j)
      !             if(NST.eq.1)eddyw02(i,J)=eddyw0(i,j)
      !             if(NST.eq.1)eddyopd2(i,J)=eddyopd(i,j)
      !       enddo
      ! enddo

      eddyqt2=0.0
      eddyqc2=0.0
      JCOLD2=0

      do i=1,MAXNZ
            do j=1,MAXNGAS
                  eddyqt2(i,j)=eddyqt1(i,j)
                  eddyqc2(i,j)=eddyqc1(i,j)
                  if(NST.eq.1)eddyqt2(i,j)=eddyqt(i,j)
                  if(NST.eq.1)eddyqc2(i,j)=eddyqc(i,j)

            enddo
      enddo
      JCOLD2=JCOLD1  

      ! eddyopd1=0.0
      ! eddyw01=0.0
      ! eddyg01=0.0
      ! do j=1,NF
      !       do i=1,MAXNZ
      !             eddyopd1(j,i)=eddyopd(j,i)
      !             eddyw01(j,i)=eddyw0(j,i)
      !             eddyg01(i,j)=eddyg0(j,i)
      !       enddo
      ! enddo

      eddyqt1=0.0
      eddyqc1=0.0
      JCOLD1=0

      do i=1,MAXNZ
            do j=1,MAXNGAS
                  eddyqt1(i,j)=eddyqt(i,j)
                  eddyqc1(i,j)=eddyqc(i,j)

            enddo
      enddo
      JCOLD1=JCOLD 

C =====================================================================================
C Load outputs into eddyblok for CLIMA. 
      !eddyopd=1./3.*(eddyopd3+eddyopd2+eddyopd1)
      !eddyg0=1./3.*(eddyg01+eddyg02+eddyg03)
      !eddyw0=1./3.*(eddyw01+eddyw02+eddyw03)

      !eddyqc=1./3.*(eddyqc1+eddyqc2+eddyqc3)
      !eddyqt=1./3.*(eddyqt1+eddyqt2+eddyqt3)

C      print*,'relhum',p_sat_relhum,qt_jeddy(1)*p_top(101)/1e6,
C     & (qt_jeddy(1)*p_top(101)/1e6)/p_sat_relhum

C      eddy_relhum=(qt_jeddy(1)*p_top(101)/1e6)/p_sat_relhum
C      eddy_nst=NST
C      do i=eJCOLD+1,101
      !FI(1,i)=qt_jeddy(1)
C      enddo


C     ndz_jeddy/V=P/(R*T)
C ========================================================================
C Eddysed Debug print statements. 
C ========================================================================

 
C      do iz=1,101
C      print*,iz,ndz(iz,1)
C      enddo 

C      print*,
C     $  grav, teff, kz_min, cloudf_min, nsub_max, 
C     $  supsat, mw_atmos, do_virtual,
C     $  nz, ze, z_top, pe, p_top, te, t_top, chf, 
C     $  engas, gas_name, gas_mmr, gas_mw, rho_p, sig, rainf,
C     $  kz, qt, qc, ndz, rg, reff, cloudf
C      print*,engas
C      print*,1



C      print*,'grav',grav
C      print*,'teff',teff
C      print*,'kz_min',kz_min
C      print*,'cloudf_min',cloudf_min
C      print*,'nsub_max',nsub_max
C      print*,'supsat',supsat
C      print*,'mw_atmos',mw_atmos
C      print*,'do_virtual',do_virtual
C      print*,'nz',nz
C      do i=1,100,25
C            print*,i,'z=    ',ze(i)
C            print*,i,'z_top=',z_top(i)
C            print*,i,'p=    ',pe(i)
C            print*,i,'p_top=',p_top(i)
C            print*,i,'t=    ',te(i)
C            print*,i,'t_top=',t_top(i)
C      enddo

C      print*,'shz ',shape(ze)
C      print*,'shzt',shape(z_top)
C      print*,'shp ',shape(pe)
C      print*,'shpt',shape(p_top)
C      print*,'sht ',shape(te)
C      print*,'shtt',shape(t_top)

C      print*,'z_top(1) =',z_top(1)
C      print*,'z_top(101)',z_top(101)

C =============================================================================
C Calculate optical properties from EddySed outputs
C =============================================================================

C
C     Need to include logic that calculates ICE and LIQUID water optical constants in 
C     init_optics. (i.e. eqscatIR_ICE, eqscatIR_LIQUID, eqextIR_ICE, eqextIR_LIQUID, 
C     ecos_qscatIR_ICE, ecos_qscatIR_LIQUID)
C

C     INIT OPTICS FOR LONGWAVE
      !PRINT*,'WAVE',WAVE
      do_IR_wav=.true.
      call init_optics(.true., .false., engas, gas_name, 
     $  NF, waveIR, nrad, eradiusIR, edrIR, eqscatIR, eqextIR,
     & ecos_qscatIR,do_IR_wav )
      !PRINT*,'WAVE',WAVE
      do_subcloud=.false.


            
      
      eopdIR=0.0
      eopd_totIR=0.0
      call calc_optics( do_subcloud, nz, engas, NF, nrad,
     $      ze, Te, gas_name, eradiusIR, edrIR, eqscatIR, eqextIR, ecos_qscatIR, !added 'e' for EddySed stuff JDW
     $  ndz, sig, rg, eopdIR, eopd_gasIR, ew0IR, eg0IR, eopd_totIR,do_IR_wav ) !Added IR for   the IR routines 
      !eg0IR=0.0
      !print*, do_subcloud, nz, engas, NF, nrad,
      ! $      ze, gas_name, eradiusIR, edrIR, eqscatIR, eqextIR, ecos_qscatIR, !added 'e' for EddySed stuff JDW
      !  $  ndz, sig, rg, eopdIR, eopd_gasIR, ew0IR, eg0IR, eopd_totIR,do_IR_wav 
      
      do_IR_wav=.false.
      call init_optics(.true., .false., engas, gas_name, 
     $  NSOL, waveSOL, nrad, eradiusSOL, edrSOL, eqscatSOL, eqextSOL,
     $  ecos_qscatSOL,do_IR_wav )
      !PRINT*,'WAVE',WAVE
      do_subcloud=.false.
      !do i =1,NF
      !      print*,'waveSOL,WAVEIR',waveSOL(i),waveIR(i)
      !enddo
      eopdSOL=0.0
      eopd_totSOL=0.0
      call calc_optics( do_subcloud, nz, engas, NSOL, nrad,
     $      ze, Te, gas_name, eradiusSOL, edrSOL, eqscatSOL, eqextSOL,
     $  ecos_qscatSOL, !added 'e' for EddySed stuff JDW
     $  ndz, sig, rg, eopdSOL, eopd_gasSOL, 
     $  ew0SOL, eg0SOL, eopd_totSOL,do_IR_wav ) !Added IR for the IR routines !Solar ASY are bad, i.e. not mathcing up in the shorter wavelegnths JDW 2021
 

      !print*,'findloc',findloc(alt_convec,MAXVAL(alt_convec),1)
      foundloc=findloc(alt_convec,MAXVAL(alt_convec),1)!findloc(relhum_eddy_july,MAXVAL(relhum_eddy_july),1)
      !print*,MAXVAL(alt_convec)
      !print*,'relhum'!qt_eddy(cloudtop)*p_top(cloudtop)/psat_top_clouds
      !end 












      !print*,relhum_eddy_july(foundloc-1)
      do_highres = .false.
      if(do_highres)then
            if(laststep)then !Only want to calculate highres SMART outputs if 
                              !climate model is on last timestep




C =============================================================================
C Calculate highres optical properties from EddySed outputs
C =============================================================================
C
C     Need to include logic that calculates ICE and LIQUID water optical constants in 
C     init_optics. (i.e. eqscatIR_ICE, eqscatIR_LIQUID, eqextIR_ICE, eqextIR_LIQUID, 
C     ecos_qscatIR_ICE, ecos_qscatIR_LIQUID)
C

C     INIT OPTICS FOR LONGWAVE
      !PRINT*,'WAVE',WAVE
      do_IR_wav=.true.
      call init_optics_highres(.true., .false., engas, gas_name, 
     $  NF_highres, waveIR_highres, nrad, eradiusIR, 
     & edrIR, 
     & eqscatIR_highres, eqextIR_highres,
     & ecos_qscatIR_highres,do_IR_wav )
      !PRINT*,'WAVE',WAVE
      do_subcloud=.false.


            
      
      eopdIR_highres=0.0
      eopd_totIR_highres=0.0
      call calc_optics_highres( do_subcloud, nz, engas, NF_highres, nrad,
     $      ze, Te, gas_name, eradiusIR, edrIR, 
     & eqscatIR_highres, eqextIR_highres, ecos_qscatIR_highres, !added 'e' for EddySed stuff JDW
     $  ndz, sig, rg, eopdIR_highres, eopd_gasIR_highres, 
     & ew0IR_highres, eg0IR_highres, eopd_totIR_highres,do_IR_wav ) !Added IR for   the IR routines 
      !eg0IR=0.0
      !print*, do_subcloud, nz, engas, NF, nrad,
      ! $      ze, gas_name, eradiusIR, edrIR, eqscatIR, eqextIR, ecos_qscatIR, !added 'e' for EddySed stuff JDW
      !  $  ndz, sig, rg, eopdIR, eopd_gasIR, ew0IR, eg0IR, eopd_totIR,do_IR_wav 





C     INIT OPTICS FOR SHORTWAVE
      do_IR_wav=.false.
      call init_optics_highres(.true., .false., engas, gas_name, 
     $  NSOL_highres, waveSOL_highres, nrad, eradiusSOL, 
     & edrSOL, eqscatSOL_highres, eqextSOL_highres,
     $  ecos_qscatSOL_highres,do_IR_wav )
      !PRINT*,'WAVE',WAVE
      do_subcloud=.false.
      !do i =1,NF
      !      print*,'waveSOL,WAVEIR',waveSOL(i),waveIR(i)
      !enddo
      eopdSOL_highres=0.0
      eopd_totSOL_highres=0.0
      call calc_optics_highres( do_subcloud, nz, engas, NSOL_highres, nrad,
     $      ze, Te, gas_name, eradiusSOL, edrSOL, 
     & eqscatSOL_highres, eqextSOL_highres,
     $  ecos_qscatSOL_highres, !added 'e' for EddySed stuff JDW
     $  ndz, sig, rg, eopdSOL_highres, eopd_gasSOL_highres, 
     $  ew0SOL_highres, eg0SOL_highres, eopd_totSOL_highres,do_IR_wav ) !Added IR for the IR routines !Solar ASY are bad, i.e. not mathcing up in the shorter wavelegnths JDW 2021

      sngas=engas
      sgas_name=gas_name
      snwave=nwave
      snrad=nrad
      sradius=eradius
      sdr=edr 
      sqscat=eqscat 
      sqext=eqext !This should be added to the IR subroutine as well.
      scos_qscat=ecos_qscat
            endif
      endif 

C =============================================================================
C End Calculate highres optical properties from EddySed outputs
C =============================================================================















C ==============================================================================
C Print optical properties to terminal / write outputs
C ==============================================================================



      print*, '   Temperature (K)      ','           Optical Depth SOL    ','     Altitude (km)'  
      do i=1,MAXNZ
            print*,t_top(i),eopdSOL(i,20),z_top(i)/1e5
            if(isnan(ndz(i,1))) stop '"NDZ is NaN", stopping'
      enddo

      print*, '   rg(um)      ','           CLIMA MMR(g/g)    ','     reff(um)'
      do i=1,MAXNZ
            print*,rg(i,1)*10000,FI(1,i+1),reff(i,1)*10000!,!ndz(i,1) 
    
      enddo

      print*, '   qc(g/g)      ','           qt(g/g)    ','     ndz(n/cm^3)'
      do i=1,MAXNZ
            !print*,qc(i,1),(qt(i,1)-qc(i,1))*(mw_atmos/18.),ndz(i,1)/(z_top(i)-z_top(i+1))
!             print*,qc(i,2),(qt(i,1)-qc(i,2))*(mw_atmos/44.01),ndz(i,2)/((r_atmos * te(i)/grav)
!      &*log( p_top(i+1)/p_top(i)))
      print*,qc(i,1),(qt(i,1)-qc(i,1))*(mw_atmos/18.),ndz(i,1)/((r_atmos* te(i)/grav)
     & *log( p_top(i+1)/p_top(i)))
      enddo
      !print*,'  Liquid water content test EddySed','        layer'
      ! do i=1,MAXNZ
      !       print*,(((qt(i,1)-qc(i,1))*(p_top(i+1)-p_top(i))/2000.0)*1e4)/abs(ze(i+1)-ze(i))/100.0, ' g/m^3', i
      ! enddo

      !print*,' Liquid water content test CLIMA', '               layer'

      ! do i = 1, MAXNZ
      !       print*,((FI(1,i)*(18./mw_atmos)*(p_top(i+1)-p_top(i))/2000.0)*1e4)/abs(ze(i+1)-ze(i))/100.0, ' g/m^3', i
      ! enddo

      ! print*,' Layer column weight CLIMA'
      ! do i=1,MAXNZ
      ! print*,((p_top(i+1)-p_top(i))/2000.0)*1e4, ' g/m^2',abs(ze(i+1)-ze(i))/100.0
      ! enddo

      !print*, '   ew0           ','           eq0(g/g)    ','     eopd_tot(/cm)'
      !do i=1,MAXNZ
      !      print*,ew0(i,1),eg0(i,1),eopd(i,1)!,!ndz(i,1)
    
      !enddo

      !PRINT*,REFF
      !Comment a few lines of code out for the time dependant stepping

C      DO i=1,size(wave)
C            print*,'wave',wave(i)*1e4
C      enddo
C      print*,'wave shape',shape(wave)
C      print*,'lam shape',shape(w)
C      do i=1,size(w)
C            print*,'LAM',LAM(i)
C      enddo
      !print*,'LAM from AV',C/AV
      write(1995,*) ' i','  Z(km)  ','    P(bar)  ','    T(K)  '
      do i=1,MAXNZ
            write(1995,666) i,z_top(i),p_top(i),t_top(i)
      enddo




      !Write cloud layer flags for SMART
      !     2 == CO2 cloud layers
      !     1 == H2O cloud layers 
      !     3 == CH4 cloud layers 
      do i =1,maxnz
            CO2_cloud_flag = 0
            H2O_cloud_flag = 0
            CH4_cloud_flag = 0


            if( qc(i,2) .gt. 0.0)then
                  CO2_cloud_flag = 2
            endif 

            if(qc(i,1) .gt. 0.0)then
                  H2O_cloud_flag = 1
            endif

            if(qc(i,3).gt.0.0)then
                  CH4_cloud_flag = 3
            endif 

            write(19991,*) CO2_cloud_flag
            write(19992,*) H2O_cloud_flag
            write(19993,*) CH4_cloud_flag

      enddo 
      rewind(19991)
      rewind(19992)
      rewind(19993)


      write(1995,*) ' i','  qc(g/g)  ',' qt(g/g)  ',' clima_qt(g/g)  '
      write(19951,*) ' i','  rg(cm)  ',' reff(cm)  ',' ndz(n/cm**3)  '
      write(19953,*) ' i','  rg(cm)  ',' reff(cm)  ',' ndz(n/cm**3)  '
      write(19952,*) ' i','  P(bar) ','chf(erg/s/cm**2)',' kzz(cm**-2)  '
      do i=1,MAXNZ
            write(1995,666) i,qc(i,1),qt(i,1),FI(1,i+1)
            write(19951,666) i,rg(i,1),reff(i,1),ndz(i,1)/((r_atmos * te(i)/grav)
     &*log( p_top(i+1)/p_top(i)))
      write(19953,666) i,rg(i,2),reff(i,2),ndz(i,2)/((r_atmos * te(i)/grav)
     &*log( p_top(i+1)/p_top(i)))
            write(19952,666) i,pf(i+1),chf(i),kz(i)
      enddo

      !write(1995,*) ' i',' qc (g/g)  ',' ew0   ',' eq0   '
      !do i=1,101
      !      write(1995,666) i,qc(i,1),ew0(i,1),eg0(i,1)
      !enddo

      write(1995,*)
      write(1997,*)
      write(1998,*)
      write(19951,*)
      write(19953,*)
      write(19952,*)
 666  FORMAT(I3,20(1x,1PE10.4),1PE9.2,2X,1PE11.4,20(1X,1PE11.4))


      !NF,MAXNZ
      write(1996,*)
      write(19961,*)

      write(1997,*)
      write(19971,*)

      write(1998,*)
      write(19981,*)


      do i=1,MAXNZ
      !write(1996,*)i
      write(1996,'(*(F14.7))')(ew0IR(i,j),j=1,NF) ! These are all zeros JDW
      write(19961,'(*(F14.7))')(ew0SOL(i,j),j=1,NSOL)

      write(1997,'(*(F14.7))')(eg0IR(i,j),j=1,NF)
      write(19971,'(*(F14.7))')(eg0SOL(i,j),j=1,NSOL)

      write(1998,'(*(F14.7))')(eopdIR(i,j),j=1,NF)
      write(19981,'(*(F14.7))')(eopdSOL(i,j),j=1,NSOL)
      !write(2, '(*(F14.7))')( real(Vec(i,j)) ,j=1,M)
      !write(1996,*)ew0(i,:)
      
      enddo
      






 667  FORMAT(/1X,"lam#",3x,"ew0",9X,"eg0",8X,"eopd")
 668  FORMAT()
C ====================================================================================
C Interpolation attempt. Not needed because init_optics utilizes CLIMA spectral grid. 
C ====================================================================================

      !Interpolated data is over (layer,wavelength) Must be done per layer! Though need (NSOL,Layer)
C      do i=1,MAXNZ

C      call pwl_value_1d( size(eopd(i,:)),wave,eopd(i,:),
C     & size(AV),C/AV,interp_eopd(i,:))


C      call pwl_value_1d( size(ew0(i,:)),wave,ew0(i,:),
C     & size(AV),C/AV,interp_ew0(i,:))
      

C      call pwl_value_1d( size(eg0(i,:)),wave,eg0(i,:),
C     & size(AV),C/AV,interp_eg0(i,:))
      !print*,'eg0',interp_eg0

C      enddo


      !interp_eopd=eopd
      !interp_eg0=eg0
      !interp_ew0=ew0

      interp_eopdIR=eopdIR
      interp_eg0IR=eg0IR
      interp_ew0IR=ew0IR

      interp_eopdSOL=eopdSOL
      interp_eg0SOL=eg0SOL
      interp_ew0SOL=ew0SOL



      !print*,'interp_eopd',shape(interp_eopd)
      !print*,'interp_ew0 ',shape(interp_ew0)
      !print*,'interp_eg0 ',shape(interp_eg0)

C ======================================================================================
C Transpose calc_optics outputs for use in CLIMA
C ======================================================================================

      !tran_int_eopd=transpose(interp_eopd)
      !print*,'tran_int_eopd',shape(tran_int_eopd)

      !tran_int_ew0=transpose(interp_ew0)
      !print*,'tran_int_ew0 ',shape(tran_int_ew0)

      !tran_int_eg0=transpose(interp_eg0)
      !print*,'tran_int_eg0 ',shape(tran_int_eg0)



      tran_int_eopdIR=transpose(interp_eopdIR)
      tran_int_ew0IR=transpose(interp_ew0IR)
      tran_int_eg0IR=transpose(interp_eg0IR)


      tran_int_eopdSOL=transpose(interp_eopdSOL)
      tran_int_ew0SOL=transpose(interp_ew0SOL)
      tran_int_eg0SOL=transpose(interp_eg0SOL)
C Save previous three cloud timseteps

!eddyopd(NF,MAXNZ), eddyw0(NF,MAXNZ),
!     & eddyg0(NF,MAXNZ),eddyqt(MAXNZ,MAXNGAS),eddyqc(MAXNZ,MAXNGAS),
!     & JCOLD











      !eddyg0=tran_int_eg0
      !eddyw0=tran_int_ew0
      !eddyopd=tran_int_eopd

      eddyg0IR=tran_int_eg0IR !Added lines for separate SOL and IR JDW
      eddyw0IR=tran_int_ew0IR
      eddyopdIR=tran_int_eopdIR
      !print*,'shape of eddyopdIR',shape(eddyopdIR),shape(eopdIR)
      !do i=1,100
      !      print*,eddyw0IR(49,i),ew0IR(i,49),i
      !enddo

      eddyg0SOL=tran_int_eg0SOL
      eddyw0SOL=tran_int_ew0SOL
      eddyopdSOL=tran_int_eopdSOL











      !do i=1,MAXNZ
      !      do j=1,NF

      !print*,'eddyg0',pack(eddyg0(j,i),eddyg0(j,i) <0.0)
      !      enddo
      !enddo
      !print*,'begin pack'
      !print*,pack(eddyg0,eddyg0 <0.0)
      !print*,'pack =',pack(ew0,ew0 .gt. 1.0)
      !print*,'wave=',wave
      !print*,'lam=',C/AV
      ! print*, eopd
      !print*,'end pack'


      !intg0=interp_linear_internal(wave,eg0,C/AV)

      ! It would make more sense here to sum over the spectral intervals and take the average values. 



C     Read in stuff for the CLIMA-SMART outputs. JDW


      sngas=engas
      sgas_name=gas_name
      snwave=nwave
      snrad=nrad
      sradius=eradius
      sdr=edr 
      sqscat=eqscat 
      sqext=eqext !This should be added to the IR subroutine as well.
      scos_qscat=ecos_qscat







      endif
      return
      end 



C ======================================================================================
C This interpolation routine is not used. JDW

      subroutine pwl_value_1d ( nd, xd, yd, ni, xi, yi ) !Interpolation for EddySed outputs 

            !*****************************************************************************80
            !
            !! PWL_VALUE_1D evaluates the piecewise linear interpolant.
            !
            !  Discussion:
            !
            !    The piecewise linear interpolant L(ND,XD,YD)(X) is the piecewise
            !    linear function which interpolates the data (XD(I),YD(I)) for I = 1
            !    to ND.
            !
            !  Licensing:
            !
            !    This code is distributed under the GNU LGPL license.
            !
            !  Modified:
            !
            !    22 September 2012
            !
            !  Author:
            !
            !    John Burkardt
            !
            !  Parameters:
            !
            !    Input, integer ( kind = 4 ) ND, the number of data points.
            !    ND must be at least 1.
            !
            !    Input, real ( kind = 8 ) XD(ND), the data points.
            !
            !    Input, real ( kind = 8 ) YD(ND), the data values.
            !
            !    Input, integer ( kind = 4 ) NI, the number of interpolation points.
            !
            !    Input, real ( kind = 8 ) XI(NI), the interpolation points.
            !
            !    Output, real ( kind = 8 ) YI(NI), the interpolated values.
            !
              implicit none
            
              integer ( kind = 4 ) nd
              integer ( kind = 4 ) ni
            
              integer ( kind = 4 ) i
              integer ( kind = 4 ) k
              real ( kind = 8 ) t
              real ( kind = 8 ) xd(nd)
              real ( kind = 8 ) yd(nd)
              real ( kind = 8 ) xi(ni)
              real ( kind = 8 ) yi(ni)
            
              yi(1:ni) = 0.0D+00
            
              if ( nd == 1 ) then
                yi(1:ni) = yd(1)
                return
              end if
            
              do i = 1, ni
            
                if ( xi(i) <= xd(1) ) then
            
                  t = ( xi(i) - xd(1) ) / ( xd(2) - xd(1) )
                  yi(i) = ( 1.0D+00 - t ) * yd(1) + t * yd(2)
            
                else if ( xd(nd) <= xi(i) ) then
            
                  t = ( xi(i) - xd(nd-1) ) / ( xd(nd) - xd(nd-1) )
                  yi(i) = ( 1.0D+00 - t ) * yd(nd-1) + t * yd(nd)
            
                else
            
                  do k = 2, nd
            
                    if ( xd(k-1) <= xi(i) .and. xi(i) <= xd(k) ) then
            
                      t = ( xi(i) - xd(k-1) ) / ( xd(k) - xd(k-1) )
                      yi(i) = ( 1.0D+00 - t ) * yd(k-1) + t * yd(k)
                      exit
            
                    end if
            
                  end do
            
                end if
            
              end do
              
              return
            end






            function interp_linear_internal(x,y,xout) result(yout)

                  implicit none
                  integer,  parameter :: dp  = kind(1.0d0)
          
                  real(dp), intent(IN)  :: x(2), y(2), xout
                  real(dp) :: yout
                  real(dp) :: alph
          
                  if ( xout .lt. x(1) .or. xout .gt. x(2) ) then
                      write(*,*) "interp1: xout < x0 or xout > x1 !"
                      write(*,*) "xout = ",xout
                      write(*,*) "x0   = ",x(1)
                      write(*,*) "x1   = ",x(2)
                      !stop
                  end if
          
                  alph = (xout - x(1)) / (x(2) - x(1))
                  yout = y(1) + alph*(y(2) - y(1))
          
                  return
          
              end function interp_linear_internal 