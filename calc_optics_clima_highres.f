cccccccccccccccccccccccccccc
      subroutine init_optics_highres( do_optics, read_mie, ngas, gas_name, 
     $  nwave, wave, nrad, radius, dr, qscat, qext, cos_qscat,do_IR_wav )

c   --------------------------------------------------------------------
c
c   Setup up a particle size grid and calculate single-particle scattering
c   and absorption efficiencies and other parameters to be used by
c   calc_optics()
c   
c   input logical:
c     do_IR_wav    .true. do IR wavelength grid.
c     do_optics    .false. means do nothing
c     read_mie     .true. means read Mie coefficients from file
c                  'gas_name'.mieff
c
c   input scalars:
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

      include 'globals_highres.h'
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


        !print*,'made it to beginning of calc_optics'
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
          filename = 'cloud_optics/' // this_gas(1:ns) // 'IR_highres.mieff' !Make use of already calculated .mie arrays JDW
          print*,'filename',filename
          else 
          filename = 'cloud_optics/' // this_gas(1:ns) // 'SOL_highres.mieff'
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


      ! print*,"Recalculating Mie"

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
          filename = 'cloud_optics/' // this_gas(1:ns) // 'IR_highres.refrind'
          else
          filename = 'cloud_optics/' // this_gas(1:ns) // 'SOL_highres.refrind'
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
        !print*,'before write extinction and scattering'
        do igas = 1, ngas

          this_gas = gas_name(igas)
          call dblank( this_gas, ns )
          IF(do_IR_wav)then
          filename =  'cloud_optics/' // this_gas(1:ns) // 'IR_highres.mieff' !ISSUE HERE WITH MIE COEFFICIENTS NEEDED TO BE KEPT SEPARATE JDW 2021
          else
          filename =  'cloud_optics/' // this_gas(1:ns) // 'SOL_highres.mieff' 
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


      subroutine calc_optics_highres( do_subcloud, nz, ngas, nwave, nrad,
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

      include 'globals_highres.h'
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

      common/smart_highres/QextsIR_highres,QExtsSOL_highres,
     & QscasIR_highres,QscasSOL_highres,g0sIR_highres,g0sSOL_highres
     &,dtauextsIR_highres,dtauextsSOL_highres,Qexts_highres,
     & Qscas_highres,g0s_highres,dtauexts_highres !JDW+TDR

      character*(*) gas_name(MAXNGAS)
      double precision bwni,wnoI,DWNI,WLNI
      double precision z(nlevel),T(nlevel)
      double precision radius(MAXNRAD,MAXNGAS)
      double precision dr(MAXNRAD,MAXNGAS)
      double precision qscat(MAXNWAVE,MAXNRAD,MAXNGAS)
      double precision qext(MAXNWAVE,MAXNRAD,MAXNGAS) ! Probably have to convert these to "_highres" JDW 2022
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

      double precision, dimension(MAXNZ,MAXNWAVE,MAXNGAS) :: QextsIR_highres
     & ,QscasIR_highres
     & ,g0sIR_highres, dtauextsIR_highres,QextsSOL_highres,QscasSOL_highres,
     & g0sSOL_highres,dtauextsSOL_highres,Qexts_highres,
     & Qscas_highres,g0s_highres,dtauexts_highres !JDW
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
      dtauextsIR_highres=dtauexts
      g0sIR_highres=g0s !JDW Hacked
      QextsIR_highres=Qexts
      QscasIR_highres=Qscas
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
      dtauextsSOL_highres=dtauexts
      g0sSOL_highres=g0s !JDW Hacked 2021
      QextsSOL_highres=Qexts
      QscasSOL_highres=Qscas

      endif 
      return
      end












