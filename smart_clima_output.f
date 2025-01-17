        subroutine smart_clima_output(nlev,p,z,t,rmix,lam,outputform,do_clouds,do_IR_wav)
        include 'globals.h'
        include 'run/freedman/params_read'
        include 'prog_params'

C Need to fix the N2 backfilling here. Instead of setting it to 100% N2 for the Smart output. JDW 2021
C ********************** BUG FOR S/E RUNS ************************** SEE COMMENT **********************
C-JDW Note that CLIMA/mixing_ratios.dat treats the condensible (CO2) and
C     non-consibles mixing ratios differently
C     e.g., if the atmosphere is 99% N2 and 1% CO2, then the CO2 fraction is
C     0.01 and N2 should be set to 1, because it is 100% of noncondensibles.
C     In practice, N2 should be = (1 - [everything but CO2]).


        common /optics/engas,gas_name,nwave,wave,nrad,radius,dr,qscat,qext,cos_qscat
C        common /smart/Qexts,Qscas,g0s,dtauexts
      
      PARAMETER (NSOL=38,NGS=8,NF=55, IK=8)  !c-rr Adding IK=8, number of sums for CO2 and H2O absorption coefficients 8/27/2012

      CHARACTER :: DIRINOUT*8
      !common/smart/Qexts,Qscas,g0s,dtauexts !Created in Calc_optics
      common/smart/QextsIR,QextsSOL,QscasIR,QscasSOL,g0sIR,g0sSOL
     &,dtauextsIR,dtauextsSOL,Qexts,Qscas,g0s,dtauexts !JDW+TDR


      
      common/optical_properties_smart/eqscatIR(NF,MAXNRAD,MAXNGAS),
     & ecos_qscatIR(NF,MAXNRAD,MAXNGAS),eqextIR(NF,MAXNRAD,MAXNGAS),
     & eqscatSOL(NSOL,MAXNRAD,MAXNGAS),
     & ecos_qscatSOL(NSOL,MAXNRAD,MAXNGAS),
     & eqextSOL(NSOL,MAXNRAD,MAXNGAS)
      common/smart_optics/sngas,sgas_name,snwave,swave,snrad,sradius,sdr   !added for the clima outputs.
     & ,sqscat,sqext,scos_qscat
      COMMON/DIR/DIRINOUT,DIRDATA
      !COMMON/WAVE/AV(NF),LAM(NF),W(NF)
C      COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4,FC2H6,FNO2, FI(NS1,ND),FH22
C      & ,FI_cloudy(NS1,ND)

C     note rmix must be created and follow the output rmix(1,j) where j is the  level and i is the species identity. 

      logical, intent(in) :: do_clouds,do_IR_wav
      character(len=200), intent(in) :: outputform
      character(len=200) :: fileout
      character*10 sgas_name(MAXNGAS)
      integer, intent(in) :: nlev
      integer :: i,j,jj,k,n,ngas,ioatm
      integer, dimension(MAXNGAS) :: iomie,iotau
      double precision, dimension(nlev) :: pbar
      double precision, dimension(nlev), intent(in) :: p,z,t
      double precision, dimension(10,nlev), intent(in) :: rmix
      double precision, dimension(NF), intent(in) :: lam
      double precision, dimension(MAXNZ,MAXNWAVE,MAXNGAS) :: QextsIR
     & ,QscasIR
     & ,g0sIR, dtauextsIR,QextsSOL,QscasSOL,g0sSOL,dtauextsSOL,Qexts,
     & Qscas,g0s,dtauexts !JDW
      double precision, dimension(3*MAXNZ) :: optout

        
c
c     Set up important atmospheric output files. Must be different from EGP. 
c
C TRIM is not working. This is strange jdw
      ngas = engas
      ioatm = 633
      fileout= TRIM(outputform)//'_smart.atm'


      !Change the cloud optical depths to be 10 for greater than 10. 
      open(unit=ioatm,file=DIRINOUT//'/'//fileout,status='replace',action='write')
      !print*,do_clouds,sgas_name
      if (do_clouds) then
        !print*,'doing clouds',fileout
        do n=1,2
          iotau(n) = ioatm + n +100
          iomie(n) = ioatm + MAXNGAS + n +100
          fileout= TRIM(outputform)//'_smart_'//TRIM(sgas_name(n))//'.tau'
          if(do_IR_wav)then
         open(unit=iotau(n),file=DIRINOUT//'/'//fileout,status='replace'
     &,action='write')
          else 
            open(unit=iotau(n),file=DIRINOUT//'/'//fileout,status='old',action='write')
          endif
          if(do_IR_wav)then
          fileout= TRIM(outputform)//'_smart_'//TRIM(sgas_name(n))//'_IR.mie'
          else 
          fileout=TRIM(outputform)//'_smart_'//TRIM(sgas_name(n))//'_SOL.mie'
          endif
          open(unit=iomie(n),file=DIRINOUT//'/'//fileout,status='replace',action='write')
        end do
      endif

c
c** write atmospheric structure
c C   FI = SPECIES MIXING RATIOS   1 = water, 2 = co2, 3 = ch4, 4 = o3, 5 = ethane


Csrmix(1,zz)=FI(1,zz) !water 1 !Numbered according to HITRAN
Cprint*,'FI(1,zz)',FI(1,zz)
Csrmix(2,zz)=FI(2,zz) !CO2 2
Csrmix(3,zz)=FI(3,zz) !O3  3
Csrmix(4,zz)=FI(4,zz) !CH4 4
Csrmix(5,zz)=FO2  !O2      5
Csrmix(6,zz)=FNO2!Nitrogen Dioxide 6
Csrmix(7,zz)=FN2  !Nitrogen 7
Csrmix(8,zz)=FI(5,zz) !C2H6 8
Csrmix(9,zz)=FH22 !H2      9
Csrmix(10,zz)=FAR  !Argon    10

      !break here??'
      !print*,nlevel,nlev
      write(ioatm,100)
     
      do j=1,nlev
        ! print*,rmix(1,j),j,'H2O'
        write(ioatm,101) p(j)*1.e5,z(j),t(j),(rmix(i,j),i=1,10)
      end do
100   format('      p          z        t      volume mixing ratios'/
     &       '     (pa)       (km)     (K)     '
     &       'H2O         CO2         O3          CH4          O2        '
     &       'NO2         N2          C2H6        H2           Ar')
101   format(1X,1P1E11.3,1X,0PF8.2,1X,0PF8.2,*(1X,1P1E11.3))

      !do zz=1,55
        !print*,eqext(zz,1,1),eqscat(zz,1,1),lam(zz)
      !enddo 
      !print*,shape(eqext),shape(eqscat)

c
c** write cloud files, if doing clouds
c   wl, qext, qsca, and g1
c
      IF(do_IR_wav)then
      if (do_clouds) then
        do j=1,nlev-1
          pbar(j) = sqrt(p(j)*p(j+1)) ! layer mid-pressure
        end do
        do n=1,ngas !replaged ngas with 2. jdw
c
          write(iomie(n),159)
159       format('  wavelength  Qext        Qsca         g0 ->')
          write(iomie(n),160) (pbar(j)*1.e5,j=1,nlev-1)
160       format( '     (um)   pbar (Pa):',*(1P1E11.3,13X,'|',11X))
          do k=1,NF !changed from nspeci
            jj = 1
            do j=1,nlev-1
              optout(jj) = QextsIR(j,k,n)   !Array of zeros !Still large array of Zeros JDW 04/1
              optout(jj+1) = QscasIR(j,k,n) !Array of Zeros 
              optout(jj+2) = g0sIR(j,k,n)   !Array of Zeros
              jj = jj + 3
            enddo
            write(iomie(n),161) lam(k), (optout(jj),jj=1,3*(nlev-1))
161         format(1X,0PF8.2,*(1X,1P1E11.3))
          end do

Crm
C Cloud differential optical depths at roughly 1 um for all pressure layers
C Note: Spectral indicies of 8 and 9 correspond to 0.99 um and 1.05 um respectivly JDW
C


!           k = 49 !Change where the iwave is set up
! c
!           write(iotau(n),167)
! 167       format('     p (Pa)   dtau at 1 um ->')
!           write(iotau(n),169) (pbar(j)*1.e5,j=1,nlev-1)
! 169       format( ' pbar:  ',*(1P1E11.3,7X,'|',5X))
!           write(iotau(n),171)    (p(j)*1.e5, 0., j=1,nlev-1)
!           write(iotau(n),171) (pbar(j)*1.e5, dtauextsIR(j,k,n), j=1,nlev-1) !Bad JDW should not be dtauextsIR but instead dtauexts, K should also be a different value
!           write(iotau(n),171)  (p(j+1)*1.e5, 0., j=1,nlev-1)
! 171       format(1X,*(1X,1P1E11.3,1X,1P1E11.3))
! c
        end do
      endif
c
c** close files
c
      close(ioatm)
      if (do_clouds) then
        do n=1,ngas
          close(iotau(n))
          close(iomie(n))
        end do
      end if

      else
              if (do_clouds) then
        do j=1,nlev-1
          pbar(j) = sqrt(p(j)*p(j+1)) ! layer mid-pressure
        end do
        do n=1,ngas !replaged ngas with 2. jdw
c
          write(iomie(n),1592)
1592       format('  wavelength  Qext        Qsca         g0 ->')
          write(iomie(n),1602) (pbar(j)*1.e5,j=1,nlev-1)
1602       format( '     (um)   pbar (Pa):',*(1P1E11.3,13X,'|',11X))
          do k=1,NSOL !changed from nspeci
            jj = 1
            do j=1,nlev-1
              optout(jj) = QextsSOL(j,k,n)   !Array of zeros !Still large array of Zeros JDW 04/1
              optout(jj+1) = QscasSOL(j,k,n) !Array of Zeros 
              optout(jj+2) = g0sSOL(j,k,n)   !Array of Zeros
              jj = jj + 3
            enddo
            write(iomie(n),1612) lam(k), (optout(jj),jj=1,3*(nlev-1))
1612         format(1X,0PF8.2,*(1X,1P1E11.3))
          end do

Crm
C Cloud differential optical depths at roughly 1 um for all pressure layers
C Note: Spectral indicies of 8 and 9 correspond to 0.99 um and 1.05 um respectivly JDW
C
          k = 20
c
          write(iotau(n),1672)
1672       format('     p (Pa)   dtau at 1 um ->')
          write(iotau(n),1692) (pbar(j)*1.e5,j=1,nlev-1)
1692       format( ' pbar:  ',*(1P1E11.3,7X,'|',5X))
          write(iotau(n),1712)    (p(j)*1.e5, 0., j=1,nlev-1)
          write(iotau(n),1712) (pbar(j)*1.e5, dtauextsSOL(j,k,n), j=1,nlev-1)
          write(iotau(n),1712)  (p(j+1)*1.e5, 0., j=1,nlev-1)
1712       format(1X,*(1X,1P1E11.3,1X,1P1E11.3))
c
        end do
      endif
c
c** close files
c
      close(ioatm)
      if (do_clouds) then
        do n=1,ngas
          close(iotau(n))
          close(iomie(n))
        end do
      end if
c
      endif
      end subroutine smart_clima_output
