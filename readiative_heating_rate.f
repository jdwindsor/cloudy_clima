      subroutine radiative_heating_rate(g,c_p,Fnet,pressure,altitude,rhr)
C =========================================================================================
C   This subroutine calculates the relative heating rate according
C   to equation (1) in DD Turner et al. 2018. DOI: https://doi.org/10.1175/JAMC-D-17-0252.1
C   Page(s): 953–968

C   RHR(z) = dT(z)/dt = (g/c_p(z)) * (Fnet_bot-Fnet_top)/(pressure_bot-pressure_top)|(z)



C   Input Vectors:
C   Fnet          = Net flux at layer midpoint   (ergs/s/cm^2)
C   pressure      = Pressure at flux coordinates (bar) 1bar=(1.e10)--100*1000*1000*100 g⋅cm−1⋅s−2
C   altitude      = altitude at flux coordinates (cm)
C   c_p           = specific heat constant at flux coordinates (erg/g/K)

C   Input Scalars:
C   g             = gravitational acceleration at flux coordinates (cm/s^2)
C ========================================================================================
      INCLUDE 'CLIMA/INCLUDE/header.inc'
      INCLUDE 'globals.h'
      INCLUDE 'prog_params' !included from EddySed JDW
      double precision g,c_p(nd),Fnet(nd),pressure(nd),altitude(nd),rhr(nd),
     & dFnet(nd),dpressure(nd)

      do i=1,nd
        c_p_layer=c_p(i)
        dFnet(i)=Fnet(i)-Fnet(i-1)
        dpressure(i)=(pressure(i)-pressure(i-1))*1000000 
        rhr(i)=(g/c_p_layer)*(dFnet(i)/dpressure(i))
      enddo
        rhr(1)=0.0
      return

      end