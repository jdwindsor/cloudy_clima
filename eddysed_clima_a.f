      double precision function advdiff( qt )
c
c   calculate divergence from advective-diffusive balance for 
c   condensate in a model layer
c
c   all units are cgs
c
c   input parameters passed through argument list:
c
c     qt              total mixing ratio of condensate + vapor (g/g)
c
c   input parameters passed through common vfall_block:
c
c     ad_qbelow       total mixing ratio of vapor in underlying layer (g/g)
c     ad_qvs          saturation mixing ratio (g/g)
c     ad_mixl         convective mixing length (cm)
c     ad_dz           layer thickness (cm) 
c     ad_rainf        rain efficiency factor 
c
c   output parameters passed through common vfall_block:
c
c     ad_qc           mixing ratio of condensed condensate (g/g)
c   
c
c   A. Ackerman Feb-2000
c
      implicit none


c   Declare common storage 

      double precision ad_qbelow, ad_qvs, ad_mixl, ad_dz, ad_rainf

      common / advdiff_block /
     $  ad_qbelow, ad_qvs, ad_mixl, ad_dz, ad_rainf



c   Declare local storage

      double precision qt, ad_qc


c   All vapor in excess of saturation condenses (supsat=0)
      ad_qc = max( 0.d0, qt - ad_qvs )
      !print*,'ad_qvs',ad_qvs
c   Difference from advective-diffusive balance 
      advdiff = 
     $   ad_qbelow*exp( - ad_rainf*ad_qc*ad_dz / ( qt*ad_mixl ) ) - qt

      return
      end

      subroutine calc_qc( gas_name, rainf_layer, rho_p, mw_cloud, 
     $           q_below, supsat, w_convect, mixl,
     $           dz_layer, grav, mw_atmos, mfp, visc, t_layer, p_layer,
     $           sig_layer, qc_layer, qt_layer, rg_layer, reff_layer, 
     $           ndz_layer, qt_top, status_r, status_q )

c   --------------------------------------------------------------------
c
c   Calculate condensate optical depth and effective radius for a layer,
c   assuming geometric scatterers.
c
c   input:
c
c     gas_name    name of condensing vapor
c     rainf_layer rain factor for layer
c     rho_p       density of condensed vapor (g/cm^3)
c     mw_cloud    molecular weight of condensing vapor (g/mol)
c     q_below     total mixing ratio (vapor+condensate) below layer (g/g)
c     supsat      fractional supersaturation persisting after condensation
c     w_convect   convective velocity scale (cm/s)
c     mixl        convective mixing length scale (cm)
c     dz_layer    thickness of layer (cm)
c     grav        gravitational acceleration (cm/s^2)
c     mw_atmos    molecular weight of atmosphere (g/mol)
c     mfp         atmospheric mean free path (cm)
c     visc        atmospheric dynamic viscosity (dyne s/cm^2)
c     t_layer     temperature of layer mid-pt (K)
c     p_layer     air pressure (dyne/cm^2)
c     sig_layer   geometric std. dev. of lognormal size distribution
c
c   output:
c
c     qc_layer    condensate mixing ratio (g/g)
c     qt_layer    gas + condensate mixing ratio (g/g)
c     rg_layer    geometric mean radius of condensate (cm)
c     reff_layer  effective (area-weighted) radius of condensate (cm)
c     ndz_layer   column of particle concentration in layer (#/cm^2)
c     qt_top      top of layer
c     status_r    error status for finding rw
c     status_q    error status for finding qt
c   
c
c   A. Ackerman Feb-2001
c
c   --------------------------------------------------------------------

      !implicit none


c   Include common data shared with eddysed()

      include 'globals.h'

      common /cvm/ readin, metallicity



c   Declare externals
      double precision metallicity
      logical readin
    
      double precision pvap_gas, vfall, advdiff
      external pvap_gas, vfall, advdiff
     

c   Declare common storage for find_root( vfall )

      double precision vf_grav, vf_mw_atmos, vf_mfp, vf_visc
      double precision vf_t, vf_p, vf_rhop

      common / vfall_block /
     $  vf_grav, vf_mw_atmos, vf_mfp, vf_visc,
     $  vf_t, vf_p, vf_rhop


c   Declare common storage for find_root( advdiff )

      double precision ad_qbelow, ad_qvs, ad_mixl, ad_dz, ad_rainf

      common / advdiff_block /
     $  ad_qbelow, ad_qvs, ad_mixl, ad_dz, ad_rainf

      COMMON/EBLOK/PG,TG,PG0,IMW,RSURF,OMEGA,POCEAN,IMOIST,
     2  BETA1,BETA2,FVDRY,PDRY

      COMMON/fh2o_ed/FH2O_e(MAXNZ+1)

c   Declare local storage

      integer status_r, status_q
      character*(*) gas_name
      double precision rainf_layer, rho_p, supsat, fs, pvap, dz_layer
      double precision w_convect, qc_layer, qt_layer, rg_layer, mfp
      double precision reff_layer, ndz_layer, qvs, mw_cloud, q_below
      double precision qlo, qhi, delta_q, t_layer, p_layer, rho_atmos
      double precision rlo, rhi, delta_v, qt_top, mw_atmos, mixl, grav
      double precision rw_layer, lnsig2, sig_layer, sig_alpha, alpha
      double precision visc
      double precision psat

c   vapor pressure of condensate (dyne/cm^2)
      pvap = pvap_gas( gas_name, t_layer , p_layer, metallicity,mw_atmos)

c   saturation factor
      fs = supsat + 1

c   atmospheric density (g/cm^3)
      rho_atmos = p_layer / ( R_GAS/mw_atmos * t_layer )

c   mass mixing ratio of saturated vapor (g/g)
      ! call satrat(t_layer,psat)
      ! if(psat .gt. PG) psat=POCEAN
      !print*,'pvap',pvap-(psat*1e6)
      !qvs = fs*pvap/ ( (R_GAS/mw_cloud) * t_layer ) / rho_atmos !Commented out this way to use the CLIMA SVP for this. 
      ! qvs = psat*1e6/p_layer*(18.01528/mw_atmos)
      !q_below=qvs + 1.0e-6 !This simple hack generates clouds. So The isssue is layers not actually saturating in EddySed. JDW
      qvs = fs*pvap/ ( (R_GAS/mw_cloud) * t_layer ) / rho_atmos
      !qvs = fs*pvap/p_layer*(mw_atmos)/(p_layer*mw_atmos)
      ! qvs=q_below
      !print*,'-qvs   ','+q_below'
      !print*,t_layer,q_below,qvs
c   --------------------------------------------------------------------
c   Layer is cloud free --For condensates to form the mixing ratio must be greaterthan the saturation vapor mixing ratio.
       !This Hack fixes the cloud-condensing problem. JDW 2022 !q_below=0.9999
! Further, why does the qvs value, not equal the q_below value? When they should be computed in the same way? 
! The scheme I am currently using should place the q_below values on the saturation vapor pressure curve. 
      if( q_below .lt. qvs )then !Windsor 2022
      !print*,'q_below vs qvs',q_below,qvs !Why does qvs explode to a large value here? 44? doesn't make sense JDW 10/2022
        qt_layer = q_below
        qt_top   = q_below
        qc_layer = 0.
        rg_layer = 0.
        reff_layer = 0.
        ndz_layer = 0.

      else

c   --------------------------------------------------------------------
c   Cloudy layer: first calculate qt and qc at top of layer,
c   then calculate layer averages

c   range of mixing ratios to search (g/g)
        qhi = q_below
        qlo = qhi / 1e8 !1e3

c   precision of advective-diffusive solution (g/g)
        delta_q = q_below /1.e16 !1e3 is general

c   load parameters into advdiff common block

        ad_qbelow = q_below
        ad_qvs = qvs
        ad_mixl = mixl
        ad_dz = dz_layer
        ad_rainf = rainf_layer

c   Find total condensate mixing ratio at top of layer

        call find_root( advdiff, ZERO, qlo, qhi, delta_q, 
     $       qt_top, status_q )

c   Use trapezoid rule (for now) to calculate layer averages
c   -- should integrate exponential
        qt_layer = 0.5*( q_below + qt_top )

c   Diagnose condensate mixing ratio
        qc_layer = max( 0.d0, qt_layer - qvs )
        !print*,'qt_layer,qvs',qt_layer,qvs
c   --------------------------------------------------------------------
c   Find <rw> corresponding to <w_convect> using function vfall()

c   load parameters into vfall common block
        vf_grav = grav
        vf_mw_atmos = mw_atmos
        vf_mfp = mfp
        vf_visc = visc
        vf_p = p_layer
        vf_t = t_layer
        vf_rhop = rho_p

c   range of particle radii to search (cm)
        rlo = 1.e-10 !When rw_layer fails this is returned
        rhi = 10.

c   precision of vfall solution (cm/s)
        delta_v = w_convect / 1000.

        call find_root( vfall, w_convect, rlo, rhi, delta_v, 
     $       rw_layer, status_r )
            !rw_layer=1.500669242387765E-002
c   geometric std dev of lognormal size distribution
        lnsig2 = 0.5*log( sig_layer )**2
        !print*,lnsig2
        !print*,'sig_layer =',sig_layer
        !print*,'rw_layer',rw_layer !major problem with rw_layer in clima


c   Compute exponent in vfall = w_convect r^alpha

c   sigma floor for the purpose of alpha calculation
        sig_alpha = max( 1.1d0, sig_layer )
        !print*,'sig_layer',sig_layer
        if(isnan(w_convect)) stop '"w_convect" is NaN, stop' !w_convect is nan. JDW
        if(isnan(sig_alpha)) stop '"sig_alpha" is NaN, stop'
        
        if( rainf_layer .gt. 1 )then

c   Bulk of precip at r > rw: exponent between rw and rw*sig
          alpha = log(
     $      vfall( rw_layer*sig_alpha ) / w_convect )
     $      / log( sig_alpha )
      !if (alpha .ne. alpha)then 
      !print*,vfall( rw_layer*sig_alpha ),rw_layer !Probably an issue with rw_layer
      !      alpha=1.09
      !endif 
        else

c   Bulk of precip at r < rw: exponent between rw/sig and rw

          alpha = log(
     $      w_convect / vfall( rw_layer/sig_alpha) )
     $      / log( sig_alpha )
            !print*,'alpha2'

        endif
        !if (alpha .ne. alpha)then 
         !   print*,'fuvk1',vfall( rw_layer*sig_alpha ) / w_convect
         !   print*,'fuvk2',w_convect / vfall( rw_layer/sig_alpha)
        !endif

c   geometric mean radius of lognormal size distribution
      if(isnan(alpha)) stop '"Alpha" is NaN, stop' !Alpha is NaN JDW
      if(isnan(rainf_layer)) stop '"rainf_layer" is NaN, stop'
      if(isnan(rw_layer)) stop '"rw_layer" is NaN, stop'
        rg_layer = rainf_layer**(1./alpha) *
     $    rw_layer * exp( -(alpha+6)*lnsig2 )
      !print*,rg_layer,w_convect,alpha
        !print*,'a1',log(
Cc     $      vfall( rw_layer*sig_alpha ) / w_convect )
      !print*,'alpha =',alpha,w_convect
      ! print*,'lnsig2 =',lnsig2 !gives infinity in CLIMA
      !print*,'sig_alpha',sig_alpha sig_alpha=2.00000
c      To force a size we go to here
c      rg_layer = 2e-4

c   droplet effective radius (cm)
        reff_layer = rg_layer*exp( 5*lnsig2 )
      !print*,'rw_layer =',rw_layer*10000 !NaN in CLIMA JDW
c   column droplet number concentration (cm^-2)
        if(isnan(rho_atmos)) stop '"rho_atmos" is NaN, stop'
        if(isnan(qc_layer)) stop '"qc_layer" is NaN, stop'
        if(isnan(dz_layer)) stop '"dz_layer" is NaN, stop'
        if(isnan(rg_layer)) stop '"rg_layer" is NaN, stop' !this is where break. 
        if(isnan(lnsig2)) stop '"lnsig2" is NaN, stop'
        ndz_layer = 3*rho_atmos*qc_layer*dz_layer /
     $    ( 4*PIE*rho_p*rg_layer**3 ) * exp( -9*lnsig2 )
      !print*,rg_layer,' = rg_layer'
      endif

      return
      end

      subroutine find_root( f, y, xlow, xhigh, delta, xnew, status )
c
c   find solution to 
c
c     f(x) - y = 0 
c
c   for x between x1 and x2 to within delta using secant method
c
c   if found, return:
c     status = 0
c     xnew   = root 
c   else if maximum number of iterations reached before convergence:
c     status = 1
c     xnew   = last estimate of root
c   else (no convergence):
c     status = -1
c     xnew   = last estimate of root
c
c   A. Ackerman Feb-2000

C    call find_root( vfall, w_convect, rlo, rhi, delta_v, 
C   $       rw_layer, status_r )
c
      implicit none


c   declare externals
    
      double precision f
      external f


c   declare local storage

      double precision x1, x2, f1, f2, slope, fnew
      double precision y, xlow, xhigh, delta, xnew
      integer status, iter


c   define maximum number of iterations
      integer MAX_ITER
      parameter( MAX_ITER = 10000 )


c   define convergence criteria for independent variable
      double precision EPSILON
      parameter( EPSILON = 1d-7 )


c   copy input range into local variables
      x1 = xlow
      x2 = xhigh


c   abort if root not bracketed by initial guesses 

      f1 = f(x1) - y
      f2 = f(x2) - y

c   default return values
      xnew = x1 
      fnew = f1 

      if( f1*f2 .gt. 0 )then 
        status = -1
        return
      endif


c   iterate until root is found or maximum number of iterations are taken

      iter = 0
      do while( ( abs(fnew) .gt. delta ) .and. ( iter .lt. MAX_ITER )
     $  .and.  ( x2/x1-1. .gt. EPSILON ) )

c   estimate root from secant between endpoints

        slope = ( f2 - f1 ) / ( x2 - x1 )
        xnew = x1 - f1/slope
        fnew = f(xnew) - y

c   pick new x1 and x2 such that f(x) crosses zero between them
      
        if( fnew*f1 .le. 0 .and. x2/xnew-1. .gt. EPSILON )then
          x2 = xnew
          f2 = fnew
        else
          x1 = xnew
          f1 = fnew
        endif

        iter = iter + 1

      enddo


c   set status flag accordingly

      if( iter .eq. MAX_ITER .and. x2/x1-1. .gt. EPSILON )then 
        status = 1
      else
        status = 0
      endif

      return
      end
      subroutine find_rootl( f, y, xlow, xhigh, delta, xnew, status )
c
c   find solution to 
c
c     log( f(x) ) - log( y ) = 0 
c
c   for x between x1 and x2 to within log( delta ) using secant method
c
c   if found, return:
c     status = 0
c     xnew   = root 
c   else if maximum number of iterations reached before convergence:
c     status = 1
c     xnew   = last estimate of root
c   else (no convergence):
c     status = -1
c     xnew   = last estimate of root
c
c   A. Ackerman Feb-2000
c
      implicit none


c   declare externals
    
      double precision f
      external f


c   declare local storage

      double precision x1, x2, f1, f2, slope, fnew
      double precision y, xlow, xhigh, delta, xnew
      integer status, iter


c   define maximum number of iterations
      integer MAX_ITER
      parameter( MAX_ITER = 50 )


c   define convergence criteria for independent variable
      double precision EPSILON
      parameter( EPSILON = 1d-6 )


c   copy input range into local variables
      x1 = xlow
      x2 = xhigh


c   abort if root not bracketed by initial guesses 

      f1 = log(f(x1)) - log(y)
      f2 = log(f(x2)) - log(y)

c   default return values
      xnew = x1 
      fnew = f1 

      if( f1*f2 .gt. 0 )then 
        status = -1
        return
      endif

c   take logarithm of precision
      delta = log( delta )

c   iterate until root is found or maximum number of iterations are taken

      iter = 0
      do while( ( abs(fnew) .gt. delta ) .and. ( iter .lt. MAX_ITER )
     $  .and. ( x2/x1-1. .gt. EPSILON ) )

c   estimate root from secant between endpoints

        slope = ( f2 - f1 ) / ( x2 - x1 )
        xnew = x1 - f1/slope
        fnew = log(f(xnew)) - log(y)

c   pick new x1 and x2 such that f(x) crosses zero between them
      
        if( fnew*f1 .le. 0 .and. x2/xnew-1. .gt. EPSILON )then
          x2 = xnew
          f2 = fnew
        else 
          x1 = xnew
          f1 = fnew
        endif

        iter = iter + 1

      enddo


c   set status flag accordingly

      if( iter .eq. MAX_ITER .and. x2/x1-1. .gt. EPSILON )then 
        status = 1
      else
        status = 0
      endif

      return
      end
      subroutine layer( nsub_max, 
     $           gas_name, grav, mw_atmos, kz_min, cloudf_min,
     $           mw_cloud, rainf, 
     $           rho_p, supsat, sig_layer, 
     $           cloudf, q_below, 
     $           t_layer, p_layer, kz, chf,cp_layer,d_molecule_layer
     &           ,eps_k_layer,z_layer,
     $           t_top, t_bot, p_top, p_bot, 
     $           qc_layer, qt_layer, rg_layer, reff_layer, 
     $           ndz_layer, 
     $           report_status_r, report_status_q )

c   --------------------------------------------------------------------
c
c   Calculate layer condensate properties by iterating on optical depth
c   in one model layer (convering on optical depth over sublayers).
c
c   input:
c
c     nsub_max    maximum number of sublayers for mesh refinement
c     gas_name    name of condensing vapor
c     kz_min      minimum eddy diffusion coefficient (cm^2/s)
c     cloudf_min  minimum cloud fractional coverage
c     mw_cloud    molecular weight of condensate (g/mol)
c     rainf       rain factor
c     rho_p       density of condensed vapor (g/cm^3)
c     supsat      fractional supersaturation persisting after condensation
c     sig_layer   geometric std deviation of lognormal size distribution
c     t_layer     temperature at layer mid-pt (K)
c     p_layer     pressure at layer mid-pt (dyne/cm^2)
c     t_top,bot   temperature at top and bottom of layer (K)
c     p_top,bot   pressure at top and bottom of layer (dyne/cm^2)
c     cp_layer    pressure specific heat capacity from CLIMA
c     d_molecule_layer  pressure specific molecular size (1/cm^2)
c     eps_k_layer layer specific viscosity
c     z_layer     altitude at layer mid-pt (cm)
c
c   output:
c
c     qc_layer    condensate mixing ratio (g/g)
c     qt_layer    gas + condensate mixing ratio (g/g)
c     rg_layer    geometric mean radius of condensate (cm)
c     reff_layer  effective (area-weighted) radius of condensate (cm)
c     ndz_layer   column of particle concentration in layer (#/cm^2)
c     opd_layer   optical depth for conservative geometric scatterers
c
c     report_status_r    report error status for finding rw
c     report_status_q    report error status for finding qt
c   
c
c   A. Ackerman Dec-2001
c
c   --------------------------------------------------------------------

      implicit none


c   Include common data shared with eddysed()

      include 'globals.h'

      !INCLUDE 'globals.h'
      common /cvm/ readin, metallicity
      !common/eddy_cp/cp_eddy !JDW 2021

      double precision metallicity
      logical readin
      double precision cp_eddy

c   Declare local storage

      logical report_status_r, report_status_q, converge
      integer status_r, status_q
      integer nsub_max, isub, nsub
      character*(*) gas_name
      double precision rainf, rho_p, supsat, sig_layer
      double precision dlnp, r_atmos, scale_h, mixl, scalef_kz, kz
      double precision p_top, t_top, chf, kz_min, cloudf, cloudf_min
      double precision p_bot, t_bot, n_atmos, q_below
      double precision qc_layer, qt_layer, rg_layer, dz_layer
      double precision reff_layer, ndz_layer, qt_below
      double precision lnsig2, dtdlnp, dp_sub, p_bot_sub, p_top_sub
      double precision dtdp, p_bar, t_bar, gradx
      double precision qc_sub, qt_sub, t_sub, p_sub, dz_sub
      double precision z_layer, d_molecule, eps_k, c_p, lapse_ratio
      double precision rg_sub, reff_sub, ndz_sub, w_convect
      double precision qt_bot_sub, qt_top
      double precision opd_layer, opd_test, dp_layer, rho_atmos
      double precision grav, mw_atmos, mfp, visc
      double precision t_layer, p_layer, mw_cloud
      double precision cp_layer,d_molecule_layer,eps_k_layer,mixl_inv
      !double precision d_molecule,eps_k

c   Set error return codes to zero
      status_r = 0
      status_q = 0

c   Number of levels of grid refinement used 
      nsub = 1

c   diameter of atmospheric molecule (cm) (Rosner, 2000)
c   (3.711e-8 for air, 3.798e-8 for N2, 2.827e-8 for H2)
      ! d_molecule = 3.711e-8 ! Generalized define global with gas mixing ratios JDW 2021
      ! Leonard Jones molcule. 
      !call molecular_size(.true.)
      d_molecule = d_molecule_layer*1.e-8
c   parameter in Lennard-Jones potential for viscosity (K) (Rosner, 2000)
c   (78.6 for air, 71.4 for N2, 59.7 for H2)
      !eps_k = 78.6!59.7 jDW 2021
      eps_k=eps_k_layer
      !print*,'d_molecule_layer',d_molecule_layer

c   specific gas constant for atmosphere (erg/K/g)
      !print*,'mw_atmos',mw_atmos
      r_atmos = R_GAS / mw_atmos

c   specific heat of atmosphere (erg/K/g)
      !c_p = 7./2. * r_atmos !Use specific heat routine from clima JDW 2021 from Ty comments
      !Need a first guess for the specific heat because the CPNT(MAXNZ) array 
      !isn't filled in CLIMA until after Eddysed is called. JDW
       c_p=cp_layer
       if(cp_layer .eq. 0.0)then
            c_p=7./2. * r_atmos
       endif

       !print*,"c_p",cp_layer/r_atmos
c   pressure thickness of layer
      dp_layer = p_bot - p_top
      dlnp = log( p_bot/p_top )

c   temperature gradient 
      dtdlnp = ( t_top - t_bot ) / dlnp
      !print*,'(c_p/r_atmos)',(c_p/r_atmos)
      !lapse_ratio = ( t_bot - t_top ) / dlnp / ( (c_p/r_atmos)*t_layer ) !For Polyatomic gasses this will change


c    ratio of lapse rate to adiabat
        if(t_top.lt. 0.0) stop '"t_top" is negative,stop'
        if(t_bot .lt. 0.0) stop '"t_bot" is negative,stop'
        if(isnan(dlnp)) stop '"dlnp" is NaN, stop'

        !print*,t_top,t_bot,dlnp
        DTDP=(DLOG( t_top ) - DLOG( t_bot ) )
     &            / dlnp
        t_bar = 0.5d0*(t_top+t_bot)
        p_bar = 0.5d0*(p_top+p_bot)/1.e6
        !call didgrad (t_bar,p_bar/1.e6,gradx)
        !call didgrad (t_bar,p_bar/1.e6,gradx) JDW
        gradx=1.0!0.17     !This is hardcoded and used as an initial hack for Eddysed JDW
C        call moistgrad (t_bar,p_bar/1.e6,gradx)
        if(isnan(gradx)) stop '"gradx" is NaN, stop'
        if(isnan(dtdp)) stop '"dtdp" is NaN, stop' !dtdp is NaN
        !lapse_ratio = dmin1(1.d0,-dtdp/gradx)
       !print *,'ratio',dtdp,gradx,lapse_ratio
cstop
       !p_layer=p_layer/1e6
c   atmospheric density (g/cm^3)
      rho_atmos = p_layer / ( r_atmos * t_layer )

c   atmospheric scale height (cm)
      scale_h = r_atmos * t_layer / grav

c   convective mixing length scale (cm): no less than 1/10 scale height
      if(isnan(lapse_ratio)) stop '"lapse_ratio" is NaN, stop' !lapse_ratio is NaN
      if(isnan(scale_h)) stop '"scale_h" is NaN, stop'
      mixl = max( 0.1d0, lapse_ratio ) * scale_h

c   mixing length = scale height matches Lunine (1989) model
c   Blackadar expression JDW 2021 Utilizing the Von Karman constant in 
C   Xu Zhang, in Uncertainties in Numerical Weather Prediction, 2021 equation 23
      mixl = scale_h
       !mixl = 0.01*scale_h
      mixl_inv = 1/0.4/z_layer + 1/mixl 
      mixl=1/mixl_inv
      mixl=max1(0.1d0*scale_h,mixl)

      !print*,'mixl/scale_h',mixl/scale_h,scale_h/(100*1000)
c   scale factor for eddy diffusion: 1/3 is baseline
      scalef_kz = 1./3.

c   vertical eddy diffusion coefficient (cm^2/s)
c   from Gierasch and Conrath (1985)
      kz = scalef_kz * scale_h * (mixl/scale_h)**(4./3.) *
     $  ( ( r_atmos*chf ) / ( rho_atmos*c_p ) )**(1./3.)

       !kz=kz/100.0

c      write (*,871) kz,scale_h,mixl,lapse_ratio,chf,p_bar,t_bar
871	format (7e15.5)

c     no less than minimum value (for radiative regions)
      kz = max( kz, kz_min )

c     convective velocity scale (cm/s)
      if(isnan(kz)) stop '"kz" is NaN, stop'
      if(isnan(mixl)) stop '"mixl" is NaN, stop' !mixl is NaN, stop.
      w_convect = kz / mixl
      !print*,'Convective Velocity Scale/Kzz = ',w_convect/1./kz
      !print*,'kz,mixl =',kz,mixl
      !print*,'pressure =',p_bar
      !print*,'dtdp,gradx',dtdp,gradx 
      !print*,'lapse_rato',lapse_ratio
       !print*,'scale_h=',scale_h
c		write (*,876) p_layer/1.d6,t_layer,w_convect,scale_h,scale_h/(w_convect*3600.d0)
876		format (F14.6,1x,f11.2,2e13.4,f12.5)

c   cloud fractional coverage
C       cloudf = cloudf_min +
C      $  max( 0.d0, min( 1.d0, 1.d0-lapse_ratio )) *
C      $  ( 1. - cloudf_min )

c   atmospheric number density (molecules/cm^3)
      n_atmos = p_layer / ( K_BOLTZ*t_layer )
      !print*,'n_atmos',n_atmos,ndz_layer/z_layer

c   atmospheric mean free path (cm)
      mfp = 1. / ( sqrt(2.)*n_atmos*PIE*d_molecule**2 )

c   atmospheric viscosity (dyne s/cm^2)
      visc = 5./16.*sqrt( PIE*K_BOLTZ*t_layer*(mw_atmos/AVOGADRO)) /
     $  ( PIE*d_molecule**2 ) /
     $  ( 1.22 * ( t_layer / eps_k )**(-0.16) )


c   --------------------------------------------------------------------
c   Top of convergence loop

      converge = .false.
      do while ( .not. converge )

c   Zero cumulative values

        qc_layer = 0.
        qt_layer = 0.
        ndz_layer = 0.
        opd_layer = 0.

c   total mixing ratio and pressure at bottom of sub-layer

        qt_bot_sub = q_below
        p_bot_sub = p_bot

c   Loop over sub-layers

        dp_sub = dp_layer / nsub
        do isub = 1, nsub

          qt_below = qt_bot_sub
          p_top_sub = p_bot_sub - dp_sub
          dz_sub = scale_h * log( p_bot_sub/p_top_sub )
          p_sub = 0.5*( p_bot_sub + p_top_sub )
          t_sub = t_bot + log( p_bot/p_sub )*dtdlnp

c   Calculate condensate mixing ratio etc for sub-layer

          call calc_qc( gas_name, rainf, rho_p, mw_cloud,
     $         qt_below, supsat, w_convect, mixl,
     $         dz_sub, grav, mw_atmos, mfp, visc, t_sub, p_sub,
     $         sig_layer, qc_sub, qt_sub, rg_sub, reff_sub, 
     $         ndz_sub, qt_top, status_r, status_q )

c   vertical sums

          qc_layer = qc_layer + qc_sub*dp_sub/grav
          qt_layer = qt_layer + qt_sub*dp_sub/grav
          ndz_layer = ndz_layer + ndz_sub
          !print*,'ndz_layer =',ndz_layer,n_atmos,ndz_sub
          if( reff_sub .gt. 0. )then
            opd_layer = opd_layer + 
     $        1.5*qc_sub*dp_sub/grav/(rho_p*reff_sub)
          endif

c   Increment values at bottom of sub-layer

          qt_bot_sub = qt_top
          p_bot_sub = p_top_sub

        enddo

c    Check convergence on optical depth

        if( nsub_max .eq. 1 )then
          converge = .true.
        elseif( nsub .eq. 1 )then
          opd_test = opd_layer
        elseif( opd_layer .eq. 0. .or. nsub .ge. nsub_max )then
          converge = .true.
        elseif( abs( 1. - opd_test/opd_layer ) .le. 1e-2 )then
          converge = .true.
        else
          opd_test = opd_layer
        endif

        nsub = nsub * 2

      enddo
c   --------------------------------------------------------------------
c   Bottom of convergence loop


c     Report problems finding root the first time it happens

      if( status_r .ne. 0 .and. report_status_r )then
        print*, 'layer():'
        write(*,'(a,i3,3a,1pe10.2)')
     $    ' find_root(vfall) status = ',status_r,
     $    ' for ',gas_name,' at p = ',p_layer/1e6
        print*,' there may be more instances not reported'
        print*,''
        print*,'status_r = ',status_r
        report_status_r = .false.
      endif

      if( status_r .ne. 0 .and. report_status_q )then
        print*, 'layer():'
        write(*,'(a,i3,3a,1pe10.2)')
     $    ' find_root(advdiff) status = ',status_q,
     $    ' for ',gas_name,' at iz,p = ',p_layer/1e6
        print*,' there may be more instances not reported'
        print*,''
        print*,'status_q = ',status_q
        report_status_q = .false.
      endif


c   Update properties at bottom of next layer
      ! print*,'old q_below',q_below
      q_below = qt_top
      ! print*,'New q_below',q_below
c   Get layer averages

      if( opd_layer .gt. 0. )then
        reff_layer = 1.5*qc_layer / (rho_p*opd_layer)
        lnsig2 = 0.5*log( sig_layer )**2
        rg_layer = reff_layer*exp( -5*lnsig2 )
      else
        reff_layer = 0.
        rg_layer = 0.
      endif

      qc_layer = qc_layer*grav / dp_layer
      qt_layer = qt_layer*grav / dp_layer

      return
      end    


      double precision function pvap_al2o3( t )
c
c     calculate saturation vapor pressure (dyne/cm^2) over Al2O3
c     Kozasa et al. Ap J. 344 325
c
c     input temperature in K
c
c     A. Ackerman Feb-2000 modifications for al2o3 M Marley Nov. 2001
c
      implicit none

      double precision t

c     pvap_al2o3 = 0.0259*exp(-73503./t + 22.005)*1e6
c     as I recall this was hand fit to give better agreement w/lodders
c     pvap_al2o3 = exp(-120000./t + 48.78)*1e6
c     as printed in Kozasa
      pvap_al2o3 = exp(-73503./t + 22.01)*1e6


      return
      end
      double precision function pvap_ch4( t )
c
c     calculate saturation vapor pressure (dyne/cm^2) over CH4
c
c     input temperature in K
c
c     Rages Jul-2000 
c
      implicit none

      integer ic
      double precision a(2),b(2),c(2)

c     quantities to calculate methane vapor pressure     
c     (WHAT ARE THE UNITS?)
c
c     AMR   -- molecular weight / ideal gas constant
c     TCRIT -- triple point temperature
c     PCRIT --    "     "   pressure
c     AS    -- specific heat at constant pressure ( gas - solid )
c     AL    --    "      "   "     "        "     ( gas - liquid )
c     ALS   -- latent heat of sublimation
c     ALV   --   "     "   "  vaporization
c
      double precision t,amr,tcrit,pcrit,as,al,als,alv
      parameter( AMR = 16.043 / 8.3143, TCRIT = 90.68, PCRIT = .11719 )
      parameter( AS = 2.213 - 2.650, AL = 2.213 - 3.370 )
      parameter( ALS = 611.10, ALV = 552.36 )
c
c     ic=1: temperature below triple point
c     ic=2: temperature above triple point
c
      ic = 1
      if (t.gt.tcrit) ic = 2

      C(1) = - AMR * AS
      C(2) = - AMR * AL
      B(1) = - AMR * ( ALS + AS * TCRIT )
      B(2) = - AMR * ( ALV + AL * TCRIT )
      A(1) = PCRIT * TCRIT ** ( -C(1) ) * EXP( -B(1) / TCRIT )
      A(2) = PCRIT * TCRIT ** ( -C(2) ) * EXP( -B(2) / TCRIT )

      pvap_ch4 = A(IC) * t**C(IC) * EXP( B(IC) / t )
   
      pvap_ch4= pvap_ch4*1e6    ! convert from bars to dyne/cm^2
c     pvap_ch4= pvap_ch4*1e1    ! convert from Pa to dyne/cm^2: T ~ 240 K
c     pvap_ch4= pvap_ch4*1e3    ! convert from mb to dyne/cm^2
c     pvap_ch4= pvap_ch4*1e2    ! nonsense, but gives cloud base T ~ 140 K

      return
      end
      double precision function pvap_fe( t )
c
c     calculate saturation vapor pressure (dyne/cm^2) over MgSiO3
c     (Lunine et al., 1988, ApJ, 338, 314, citing Barshay and Lewis 1976)
c
c     input temperature in K
c
c     A. Ackerman Feb-2000 modifications for iron by M. Kress March-2000
c
      implicit none

      double precision t

      goto 9999   ! use new Fe expression

      if( t .gt. 1800. )then
            pvap_fe = dexp(-37120./t + 9.86)
      else 
            pvap_fe = dexp(-47664./t + 15.71)
      endif
   
      pvap_fe = pvap_fe*1e6    ! convert from bars to dyne/cm^2

9999  continue

c    OLD EXPRESSION BELOW
c     fit to Fe vapor pressure from email from Frank Ferguson
c     dated 3/25/03
c       pvap_fe = 10.d0**(-19140./t + 8.641)
c       pvap_fe = pvap_fe*1.332*1e3  ! convert from Torr to dyne/cm^2

c    NEW EXPRESSION from Channon Visscher, correspondance on 6/3/11, added 7/27/11 (cvm)

       pvap_fe = 10.d0**(7.09-20833./t)
       pvap_fe = pvap_fe * 1e6   ! convert from bars to dyne/cm^2
 
      return
      end


      double precision function pvap_gas( gas_name, t_layer , p_layer, metallicity,mw_atmos)
c
c   calculate vapor pressure for a gas
c
c   A. Ackerman Nov-2001
c
      !implicit none

c  Declare externals
      double precision mw_atmos
      double precision pvap_ch4, pvap_nh3, pvap_h2o
      double precision pvap_fe, pvap_kcl, pvap_mgsio3, pvap_al2o3, pvap_mg2sio4
      double precision pvap_mns, pvap_zns, pvap_na2s, pvap_cr, pvap_tholin, pvap_soot

      external pvap_ch4, pvap_nh3, pvap_h2o
      external pvap_fe, pvap_kcl, pvap_mgsio3, pvap_al2o3, pvap_mg2sio4
      external pvap_mns, pvap_zns, pvap_na2s, pvap_cr, pvap_tholin, pvap_soot
      COMMON/EBLOK/PG,TG,PG0,IMW,RSURF,OMEGA,POCEAN,IMOIST,
     2  BETA1,BETA2,FVDRY,PDRY

c  Declare passed arguments

      double precision t_layer, p_layer, metallicity
      character*(*) gas_name
      double precision cpsat,co2_psat
c  Evaluate the vapor pressure using the appropriate function

      if( gas_name .eq. 'CH4' )then
        pvap_gas = pvap_ch4( t_layer ) 
      elseif( gas_name .eq. 'NH3' )then
        pvap_gas = pvap_nh3( t_layer )
      elseif( gas_name .eq. 'H2O' )then
            call satrat(t_layer,cpsat)
            if(cpsat .gt. PG) cpsat=POCEAN
            pvap_gas=pvap_h2o( t_layer )
            !print*,'pvap_gas 1',pvap_gas
        !pvap_gas = cpsat*1e6!pvap_h2o( t_layer ) !decreasing the saturation vapor pressure makes clouds form readily in the runaway greenhouse limit JDW 10/2022
        !print*,'pvap_gas 2',pvap_gas
      elseif( gas_name .eq. 'Fe' )then
        pvap_gas = pvap_fe( t_layer )
      elseif( gas_name .eq. 'KCl' )then
        pvap_gas = pvap_kcl( t_layer, metallicity )
      elseif( gas_name .eq. 'MgSiO3' )then
        pvap_gas = pvap_mgsio3( t_layer )
      elseif( gas_name .eq. 'Mg2SiO4' )then
        pvap_gas = pvap_mg2sio4( t_layer , p_layer )
      elseif( gas_name .eq. 'Al2O3' )then
        pvap_gas = pvap_al2o3( t_layer )
      elseif( gas_name .eq. 'MnS' )then
	        pvap_gas = pvap_mns( t_layer, metallicity )
	  elseif( gas_name .eq. 'Na2S' )then
	        pvap_gas = pvap_na2s( t_layer , metallicity)
	  elseif( gas_name .eq. 'ZnS' )then
            pvap_gas = pvap_zns( t_layer , metallicity)
      elseif( gas_name .eq. 'Cr' )then
		    pvap_gas = pvap_cr( t_layer, metallicity )
      elseif( gas_name .eq. 'tholin' )then
		    pvap_gas = pvap_tholin( t_layer )		
      elseif( gas_name .eq. 'soot' )then
                pvap_gas = pvap_soot( t_layer )
      elseif( gas_name .eq. 'CO2')then
            call satco2(t_layer,co2_psat)
            pvap_gas = co2_psat*1e6
            !pvap_gas = pvap_co2(t_layer)
      else
        print*,'stop in pvap_gas(), bad gas_name = ',gas_name
      endif
      
      return
      end


      double precision function pvap_co2( t )

c   --------------------------------------------------------------------
c
c   calculate saturation vapor pressure (dyne/cm^2) over carbon dioxide:
c   liquid when T > 100 K 
c   ice when colder
c
c
c   input temperature in K 
c     
c
c     
c   J. Windsor Oct-2020
c   ---------------------------------------------------------------------

      implicit none 
      double precision A,B,C 
      double precision t
      parameter(A=6.81228)
      parameter(B=1301.679)
      parameter(C=-3.494)
      
      pvap_co2=10**(A-(B/(t+C)))
c     log10(P) = A − (B / (T + C))     

      if( t .lt. 145.26 ) pvap_co2=0.00000001
      if(t .gt. 195.89) pvap_co2=0.00000001

      

      return 
      end 


      



      double precision function pvap_h2o( t )

c   --------------------------------------------------------------------
c
c   calculate saturation vapor pressure (dyne/cm^2) over water:
c   liquid when T > 273 K
c   ice when colder
c
c   input temperature in K
c
c   A. Ackerman Feb-2000
c
c   --------------------------------------------------------------------

      implicit none

c   --------------------------------------------------------------------
c   define constants used in Buck's expressions
c   Buck, 1981 (J. Atmos. Sci., 20, p. 1527)

      double precision BAL, BBL, BCL, BDL
      double precision BAI, BBI, BCI, BDI

      parameter( BAL = 6.1121e3 )
      parameter( BBL = 18.729 )
      parameter( BCL = 257.87 )
      parameter( BDL = 227.3 )

      parameter( BAI = 6.1115e3 )
      parameter( BBI = 23.036 )
      parameter( BCI = 279.82 )
      parameter( BDI = 333.7 )

c   --------------------------------------------------------------------
c   define constants used in Wexler formulas
c   (see Flatau et al., 1992, J. Appl. Meteor. p. 1507)

      double precision GG0, GG1, GG2, GG3, GG4, GG5, GG6, GG7
      double precision HH0, HH1, HH2, HH3, HH4, HH5

      parameter( GG0 =-0.29912729e+4 )
      parameter( GG1 =-0.60170128e+4 )
      parameter( GG2 = 0.1887643854e+2 )
      parameter( GG3 =-0.28354721e-1 )
      parameter( GG4 = 0.17838301e-4 )
      parameter( GG5 =-0.84150417e-9 )
      parameter( GG6 = 0.44412543e-12 )
      parameter( GG7 = 0.28584870e+1 )

      parameter( HH0 = -0.58653696e+4 )
      parameter( HH1 =  0.2224103300e+2 )
      parameter( HH2 =  0.13749042e-1 )
      parameter( HH3 = -0.34031775e-4 )
      parameter( HH4 =  0.26967687e-7 )
      parameter( HH5 =  0.6918651 )

c   --------------------------------------------------------------------
c   define DO_BUCK: .true.  means use Buck's expression,
c                   .false. means use Wexler

      logical DO_BUCK
      parameter( DO_BUCK = .true. )


c   declare local storage

      double precision t, tc
      !print*,'DO_BUCK',DO_BUCK
c   --------------------------------------------------------------------
c   branch on temperature for liquid or ice
c   --------------------------------------------------------------------

      if( t .lt. 273.16 )then

c   --------------------------------------------------------------------
c   saturation vapor pressure over ice

        if( DO_BUCK )then 
          tc = t - 273.16
          pvap_h2o = BAI * exp( (BBI - tc/BDI)*tc / (tc + BCI) )
        else 
          pvap_h2o = 10*exp( 1.0/t* 
     $        ( HH0+(HH1+HH5*log(t)+
     $        ( HH2+(HH3+HH4*t)*t)*t)*t )  )
        endif

      else
    
c   --------------------------------------------------------------------
c   saturation vapor pressure over water
c   for T > 1050 K, fix at 600 bars 

        if( DO_BUCK )then 
          if( t .lt. 1048. )then 
            tc = t - 273.16
            pvap_h2o = BAL * exp( (BBL - tc/BDL)*tc / (tc + BCL) )
          else
            pvap_h2o = 600.e6
          endif
        else 
          pvap_h2o = 10*exp( (1.0/(t*t))* 
     $        ( GG0+(GG1+(GG2+GG7*log(t)+
     $        ( GG3+(GG4+(GG5+GG6*t)*t)*t)*t)*t)*t ) )

        endif

      endif

      return
      end


c
c       double precision function pvap_kcl( t )
c
c       implicit none
c
c       double precision t, pvaplog
c       double precision a1,b1,c1
c       double precision a2,b2,c2
c
c* data from NIST web site (http://webbook.nist.gov)
c*  A B and C for T = 1170-1466 K
c       data a1,b1,c1 / 4.61668, 6910.833, -176.083 /
c*  A B and C for T = 1094-1680 K
c       data a2,b2,c2 / 4.78236, 7440.691, -122.709 /
c
c* This is the first way to calculate the KCl vapor pressure from the NIST
c* database. Valid from T = 1094 - 1680 K.
cc         if (t .lt. 1094.) pvap_kcl = 0.0
cc         if (t .ge. 1094. .and. t .le. 1680.) then
cc         if (t .ge. 1094)  then
c            pvaplog = a2 - b2/(t+c2)
c            pvap_kcl = 10**(pvaplog) 
cc       Gives answer in bars -- need to convert to dynes/cm^2 to be consistent?????? (cvm)
cc         endif
c
c       return 
c       end

cc     cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc     ALTERNATE EXPRESSION TO CALCULATE PVAP_KCL
      double precision function pvap_kcl( t, metallicity )

c     calculate saturation vapor pressure (dyne/cm^2) of KCl
c     (Fit from Lodders data)
c     input temperature in K
c     (cvm)
      implicit none
      double precision t, pvap_kcl_bars, metallicity

      pvap_kcl_bars = 10.d0**(7.6106 - 11382./t)
c     Then convert from bars to dynes/cm^2    
      pvap_kcl = pvap_kcl_bars*1e6   

      return
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function pvap_mg2sio4( t , p )
c
c     calculate saturation vapor pressure (dyne/cm^2) over Mg2SiO4
c     Kozasa et al. Ap J. 344 325
c
c     input temperature in K
c
c     A. Ackerman Feb-2000 modifications for silicates by M. Kress March-2000
c
      implicit none

      double precision t, p, metallicityMH

c     OLD EXPRESSION BELOW 
c       pvap_mg2sio4 = exp(-62279./(t) + 20.94)
c       pvap_mg2sio4 = pvap_mg2sio4*1e6    ! convert from bars to dyne/cm^2

c    NEW EXPRESSION from Channon Visscher, correspondance on 7/27/11, added 8/1/11 (cvm)
c       metallicityMH = 0.0
c       pvap_mg2sio4 = 10.d0**(20.04 - 40540./t - 2*metallicityMH)
c       pvap_mg2sio4 = pvap_mg2sio4 * 1e6  ! convert from bars to dyne/cm^2

c Another new expression from Channon Visscher, correspondance on 10/6/11, includes total pressure dependence and met dep. 
        metallicityMH = 0.0
        pvap_mg2sio4 = 10.d0**(-32488./t + 14.88 - 0.2*log10(p*1e6) - 1.4*metallicityMH) * 1e6 !convered from bars to dynes/cm2

      return
      end



      double precision function pvap_mgsio3( t )
c
c     calculate saturation vapor pressure (dyne/cm^2) over MgSiO3
c     (Lunine et al., 1986, ApJ, 338, 314, citing Barshay and Lewis 1976)
c
c     input temperature in K
c
c     A. Ackerman Feb-2000 modifications for silicates by M. Kress March-2000
c
      implicit none

      double precision t, metallicityMH

c      pvap_mgsio3 = exp(-58663./t + 25.37)   
c      pvap_mgsio3 = pvap_mgsio3*1e6    ! convert from bars to dyne/cm^2
c      NEW EXPRESSIONS FROM CV
c       metallicityMH=0.0
       pvap_mgsio3 = 10.d0**(11.83 - 27250./t - metallicityMH)
       pvap_mgsio3 = 1e6 * pvap_mgsio3 !convert bars -> dynes/cm^2

c Another new expression from Channon Visscher, correspondance on 10/6/11, includes total pressure dependence and met dep. 
c This is the saturation vapor pressure of SiO (which is now the limiting factor, instead of Mg), assuming that the forsterite has 
c already formed underneath. You would need a different expression if only the enstatite cloud existed. 
c        pvap_mgsio3 = 10.d0**(-28665./t + 13.43) * 1e6
       
      return
      end


c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function pvap_nh3( t )
c
c     calculate saturation vapor pressure (dyne/cm^2) over NH3
c
c     input temperature in K
c
c     Marley/Beckmann Jul-2000 
c

      implicit none

      double precision t

      pvap_nh3 = dexp(-86596./t**2 - 2161./t + 10.53)
   
      pvap_nh3 = pvap_nh3*1e6    ! convert from bars to dyne/cm^2

      return
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      double precision function pvap_mns( t, metallicity )

c     calculate saturation vapor pressure (dyne/cm^2) of MnS
c     (Channon's Email, 2011/06/21)
c     input temperature in K
c     (cvm)
      implicit none
      double precision t, pvap_mns_bars, metallicityMH, metallicity
      metallicityMH = metallicity
c     Mn vapor pressure above cloud 
      pvap_mns_bars = 10.d0**(11.5315-23810./t - metallicityMH)
c     Then convert from bars to dynes/cm^2    
      pvap_mns = pvap_mns_bars*1e6   

      return
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      double precision function pvap_na2s( t, metallicity )

c     calculate saturation vapor pressure (dyne/cm^2) of Na2S
c     (Channon's Email, 2011/06/03
c     input temperature in K
c     (cvm)
      implicit none
      double precision t, pvap_na2s_bars, metallicityMH, metallicity
      metallicityMH = metallicity
c     Na vapor pressure above cloud 
      pvap_na2s_bars = 10.d0**(8.5497-13889./t-0.5*metallicityMH)
c     Then convert from bars to dynes/cm^2    
      pvap_na2s = pvap_na2s_bars*1e6   

      return
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      double precision function pvap_zns( t, metallicity )

c     calculate saturation vapor pressure (dyne/cm^2) of ZnS
c     (Channon's Email, 2011/06/21)
c     input temperature in K
c     (cvm)
      implicit none
      double precision t, pvap_zns_bars, metallicityMH, metallicity
      metallicityMH = metallicity
c     Zn vapor pressure above cloud 
      pvap_zns_bars = 10.d0**(12.8117-15873./t - metallicityMH)
c     Then convert from bars to dynes/cm^2    
      pvap_zns = pvap_zns_bars*1e6   

      return
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      double precision function pvap_cr( t, metallicity )

c     calculate saturation vapor pressure (dyne/cm^2) of Cr
c     (Fit from data in email - Fortney 4/22/11) (cvm)
c     input temperature in K
c     (cvm)
      implicit none
      double precision t, pvap_cr_bars, metallicity

c     Cr vapor pressure above cloud 
      pvap_cr_bars = 10.d0**(7.2688-20353./t)
c     Then convert from bars to dynes/cm^2    
      pvap_cr = pvap_cr_bars*1e6   

      return
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      double precision function pvap_tholin( t )

c     calculate saturation vapor pressure (dyne/cm^2) of tholin. 
c 	  right now I'm just going to guess
c     so this is a PLACEHOLDER and designed to be TEMPORARY
c     (cvm, 7/26/12)
      implicit none
      double precision t, pvap_tholin_bars

c     tholin vapor pressure above cloud 
      pvap_tholin_bars = 10.d0**(7.6106 - 11382./t) !this is actually the KCl vapor pressure, not a real tholin one. 
c     Then convert from bars to dynes/cm^2    
      pvap_tholin = pvap_tholin_bars*1e6   

      return
      end
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function pvap_soot( t )

c     calculate saturation vapor pressure (dyne/cm^2) of soot. 
c 	  right now I'm just going to guess
c     so this is a PLACEHOLDER and designed to be TEMPORARY
c     (cvm, 7/26/12)
c     (not actually used in real calculation, which is done in ackerman3.1.f )
      implicit none
      double precision t, pvap_soot_bars

c     tholin vapor pressure above cloud 
      pvap_soot_bars = 10.d0**(7.6106 - 11382./t) !this is actually the KCl vapor pressure, not a real soot one. 
c     Then convert from bars to dynes/cm^2    
      pvap_soot = pvap_soot_bars*1e6   

      return
      end


      double precision function qvs_below( p_test )
c
c   calculate saturation mixing ratio for a gas, 
c   extrapolated below model domain
c
c   A. Ackerman Nov-2001
c
      implicit none

c   Declare externals

      external pvap_gas
      double precision pvap_gas

c   Declare common storage 
      common /cvm/ readin, metallicity
      double precision metallicity
      logical readin

      double precision qv_dtdlnp, qv_p, qv_t, qv_factor
      character*10 qv_gas_name

      common / qvs_below_block /
     $  qv_dtdlnp, qv_p, qv_t, qv_factor, qv_gas_name

c  Declare passed arguments and local storage

      double precision p_test, t_test, pvap_test


c  Extrapolate temperature lapse rate to test pressure

      t_test = qv_t + log( qv_p / p_test )*qv_dtdlnp
 
c  Compute saturation mixing ratio

      pvap_test = pvap_gas( qv_gas_name, t_test, p_test, metallicity )
      qvs_below = qv_factor * pvap_test / p_test

      return
      end
      double precision function vfall( r )
c
c   calculate fallspeed for a spherical particle at one layer in an
c   atmosphere, depending on Reynolds number for Stokes flow.
c
c   For Re_Stokes < 1, use Stokes velocity with slip correction
c   For Re_Stokes > 1, use fit to Re = exp( b1*x + b2*x^2 )
c     where x = log( Cd Re^2 / 24 )
c     where b2 = -0.1 (curvature term) and 
c     b1 from fit between Stokes at Re=1, Cd=24 and Re=1e3, Cd=0.45
c
c   and Precipitation, Reidel, Holland, 1978) and Carlson, Rossow, and
c   Orton (J. Atmos. Sci. 45, p. 2066, 1988)
c
c   all units are cgs
c
c   input parameters passed through argument list:
c
c     r         particle radius (cm)
c
c   input parameters passed through common vfall_block
c
c     vf_grav        acceleration of gravity (cm/s^2)
c     vf_mw_atmos    atmospheric molecular weight (g/mol)
c     vf_mfp         atmospheric molecular mean free path (cm)
c     vf_visc        atmospheric dynamic viscosity (dyne s/cm^2)
c     vf_t           atmospheric temperature (K)
c     vf_p           atmospheric pressure (dyne/cm^2)
c     vf_rhop        density of particle (g/cm^3)
c
c   A. Ackerman Feb-2000
c
      implicit none


c   declare common storage 

      double precision vf_grav, vf_mw_atmos, vf_mfp, vf_visc
      double precision vf_t, vf_p, vf_rhop

      common / vfall_block /
     $  vf_grav, vf_mw_atmos, vf_mfp, vf_visc, 
     $  vf_t, vf_p, vf_rhop

c   declare and define local storage

      double precision r, knudsen, rho_atmos, drho
      double precision reynolds, slip, x, y, b1, b2, cdrag

      parameter( b1 = 0.8 )            ! Ackerman
c     parameter( b1 = 0.86 )           ! Rossow
c     parameter( b1 = 0.72 )           ! Carlson
      parameter( b2 = -0.01 ) 
      parameter( cdrag = 0.45 )        ! Ackerman
c     parameter( cdrag = 0.2 )         ! Rossow
c     parameter( cdrag = 2.0 )         ! Carlson


c   universal gas constant (erg/mol/K)
      double precision R_GAS
      parameter( R_GAS = 8.3143e7 )


c   calculate vfall based on Knudsen and Reynolds numbers

      knudsen = vf_mfp / r
      rho_atmos = vf_p / ( (R_GAS/vf_mw_atmos) * vf_t )
      drho = vf_rhop - rho_atmos

c   Cunningham correction (slip factor for gas kinetic effects)
      slip = 1. + 1.26*knudsen

c   Stokes terminal velocity (low Reynolds number)
      vfall = slip*(2./9.)*drho*vf_grav*r**2 / vf_visc
      reynolds = 2.*r*rho_atmos*vfall / vf_visc

      if( reynolds .gt. 1. )then

c   correct drag coefficient for turbulence (Re = Cd Re^2 / 24)

        x = log( reynolds )
        y = b1*x + b2*x**2
        reynolds = exp(y)
        vfall = vf_visc*reynolds / (2.*r*rho_atmos)

        if( reynolds .gt. 1e3 )then

c   drag coefficient independent of Reynolds number

          vfall = slip*sqrt( 8.*drho*r*vf_grav / (3.*cdrag*rho_atmos) )

        endif

      endif

      return
      end

      subroutine eddysed( 
     $  grav, kz_min, cloudf_min, nsub_max, 
     $  supsat, mw_atmos, do_virtual,
     $  nz, z, z_top, p, p_top, t, t_top, chf, 
     $  ngas, gas_name, gas_mmr, gas_mw, rho_p, sig, rainf, 
     $  kz,cp_layer, qt, qc, ndz, rg, reff, cloudf )

c   --------------------------------------------------------------------
c
c   Given an atmosphere and condensates, calculate size and concentration
c   of condensates in balance between eddy diffusion and sedimentation.
c
c   input scalars:
c
c     grav         gravitational acceleration (cm/s^2)
c     teff         effective temperature (K) !Removed because not used JDW August 2021
c     kz_min       minimum eddy diffusion coefficient (cm^2/s)
c     cloudf_min   minimum cloud fractional coverage
c     nsub_max     maximum number of sublayers for adaptive mesh refinement
c     supsat       supersaturation after condensation (fraction)
c     mw_atmos     molecular weight of atmosphere (g/mol)
c     do_virtual   reduce mixing ratio due to decrease below cloud base
c     nz           number of layers
c     ngas         number of condensing gases
c
c   input vectors:
c
c     z            altitude at layer mid-pt (by pressure) (cm)
c     z_top        altitude at top of layer (cm)
c     p            pressure at layer mid-pt (by pressure) (dyne/cm^2)
c     p_top        pressure at top of layer (dyne/cm^2)
c     t            temperature at layer mid-pt (by pressure) (K)
c     t_top        temperature at top of layer (K)
c     chf          layer convective heat flux (erg/s/cm^2)
c
c     gas_name     names of condensing gases
c     gas_mmr      mass mixing ratio of gas below cloud base (g/g)
c     gas_mw       molecular weight of gas (g/mol)
c     rho_p        density of condensed vapor (g/cm^3)
c     sig          geometric std deviation of lognormal size distribution
c     rainf        ratio of microphysical sed flux to eddy sed flux
c
c   output vectors:
c
c     kz           eddy diffusion coefficient (cm^2/s)
c     qt           total (gas+condensed) mixing ratio of condensate (g/g)
c     qc           mixing ratio of condensed condensate (g/g)
c     ndz          number column density of condensate (cm^-3) !JDW the units should be (#/cm^-2)
c     rg           geometric mean radius of lognormal size distribution
c     reff         droplet effective radius (second moment of size distrib, cm)
c     cloudf       cloud fractional coverage (not applied to clouds)
c   
c
c   A. Ackerman Feb-2000
c
c   --------------------------------------------------------------------

      implicit none


c   Include common data shared with eddysed()

      include 'globals.h'
c   Declare externals

      double precision qvs_below, pvap_gas
      external qvs_below, pvap_gas
      integer j
      integer eJCOLD

c   Declare common storage for qvs_below

      double precision qv_dtdlnp, qv_p, qv_t, qv_factor
      character*10 qv_gas_name

      common / qvs_below_block /
     $  qv_dtdlnp, qv_p, qv_t, qv_factor, qv_gas_name

      common / cloudbase / 
     $  p_cloudbase 

      common/colblok/eJCOLD,upatm_mix
      common/coldtrap/gas_mmr_cold


      !common/eddyblok/eddyopd(NF,MAXNZ), eddyw0(NF,MAXNZ),
      !& eddyg0(NF,MAXNZ),eddyqt(MAXNZ,MAXNGAS),eddyqc(MAXNZ,MAXNGAS),
      !& JCOLD,eJCOLD

      common /cvm/ readin, metallicity
      COMMON/ROSNER/cumul_molec_viscosity(MAXNZ),cumul_molec_size(MAXNZ)
      double precision cumul_molec_viscosity,cumul_molec_size
      double precision upatm_mix
      double precision p_cloudbase(MAXNGAS)
      double precision metallicity
      logical readin
      double precision gas_mmr_cold(MAXNGAS)
c   Declare local storage

      logical report_status_r(MAXNGAS)
      logical report_status_q(MAXNGAS)
      logical do_virtual,converge
      integer cloud_base(MAXNGAS)
      integer status_p
      integer nsub_max
      integer nz, ngas, iz, igas
      integer itop, ibot, incr
      character*(*) gas_name(ngas)
      double precision z(MAXNZ)
      double precision cp_layer(MAXNZ)
      double precision p(MAXNZ)
      double precision t(MAXNZ)
      double precision z_top(MAXNZ+1)
      double precision p_top(MAXNZ+1)
      double precision t_top(MAXNZ+1)
      double precision chf(MAXNZ)
      double precision kz(MAXNZ)
      !double precision gas_mmr(MAXNGAS,MAXNZ)
      double precision gas_mmr(MAXNGAS)
      double precision gas_mw(MAXNGAS)
      double precision rho_p(MAXNGAS)
      double precision qt(MAXNZ,MAXNGAS)
      double precision qc(MAXNZ,MAXNGAS)
      double precision ndz(MAXNZ,MAXNGAS)
      double precision sig(MAXNZ,MAXNGAS)
      double precision rg(MAXNZ,MAXNGAS)
      double precision rainf(MAXNZ,MAXNGAS)
      double precision reff(MAXNZ,MAXNGAS)
      double precision cloudf(MAXNZ)
      double precision rw(MAXNZ,MAXNGAS)
      double precision qc_path(MAXNGAS)
      double precision qt_path(MAXNGAS)
      double precision q_below(MAXNGAS)
      double precision rho_atmos, n_atmos, mw_atmos
      double precision kz_min
      double precision lnsig2, cloudf_min
      double precision c_p, d_molecule
      double precision qlo, qhi, delta_q, eps_k
      double precision scalef_kz, supsat, qvs_factor
      double precision qc_layer, qt_layer, rg_layer, reff_layer
      double precision ndz_layer, kz_layer, dtdlnp
      double precision pvap, p_lo, p_hi, p_base, t_base
      double precision t_bot, p_bot, qvs, t_layer, p_layer, grav
      double precision opd_below
      double precision eps_k_layer,d_molecule_layer
      double precision cpsat
      do igas=1,size(t_top)
            !print*,'t_top',t_top(igas),igas
            if(t_top(igas).lt. 0.0) stop '"t_top is negative, stop' !t_top has a negative value somewhere. 
      enddo
c   Initialize cloud base indices and logical flags that shut off 
c   reporting root problems after the first time

      do igas = 1, ngas
        cloud_base(igas) = 0
        report_status_r(igas) = .true.
        report_status_q(igas) = .true.
      enddo

c   Calculate indices of bottom and top of domain

      if( z(2) .gt. z(1) )then
        ibot = 1
        itop = nz
        incr = 1
      else
        ibot = nz
        itop = 1 
        incr = -1 
      endif

c   --------------------------------------------------------------------
c   Loop over condensates

      do igas = 1, ngas 

c   Start at bottom of domain: get p, T, and qt 

        t_bot = t_top(ibot-incr)
        p_bot = p_top(ibot-incr)
            !print*,'ibot-incr',ibot-incr
        q_below(igas) =gas_mmr(igas)!,ibot)!-incr) !This is where I could make it more physical maybe. JDW
        !print*,'q_below(igas)',q_below(igas)
        
c   Initialize vertical path of condensate <qc_path> to zero
        qc_path(igas) = 0.
        qt_path(igas) = 0.

c   Adjust mixing ratio at bottom of model domain if bottom layer 
c   is saturated

!         if ( do_virtual ) then
!         !print*,'mw_atmos',mw_atmos
!         mw_atmos=2.20
!         qvs_factor = (supsat+1)*gas_mw(igas)/mw_atmos
!         pvap = pvap_gas( gas_name(igas), t_bot , p_bot, metallicity)
!         qvs = qvs_factor*pvap/p_bot

!         if( qvs .le. q_below(igas) )then
! c   Find pressure at cloud base
! c   parameters for finding root 
!           p_lo = p_bot
!           p_hi = p_bot * 1e2
!           delta_q = q_below(igas) / 1e2

! c   temperature gradient
!           dtdlnp = ( t_top(ibot) - t_bot ) / log( p_bot/p_top(ibot) )

! c   load parameters into qvs_below common block
              
!           qv_dtdlnp = dtdlnp
!           qv_p = p_bot
!           qv_t = t_bot
!           qv_gas_name = gas_name(igas) 
!           qv_factor = qvs_factor

!           call find_rootl( qvs_below, q_below(igas), p_lo, p_hi,
!      $         delta_q, p_base, status_p )

 
              
!           if( status_p .ne. 0 )then
!             print*, ''
!             print*, 'unable to find cloud base pressure in eddysed():'
!             print*, ' find_rootl(qvs_below) status = ',status_p,
!      $        ' for ', gas_name(igas)
!             print*,''
!           endif

!           t_base = t_bot + log( p_bot/p_base )*dtdlnp
          

! c   Calculate temperature and pressure below bottom layer
! c   by adding a virtual layer 

!           p_layer = 0.5*( p_bot + p_base )
!           t_layer = t_bot + log( p_bot/p_layer )*dtdlnp


! c   Calculate qc, qt, rg, reff, and ndz for virtual layer
!           !call molecular_size(.false.)
!           !call molecular_viscosity(.false.)
!           call layer( nsub_max, 
!      $         gas_name(igas), grav, mw_atmos, kz_min, cloudf_min,
!      $         gas_mw(igas), rainf(ibot,igas),
!      $         rho_p(igas), supsat, sig(ibot,igas), 
!      $         cloudf(ibot), q_below(igas), 
!      $         t_layer, p_layer, kz_layer, chf(ibot),cp_layer(ibot),
!      $         cumul_molec_size(ibot),cumul_molec_viscosity(ibot),z(iz),
!      $         t_bot, t_base, p_bot, p_base,
!      $         qc_layer, qt_layer, rg_layer, reff_layer, 
!      $         ndz_layer, 
!      $         report_status_r(igas), report_status_q(igas) )
!            !print*,'sig ='!,sig(ibot,igas)

! c   Report optical depth below domain

! c          print*, ''
! c          print*, ' eddysed(): condensing gas = ', gas_name(igas)
! c          print*, ' cloud base at p, T = ', p_base/1e6, t_base
! c          print*, ' optical depth below domain = ',
! c     $        1.5*qc_layer*( p_base - p_bot )/grav /
! c     $        ( rho_p(igas)*reff_layer )

! c          print*, ''

!         endif
!       endif
 
c   --------------------------------------------------------------------
c   Loop over atmospheric layers from the bottom up
        do iz = ibot, itop, incr

            !print*,'q_below 1',q_below(igas)

            !do j=70,101
            !      gas_mmr(1,j)=12./1000.
            !enddo 

            !q_below(igas)=gas_mmr(1,iz)
            !print*,'gas_mmr',iz,gas_mmr(1,ibot-iz:iz) !This makes sure the Cold trap layer from Clima is followed. 
                  !PRINT*,iz,ejcold
             if (iz .lt. eJCOLD)then
             q_below(1)=gas_mmr_cold(1) !upatm_mix!!gas_mmr(1,iz) JWD 2022 testing uncommented
             !q_below(2)=gas_mmr_cold(2)
            !       !print*,'changing',q_below(1)
             endif !Uncomment if broken jdw 2021   Uncomment after testing for CO2 done JDW 2022


             !if(q_below(1) .lt. gas_mmr_cold) q_below(1)= gas_mmr_cold
            if(q_below(1).lt.upatm_mix) q_below(1)=upatm_mix !Uncomment after CO2 testing JDW 2022
            !print*,q_below(igas),z_top(iz),iz



c   Calculate layer qc, qt, rg, reff, and ndz for each layer
            !print*,'cumul_molec_size',cumul_molec_size
          call layer( nsub_max, 
     $         gas_name(igas), grav, mw_atmos, kz_min, cloudf_min,
     $         gas_mw(igas), rainf(iz,igas),
     $         rho_p(igas), supsat, sig(iz,igas), 
     $         cloudf(iz), q_below(igas), 
     $         t(iz), p(iz), kz(iz), chf(iz),cp_layer(iz), 
     $         cumul_molec_size(iz),cumul_molec_viscosity(iz),z(iz),
     $         t_top(iz), t_top(iz-incr), p_top(iz), p_top(iz-incr),
     $         qc(iz,igas), qt(iz,igas), rg(iz,igas), reff(iz,igas), 
     $         ndz(iz,igas), 
     $         report_status_r(igas), report_status_q(igas) )

            !print*,'sig =', iz,igas,sig(iz,igas)
c   Accumulate vertical path of condensate 
      !print*,q_below(igas),z_top(iz)
      !print*,qt(iz,1)-qc(iz,1),t_top(iz)
      !print*,'qvs',qvs
      !qc(iz,1)=gas_mmr(1,iz)
          qc_path(igas) = qc_path(igas) + qc(iz,igas)*
     $       ( p_top(iz-incr) - p_top(iz) )/grav
          qt_path(igas) = qt_path(igas)  + (qt(iz,igas)-qc(iz,igas))*
     $       ( p_top(iz-incr) - p_top(iz) )/grav

!       print*,'layer liquid water path',((qt(iz,igas)-qc(iz,igas))*
!      $       ( p_top(iz-incr) - p_top(iz) )/grav)*1e1, ' kg/m^2',iz

        enddo                          ! nz

c   Print some diagnostics

      print*, ''
      print*, ' eddysed(): condensing gas = ', gas_name(igas)
      print*, ' condensate path = ', qc_path(igas)*1e4, ' g/m^2'
      print*, ' liquid water path = ',qt_path(igas)*1e4,' g/m^2'
      print*, ''

      enddo                            ! ngas
C      do iz=1,101
C           C        do i=1,6
C                                print*,iz,ndz(iz,1),z_top(iz)/1e5
C            C        enddo
C                              enddo

C     Debug check for EddySed Relative Humidity Calculations. 

C      Do iz=1,maxnz
C            call satrat(t(iz),cpsat)
C            pvap = pvap_gas( gas_name(1), t(iz) , p(iz), metallicity)
C            rho_atmos=p(iz)/(R_GAS/mw_atmos*t(iz))
            !pvap/((R_GAS/mw_atmos)*t(iz))/rho_atmos


            !print*,'layer MMR',(qt(iz,1)-qc(iz,1))*(29/18.0),pvap,cpsat*1e6
C             print*,'layer MMR',(qt(iz,1)-qc(iz,1))*(29/18.0),
C      & pvap/P(iz)


C            print*,((qt(iz,1)-qc(iz,1))*(29/18.0))*P(iz)/
C     & (cpsat*1e6),((qt(iz,1)-qc(iz,1))*(29/18.0))*P(iz)/
C     & (pvap),iz
C      enddo



c   vapor pressure of condensate (dyne/cm^2)
c      pvap = pvap_gas( gas_name, t_layer , p_layer, metallicity)
c   saturation factor
c      fs = supsat + 1
c   atmospheric density (g/cm^3)
c      rho_atmos = p_layer / ( R_GAS/mw_atmos * t_layer )
c   mass mixing ratio of saturated vapor (g/g)
      !call satrat(t_layer,psat)
      !print*,'pvap',pvap-(psat*1e6)
c      qvs = fs*pvap/ ( (R_GAS/mw_cloud) * t_layer ) / rho_atmos


      return
      end






         subroutine didgrad(t,p,gradx)
         implicit double precision (a-h,o-z)
         common /cps/ tlog(53),plog(26),cp(53,26),grad(53,26)
         !print*,'grad',grad
         tl = dlog10(t)
         pl = dlog10(p)

       CALL LOCATE (tlog,53,tl,kt)  
       CALL LOCATE (plog,26,pl,kp)  

		ipflag = 0
        if (kp.eq.0) then
c       we are at low pressure, use the lowest pressure point
          factkp = 0.0
          kp = 1
		  ipflag = 1
        endif

      if (kp.ge.26) then
c       we are at high pressure, use the highest pressure point
          factkp = 1.0
          kp = 25
		  ipflag = 1
        endif

		itflag = 0
        if (kt.ge.53) then
c       we are at high temp, use the highest temp point
          factkt = 1.0
          kt = 52
		  itflag = 1
        endif

        if (kt.eq.0) then
c       we are at temp, use the highest temp point
          factkt = 0.0
          kt = 1
		  itflag = 1
        endif


      if (((kp.gt.0).and.(kp.lt.26)).and.(ipflag.eq.0)) then
c        print *,'kp',kp,plog(kp),pl,plog(kp+1)
         FACTkp= (-Plog(Kp)+Pl)/(Plog(Kp+1)-Plog(Kp))
      endif

      if (((kt.gt.0).and.(kt.lt.53)).and.(itflag.eq.0)) then
         FACTkt= (-Tlog(Kt)+Tl)/(Tlog(Kt+1)-Tlog(Kt))
      endif

      gp1 = grad(kt,kp)
      gp2 = grad(kt+1,kp)
      gp3 = grad(kt+1,kp+1)
      gp4 = grad(kt,kp+1)

        gradx = (1.d0-factkt)*(1.d0-factkp)*gp1 +
     1   factkt*(1.d0-factkp)*gp2 + factkt*factkp*gp3 + 
     2   (1.d0-factkt)*factkp*gp4
      !print*,'gradx =',gradx
       !gradx=0.17
      RETURN
      END


C
C  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      SUBROUTINE LOCATE(XX,N,X,J)
        implicit double precision (a-h,o-z)
C	Table searching routine from Numerical Recipes.  For
C       N=14 it is about twice as fast as the previous method.
      DIMENSION XX(N)
	
      JL=0
      JU=N+1
10    IF (JU-JL.GT.1) THEN
      JM=(JU+JL)/2
      IF ((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM))) THEN
      JL=JM
      ELSE
      JU=JM
      ENDIF
      GOTO 10
      ENDIF
      J=JL
      RETURN
      END