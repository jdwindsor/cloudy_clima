c     sum over N gauss points for each window at some P:T
c     Allow for both single and double regions [use gfrac to choose]
         common /gauss_pts_wts/
     +   gauss_pts(ngauss*2), gauss_wts(ngauss*2), gfrac,
     +   ngauss1, ngauss2, igtot

c     allow for variable # of coarse pressures at each coarse T.
c     All inputs should be in ascending order of P:T for files & interpolations
      real log_pc
      common /coarse_grid/
     +   nc_p(max_tc),
     +   pc(max_pc,max_tc),    ! note order
     +   tc(max_tc),
     +   tc_inv(max_tc),       ! for interpolation of abundances/Kcoeffs
     +   log_pc(max_tc,max_pc),plogx

      real kappa_gas
c     for this version, preserve individual gauss points
      common /results/ kappa_gas(max_windows,ngauss*2,max_pc,max_tc)

c number of elements for abunds has to be exact
      common /read_dat/ abunds(max_pc,max_tc,max_elements), n_window,     ! MY mixing ratios - read 
                                                                          ! from Kcoeff files
     + n_species,               ! check use of n_mols_used in other versions
     + window_c(max_windows),   ! # of windows must be = for all files
     + delta_nu(max_windows)    ! # of windows must be = for all files

      character*8 molid_out    ! derive from input

      integer off_setp             ! each T subdirectory has an off_set
      integer off_setT             ! for skipping lowest T directories

c     remove mols_in : read this from file
c     include two offsets: for P and T 
c     assume that only a single value of off_setT is needed:
      common /T1/ nc_t, molid_out(max_elements), off_setp(max_tc),off_setT,ilayers ! add ilayers : Thu Feb 28 2002

