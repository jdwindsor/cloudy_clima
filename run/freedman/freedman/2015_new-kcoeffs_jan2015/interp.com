c                       These are the interpolated abunds at a single P:T point
c                       used to calculate the mean mol wt
      common /interp_A/ output_abunds(max_mol_used),jlowP(2),jlowT ! need 2 for two temperatures

c                       these are the interpolated Kcoeffs for a single P:T point
      real kappa_gas_int_pt 
      common /interp_K/ kappa_gas_int_pt(max_windows,ngauss*2)

