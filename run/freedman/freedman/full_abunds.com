      character *5 abunds_in(5)

c     f_abunds are the abunds used to calculate the mean mol wt
      common /full_abunds/ f_abunds(max_tc,max_pc,max_mol_used),
     +                     f_log_abunds(max_tc,max_pc,max_mol_used)
      real itc(max_tc,max_pc),ipc(max_tc,max_pc)
      dimension ncp(max_tc) 
c  this should be the same as nc_p [read in from Kcoeff file]
