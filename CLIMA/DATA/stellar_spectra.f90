Program stellar_spectra
! Program that creates the STELLAR_SPECTRA.pdat data table. Written by Ramses Ramirez
Implicit none

real :: Sun(38), GJ581(38), ADLEO(38), GJ644(38), M2600(38),M2700(38),& 
  M2800(38), M3000(38), M3200(38), M3400(38), M3600(38), M3800(38)
integer :: I, WL1(38), WL2(38)

data WL1/2376, 2750, 2850, 3071, 3292, 3412, 3900, 4500, 5400, 5495, 5666, 6050, 6250, 6667, 6910, 7520, 7840, & 
8420, 8910, 9620, 10360, 10700, 11300, 12030, 13070, 14310, 15650, 16880, 18620, 20200, 22030, 24810, 26600, 29200, & 
32390, 35770, 40100, 41720/

data WL2/2750, 2850, 3071, 3292, 3412, 3900, 4500, 5400, 5495, 5666, 6050, 6250, 6667, 6910, 7520, 7840, 8420, &
 8910, 9620, 10360, 10700, 11300, 12030, 13070, 14310, 15650, 16880, 18620, 20200, 22030, 24810, 26600, 29200, 32390, &  
 35770, 40100, 41720, 45450/

data Sun/4.19622E+03,2.52381E+03,1.21123E+04,1.88685E+04,1.28827E+04,5.36868E+04,9.43845E+04,1.72820E+05,1.69032E+04, & 
2.97693E+04,6.52759E+04,3.27690E+04,	6.49410E+04, 3.48310E+04, 8.12603E+04, 3.85621E+04, 6.33436E+04, 4.75263E+04, &
6.19057E+04, 5.64111E+04,2.28957E+04, 3.64219E+04, 3.85110E+04, 4.60209E+04, 	4.47818E+04, 3.93782E+04, 2.90537E+04, &
 3.04737E+04, 1.91326E+04, 1.66074E+04,1.89592E+04, 9.15318E+03, 1.04358E+04, 8.96825E+03, 5.99994E+03, 5.08508E+03, &
 1.43445E+03, 1.17086E+04/


data GJ581/1.e-10 ,6.150852e+00 ,3.783172e+00 ,3.930578e+01 ,9.897419e+01 ,6.936504e+02 ,3.831481e+03 ,1.796130e+04 &
,3.497569e+03 ,5.408555e+03 ,1.325037e+04 ,6.701170e+03 ,2.151627e+04 ,1.041745e+04 ,5.111270e+04 ,2.992022e+04 ,6.813546e+04 &
 ,5.592080e+04 ,7.873554e+04 ,8.720684e+04 ,4.095957e+04 ,7.146951e+04 ,7.897122e+04 ,1.015006e+05 ,8.868386e+04 ,8.837347e+04 &
 ,7.557898e+04 ,8.022509e+04 ,4.984924e+04 ,4.981707e+04 ,5.611716e+04 ,2.097524e+04 ,2.286106e+04 ,2.361345e+04 ,2.143292e+04 &
 ,2.029900e+04 ,5.649430e+03 ,9.165488e+03/

data ADLEO/2.147597e+02 ,1.140012e+02 ,1.098568e+02 ,1.446477e+02 ,1.063284e+02 ,7.325340e+02 ,4.003372e+03 ,1.537698e+04 &
,2.420347e+03 ,4.127380e+03 ,1.028931e+04 ,5.401970e+03 ,1.777807e+04 ,8.649839e+03 ,3.984030e+04 ,2.869007e+04 ,6.021743e+04 &
 ,5.285169e+04 ,8.232664e+04 ,8.548414e+04 ,3.685956e+04 ,6.150306e+04 ,7.074994e+04 ,9.157912e+04 ,9.875700e+04 ,9.691489e+04 &
 ,8.242195e+04 ,9.213580e+04 ,5.922723e+04 ,5.508074e+04 ,6.174431e+04 ,2.337758e+04 ,2.648419e+04 ,2.582309e+04 ,2.209050e+04 &
 ,2.078379e+04 ,5.846160e+03 ,9.741422e+03/

data GJ644/2.665137e+01 ,2.452364e+01 ,4.140269e+01 ,1.199347e+02 ,9.211189e+01 ,3.794875e+02 ,3.158849e+03 ,1.467937e+04 &
 ,2.309535e+03 ,4.364567e+03 ,1.029097e+04 ,5.792553e+03 ,1.595596e+04 ,8.126202e+03 ,3.952357e+04 ,2.403848e+04 ,5.663189e+04 &
 ,4.894066e+04 ,7.229195e+04 ,8.157967e+04 ,3.887763e+04 ,6.844841e+04 ,7.706916e+04 ,1.010653e+05 ,9.041903e+04 ,9.145840e+04 &
 ,8.031340e+04 ,8.791293e+04 ,5.587318e+04 ,5.738482e+04 ,6.710223e+04 ,2.567152e+04 ,2.818554e+04 ,2.941747e+04 ,2.702522e+04 &
 ,2.596668e+04 ,7.314782e+03 ,1.212600e+04/ 

data M2600/2.120159e-03 ,2.378086e-03 ,5.801793e-02 ,7.598821e+00 ,7.535076e+00 ,1.766032e+01 ,2.197150e+02 ,1.673210e+03 &
,2.405589e+02 ,3.850012e+02 ,8.472106e+02 ,6.943047e+02 ,2.209915e+03 ,1.045246e+03 ,9.771773e+03 ,6.989924e+03 ,3.372571e+04 &
 ,2.210422e+04 ,6.297886e+04 ,7.698870e+04 ,3.984710e+04 ,7.476363e+04 ,8.759767e+04 ,1.298871e+05 ,1.082351e+05 ,9.647889e+04 &
 ,1.044576e+05 ,1.011795e+05 ,6.055138e+04 ,6.867367e+04 ,8.287106e+04 ,2.897657e+04 ,3.284950e+04 ,3.214451e+04 ,3.144138e+04 &
 ,3.465622e+04 ,9.856537e+03 ,1.562530e+04/

data M2700/9.241440e-03 ,8.499228e-03 ,2.150301e-01 ,1.588935e+01 ,1.610090e+01 ,4.573786e+01 ,4.590903e+02 ,2.829946e+03 &
 ,4.157019e+02 ,6.889884e+02 ,1.535109e+03 ,1.174574e+03 ,3.502421e+03 ,1.629817e+03 ,1.415113e+04 ,1.019521e+04 ,3.938495e+04 &
 ,2.764784e+04 ,6.728791e+04 ,7.880205e+04 ,3.965841e+04 ,7.279897e+04 ,8.509221e+04 ,1.233203e+05 ,1.078435e+05 ,9.684564e+04 &
 ,9.930748e+04 ,9.916459e+04 ,6.074258e+04 ,6.644722e+04 ,8.011971e+04 ,2.898854e+04 ,3.276601e+04 ,3.139289e+04 ,2.974120e+04 &
 ,3.200132e+04 ,9.141643e+03 ,1.484515e+04/

data M2800/3.337671e-02 ,2.751169e-02 ,5.604777e-01 ,2.945027e+01 ,3.019024e+01 ,9.185499e+01 ,8.091692e+02 ,4.339080e+03 &
,6.349799e+02 ,1.074737e+03 ,2.406201e+03 ,1.761073e+03 ,4.992919e+03 ,2.296092e+03 ,1.862050e+04 ,1.340710e+04 ,4.450646e+04 &
 ,3.259496e+04 ,7.101577e+04 ,8.033179e+04 ,3.947748e+04 ,7.125612e+04 ,8.304688e+04 ,1.180604e+05 ,1.070657e+05 ,9.681274e+04 &
 ,9.524914e+04 ,9.708461e+04 ,6.033156e+04 ,6.431643e+04 ,7.721069e+04 ,2.847650e+04 ,3.210574e+04 ,3.026682e+04 ,2.807782e+04 &
 ,2.970492e+04 ,8.504778e+03 ,1.400876e+04/ 

data M3000/3.712237e-01 ,2.923512e-01 ,3.485603e+00 ,7.985349e+01 ,8.358737e+01 ,3.098664e+02 ,2.088471e+03 ,9.111786e+03 &
 ,1.329167e+03 ,2.389739e+03 ,5.683290e+03 ,3.711449e+03 ,9.610242e+03 ,4.543683e+03 ,2.882376e+04 ,2.034515e+04 ,5.321846e+04 &
 ,4.154498e+04 ,7.613668e+04 ,8.146363e+04 ,3.853404e+04 ,6.784755e+04 ,7.859384e+04 ,1.084594e+05 ,1.040382e+05 ,9.524351e+04 &
 ,8.777449e+04 ,9.239750e+04 ,5.878750e+04 ,5.995558e+04 ,7.101206e+04 ,2.729420e+04 ,3.081530e+04 ,2.824774e+04 ,2.505237e+04 &
 ,2.564432e+04 ,7.363577e+03 ,1.246092e+04/

data M3200/1.799251e+00 ,1.547190e+00 ,1.333938e+01 ,1.662205e+02 ,1.749443e+02 ,7.434585e+02 ,4.123831e+03 ,1.607626e+04 &
 ,2.289678e+03 ,4.351798e+03 ,1.105653e+04 ,6.553020e+03 ,1.581603e+04 ,7.954888e+03 ,3.903769e+04 ,2.644970e+04 ,5.941905e+04 &
 ,4.790352e+04 ,7.866632e+04 ,8.095183e+04 ,3.728536e+04 ,6.456901e+04 ,7.430080e+04 ,1.003812e+05 ,9.979576e+04 ,9.230964e+04 &
 ,8.145443e+04 ,8.766116e+04 ,5.660234e+04 ,5.584065e+04 ,6.491281e+04 ,2.570168e+04 ,2.923041e+04 ,2.622702e+04 ,2.234748e+04 &
 ,2.223891e+04 ,6.388678e+03 ,1.100123e+04/

data M3400/5.005910e+00 ,4.349718e+00 ,3.295816e+01 ,2.870943e+02 ,2.997067e+02 ,1.361966e+03 ,6.722434e+03 ,2.445354e+04 &
 ,3.354734e+03 ,6.631291e+03 ,1.754249e+04 ,9.668833e+03 ,2.232631e+04 ,1.168750e+04 ,4.793237e+04 ,3.105195e+04 ,6.356395e+04 &
 ,5.166866e+04 ,7.950505e+04 ,7.960856e+04 ,3.597679e+04 ,6.152284e+04 ,7.030597e+04 ,9.346217e+04 ,9.513264e+04 ,8.902870e+04 &
 ,7.630857e+04 ,8.359724e+04 ,5.461319e+04 ,5.251905e+04 ,5.944058e+04 ,2.389686e+04 ,2.742112e+04 ,2.426813e+04 ,2.010252e+04 &
 ,1.945964e+04 ,5.576802e+03 ,9.658441e+03/

data M3600/9.729726e+00 ,8.387493e+00 ,5.936062e+01 ,4.261795e+02 ,4.381047e+02 ,2.034149e+03 ,9.488307e+03 ,3.331769e+04 &
 ,4.398841e+03 ,9.029352e+03 ,2.433230e+04 ,1.282618e+04 ,2.871982e+04 ,1.563198e+04 ,5.554282e+04 ,3.434408e+04 ,6.612400e+04 &
 ,5.363900e+04 ,7.919148e+04 ,7.775020e+04 ,3.462272e+04 ,5.865452e+04 ,6.657342e+04 ,8.752331e+04 ,9.043011e+04 ,8.543118e+04 &
 ,7.176290e+04 ,7.993545e+04 ,5.292417e+04 ,4.975158e+04 ,5.470829e+04 ,2.246636e+04 ,2.598860e+04 ,2.282101e+04 ,1.839566e+04 &
 ,1.725099e+04 ,4.922755e+03 ,8.525020e+03/ 

data M3800/1.811244e+01 ,1.525373e+01 ,9.763815e+01 ,5.939752e+02 ,5.962422e+02 ,2.783160e+03 ,1.252018e+04 ,4.338422e+04 &
,5.517728e+03 ,1.176215e+04 ,3.175160e+04 ,1.636564e+04 ,3.573364e+04 ,2.049394e+04 ,6.304434e+04 ,3.687125e+04 ,6.764207e+04 &
 ,5.456927e+04 ,7.805623e+04 ,7.543941e+04 ,3.317765e+04 ,5.573550e+04 ,6.285575e+04 ,8.197085e+04 ,8.544390e+04 ,8.119084e+04 &
 ,6.729879e+04 ,7.585331e+04 ,5.073087e+04 ,4.668570e+04 ,5.014345e+04 ,2.133125e+04 ,2.472708e+04 ,2.162700e+04 ,1.684774e+04 &
 ,1.527797e+04 ,4.338115e+03 ,7.508182e+03/
   

open(2,file='STELLAR_SPECTRA.pdat')
!Write labels
write(2,200)

do I =1,38
   write(2,100), I, WL1(I), WL2(I), Sun(I), GJ581(I), ADLEO(I), GJ644(I), M2600(I),M2700(I), M2800(I), M3000(I), & 
               M3200(I), M3400(I), M3600(I),M3800(I) 
enddo

write(2,*),'!Spectral fluxes are normalized to solar luminosity (1.36e6 ergs/cm^2/sec). Wavelengths(WL) are in Angstroms'


 close(2)

100 format (i2, 2(3x, i5), 12(3x, 1pe12.5))
200 format (1x,'I',4x,'WL1',5x,'WL2',8x,'Sun',11x,'GJ581', 11x,'ADLEO', 10x,'GJ644', 10x, 'M2600', 10x, &
           'M2700', 10x,'M2800', 10x, 'M3000', 10x, 'M3200', 10x, 'M3400', 10x, 'M3600', 10x, 'M3800')
                 
end
































































