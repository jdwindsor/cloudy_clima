readcol,'NH3.refrind',n,w1,n1,k1
readcol,'H2O.refrind',n,w2,n2,k2
readcol,'MgSiO3.refrind',n,w3,n3,k3
readcol,'Mg2SiO4.refrind',n,w4,n4,k4
readcol,'Fe.refrind',n,w5,n5,k5
readcol,'Al2O3.refrind',n,w6,n6,k6

!p.thick=3
!p.charthick=3

set_plot,'ps'
device,filename='opt2011.ps',/color
tek_color

plot,w1,n1,/xlog,/xsty,xthi=3,ythi=3,charthi=3,/nodata,yran=[.1,100],/ysty,/ylog,$
  xtit='Wavelength (microns)', ytit='n (solid)    k (dotted)',charsi=1.4,$
   tit='Fortney'

oplot,w1,n1,col=2
oplot,w1,k1,col=2,line=1

oplot,w1,n2,col=3
oplot,w1,k2,col=3,line=1

oplot,w1,n3,col=4
oplot,w1,k3,col=4,line=1

oplot,w1,n4,col=5
oplot,w1,k4,col=5,line=1

oplot,w1,n5,col=6
oplot,w1,k5,col=6,line=1

oplot,w1,n6,col=8
oplot,w1,k6,col=8,line=1


xyouts,.4,30,'NH3',col=2
xyouts,.7,30,'H2O',col=3
xyouts,1.2,30,'MgSiO3',col=4
xyouts,3,30,'Mg2SiO4',col=5
xyouts,8,30,'Fe',col=6
xyouts,20,30,'Al2O3',col=8


device,/close
set_plot,'x'



spawn, 'open opt2011.ps'

end
