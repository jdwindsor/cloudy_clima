PRO INTERPOLATOR

  ; read in CIA
  Ntemp = 115
  Nwn   = 1000
  Nhdr  = 1
  k     = FLTARR(Nwn,Ntemp)
  temp  = FLTARR(Ntemp)
  fn    = 'final1_abel_CIA.dat'
  FOR i=0L,Ntemp-1 DO BEGIN
    READCOL, fn, tempIN, NUMLINE=1, SKIPLINE=Nhdr + (Nwn+1)*i, /SILENT
    READCOL, fn, wnIN, kIN, NUMLINE=Nwn, SKIPLINE=Nhdr + (Nwn+1)*i + 1, /SILENT
    PRINT, tempIN
    temp[i] = tempIN
    k[*,i]  = kIN
  ENDFOR
  wn = wnIN

  ; output temperature grid
  Tout = [60.0000,100.000,150.000,200.000,250.000,$
          300.000,350.000,400.000,450.000,500.000,$
          600.000,700.000,800.000,900.000,1000.00,$
          1250.00,1500.00,1750.00,2000.00,2500.00,$
          3000.00]
  Nout = N_ELEMENTS(Tout)

  ; loop and print
  iwn = WHERE(wn GT 10000. AND wn LE 19980.0) ; only output for wavenumbers greater than this
  FOR i=0,N_ELEMENTS(iwn)-1 DO BEGIN
    kHR   = k[iwn[i],*]
    kout  = 10.^(INTERPOL(kHR,temp,Tout))
    izero = WHERE(kout lt 2.e-33)
    IF( izero[0] NE -1 ) THEN kout[izero] = 0.
    PRINT, wn[iwn[i]],kout[0],kout[1],kout[2],kout[3],kout[4],kout[5],kout[6],kout[7],kout[8],kout[9],kout[10],kout[11],kout[12],kout[13],kout[14],kout[15],kout[16],kout[17],kout[18],kout[19],kout[20]
  ENDFOR

END
