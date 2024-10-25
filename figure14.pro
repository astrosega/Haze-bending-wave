pro figure14

  ;Program that plots different viscosities from the Mimas 5:3 beding.
  ;
  !p.multi=0

  wi, 2, wsize=[1350,1250]
  ;Gresh
  cgplot, [131.902], [340], psym='filledcircle',xrange=[131.3,132.8],yrange=[0,700],err_yhigh=200,err_ylow=100,$
    err_color=cgcolor('blue'),color='blue', symsize = 4,ytitle='Ring viscosity (cm$\up2$ s$\up-1$)',$
    xtitle='Radius [Mm]',charsize=4.2,thick=2,err_thick=2,xthick=3,ythick=3,charthick=2



  ;Lissauer
  cgoplot, [131.902], [260], psym='filleddiamond',xrange=[131400,133150], symsize=4,thick=2,color='cyan'
  ;,/window

  ;Tiscareno
  cgoplot, [131.5903], [280], psym='filledsquare',xrange=[131400,133100], symsize=4,thick=2
  cgtext, [131.3903], [250], 'Prometheus 11:12 ILR',charsize=3
  cgoplot, [131.474], [150], psym='filledsquare',xrange=[131400,133100],symsize=4,thick=2
cgtext, [131.32], [120], 'Pandora 9:5 ILR',charsize=3
  ;,/window

  ;Esposito
  cgoplot, [131.902], [280], psym='opensquare',xrange=[131400,133100],yrange=[50,850],err_yhigh=260,err_ylow=180,symsize=4, err_color=cgcolor('red'),color='red',thick=2,err_width=0.015,err_thick=2;,/window
  cgtext, [131.77], [59], 'Mimas 5:3 IVR',charsize=3,charthick=1
  ;cgoplot, [131.474], [540], psym='opensquare',xrange=[131400,133100],yrange=[50,850],symsize=4,color='red',thick=2
  ;,/window
  cgoplot, [132.7], [60], psym='opensquare',xrange=[131400,133100],yrange=[50,850],symsize=4,color='red',thick=2
  cgtext, [132.39], [20], 'Pandora 11:10 ILR',charsize=3,charthick=1
  ;,/window

  ;; Isotropic model
  cgoplot, [131.902], [580], psym=17,xrange=[131400,133100],symsize=4, err_color=cgcolor('magenta'),color='magenta',thick=2
  ;,/window
  ;cgoplot, [131902], [590], psym=17,xrange=[131400,133100],yrange=[50,850],symsize=2, err_color=cgcolor('dark green'),color='dark green';,/window
  ;cgoplot, [131902], [500], psym=17,xrange=[131400,133100],yrange=[50,850],symsize=2, err_color=cgcolor('brown'),color='brown';,/window


  cgoplot, [131.902], [340], psym='filledcircle',xrange=[131.4,133.15],yrange=[0,850],err_yhigh=200,err_ylow=100,$
    err_color=cgcolor('blue'),color='blue', symsize = 4,ytitle='Ring viscosity (cm$\up2$ s$\up-1$)',$
    xtitle='Radius (Mm)',charsize=4.2,thick=2,err_thick=2

  ; cgoplot, [131902], [528], psym=18,xrange=[131400,133100],yrange=[50,850],symsize=2
  ; cgoplot, [131902], [518], psym=18,xrange=[131400,133100],yrange=[50,850],symsize=2
  ; cgoplot, [131902], [420], psym=18,xrange=[131400,133100],yrange=[50,850],symsize=2



  ;cgLegend, Psym=[16, 14, 15,6 ,17,18], Location=[0.42, 0.85],  Titles=['Gresh et al.', 'Lissauer (1984)','Tiscareno et al. (2007)', 'Esposito et. al (1983)', 'Isotropic Model (this work)', 'Anisotropic Model (this work)'],$
  ;   Length=0, VSpace=2.75,/Background, bg_color='white',charsize=3.5,symthick=2,symsize=3,color=['blue','red','black','orange','cyan','black'];,/window

  colors=['blue','cyan','black','red','magenta'];
  items =['Gresh et al. (1986)', 'Lissauer (1984)','Tiscareno et al. (2007)', 'Esposito et. al (1983)', 'This work']
  Psym=[16, 14, 15,6 ,17]

  al_legend,items,psym=psym,color=colors,charsize=3.5,thick=2.3,symsize=3.3,position=[132, 600];,/top,/right

  ;cgLegend, Psym=[16, 14, 15,6 ,18], Location=[0.65, 0.85],  Titles=['Gresh et al.', 'Lissauer','Tiscareno et al.', 'Esposito et. al', 'This Work'],$
  ;  Length=0, /Box, VSpace=2.75,/Background, bg_color='white',charsize=2.5,symthick=2;,/window

  write_png,'visc1'+ '.jpg',TVRD(/TRUE)



end