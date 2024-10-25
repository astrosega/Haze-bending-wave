pro nowakepotential_torques

  ;restore,'nwpl236dt97z.sav' ;original plot

;  restore,'nwpl236dt97fix.sav'
  
  restore,'/home/sega/IDLWorkspace/Default/nwpl236dt97ultimate2.sav'

  slope = .25
  Ts=3
  dhdot  = -omega*(slope*sin(omega*seconds))
  ;if flag then dhdot  = omega*(slope*cos(omega*seconds))

  deltav = (-omegax + dhdot*y0s*z2s + vkdr *y0s*z1s) * (l/2.)  + (omegay + dhdot*x0s*z2s + vkdr * x0s*z1s) * (W/2.)

  deltavx =( omegaz + dhdot*y0s*x2s + vkdr *y0s*x1s) * (l/2.)

  ;omegay[i] =  omegay[i-1] -2*cos(omega*i*dt)*slope*(omega^2)*(rho*h*l*w*w*w/24)*(cos(phi)- sin(phi)*thetaz[i-1])/(Iyy)*dt
  ;omegax[i] =  omegax[i-1] +2*cos(omega*i*dt)*slope*(omega^2)*(rho*h*w*l*l*l/24)*(sin(phi) +cos(phi)*thetaz[i-1])/(Ixx)*dt + (Iyy - Izz)/Ixx*omegay[i-1]*omegaz[i-1]
  ;plot, seconds*(2*!dpi/omega)^(-1), thetay,xrange=[0,25],yrange=[-3,3]

  deltavk = - omegax * z2s * l/2 + omegay * W/2* z2s + omegaz * l/2. * x2s +  dhdot * y0s * l/2 + dhdot * x0s * W/2 ;+ vzrecord

  deltavj = - omegax * z1s * l/2 + omegay * W/2* z1s + omegaz * l/2. * x1s +  vkdr *y0s* l/2 + vkdr  * x0s * W/2

  ; deltavi = (omegaz*cos(thetay)*sin(thetaz + phi)  - vkdr *cos(thetax)*sin(phi + thetaz])) * (l/2.) + (omegay*cos(thetay)*sin(thetax)  - vkdr* cos(thetay)*cos(thetaz + phi))*W/2 + (-omegax*cos(thetay)*sin(thetax)   - vkdr*sin(thetaz + phi)*cos(thetay))*l/2



  ; plot,  seconds*(2*!dpi/omega)^(-1), deltav;,xrange=[0,2]

  ; plot,  seconds*(2*!dpi/omega)^(-1), sqrt(deltavk^2 + deltavj^2);,xrange=[0,2]
  ; oplot,  seconds*(2*!dpi/omega)^(-1),  sqrt(deltav^2 + deltavx^2), color=1000
  ; oplot, seconds*(2*!dpi/omega)^(-1) ,1e-2*sin(omega*seconds),color=50000
  !p.multi=0
  ; wi,1,wsize = [1180.,680]
  wi,1,wsize = [1880.,980]
  ;cgwindow,'cgPlot',wxsize = [1280.],wysize=[980]


  plot,[1,2]
  ;plot,[1,2]
  deltavk= deltavk[0:-2]
  phis2 = phis

  i = where(abs((shift(phis,-1) - phis)) gt 3)

  ;phis[i[0]+1:-1] =phis[i[0]+1:-1]+2*!dpi
  ;phis[i[1]+1:i[2]] = phis[i[1]+1:i[2]] + 2*!dpi
  ;phis[i[2]+1:i[3]] = phis[i[2]+1:i[3]] + 4*!dpi
  ;phis[i[3]+1:i[4]] = phis[i[3]+1:i[4]] + 2*!dpi
  ;phis[i[4]+1:i[5]] = phis[i[4]+1:i[5]] + 4*!dpi
  ;phis[i[5]+1:-1] = phis[i[5]+1:-1] + 2*!dpi
  ;phis[i[6]+1:i[7]] = phis[i[6]+1:i[7]] - 2*!dpi
  ;phis[i[7]+1:i[8]] = phis[i[7]+1:i[8]] + 2*!dpi
  ;phis[i[8]+1:i[9]] = phis[i[8]+1:i[9]] - !dpi
  ;for j=0, 2 do begin
  phis0 =phis[i[0]+1:i[1]]
  phis1 = phis[i[1]+1:i[2]]
  phis2 = phis[i[2]+1:i[3]]
  phis3 = phis[i[3]+1:i[4]]
  phis4 = phis[i[4]+1:i[5]]
  phis5 = phis[i[5]+1:-1]
  ;phis6 = phis[i[6]+1:i[7]]

  ;   if j eq 0 or j eq 1 then phis2[i[j]+1:-1] = phis2[i[j]+1:-1] + 2*!dpi

  ;  if j mod 2 le 0 then phis2[i[j]+1:-1] = phis[i[j]+1:-1] +2*!dpi else phis[i[j]+1:-1] = phis[i[j]+1:-1] + 2*!dpi
  phis00=make_array(n_elements(seconds),value=25*!dpi/180.)
  ;  print,i[j]
  ;endfor
  ;phis[i[0]+1:-1] =phis[i[0]+1:-1]+2*!dpi

  cgSet_TTFont,'dejavusans'
  pos = cgLayout([1,3], OXMargin=[24, 5], OYMargin=[20, 10], XGap=5, YGap=2);,iymargin=[2,2]
  letra=3

  p=pos[*,0]
  cgplot,seconds*(2*!dpi/omega)^(-1),phis00*180/!dpi,color=cgcolor('black'),linestyle=2,charsize=letra, ytitle='Orientation [deg]',yrange=[-180,180],xrange=[0,Ts],thick=2,xmargin=[2,2],XTICKFORMAT="(A1)",font=1,position=[p[0],p[1],p[2],p[3]+.05]
  cgoplot,seconds*(2*!dpi/omega)^(-1),     thetaxs*180/!dpi ,color=cgcolor('blue'),thick=2,nsum=4
  cgoplot,seconds*(2*!dpi/omega)^(-1),    thetays*180/!dpi, color = cgcolor('red'),thick=2,aspect=[1./3.],nsum=4
  cgoplot,seconds[0:i[0]]*(2*!dpi/omega)^(-1),    phis[0:i[0]+1]*180/!dpi,color= cgcolor('green'),thick=2,aspect=[1./3.]
  cgoplot,seconds[i[0]+1:i[1]]*(2*!dpi/omega)^(-1),    phis0*180/!dpi,color= cgcolor('green'),thick=2,aspect=[1./3.]
  cgoplot,seconds[i[1]+1:i[2]]*(2*!dpi/omega)^(-1),    phis1*180/!dpi,color= cgcolor('green'),thick=2,aspect=[1./3.]
  cgoplot,seconds[i[2]+1:i[3]]*(2*!dpi/omega)^(-1),    phis2*180/!dpi,color= cgcolor('green'),thick=2,aspect=[1./3.]
  cgoplot,seconds[i[3]+1:i[4]]*(2*!dpi/omega)^(-1),    phis3*180/!dpi,color= cgcolor('green'),thick=2,aspect=[1./3.]
  cgoplot,seconds[i[4]+1:i[5]]*(2*!dpi/omega)^(-1),    phis4*180/!dpi,color= cgcolor('green'),thick=2,aspect=[1./3.]
  cgoplot,seconds[i[5]+1:-1]*(2*!dpi/omega)^(-1),    phis5*180/!dpi,color= cgcolor('green'),thick=2,aspect=[1./3.]

  ;cgoplot,seconds*(2*!dpi/omega)^(-1),phis00*180/!dpi,color=cgcolor('black'),thick=2,aspect=[1/.3.]
  ;cgoplot,seconds*(2*!dpi/omega)^(-1),    phis6*180/!dpi,color= cgcolor('green'),thick=2,aspect=[1./3.]
  cgoplot,seconds*(2*!dpi/omega)^(-1),phis00*180/!dpi,color=cgcolor('black'),thick=2,aspect=[1./3.],linestyle=2

  colors = ['Blue','Red', 'Green','white','white','Black']
  items = ['!4h!X$\downx$', '!4h!X$\downy$','!4h!X$\downz$','','','!4h!X$\downz$ (flat ring)']
  ;al_legend,items,linestyle=[0,0,0],color=colors,charsize=3.5,linsize=0.2,thick=2,position=[0.5, 0],box=0
  al_legend,items,linestyle=[0,0,0,0,0,2],color=colors,charsize=3.5,linsize=0.2,thick=2,position=[.8, 195],box=0

  p=pos[*,1]

  cgplot,seconds*(2*!dpi/omega)^(-1),   alog10(sqrt(kepx^2+kepy^2+kepz^2)) ,color=cgcolor('blue'),charsize=letra, ytitle='Torques [1/s$\up2$]',xrange=[0,3],thick=2,xmargin=[2,2],position=[p[0],p[1]-.14,p[2],p[3]],XTICKFORMAT="(A1)",nsum=5,/noerase,font=1;,yrange=[-2,2]
  cgoplot,seconds*(2*!dpi/omega)^(-1),  alog10(sqrt(gravx^2+gravy^2+gravz^2)), color = cgcolor('red'),thick=2,nsum=5
  cgoplot,seconds*(2*!dpi/omega)^(-1),  alog10(sqrt(tidalx^2+tidaly^2+tidalz^2)),color= cgcolor('green'),thick=2,nsum=5
  cgoplot,seconds*(2*!dpi/omega)^(-1),  alog10(sqrt(couplex^2+coupley^2)),color= cgcolor('magenta'),nsum=5,thick=2

  ;cgplot,seconds*(2*!dpi/omega)^(-1),   kepx ,color=cgcolor('blue'),charsize=letra, ytitle='Torques [1/s$\up2$]',xrange=[0,3],thick=2,xmargin=[2,2],position=[p[0],p[1]-.14,p[2],p[3]],XTICKFORMAT="(A1)",/noerase,font=1;,yrange=[-2,2]
  ;cgoplot,seconds*(2*!dpi/omega)^(-1),  kepy , color = cgcolor('red'),thick=2
  ;cgoplot,seconds*(2*!dpi/omega)^(-1),  kepz ,color= cgcolor('green'),thick=2





  ;colors = ['Blue','Red', 'Green','Black']
  ;items = ['!4h!X$\downx$', '!4h!X$\downy$','!4h!X$\downz$','!4h!X$\downz$ flat ring']
  ;;al_legend,items,linestyle=[0,0,0],color=colors,charsize=3.5,linsize=0.2,thick=2,position=[0.5, 0],box=0
  ;al_legend,items,linestyle=[0,0,0,0],color=colors,charsize=4,linsize=0.2,thick=2,position=[.8, 190],box=0

  ;  cgplot,seconds*(2*!dpi/omega)^(-1),   alphax ,color=cgcolor('blue'),charsize=2.8, xtitle='Time (Orbital Periods)', ytitle='Angular acceleraion [1/s^2]',xrange=[0,2];,yrange=[-2,2]
  ;  cgoplot,seconds*(2*!dpi/omega)^(-1),  alphay, color = cgcolor('red')
  ;  cgoplot,seconds*(2*!dpi/omega)^(-1),  alphaz,color= cgcolor('green')
  ;    cgplot, seconds*(2*!dpi/omega)^(-1),    omegax, color=cgcolor('blue'),charsize=2.8, xtitle='Time (Orbital Periods)', ytitle='Angular velocity [1/s]',xrange=[0,30],yrange=[-6e-4,15e-4]
  ;    cgoplot, seconds*(2*!dpi/omega)^(-1),   omegay ,color= cgcolor('red')
  ;    cgoplot, seconds*(2*!dpi/omega)^(-1),   omegaz,color= cgcolor('green')
  ;  cgAxis, YAxis=1, YRange=[0, 0.01],charsize=2.8, /Save, ytitle='Relative speed [m/s]'
  p=pos[*,2]
  cgplot,seconds*(2*!dpi/omega)^(-1),  deltavk*100, nsum=3, xtitle='Time (Orbital Periods)',  ytitle='Rel. speeds [cm/s]',charsize=letra,xrange=[0,Ts],yrange=[-4,4],thick=2,position=[p[0],p[1]-.05,p[2],p[3]-.14],/noerase,font=1
  ;cgLegend, Color=['Blue','Red', 'Green'], Location=[0.11, 0.88],  Titles=['!4h!Xx', '!4h!Xy','!4h!Xz'], Length=0.075, /Box, VSpace=4.75,/Background, bg_color='white',charsize=4.2

  ;   cgLegend, Color=['Blue','Red', 'Green'], Location=[0.2, 0.48],  Titles=['!4a!Xx', '!4a!Xy','!4a!Xz'], Length=0.075, /Box, VSpace=2.75,/Background, bg_color='white',charsize=2.2
  ;cgLegend, Color=['Black'], Location=[0.11, 0.51],  Titles=['!4D!Xv$\downk$'], Length=0.075, VSpace=2.75,/Background, bg_color='white',charsize=3.2
  ;colors = ['black','red']
  items = ['!4D!Xv$\downz$']
  al_legend,items,linestyle=[0],color=['black'],charsize=4,linsize=0.2,thick=2,/top,/left,position=[0.1, -1],box=0

  print,mean(phis[1000:-1])
  write_png,'nowakepotential.png',tvrd(/true)

  ;cgcontrol,output='nowakepotential.png'

end