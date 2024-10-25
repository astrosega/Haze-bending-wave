
Pro figure9
;Description
;Program that plors Figure 9 from Sega et al 2024. It rans two simulations wpl232dt97 (standing for "wake potential simulation with L=232 and a timestep of 97 seconds") and wpl77dt97  (standing for "wake potential simulation with L=77 and a timestep of 97 seconds"). As described in Sega et al 2024

   wpl232dt97, seconds, omega, thetaxs, thetays, phis, deltav, deltavk, vkdr, omegax, omegay, omegaz, w, l, h, x0s, y0s, x2s, x1s, z1s, z2s, dhdot
   wpl77dt97, seconds, omega, thetaxs36, thetays36, phis36, deltav36, deltavk36, vkdr, omegax36, omegay36, omegaz36, w, l36, h, x0s36, y0s36, x2s36, x1s36, z1s36, z2s36, dhdot


  slope = .25
  Ts=3
  dhdot  = -omega*(slope*sin(omega*seconds))
  wi,1,wsize = [1880.,980]

     !p.multi=[0,1,2]
 
  i = where(abs((shift(phis,-1) - phis)) gt 3)
 
   phis0 =phis[0:i[0]]
  phis1 =phis[i[0]+1:i[1]] ;- 2*!dpi*95/100.
  phis2 = phis[i[1]+1:i[2]]
  phis3 = phis[i[2]+1:-1] ;-2*!dpi

  i = where(abs((shift(phis,-1) - phis)) gt 3)
 pos = cgLayout([1,2], OXMargin=[24, 5], OYMargin=[20, 5], XGap=5, YGap=5)

  phis[i[0]+1:-1] =phis[i[0]+1:-1]+2*!dpi

cgplot,seconds[0:i[0]]*(2*!dpi/omega)^(-1),    phis0*180/!dpi,color= cgcolor('black'), ytitle='Orientation [deg]',yrange=[-180,180],xmargin=[2,2],xrange=[0,Ts],thick=2,charsize=3.5,font=1,position=pos[*,0],XTICKFORMAT="(A1)"
cgoplot,seconds[i[0]+2:i[1]+3]*(2*!dpi/omega)^(-1),    phis1*180/!dpi,color= cgcolor('black'), ytitle='Degrees',yrange=[-180,180],xmargin=[2,2],xrange=[0,Ts],thick=2,charsize=3.5,font=1,position=pos[*,0],XTICKFORMAT="(A1)"
cgoplot,seconds[i[1]:i[2]]*(2*!dpi/omega)^(-1),    phis2*180/!dpi,color= cgcolor('black'), ytitle='Degrees',yrange=[-180,180],xmargin=[2,2],xrange=[0,Ts],thick=2,charsize=3.5,font=1,position=pos[*,0],XTICKFORMAT="(A1)"
cgoplot,seconds[i[2]:-1]*(2*!dpi/omega)^(-1),    phis3*180/!dpi,color= cgcolor('black'), ytitle='Degrees',yrange=[-180,180],xmargin=[2,2],xrange=[0,Ts],thick=2,charsize=3.5,font=1,position=pos[*,0],XTICKFORMAT="(A1)"

omegax=omegaz[0:-2]
 thetax = tsumcumul(seconds,omegax)

 
 cgoplot,seconds*(2*!dpi/omega)^(-1),    phis36*180/!dpi,color= cgcolor('red'),thick=2,aspect=[1./3.],linestyle=2
 
 
 lineitems = ['solid','DASHED']
 colors = ['black','red']
 items = ['!4h!X$\downz$ for L = 232m','!4h!X$\downz$ for L = 77m']
 al_legend,items,linestyle=[0,2],color=colors,position=[1.3,-80],charsize=3.5,linsize=0.2,thick=2,box=0
 
 
  cgplot,seconds*(2*!dpi/omega)^(-1),  abs(deltavk*100), nsum=2, xtitle='Time [Orbital Periods]',  ytitle='Relative speed [cm/s]',xrange=[0,Ts],thick=2,charsize=3.5,font=1,yrange=[0,6] ,position=pos[*,1],/noerase;,yrange=[-8,8];
  cgoplot,seconds*(2*!dpi/omega)^(-1), abs(deltavk36*100), nsum=1, xtitle='Time [Orbital Periods]',  ytitle='Relative speed [cm/s]',charsize=3.5,xrange=[0,Ts],yrange=[-1,1],color=cgcolor('red'),linestyle=2,thick=2;

lineitems = ['solid','DASHED']
colors = ['black','red']
items = ['|!4D!Xv$\downz$| for L = 232m','|!4D!Xv$\downz$| for L = 77m']

al_legend,items,linestyle=[0,2],color=colors,/bottom,/left,charsize=3.5,linsize=0.2,thick=2,position=[0.05,4],box=0

  write_png,'wakepotential.png',tvrd(/true)
 
 print,mean(phis[1000:-1])
  
 print,min(omegax)/max(dhdot)



end
