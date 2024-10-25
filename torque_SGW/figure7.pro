pro figure7

;Purpose
;To produce figure 7 in Sega et al 2024. It runs equilibrium_h.pro to compute the equilibrium value for the self-gravity wake pitch angle considering all the torques mentioned in Sega et al 2024.

  Mm = 3.7493e19
  Ms = 5.6e26
  Me = 5.972e24
  ;Ms = 1.989e30

  G  = 6.67408e-11      ;SI
  kappa= 1.2822e-4
  sigma = 363.

  Q = 1

  lambda_c = (!dpi*Q/kappa)^2*G*sigma ;this is half the crit wavelengh thing'

  rho = 492.
  ;rho = 900.    ;             330
  epsilons = [.5]         ;.25    5       5
  ;epsilons = [.1,.2,.3,.4,.5,.6,.7]
  ;  rho_mean =  581 ;674.756;350d;DDA;350d    ;440   370   330
  rho_meanDDA = 360d
  l = 232d
  h   = 4d
  W   = 18d

  ;l = 4*lambda_c
  ;W   = lambda_c/3.
  ;h   = W/4.

  print,h,w,l,' dimensions'
  ;First I draw the wave

  G  = 6.67408e-11      ;SI
  cD = 1.50217913e-7    ;this is the curly D in SCL, (1/s^2)
  mu = 1.2957087e-4     ;(1/s) vertical frequency of particles
  muM= .771634e-4       ;(1/s) vertical frequency of mimas
  aM = 185539.          ;Km  mimas semi-major axis
  iM = 1.574 * !Dpi/180 ;rads, mimas inclination

  ;orbital paramenters from Mimas taken from Nasa webpage


  j=complex(0,1)
  ;initialize your variable array. Yes we have to do it again because each occ is optimized for a particular t-vector

  reso  =   0.001 ;km ;This is the resolution of the x axis, not of the ray-tracing
  km    = 300
  start = -21
  x     = dindgen(round(km*(1/reso)), Increment=reso, start=start) ;radial axis, in Km. t=0 at resonance

  rv    = 128000000.
  d     = .015; d is the vertical thinckness of the ring in Km
  d_nor = d
  v     = .05 ;viscosity (m^2/s)
  sigma = 363.
  phi = 25./180. *!dpi

  sigma = 492. * 0.5 * sqrt(!dpi)
  k = 2*!dpi/lambda_c
  n0    = 9.5
  a0    = 10e-2
  amax  = 9.5
  amin  = 11e-3
  q     = -3


  rho_p = 400.

  as = dindgen(10000, start = amin, increment = amax/10000)

  n      = n0*(as/a0)^q
  a_mean = total(n*as)/total(as)

  rho_dist = rho_p * !dpi* (4./3.) * as^3 *  n0*(as/a0)^q

  ;rho_mean = 7*mean(rho_dist)

  ;  rho = 900.
  M   = rho*h*w*l
  Ms=5.6834e26
  omega = 1.28897e-4
  ;omega = 0.000134859 ; a=12.8Mm

  GM = G*M

  rho = 492.


  x = linspace(0.1,!dpi/2,1000)

  ds = .05

  thetas = make_array(2./ds,n_elements(epsilons))
  thetasDDA = make_array(2./ds,n_elements(epsilons))

  ss = make_array(2./ds)
  junk  = make_array(2./ds)
  j=0
  num=.8

  foreach epsilon, epsilons do begin
    ; bmax= Sqrt((1/16 * L^2 * Sin(x)^2) + ((4 * 2^(1/3) * GM * (L + L * Cos(x))) / (27d * L^2 * s^4 * GM * omega^4 * (L + L * Cos(x)) * Sin(x)^2 +$
    ;    Sqrt(55296 * s^6 * GM^3 * omega^6 * (L + L * Cos(x))^3 + 729 * L^4 * s^8 * GM^2 * omega^8 * (L + L * Cos(x))^2 * Sin(x)^4))^(1/3) - (27 * L^2 * s^4 * GM * omega^4 * (L + L * Cos(x)) * Sin(x)^2 + $
    ;    Sqrt(55296 * s^6 * GM^3 * omega^6 * (L + L * Cos(x))^3 + 729 * L^4 * s^8 * GM^2 * omega^8 * (L + L * Cos(x))^2 * Sin(x)^4))^(1/3)) / (6 * 2^(1/3) * s^2 * omega^2))

    ;   A= (27 * L^2 * s^4 * GM * omega^4 * (L + L * Cos(x)) * Sin(x)^2 + $
    ;    Sqrt(55296 * s^6 * GM^3 * omega^6 * (L + L * Cos(x))^3 + 729 * L^4 * s^8 * GM^2 * omega^8 * (L + L * Cos(x))^2 * Sin(x)^4))^(1/3)) / (6 * 2^(1/3) * s^2 * omega^2))
    ;   C=((4 * 2^(1/3) * GM * (L + L * Cos(x))) / (27d * L^2 * s^4 * GM * omega^4 * (L + L * Cos(x)) * Sin(x)^2 +$
    ;    Sqrt(55296 * s^6 * GM^3 * omega^6 * (L + L * Cos(x))^3 + 729 * L^4 * s^8 * GM^2 * omega^8 * (L + L * Cos(x))^2 * Sin(x)^4))^(1/3)
    ;  D=0.25*L*sin(x)

    COMMON SHAREs, s
    COMMON SHARErho, rho_mean
    COMMON SHAREbmin, bmin

    bmin=0


    rho_mean =  826.547 ;bmax approx
    rho_mean=550.547    ;hill

    rho_mean=458 ;shear2 THIS IS THE VALUEEEE FOR THE SHEAR PLOTT

    ;rho_mean=445.547 ;BW paper parameters for l and w
    ;rho_mean=500.547 ;BW paper parameters for l and w




    for s = .01d, 2d, ds do begin
      COMMON SHAREs, s
      COMMON SHARErho, rho_mean

      ;s=1.5

      ;     s=1.5
      ;     x=.4
      A= (27d * L^2 * s^4 * GM * omega^4 * (L + L * Cos(x)) * Sin(x)^2 + $
        Sqrt(55296d * s^6 * GM^3 * omega^6 * (L + L * Cos(x))^3 + 729d * L^4 * s^8 * GM^2 * omega^8 * (L + L * Cos(x))^2 * Sin(x)^4))^(1/3d) / (6d * 2^(1/3d) * s^2 * omega^2)
      C=(4 * 2^(1/3d) * GM * (L + L * Cos(x))) / (27d * L^2 * s^4 * GM * omega^4 * (L + L * Cos(x)) * Sin(x)^2 +$
        Sqrt(55296d * s^6 * GM^3 * omega^6 * (L + L * Cos(x))^3 + 729d * L^4 * s^8 * GM^2 * omega^8 * (L + L * Cos(x))^2 * Sin(x)^4))^(1/3d)
      D=0.25*L*sin(x)
      ;   print,'A= ',A
      ;   print, 'B= ',C

      bmax = 1/2d*D + 1/2.*sqrt(D^2+C-A) + 1/2.*sqrt(2*D^2+A-C+2*D^3/sqrt(D^2+C-A))


      bmax0 = (1/8d) * L * Sin(x) + (1/2d) * Sqrt((1/16d * L^2 * Sin(x)^2) + ((4 * 2^(1/3d) * GM * (L + L * Cos(x))) / (27d * L^2 * s^4 * GM * omega^4 * (L + L * Cos(x)) * Sin(x)^2 +$
        Sqrt(55296d * s^6 * GM^3 * omega^6 * (L + L * Cos(x))^3 + 729 * L^4 * s^8 * GM^2 * omega^8 * (L + L * Cos(x))^2 * Sin(x)^4))^(1/3d) - (27d * L^2 * s^4 * GM * omega^4 * (L + L * Cos(x)) * Sin(x)^2 + $
        Sqrt(55296d * s^6 * GM^3 * omega^6 * (L + L * Cos(x))^3 + 729 * L^4 * s^8 * GM^2 * omega^8 * (L + L * Cos(x))^2 * Sin(x)^4))^(1/3d)) / (6d * 2^(1/3d) * s^2 * omega^2)) + (1/2d) * $
        Sqrt((1/8d * L^2 * Sin(x)^2) - ((4 * 2^(1/3d) * GM * (L + L * Cos(x))) / (27 * L^2 * s^4 * GM * omega^4 * (L + L * Cos(x)) * Sin(x)^2 + $
        Sqrt(55296d * s^6 * GM^3 * omega^6 * (L + L * Cos(x))^3 + 729d * L^4 * s^8 * GM^2 * omega^8 * (L + L * Cos(x))^2 * Sin(x)^4))^(1/3d) + (27d * L^2 * s^4 * GM * omega^4 * (L + L * Cos(x)))))

      ; print,bmax

      ; print,s

      junk[(s/ds)] = min(abs(G*Ms/(rv^3) *cos(x)*(l/2)^3*rho *H * W - 2*(1+epsilon)*GM*H*rho_meanDDA*((l/2)*(1 + cos(x)) + (s*omega)^2*((l/2)*sin(x))^4/(4*GM))) , i); DDA

      f1=2*G*Ms/(rv^3)*cos(x)*sin(x)*(l/2)^3*rho *H*W
      f2=2*(1+epsilon)*GM*H*rho_mean*(bmax*(1 + cos(x)) + sin(x)*(s*omega)^2/(4*GM)*(bmax)^4)

      func = Lambda(x: abs(kep(x)-tidal(x)))
      ;root = FX_ROOT([x[i]+.1,x[i],x[i]-.1], func,tol=100.,double=1)
      ;root = FX_ROOT([.9,.8,.7], func,tol=1e-3,double=1)
      root = equilibrium_h(s,rho_mean,l,h,w,omega,phi=x[i])
      print,root*180/!dpi
      ;print,'diference',abs(kep(root)-tidal(root))


      ; print,root

      ; junk[(s/ds)] = min(abs(2*G*Ms/(rv^3)*cos(x)*sin(x)*(l/2)^3*rho *H*W - 2*(1+epsilon)*GM*H*rho_mean*(bmax*(1 + cos(x)) + sin(x)*(s*omega)^2/(4*GM)*(bmax)^4)) , i)

      ; print,junk
      thetas[(s/ds),j] = root
      thetasDDA[(s/ds),j] = x[i]

      ss[(s/ds)] = s



      ;plot,2*G*Ms/(rv^3)*cos(x)*sin(x)*(l/2)^3*rho *H*W
      ;oplot,2*(1+epsilon)*GM*H*rho_mean*(bmax*(1 + cos(x)) + sin(x)*(s*omega)^2/(4*GM)*(bmax)^4d),color=1000

      ;      if s gt 1.46 then begin
      ;  print,'f1= ',f1[20]
      ; print,'f2= ',f2[20]
      ; print,'x= ',x[i]
      ; print,'root= ',root
      ; print,'stop'

      ;     endif

    endfor
    j = j + 1
    num=num-.1

  endforeach

  ; !p.multi = [1,1,1]

  ;WINDOW, 0, XSIZE=1366, YSIZE=768, TITLE='Haze'

  ;restore, 'eqVW.sav'
  wi,2, wsize = [1366, 768]

  cgplot,ss,thetas[*,0]*180/!dpi, yrange =[0,90],ytitle='Pitch-angle (!4h!X)',xtitle='Shear rate (!4C!X)', background = cgcolor('white'),charsize=2.5,xthick=2.4,ythick=2,thick=2

  ;cgoplot,ss,64.35+2.87 - ss*(36.62-5.53), yrange =[0,90],color=cgcolor('red'),linestyle = 0,thick=1
  ;cgoplot,ss,64.35-2.87 - ss*(36.62+5.53), yrange =[0,90],color=cgcolor('red'),linestyle = 0,thick=1
  cgoplot,ss,64.35 - ss*36.62, yrange =[0,90],color=cgcolor('red'),linestyle = 2,thick=2

  cgoplot,ss,atan(2./7 * sqrt(4-2*ss)/ss)*180/!dpi, yrange =[0,90],color=cgcolor('blue'),linestyle = 2,thick = 2
  ;  cgoplot,ss,  atan(1.932-5.186*.5*ss+4.704*(.5*ss)^2)*180/!dpi, yrange =[0,90],color=cgcolor('green')
  cgoplot,[1.5],[23.5], yrange =[0,90],psym=4,symsize=2,err_yhigh=2,err_ylow=2
  cgoplot,[1.5,1,.5],[25,34,48], yrange =[0,90],psym=5,symsize=2,color=cgcolor('Black')
  ;cgoplot,[1.5,1,.5],[25.4,37,47.5], yrange =[0,90],psym=5,symsize=2,color=cgcolor('magenta')


  cgLegend, Color=['Blue','Red','Black','White','White'],Symcolor=['Blue','Red', 'Green','Black','Black'], Psym=[3,3,3,4,5], linestyle=[2,2,0,0,0], Location=[0.3, 0.8],  Titles=['Michikoshi & Kokubo 2014 (Galactic/Simulations)','Seigar et al. 2006 (Galactic/Observations)', 'This Work','Jerousek+ 2016 (Cassini occultations)','Salo+ 2018 (Numerical simulations)'], Length=0.075, /Box, VSpace=2.75,/Background, bg_color='white',charsize=2.2,symsize=2

  if n_elements(epsilons) gt 2 then cgplot,epsilons,thetas[where(ss eq 1.5),*]*180/!dpi, yrange =[0,90],xtitle='Coefficient of rerstitution',ytitle='!4h!X', background = cgcolor('white'),charsize=2.5

  wi,3, wsize = [1366, 768]
  cgplot,ss,thetasDDA[*,0]*180/!dpi, yrange =[0,90],ytitle='Pitch-angle (!4h!X)',xtitle='Shear rate (!4C!X)', background = cgcolor('white'),charsize=2.5,xthick=2,ythick=2,thick=2,title='DDA'

  ;cgoplot,ss,64.35+2.87 - ss*(36.62-5.53), yrange =[0,90],color=cgcolor('red'),linestyle = 0,thick=1
  ;cgoplot,ss,64.35-2.87 - ss*(36.62+5.53), yrange =[0,90],color=cgcolor('red'),linestyle = 0,thick=1
  cgoplot,ss,64.35 - ss*36.62, yrange =[0,90],color=cgcolor('red'),linestyle = 2,thick=2

  cgoplot,ss,atan(2./7 * sqrt(4-2*ss)/ss)*180/!dpi, yrange =[0,90],color=cgcolor('blue'),linestyle = 2,thick = 2
  ;  cgoplot,ss,  atan(1.932-5.186*.5*ss+4.704*(.5*ss)^2)*180/!dpi, yrange =[0,90],color=cgcolor('green')
  cgoplot,[1.5],[23.5], yrange =[0,90],psym=4,symsize=2,err_yhigh=2,err_ylow=2
  cgoplot,[1.5,1,.5],[25,34,48], yrange =[0,90],psym=5,symsize=2,color=cgcolor('Black')

  cgLegend, Color=['Blue','Red','Black','White','White'],Symcolor=['Blue','Red', 'Green','Black','Black'], Psym=[3,3,3,4,5], linestyle=[2,2,0,0,0], Location=[0.3, 0.8],  Titles=['Michikoshi & Kokubo 2014 (Galactic/Simulations)','Seigar et al. 2006 (Galactic/Observations)', 'This Work','Jerousek+ 2016 (Cassini occultations)','Salo+ 2018 (Numerical simulations)'], Length=0.075, /Box, VSpace=2.75,/Background, bg_color='white',charsize=2.2,symsize=2


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;restore, 'eqHill.sav'
  ;save, thetas, ss, thetasDDA,filename='eqVW.sav'
  restore,'eqVW.sav'

  thisletter="161B

  wi,3, wsize = [1266, 900]
  cgplot,ss[0:-8],thetas[0:-8]*180/!dpi, yrange =[20,52],ytitle='Pitch-angle (!9' + string(thisletter) +'!X$\downw$) [degrees]',xtitle='Shear rate (!18q!X) [-]', background = cgcolor('white'),charsize=3.8,xthick=2.4,ythick=2,thick=3,xrange=[0,1.6],font=1

  cgoplot,[1.5],[23.5], yrange =[20,52],psym=4,symsize=3,err_yhigh=2,err_ylow=2,thick=3
  cgoplot,[1.5,1,.5],[25,34,48], yrange =[0,90],psym=5,symsize=3,color=cgcolor('Blue'),thick=3


  ;cgLegend, Color=['Black','White','White'],Symcolor=['Green','Black','Blue'], Psym=[3,4,5], linestyle=[0,0,0], Location=[0.15, 0.32],  Titles=['This Work','Jerousek+ 2016 (occultations)','Salo+ 2018 (simulations)'], Length=0.075, VSpace=2.75, bg_color='white',charsize=3,symsize=3,symthick=2,thick=3,charthick=1;,/box
  colors=['Black','Black','Blue'];
  items =['This Work','Jerousek+ 2016 (occultations)','Salo+ 2018 (simulations)']
  Psym=[0, 4, 5]
  linestyle = [0,0,0]

  al_legend,items,psym=psym,color=colors,charsize=3.2,thick=3,symsize=3.3,/bottom,/left,linsize=[0.1],box=0


  write_png,'shear4.png',TVRD(/TRUE)
  print,'hell', 4*GM*H*rho_mean*sin(phi)*((l/2)*(1 + cos(phi)) + (1.5*omega)^2*((l/2)*sin(phi))^4/(4*GM))
  print,  G*Ms/(rv^3) *cos(phi)*sin(phi)*(l/2.)^3 *rho *H * W

  print,thetas[where(s eq 1.5)]*180/!dpi
  print,lambda_c

  ;  save, thetas, ss, thetasDDA,filename='eqVW.sav'




end