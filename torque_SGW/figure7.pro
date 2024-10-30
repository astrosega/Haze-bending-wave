pro figure7

;Purpose
;To produce figure 7 in Sega et al 2024. It runs equilibrium_h.pro to compute the equilibrium value for the self-gravity wake pitch angle considering all the torques mentioned in Sega et al 2024. (reduced version of equilibrium_bwpaper.pro)
;Since equilibitum.h runs something aking to wake_rot24.pro multiple times, this program takes a while to run in 2024 (in the hour range)

;Requirements (all in the Github directory).
;equilibrium_h.pro
;tsum.pro
;linspace.pro
;HILLCOL4DEQ2.pro
;TSUM.pro
;WHERE_XYZ.pro


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

      root = equilibrium_h(s,rho_mean,l,h,w,omega,phi=x[i])
      print,root*180/!dpi

      thetas[(s/ds),j] = root
      thetasDDA[(s/ds),j] = x[i]

      ss[(s/ds)] = s


    endfor
    j = j + 1
    num=num-.1

  endforeach


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