pro figure11

;Produces the vizualizaton of the haze (figure 11) shown in Sega et al 2024. The code uses the routine haze_creation.pro which is the core simulation that computes the trajectories of particles after their collisions with the self-gravity wakes
;
;Change Log
;
;Date created 09-09-2020

  G  = 6.67408e-11      ;SI
  cD = 1.50217913e-7    ;this is the curly D in SCL, (1/s^2)
  mu = 1.2957087e-4     ;(1/s) vertical frequency of particles
  muM= .771634e-4       ;(1/s) vertical frequency of mimas
  aM = 185539.          ;Km  mimas semi-major axis
  iM = 1.574 * !Dpi/180 ;rads, mimas inclination
  rv    = 131902.0


  J4 = -935.83e-6;

  J2 = 16290.71e-6;
  M = 5.6834e26;
  R = 60330000;

  omegaM = 0.0000771634
  muM    = 0.000077365

  ; omegap = ((G*M/rv^3) (1 + (J2)*(3/2)*(R/rv)^2 - (J4) (15/8)*(R/rv)^4))^.5


  phis = linspace(0,2*!dpi,6)
  sigma  =  363
  Hsigma = (sigma/1000.)
  A      = (476.92d/(Sqrt(Hsigma)))/1000d

  j=complex(0,1)
  ;initialize your variable array. Yes we have to do it again because each occ is optimized for a particular t-vector

  reso  =   0.001 ;km ;This is the resolution of the x axis, not of the ray-tracing
  km    = 300
  start = -21
  x     = dindgen(round(km*(1/reso)), Increment=reso, start=start) ;radial axis, in Km. t=0 at resonance
  i = -start/reso + 100/reso
  epsilon = 0.1
  S       = .000025

  omega = 4*((G*M/((rv-x[i])*1000)^3)*(1 + J2*(3./2.)*(R/((rv-x[i])*1000))^2 - J4*(15./8.)*(R/((rv-x[i])*1000))^4))^.5 - 4*OmegaM - muM



  d     = .015; d is the vertical thinckness of the ring in Km
  v     = .05 ;viscosity (m^2/s)
  sigma = 363.
  dampingfactor = (1./3.)*(cD^2)*mu/(2*!Dpi*G*sigma)^3             ;this is in SI (meters. remember t and rv are in km)
  damp          =  Exp(-((x^3/rv^2)*1000d*(dampingfactor)*(v)))    ;the 1000 is to transform x into meters since the viscosity is in meters

  Hsigma=(sigma/1000.)


  A=(476.92d/(Sqrt(Hsigma)))/1000d ;amplitud of the wave
  e=(2d*!DPI*G*sigma/(rv*1000d*cD)) ;the 1000 is to transform rv into m
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  upper = (rv/(rv+x))^0.5 * (A/(sqrt(!DPI))) * damp * (Exp(j*x^2/(2*e*rv*rv))*(0.5*sqrt(!DPI/2)-0.5*j*sqrt(!DPI/2) + fresnel_complex(x/(sqrt(2*e)*rv))))
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  slope     = -(shift(upper, 1) - upper)/reso
  slope[0]  = 0
  slope=abs(slope)
  Amp = Abs(upper)
  k = slope/amp



  theta = linspace(0, 3*!dpi)

  !p.multi = [0,0,0]

  ; WINDOW, 0, XSIZE=1920, YSIZE=1200, TITLE='Haze'

  graph =  plot(theta/(2*!dpi),A*damp[i]*cos(theta), color='red', yrange = [-.5,.5], xtitle = 'Time in vertical periods!X', font_size=32,thick =4,dim=[1820,720],Name='Wave motion',layout=[2,1,1],margin=[.2, .15, .02, .05],ytitle = 'z [km]')

  ;cgLegend, Color=['Red','Blue', 'Orange','Cyan'], Location=[0.2, 0.88],  Titles=['Ring', 'Trajectories after collision with SGW','Trajectories after 1st pass','Trajectories after 2nd pass'], Length=0.075, /Box, VSpace=2.75,/Background, bg_color = 'white', charsize=2.2
  ;  cgLegend, Color=['Red','Blue'], Location=[0.2, 0.88],  Titles=['Ring', 'Trajectories after collision with SGW'], Length=0.075, /Box, VSpace=2.75,/Background, bg_color = 'white', charsize=3


  ind2 = 0
  zs = make_array(n_elements(phis))

  foreach phi, phis do begin

    haze_creation, phi, theta, omega, i, sigma, 3, epsilon, x, upper, slope, damp, k, A = A, koeffs = koeffs,s=S, plotf=1, landr=1,c=c

    graph=   plot(theta/(2*!dpi), c[0,0]*cos(theta) + c[1,0]*sin(theta),color=cgcolor('red'),overplot=1,Name=['Motion after collision'])
    if phi eq phis[0] then leg  = legend(position=[.34,.93], font_size = 21,linestyle=6)
    graph1=  plot(theta/(2*!dpi),c[0,1]*cos(theta) + c[1,1]*sin(theta),color=cgcolor('red'),overplot=1)

  endforeach
  tic

  indi = 0

  G  = 6.67408e-11      ;SI
  cD = 1.50217913e-7    ;this is the curly D in SCL, (1/s^2)
  mu = 1.2957087e-4     ;(1/s) vertical frequency of particles
  muM= .771634e-4       ;(1/s) vertical frequency of mimas
  aM = 185539.          ;Km  mimas semi-major axis
  iM = 1.574 * !Dpi/180 ;rads, mimas inclination
  rv    = 131902.0


  J4 = -935.83e-6;

  J2 = 16290.71e-6;
  M = 5.6834e26;
  r = 185539000;
  R = 60330000;

  omegaM = 0.0000771634
  muM    = 0.000077365

  epsilon = .1


  ; omegap = ((G*M/rv^3) (1 + (J2)*(3/2)*(R/rv)^2 - (J4) (15/8)*(R/rv)^4))^.5


  phis = linspace(0,2*!dpi,12)
  sigma  =  363.
  Hsigma = (sigma/1000.)
  A      = (476.92/(Sqrt(Hsigma)))/1000.
  j=complex(0,1)
  ;initialize your variable array. Yes we have to do it again because each occ is optimized for a particular t-vector

  reso  =   0.001 ;km ;This is the resolution of the x axis, not of the ray-tracing
  km    = 300
  start = -21
  x     = dindgen(round(km*(1/reso)), Increment=reso, start=start) ;radial axis, in Km. t=0 at resonance

  ; ri = 1d
  ;rf = 150d
  ri = double(start) +.5
  rf = double(km + start)


  res = .4d

  rad = make_array((rf - ri)/res + 1)
  zmax = make_array((rf - ri)/res+ 1)
  zmin = make_array((rf - ri)/res+ 1)

  d     = .015; d is the vertical thinckness of the ring in Km
  v     = .05 ;viscosity (m^2/s)
  sigma = 363.
  dampingfactor = (1./3.)*(cD^2)*mu/(2*!Dpi*G*sigma)^3             ;this is in SI (meters. remember t and rv are in km)
  damp          =  Exp(-((x^3/rv^2)*1000.*(dampingfactor)*(v)))    ;the 1000 is to transform x into meters since the viscosity is in meters

  Hsigma=(sigma/1000.)


  A=(476.92/(Sqrt(Hsigma)))/1000. ;amplitude of the wave

  ;The 476.92 is calculated by estimating the amplitude of the driving force (MIMAS g field). It depends on the driving frequency and
  ;natural frequency as well as some parameters of MIMAS orbit. This value is computed in SCL. A is in Kilometers. Note that in Gresh et al.
  ;(Radio occultations) the best fitted amplitude differed by this one by a factor of 4 while in SCL (optical) it differed by 1.15


  ;e is a parameter that appears in the solution of the DE and makes the exponent, which is part of the solution, dimensionless. It is related to the wave
  ;number of the bending wave and it depends on G, sigma, and the difference in the driving and natural frequencies. It also depends on the distance
  ;from saturn at resonance. This value of e was calculated by me using the Keplerian approximation (vertical frequency=orbital frequency).
  ; e is defined in SCL
  
  e=(2.*!DPI*G*sigma/(rv*1000d*cD)) ;the 1000 is to transform rv into m

  ;Here comes the wave profile at maximum amplitude at resonance (this sets the phase to -Pi/4). The sqrt of pi in the amplitude comess from the fresnel integral, don't panick.
  ;it is not part of the amplitude
  fase = 1.5
  upper = (rv/(rv+x))^0.5 * (A/(sqrt(!DPI))) * exp(j*(fase + !dpi/4.)) * damp * (Exp(j*x^2/(2*e*rv*rv))*(0.5*sqrt(!DPI/2)-0.5*j*sqrt(!DPI/2) + fresnel_complex(x/(sqrt(2*e)*rv))))

  slope     = -(shift(upper, 1) - upper)/reso
  slope[0]  = 0

  s = .0000242d
  ;nhits = 3
  ;slope1 = -(shift(slope, 1) - slope)/reso
  slope=abs(slope)

  Amp = Abs(upper)
  k = slope/amp


  for ind = ri, rf, res do begin ;This loop goes to all the considered radial location between ri and rf and computes the results of the collision there.
    i = -start/reso + ind/reso
    omega = 4*((G*M/((rv-x[i])*1000)^3)*(1 + J2*(3./2.)*(R/((rv-x[i])*1000))^2 - J4*(15./8.)*(R/((rv-x[i])*1000))^4))^.5 - 4*OmegaM - muM

    !p.multi = [0,0,0]
    nhits = 2
    ind2 = 0
    zs = make_array(n_elements(phis),nhits)
    koeffs = 0


    foreach phi, phis do begin
      haze_creation, phi, theta, omega, i, sigma, nhits, epsilon, x, upper, slope, damp, k, A = A, koeffs = koeffs, s=s, fase = fase,landr=1


      Zs[ind2,*] = koeffs
      ind2+=1
    endforeach

    ; print, indi

    Zmax[indi] = max(zs)

    Zmin[indi] = min(zs)

    rad[indi] = ind

    indi=indi+1

  endfor

  toc

  upper = (rv/(rv+x))^0.5 * (A/(sqrt(!DPI))) * exp(j*(fase+!dpi/2)) * damp * (Exp(j*x^2/(2*e*rv*rv))*(0.5*sqrt(!DPI/2)-0.5*j*sqrt(!DPI/2) + fresnel_complex(x/(sqrt(2*e)*rv))))

  graf = plot(-x,A*damp*cos(tsumcumul(x,k) - !dpi/4 + fase),yrange = [-1,1], xrange = [-150,0], xtitle = 'Distance from resonance [km]',thick =3,color ='red', Name = 'Ring',dim=[1820,980],linestyle=2,/current,layout=[2,1,2],margin=[.12, .15, .05, .05])
  graf.font_size = 32
  ; graf =plot(rad,zmax, overplot = 1, transparency=70.,fill_background=1,fill_color="blue",fill_transparency=90., fill_level=0,color='cornflower', Name = 'Haze',thick=4)
  graf2 =fillplot(-rad,[[zmax],[zmin]], overplot = 1,thick=1,fill_color="blue",fill_transparency=80.,color ='cornflower',name='Haze')
  leg  = legend(target =[graf,graf2],position=[.973,.93], font_size = 21,linestyle=6);,position=[.66,.93]

  graf = plot(-x,A*damp*cos(tsumcumul(x,k) - !dpi/4 + fase),xrange = [-73.75,-71.75], yrange = [-.62,-0.42], thick =3,color ='red',dim=[1820,980],layout=[2,1,2],linestyle=2,/current,margin=[0.8],position=[0.565,.20,0.7,.39],XTICKFONT_SIZE=0,YTICKFONT_SIZE=0)
  ;[0.26,0.32,0.5,0.5] xrange = [-71.75,-73.75]yrange = [-.62,-0.42], xrange = [-59.5,-61.5];,yrange = [-.75,-0.55];xrange = [-71.75,-73.75]yrange = [-.62,-0.42]

  yaxis = AXIS('Y', LOCATION='right'$
    ,target=graf,TICKFONT_SIZE=12.5)


  graf2 =fillplot(-rad,[[zmax],[zmin]], overplot = 1,thick=1,fill_color="blue",fill_transparency=80.,color ='cornflower',name='Haze')

  graf = plot(-x,A*damp*cos(tsumcumul(x,k) - !dpi/4 + fase),yrange = [-.1,0.1], xrange = [-77,-75], thick =3,color ='red',dim=[1820,980],layout=[2,1,2],linestyle=2,/current,margin=[0.8],position=[0.565,.73,0.7,.93],XTICKFONT_SIZE=0,YTICKFONT_SIZE=0,name='inlet1');[0.65,.74,0.770,.94][80,82]

  yaxis = AXIS('Y', LOCATION='right' $
    ,target=graf,TICKFONT_SIZE=12.5)

  graf2 =fillplot(-rad,[[zmax],[zmin]], overplot = 1,thick=1,fill_color="blue",fill_transparency=80.,color ='cornflower',name='inlet2')

  a=max(zmax(where(rad gt 40)),i)
  b=min(zmin(where(rad gt 40)),j)

  print,(a-b)/(rad[where(zmax eq a)]-rad[where(zmin eq b)])

  print,max(zmax-zmin)

  ; PLOT,rad,zmax-zmon


  der = ((shift(zmin,-1)-zmax)/(rad[10]-rad[11]))
  print,max(abs(der[0:-2]))


end

