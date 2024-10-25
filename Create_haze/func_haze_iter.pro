function func_haze_iter, upper


 
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
  
  epsilon = 1
 

  ; omegap = ((G*M/rv^3) (1 + (J2)*(3/2)*(R/rv)^2 - (J4) (15/8)*(R/rv)^4))^.5


  phis = linspace(0,2*!dpi,20)
  ;phis = [0,!dpi/12.,!dpi/10, !dpi/8, !dpi/3. ,!dpi/2.5, !dpi/2]
  ;phis = [!dpi/3. ,!dpi/2.5, !dpi/2]
  ;phis = [!dpi/12.,!dpi/2.1]

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
  
  
  res = .2d
  
  rad = make_array((rf - ri)/res + 1)
   zmax = make_array((rf - ri)/res+ 1)
    zmin = make_array((rf - ri)/res+ 1)
    
    epsilon = 1
    

    d     = .015; d is the vertical thinckness of the ring in Km
    v     = .04 ;viscosity (m^2/s)
    sigma = 363.
    dampingfactor = (1./3.)*(cD^2)*mu/(2*!Dpi*G*sigma)^3             ;this is in SI (meters. remember t and rv are in km)
    damp          =  Exp(-((x^3/rv^2)*1000.*(dampingfactor)*(v)))    ;the 1000 is to transform x into meters since the viscosity is in meters

    Hsigma=(sigma/1000.)


    A=(476.92/(Sqrt(Hsigma)))/1000. ;amplitud of the wave

    ;The 476.92 is calculated by estimating the amplitude of the driving force (MIMAS g field). It depends on the driving frequency and
    ;natural frequency as well as some parameters of MIMAS orbit. This value is computed in SCL. A is in Kilometers. Note that in Gresh et al.
    ;(Radio occultations) the best fitted amplitude differed by this one by a factor of 4 while in SCL (optical) it differed by 1.15


    ;e is a parameter that appears in the solution of the DE and makes the exponent, which is part of the solution, dimensionless. It is related to the wave
    ;number of the bending wave and it depends on G, sigma, and the difference in the driving and natural frequencies. It also depends on the distance
    ;from saturn at resonance. This value of e was calculated by me using the Keplerian approximation (vertical frequency=orbital frequency).
    ; e is defined in SCL

    e=1
    e=double(e)
    e=(2.*!DPI*G*sigma/(rv*1000.*cD)) ;the 1000 is to transform rv into m

    ;Here comes the wave profile at maximum amplitude at resonance (this sets the phase to -Pi/4). The sqrt of pi in the amplitude comess from the fresnel integral, don't panick.
    ;it is not part of the amplitude

    upper = (rv/(rv+x))^0.5 * (A/(sqrt(!DPI))) * damp * (Exp(j*x^2/(2*e*rv*rv))*(0.5*sqrt(!DPI/2)-0.5*j*sqrt(!DPI/2) + fresnel_complex(x/(sqrt(2*e)*rv))))

    slope     = -(shift(upper, 1) - upper)/reso
    slope[0]  = 0
    
    s = .0003d
    ;nhits = 3
    ;slope1 = -(shift(slope, 1) - slope)/reso
    slope=abs(slope)

    Amp = Abs(upper)
    k = slope/amp
 
  
  for ind = ri, rf, res do begin
    i = -start/reso + ind/reso
    
    

 ; omega = 4*((G*M/((rv-x[i])*1000)^3)*(1 + J2*(3./2.)*(Rs/((rv-x[i])*1000))^2 - J4*(15./8.)*(Rs/((rv-x[i])*1000))^4))^.5 - 4*OmegaM - muM
   omega = 4*((G*M/((rv-x[i])*1000)^3)*(1 + J2*(3./2.)*(R/((rv-x[i])*1000))^2 - J4*(15./8.)*(R/((rv-x[i])*1000))^4))^.5 - 4*OmegaM - muM

 ; 
 ; omega = 1.2957087e-4 



  !p.multi = [0,0,0]

 ; WINDOW, 0, XSIZE=1920, YSIZE=1200, TITLE='Haze'
 ; cgplot,theta/(2*!dpi),A*damp[i]*cos(theta), color=cgcolor('red'), yrange = [-3,3], xtitle = 'Time in vertical periods', ytitle = 'z (Km)', title = 'Haze particle vertical motion relative to the ring', charsize=2.2, background = cgcolor('white'),thick =4

 ; cgLegend, Color=['Red','Blue', 'Orange','Cyan'], Location=[0.2, 0.88],  Titles=['Ring', 'Trajectories after collision with SGW','Trajectories after 1st pass','Trajectories after 2nd pass'], Length=0.075, /Box, VSpace=2.75,/Background, bg_color = 'white', charsize=2.2

  ind2 = 0
  zs = make_array(n_elements(phis),1)
  coeffs = 0


  foreach phi, phis do begin
   
                ;       UNDEFINE, R
    ;;;;;;;;;;;This starts me pasting what used to be a function
                 ;       C = make_array(2)
                       ;SCL theory damping
                    
                  ;      mu = mu + epsilon*k[i] * slope[i] * s
                   
                   ;     B = [[cos(phi), sin(phi)],[-mu*sin(phi), mu*cos(phi)]]
                    ;    R = [A*damp[i]*cos(phi), -omega[0]*A*damp[i]*sin(phi) - epsilon*slope[i]*s]
                       
                     ;   C = Invert(B) ## R ; I verified this is the right operator for this order of opertations
                    
                        
                        ;phase  = int_tabulated(x[0:i],k[0:i])
                       ; phase  = tsumcumul(x[0:i],k[0:i])
                      ;  
                       ; coeffs = [c[0]*cos(phase[-1]) + c[1]*sin(phase[-1])]
   ;;;;;;;;;;;;;;;;;;This ends it
     haze_creation, phi, theta, omega, i, sigma, 1, epsilon, x, upper, slope, damp, k, A = A, coeffs = coeffs, s=s
           

    Zs[ind2] = coeffs
    ind2+=1
  endforeach
  
 ; print, indi

  Zmax[indi] = max(zs)

  Zmin[indi] = min(zs)
  
  rad[indi] = ind
  
  indi=indi+1
  
  endfor

  ;  WINDOW, 0, XSIZE=1820, YSIZE=1100, TITLE='Haze'
  ;  cgplot,x,A*damp*cos(tsumcumul(x,k)),yrange = [-3,3], xrange = [0,150], xtitle = 'Distance from resonance (Km)', ytitle = 'z (Km)', title = 'Haze particle motion relative to the ring', charsize=2.2, background = cgcolor('white'),thick =4,color = cgcolor('red')
  ;  cgoplot,rad,zmax;,yrange = [-3,3], xrange = [75,125], xtitle = 'Distance from resonance (Km)', ytitle = 'z (Km)', title = 'Haze particle motion relative to the ring', charsize=2.2, background = cgcolor('white'),thick =4,color = cgcolor('red')
  ;  cgoplot,rad,zmin;,yrange = [-3,3], xrange = [75,125], xtitle = 'Distance from resonance (Km)', ytitle = 'z (Km)', title = 'Haze particle motion relative to the ring', charsize=2.2, background = cgcolor('white'),thick =4,color = cgcolor('red')

  ;   graf = plot(x,A*damp*cos(tsumcumul(x,k)),yrange = [-3,3], xrange = [0,150], xtitle = 'Distance from resonance (Km)', ytitle = 'z (Km)', title = 'Haze particle motion relative to the ring',thick =4,color ='red', Name = 'Ring')
  ;   graf =plot(rad,zmax, overplot = 1, transparency=70.,fill_background=1,fill_color="blue",fill_transparency=90., fill_level=0,color='cornflower', Name = 'Haze',thick=4)
  ;      leg  = legend(position=[.8,.8], font_size = 16)
  ;      graf =plot(rad,zmax, overplot = 1, fill_level=0,color='white', Name = 'Haze',thick=4)
  ;   graf =plot(rad,zmin, overplot = 1, transparency=70.,fill_background=1,fill_color="blue",fill_transparency=90., fill_level=0,color='white')


  y_haze[0] = zmax
  y_haze[1] = zmin

  return, y_haze



end
