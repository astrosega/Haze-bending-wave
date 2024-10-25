function haze_func , upper, damp, res, sigma, fase, x,epsilon,reso


 
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
  
 ; epsilon = 1
 

  ; omegap = ((G*M/rv^3) (1 + (J2)*(3/2)*(R/rv)^2 - (J4) (15/8)*(R/rv)^4))^.5


  phis = linspace(0,2*!dpi,12)
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
 ; x     = dindgen(round(km*(1/reso)), Increment=reso, start=start) ;radial axis, in Km. t=0 at resonance
  
 ; ri = 1d
  ;rf = 150d
  ri = double(start) +.5
  rf = double(km + start)
  
  
  ;res = .2d
  
  rad = make_array((rf - ri)/res + 1)
   zmax = make_array((rf - ri)/res+ 1)
    zmin = make_array((rf - ri)/res+ 1)
    
   ; epsilon = 0.5

    slope     = -(shift(upper, 1) - upper)/reso
    slope[0]  = 0
    
    s = .0000242d

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

  nhits = 2
  ind2 = 0
  zs = make_array(n_elements(phis),nhits)
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
     haze_creation, phi, theta, omega, i, sigma, nhits, epsilon, x, upper, slope, damp, k, A = A, koeffs = coeffs, s=s, fase = fase,landr=1
           

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


  zmax = interpol(zmax,rad,x)
  zmin = interpol(zmin,rad,x)
  
   y_haze = make_array(n_elements(zmax),2)
  
  y_haze[*,0] = zmax
  y_haze[*,1] = zmin

  return, y_haze



end