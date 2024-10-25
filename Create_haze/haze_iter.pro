pro haze_iter 

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
 ; r = 185539000;
 
  R = 60330000;
  
  omegaM = 0.0000771634
  muM    = 0.000077365
  
 ; omegap = ((G*M/rv^3) (1 + (J2)*(3/2)*(R/rv)^2 - (J4) (15/8)*(R/rv)^4))^.5


phis = linspace(0,2*!dpi,6)
;phis = [0,!dpi/12.,!dpi/10, !dpi/8, !dpi/3. ,!dpi/2.5, !dpi/2]
;phis = [!dpi/3. ,!dpi/2.5, !dpi/2]
;phis = [!dpi/12.,!dpi/2.1]

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
epsilon = .5
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

WINDOW, 0, XSIZE=1920, YSIZE=1200, TITLE='Haze'
cgplot,theta/(2*!dpi),A*damp[i]*cos(theta), color=cgcolor('red'), yrange = [-1,1], xtitle = 'Time in vertical periods', ytitle = 'z (Km)', title = 'Haze particle vertical motion relative to the ring', charsize=2.2, background = cgcolor('white'),thick =4

;cgLegend, Color=['Red','Blue', 'Orange','Cyan'], Location=[0.2, 0.88],  Titles=['Ring', 'Trajectories after collision with SGW','Trajectories after 1st pass','Trajectories after 2nd pass'], Length=0.075, /Box, VSpace=2.75,/Background, bg_color = 'white', charsize=2.2
 cgLegend, Color=['Red','Blue'], Location=[0.2, 0.88],  Titles=['Ring', 'Trajectories after collision with SGW'], Length=0.075, /Box, VSpace=2.75,/Background, bg_color = 'white', charsize=2.2


ind2 = 0
zs = make_array(n_elements(phis))

foreach phi, phis do begin

haze_creation, phi, theta, omega, i, sigma, 2, epsilon, x, upper, slope, damp, k, A = A, coeffs = coeffs,s=S, plotf=1, landr=1

endforeach




end
