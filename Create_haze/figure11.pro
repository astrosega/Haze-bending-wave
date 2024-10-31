pro figure11

;Produces the vizualizaton of the haze (figure 11) shown in Sega et al 2024. Requires IDL 8.7.2 or above due to the use of fillplot

;Requirements
;haze_creation.pro   -> (included in this directory) which is the core simulation that computes the trajectories of particles after their collisions with the self-gravity wakes
;Fresnel_complex.pro -> (in the root of the repositorty). Computes the complex fresnel integral nescessary to draw the SCL wave.
;linspace.pro        -> (in the root of the repositorty).
;tsumcumul.pro        -> (in the root of the repositorty).
;
;Change Log
;
;Date created 09-09-2022


;;;;;Constants;;;;;;

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

  epsilon = .5         ;coefficient of resitution

  phis = linspace(0,2*!dpi,6) ;this corresponds to the phases (in a vertical period of a self-gravity wake) in which the collision happen (for instance, a collision at a phase of 0 occurs when the wave is flat and the elevation is maximum). We model six collisions per period.
  sigma  =  363
  Hsigma = (sigma/1000.)
  A      = (476.92d/(Sqrt(Hsigma)))/1000d

  j=complex(0,1)
  
  ;initialize your radial axis.

  reso  =   0.001 ;km ;This is the resolution of the x axis.
  km    = 300
  start = -21
  x     = dindgen(round(km*(1/reso)), Increment=reso, start=start) ;radial axis, in Km. t=0 at resonance
  i = -start/reso + 100/reso
 
  
  s       = .000025 ;constant 'a' in Sega et al 2024. ICARUS

  omega = 4*((G*M/((rv-x[i])*1000)^3)*(1 + J2*(3./2.)*(R/((rv-x[i])*1000))^2 - J4*(15./8.)*(R/((rv-x[i])*1000))^4))^.5 - 4*OmegaM - muM ;This is the vertical driven frequency.



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



  theta = linspace(0, 3*!dpi) ;This corrrespoinds to for how many periods are the collisions computed
  COMMON SHAREname, nhits
  nhits = 1              ; this is the number of times a particle is hit. The particle can get kicked up while returning to the rings, mantaining its momemtum, or just be reabsorbed after a vertical period

  !p.multi = [0,0,0]

  ; WINDOW, 0, XSIZE=1920, YSIZE=1200, TITLE='Haze'
  graph =  plot(theta/(2*!dpi),A*damp[i]*cos(theta), color='red', yrange = [-.5,.5], xtitle = 'Time in vertical periods!X', font_size=32,thick =4,dim=[1820,720],Name='Wave motion',layout=[2,1,1],margin=[.2, .15, .02, .05],ytitle = 'z [km]')

  ind2 = 0
  zs = make_array(n_elements(phis))

  foreach phi, phis do begin

    haze_creation, phi, theta, omega, i, sigma, nhits, epsilon, x, upper, slope, damp, k, A = A, koeffs = koeffs,s=S, plotf=0,c=c , landr=1

    graph=   plot(theta/(2*!dpi), c[0,0]*cos(theta) + c[1,0]*sin(theta),color=cgcolor('red'),overplot=1,Name=['Motion after collision'])
    if phi eq phis[0] then leg  = legend(position=[.34,.93], font_size = 21,linestyle=6)
    ;graph1=  plot(theta/(2*!dpi),c[0,1]*cos(theta) + c[1,1]*sin(theta),color=cgcolor('red'),overplot=1)

  endforeach
  tic


  indi = 0



  phis = linspace(0,2*!dpi,12)  ; see other definition of phis above
  fase = 1.5                    ;phase of the wave modeled

  ri = double(start) +.5
  rf = double(km + start)


  res = .4d

  rad = make_array((rf - ri)/res + 1) ;Creates array to record the haze envelop
  zmax = make_array((rf - ri)/res+ 1) ;This will record the position of the highest haze particle
  zmin = make_array((rf - ri)/res+ 1) ;This will record the position of the lowest haze particle


  for ind = ri, rf, res do begin ;This loop goes to all the considered radial location between ri and rf and computes the results of the collision there.
    i = -start/reso + ind/reso
    omega = 4*((G*M/((rv-x[i])*1000)^3)*(1 + J2*(3./2.)*(R/((rv-x[i])*1000))^2 - J4*(15./8.)*(R/((rv-x[i])*1000))^4))^.5 - 4*OmegaM - muM ;This is the driven vertical freq. inside of the wave (it is the same as the natural vertical freq. at resonance).

    !p.multi = [0,0,0]
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


  upper = (rv/(rv+x))^0.5 * (A/(sqrt(!DPI))) * exp(j*(fase+!dpi/2)) * damp * (Exp(j*x^2/(2*e*rv*rv))*(0.5*sqrt(!DPI/2)-0.5*j*sqrt(!DPI/2) + fresnel_complex(x/(sqrt(2*e)*rv))))

  graf = plot(-x,A*damp*cos(tsumcumul(x,k) - !dpi/4 + fase),yrange = [-1,1], xrange = [-150,0], xtitle = 'Distance from resonance [km]',thick =3,color ='red', Name = 'Ring',dim=[1820,980],linestyle=2,/current,layout=[2,1,2],margin=[.12, .15, .05, .05])
  graf.font_size = 32
  ; graf =plot(rad,zmax, overplot = 1, transparency=70.,fill_background=1,fill_color="blue",fill_transparency=90., fill_level=0,color='cornflower', Name = 'Haze',thick=4)
  graf2 =fillplot(-rad,[[zmax],[zmin]], overplot = 1,thick=1,fill_color="blue",fill_transparency=80.,color ='cornflower',name='Haze')
  leg  = legend(target =[graf,graf2],position=[.973,.93], font_size = 21,linestyle=6);,position=[.66,.93]

  graf = plot(-x,A*damp*cos(tsumcumul(x,k) - !dpi/4 + fase),xrange = [-73.75,-71.75], yrange = [-.67,-0.42], thick =3,color ='red',dim=[1820,980],layout=[2,1,2],linestyle=2,/current,margin=[0.8],position=[0.565,.18,0.7,.37],XTICKFONT_SIZE=0,YTICKFONT_SIZE=0)
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
  junk = min(abs(rad - 25),i)
  junk = min(abs(rad - 125),j)
  
 

 ; print,(a-b)/(rad[where(zmax eq a)]-rad[where(zmin eq b)])

  print,'maximum thickness of the haze [km] =', max(zmax-zmin) ;For epsilon=0,1 I get 0.052 km. For epsilon = 1, 0.32 km. For epsilon=0,5 I get 0.164713 km
  print,'minimum thickness of the haze [km] within the wave =', min(zmax[i:j]-zmin[i:j],k) ;For epsilon=0.5 I get 0.0704 km.
 ; rad2 =rad[i:j] 
 ; print,rad2[k]
 ; fig= PLOT(rad,(zmax-zmin)*1000.,thick =3,xtitle = 'Distance from resonance [km]',ytitle = 'Haze vertical thickness [m]')
 ; fig.font_size = 32

  der = ((shift(zmax,-1)-zmax)/(rad[10]-rad[11]))
  print, 'maximum slope of the haze =',max(abs(der[0:-2]))


end

