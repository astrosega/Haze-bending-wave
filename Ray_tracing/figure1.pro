;Outline
;A version of saturn_star that is written to generate Figure 1 from Sega et al 2024. with minimum setup requirements. It demostrates the geometry and optics of an UVIS occultation
;
;Requirements
;
;stars_robust.sav                    ---- contains the geometrical parameters of the diffent stars occultated by the rings used in this study
;GamPeg032I_Mimas53BW_400m_res.sav   ---- the occultation to be simulated binned to 400m of radial resolution. These and all the other stellar occultation files used are under a folder named "magnus" in the GitHub repository
;Fresnel_complex.pro -> (in the root of the repositorty). Computes the complex fresnel integral nescessary to draw the SCL wave.
;oplotmimas.pro      -> In this directory. Plot warper.
;func_shugampeg.pro  -> in this directory,
;The Coyote library
;
;Authors
;Daniel Sega

;Change Log:

;9 / 26 / 2024 -> Created

pro figure1

  !p.multi=[0,1,2]
  beta =  0.
  v    =  0.0260
  ; v    =  0.0220

  sigma=  365.7
  W    = 58.5

  name='GamPeg032*'

  restore, 'stars_robust.sav' ;this is the list of occultations
  noccs = n_elements(stars)
  chis  = fltarr(noccs) ; this is an array to store the chi squared for each occultation

  if keyword_set(name) then begin
    stars = stars(where(strmatch(stars, name) EQ 1))
    suba = where(strmatch(stars, name) EQ 1)
    noccs = n_elements(stars)
    chis  = fltarr(noccs) ; this is an array to store the chi squared for each occultation

  endif




  ;Definitions of some constants

  G  = 6.67408e-11      ;SI
  cD = 1.50217913e-7    ;this is the curly D in SCL, (1/s^2)
  mu = 1.2957087e-4     ;(1/s) vertical frequency of particles at resonance
  muM= .771634e-4       ;(1/s) vertical frequency of mimas
  aM = 185539.          ;Km  mimas semi-major axis
  iM = 1.574 * !Dpi/180 ;rads, mimas inclination
  Omega= 1.28897e-4       ;keplerian frequency of ring particles at resonance
  kappa= 1.2822e0-4        ; Linmbad frequency of ring particles at resonance
  ;orbital paramenters from Mimas taken from Nasa webpage


  j=complex(0,1)
  ;initialize your variable array. Yes we have to do it again because each occ is optimized for a particular t-vector

  reso  =   0.001d ;km ;This is the resolution of the x axis, not of the ray-tracing
  ;km    = 300
  ;start = -21
  km=500
  start=-221
  x     = dindgen(round(km*(1/reso)), Increment=reso, start=start) ;radial axis, in Km. t=0 at resonance


  d=0.015; d is the vertical thinckness of the ring in Km

  stars_format = stars
  star_file = file_search(filepath(stars_format+ '*', subdir = ['magnus']))

  ;print,star_file

  ;Open a file with the data already converted into optical depth (tao), and a radius (radius) vector
  Restore, star_file

  ;I'm reading off the parameter of Mimia's orbit

  ; print,star_rev[i_mimas]
  ;print,rev[i_geo]

  phi   =  2.4404356017937312
  B = -0.35403773912744846

  B   = B[0]
  phi = phi[0]

  phase = -2.2164771582762475 ;phase used for the occultation shown in Fig. 1, Sega et al 2024.
  phase = phase[0]
  phase_print = phase
  while (phase_print LT 0) do begin
    phase_print = phase_print + 2*!dpi
  endwhile
  phase_print = phase_print*180/!dpi




  ;print,'phase =',phase_print, '   star = ', star, '  phi = ', phi*180/!dpi


  B = abs(B)
  Beff=abs(atan(tan(B)/cos(phi))) ;this is the effective angle (alpha in Jerousek 2012)

  d_nor = d
  res_int = round(1/reso)

  if Beff GT 0.85 then begin ;0.872          ; for occultations of phi GT 47 degrees, the thickness of the ring changes the optical depth very little
    d_nor = 1.                  ; (smaller than the noise in the data of AlpVir008E). After this the effect of the thickness in the model is that of
    Res_int = 100.                 ; increasing the resolution (see Gresh 1986 et al., their model has this feature as well),
    ;Resolution of intercept            ; so we take a computational shortcut and use low resolution (0.001 Km is low for normal occultations)
  endif
  ;Hamp is the normal optical depth outsiede of the wave for this occultation. We get it with the granola bar model


  if ~keyword_set(sigma) then sigma = 363.                        ;sigma is the surface mass density in  Kg/m^2 ;Also, this is a best fit for gampeg32
  rv    = 131902.0 ;resonant radius (km)


  v=v ;viscosity (m^2/s)
  ; alpha = G*sigma*rv*.4
  ;SCL theory damping

  dampingfactor = (1./3.)*(cD^2)*mu/(2*!Dpi*G*sigma)^3             ;this is in SI (meters. remember t and rv are in km)
  damp          =  Exp(-((x^3/rv^2)*1000.*(dampingfactor)*(v)))    ;the 1000 is to transform x into meters since the viscosity is in meters

  ;This is the slope of the ring, and some of the processes of the simulation depend on it.


  ;The slope of the ray is related to the opening angle by
  s = Tan(beff)

  ;Hsigma in the surface mass density measured in units of 100 g/cm^3 (Hectogram/cm^3) We are a best fit value.
  Hsigma=(sigma/1000.)


  A=(476.92/(Sqrt(Hsigma)))/1000. ;amplitud of the wave, Km

  ;The 476.92 is calculated by estimating the amplitude of the driving force (MIMAS g field). It depends on the driving frequency and
  ;natural frequency as well as some parameters of MIMAS orbit. This value is computed in SCL, but I recomputed it with more recent values. A is in Kilometers. Note that in Gresh et al.
  ;(Radio occultations) the best fitted amplitude differed by this one by a factor of 4 while in SCL (optical) it differed by 1.15

  ;e is a parameter that appears in the solution of the DE and makes the exponent, which is part of the solution, dimensionless. It is related to the wave
  ;number of the bending wave and it depends on G, sigma, and the difference in the driving and natural frequencies. It also depends on the distance
  ;from saturn at resonance. This value of e was calculated by me using the Keplerian approximation (vertical frequency=orbital frequency).
  ; e is defined in SCL

  e=1
  e=double(e)
  e=(2.*!DPI*G*sigma/(rv*1000.*cD)) ;the 1000 is to transform rv into m




  ;Here comes the wave profile at maximum amplitude at resonance (this sets the phase to -Pi/4). The sqrt of pi in the amplitude comes from the fresnel integral, don't panick;
  ;it is not part of the amplitude
  upper= A * damp * (1/(sqrt(!DPI))) * Exp(j*(phase + !DPI/4 +!DPI/2+!DPI/6)) * (Exp(j*x^2/(2*e*rv*rv))*(0.5*sqrt(!DPI/2) - 0.5*j*sqrt(!DPI/2) + fresnel_complex(x/(sqrt(2*e)*rv))))
  lower= A * damp * (1/(sqrt(!DPI))) * Exp(j*(phase + !DPI/4 + !DPI/2+ !DPI/6)) * (Exp(j*x^2/(2*e*rv*rv))*(0.5*sqrt(!DPI/2) - 0.5*j*sqrt(!DPI/2) + fresnel_complex(x/(sqrt(2*e)*rv)))) -d_nor/2.


  slope     = -(shift(upper, 1) - upper)/reso
  slope[0]  = 0
  slope=abs(slope)
  wi,1,wsize=[1900,1080]
  cgplot,-x,upper,xrange=[-200,50],yrange=[-1,1],ytitle="Elevation [km]",thick=4, background=cgcolor('white'),color=0, xthick=2,ythick=2,charsize=3.5,charthick=2, XTICKFORMAT="(A1)",position = [.1,.55,.95,.95]
  ;  oplot,-x,lower,thick=4

  colors = ['Black','Red']
  items = ['Ring', 'Lightrays']
  al_legend,items,linestyle=[0,0],color=colors,charsize=3.5,linsize=0.2,thick=2,position=[-15, -.2],box=1,charthick=2


  ;  y_int = real_part(start * s)
  ray= s*(x + 36); + y_int
  oplot,x,ray,linestyle=0,thick=2,color=1000
  ray= s*(x + 74.5)
  oplot,x,ray,linestyle=0,thick=2,color=1000
  ray= s*(x + 105)
  oplot,x,ray,linestyle=0,thick=2,color=1000
  ray= s*(x + 138.1)
  oplot,x,ray,linestyle=0,thick=2,color=1000

  junk = func_shugampeg()
  rscl = junk[*,0]
  optd = junk[*,1]

  cgplot, rscl-rv, optd,thick=2,xtitle="Position from resonance [km]",ytitle="Optical Depth",background=cgcolor('white'),color=0, xthick=2,ythick=2,charsize=3.5,charthick=2,ymargin=[-.1,-10],position = [.1,.12,.95,.52]

  oplot,x, -(120*s)*(x+37),  color=cgcolor('red'),linestyle=5,thick=2
  oplot,x, -(120*s)*(x+74.3),color=cgcolor('red'),linestyle=5,thick=2
  oplot,x, -(120*s)*(x+105), color=cgcolor('red'),linestyle=5,thick=2
  oplot,x, -(120*s)*(x+137.6), color=cgcolor('red'),linestyle=5,thick=2

  colors = ['Black']
  items = ['Optical depth']
  al_legend,items,linestyle=[0],color=colors,charsize=3.5,linsize=0.2,thick=2,position=[-30, 4],box=1,charthick=2


 ; write_png,'bendingwave.png',tvrd(/true)
end