;Outline
;Version of saturn_haze_new made to take in a string with the star name and revolution number, and to be called as a procedure in the simultanious fit
;routine saturn_magnus. It takes the 7 relevant parameters for the model: Amplitud of hump, position of peak of hump
; standard deviation of hump, the hight of the haze (haze_amp) and the viscosity (v), thickness of ring (d).
;It returns Chi Square for the simulation


;Change Log
;7/20/17 -> Created

;Authors
;Daniel Sega

Pro SATURN_STAR, star, Amp, Peak, Dev, Ratiorad, Haze_amp, v, d, chi, Plot=plot, sigma=sigma


  ;Definitions of some constants

  G   = 6.67408e-11      ;SI
  cD  = 1.50217913e-7    ;this is the curly D in SCL, (1/s^2)
  mu  = 1.2957087e-4     ;(1/s) vertical frequency of particles
  muM = .771634e-4       ;(1/s) vertical frequency of mimas
  aM  = 185539.          ;Km  mimas semi-major axis
  iM  = 1.574 * !Dpi/180 ;rads, mimas inclination
  ;Msat= 5.683 * 10 ^ 26 ;Kg

  ;orbital paramenters from Mimas taken from Nasa webpage


  j=complex(0,1)
  ;initialize your variable array. Yes we have to do it again because each occ is optimized for a particular t-vector

  reso  =   0.001 ;km ;This is the resolution of the x axis, not of the ray-tracing
  km    = 300
  start = -21
  x     = dindgen(round(km*(1/reso)), Increment=reso, start=start) ;radial axis, in Km. t=0 at resonance


  d=d; d is the vertical thinckness of the ring in Km
  
  
  stars_format = Rename_fn(star)
 ; print,stars_format
  star_file = file_search(filepath(stars_format+ '*', subdir = ['magnus']))
  
  ;print,star_file

  ;Open a file with the data already converted into optical depth (tao), and a radius (radius) vector
  Restore, star_file
 
 ;I'm reading off the parameter of Mimia's orbit
  Restore, Filepath('Mimas_lon_z_danielOX.sav', subdir = ['Mimas_data'])
  
  ;I'm reading off the geometry of the occultation
  Restore, 'geometries.sav'
  

  
  i_geo  = where(star EQ rev)
  i_mimas= i_geo + 1
  
 ; print,star_rev[i_mimas]
  ;print,rev[i_geo] 
 
  flip1 = 0
  flip2 = 0 
  if B[i_geo[0]] LT 0 then flip1 = !dpi
  if ((phi[i_geo[0]] LE 90 AND phi[i_geo[0]] GE 0) OR ((phi[i_geo[0]] LE 360 AND phi[i_geo[0]] GE 270))) then flip2 = !dpi 
  
  B   = ( B[i_geo] /  180d)   * (!Dpi)
  phi = (phi[i_geo]/  180d) * (!Dpi)

  B   = B[0]
  phi = phi[0]
    

  
  ;For some occs on this catalog of geometries, B is given in a negative number, which means the ring is fliped upsidedown
  ;for the values used to recreate that occultation. Because of how the code is written, we rather not change the B angle but
  ;add a phase of Pi to the wave (it produces the same result)
    fase = 0
    sign = 1
  if (UP_DOWN[i_mimas[0]] EQ '-') then begin
    sign = -1
    fase = !dpi
  endif                                  ;all of this is because the sign of the velocity of mimas comes into play
  
  phase = -4*(lon_mimas[i_mimas[0]])*!dpi/180d + sign*asin(z_mimas[i_mimas[0]]/(sin(im)*am)) + fase - !dpi/4 + flip1 + flip2;eqn 47a in SCL\
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
               Res_int = 10.                 ; increasing the resolution (see Gresh 1986 et al., their model has this feature as well), 
         ;Resolution of intercept            ; so we take a computational shortcut and use low resolution (0.001 Km is low for normal occultations)
endif
  ;Hamp is the normal optical depth outsiede of the wave for this occultation. We get it with the granola bar model
  


  Hamp  = granolabar3(B,phi)*Sin(B)
  if ~keyword_set(sigma) then sigma = 363.                        ;sigma is the surface mass density in  Kg/m^2 ;Also, this is a best fit for gampeg32
  rv    = 131902.0 ;resonant radius (km)
  

  v=v ;viscosity (m^2/s)

  ;SCL theory damping

  dampingfactor = (1./3.)*(cD^2)*mu/(2*!Dpi*G*sigma)^3             ;this is in SI (meters. remember x and rv are in km)
  damp          =  Exp(-((x^3/rv^2)*1000.*(dampingfactor)*(v)))    ;the 1000 is to transform x into meters


  ;The slope of the ray is related to the opening angle by
  s=Tan(beff)

  ;Hsigma in the surface mass density measured in units of 100 g/cm^3 (Hectogram/cm^3) We are using the constant value used in Gresh et al.
  Hsigma=(sigma/1000.)


  A=(476.92/(Sqrt(Hsigma)))/1000. ;amplitud of the wave

  ;The 476.92 is calculated by estimating the amplitude of the driving force (MIMAS g field). It depends on the driving frequency and
  ;natural frequency as well as some parameters of MIMAS orbit. This value is computed in SCL. A is in Kilometers. Note that in Gresh et al.
  ;(Radio occultations) the best fitted amplitude differed by this one by a factor of 4 while in SCL (optical) it differed by 1.15


  ;e is a parameter that appears in the solution of the DE and is is dimensionless. It is related to the wave
  ;number of the bending wave and it depends on G, sigma, and the difference in the driving and natural frequencies. It also depends on the distance
  ; e is defined in SCL eqn 25

  e=1
  e=double(e)
  e=(2.*!DPI*G*sigma/(rv*1000.*cD)) ;the 1000 is to transform rv into m. e is dimensionless

  ;Here comes the wave profile at maximum amplitude at resonance (this sets the phase to -Pi/4). The sqrt of pi in the amplitude comess from the fresnel integral, don't panick.
  ;it is not part of the amplitude

  upper=(rv/(rv+x))^0.5 * (A/(sqrt(!DPI))) * damp *Exp(j*(phase + !DPI/4))* (Exp(j*x^2/(2*e*rv*rv))*(0.5*sqrt(!DPI/2)-0.5*j*sqrt(!DPI/2) + fresnel_complex(x/(sqrt(2*e)*rv)))) +d_nor/2
  lower=(rv/(rv+x))^0.5 * (A/(sqrt(!DPI))) * damp *Exp(j*(phase + !DPI/4))* (Exp(j*x^2/(2*e*rv*rv))*(0.5*sqrt(!DPI/2)-0.5*j*sqrt(!DPI/2) + fresnel_complex(x/(sqrt(2*e)*rv)))) -d_nor/2


  ;Initialize the arrays in which you're going to put your results for the optical depth.
  ;How long is the array depends on the geometry, assuming we always simulate light rays covertin 200 km from resonance

 ; plot,x,upper,xrange=[200,-10],yrange=[-10,10],xtitle="Position from resonance (Km)",ytitle="Elevation from the ring's plane (Km)"
 ; oplot,x,lower


  n       = 200
  step    = 0.4
  points_in_model = round(n/step) ;number of rays in the model
  optd=FLTARR(round(n * (1/step)) + 1)

  ;initialize the array where you are going to dump your results
  ;and the r value (radius from Saturn's center) for each result (we take the x(radious - rv) where the light intersects the upper ring
  r= -dindgen(points_in_model, increment=step, start=start) +rv
  ;initialize all the indices

  i= 0d ;this is so that the programs finds where the ring intercepts the light ray quickly

  m = 0
  m = long(m)
  p = 0
  p = long(p)

  ;This are for security
  Noti=0l
  Notm=0l
  r1=0

  y_int =  start * s

  FOR p=0.0, n, step do begin


    ;next, plot the ray
    ray= -s*(x -p) + y_int
       ; oplot,x,ray


    ;set all the while loops to false
    upperint=0
    lowerint=0
    uppermiss=0
    ;this is to exit the loop
    Noti=i
   

    while (upperint EQ 0) do begin
      ;print,i
      ;   print,upper[i]
      ;   print,ray[i]
      IF (round(upper[i]*res_int) EQ round(ray[i]*res_int)) THEN BEGIN

        ;print,x[i]

        upperint=1
        ;in case you miss the next one
        rmiss=r1

        r1=x[i]

      endif  Else begin

        ;If you don't find the ring get out of the loop
        IF (i GT (Noti+round((1/reso))*10)) then begin ;if you go 10 km and you still don't find the ring, then

          upperint=1
          uppermiss=1
          i=Noti
          ;print,"I missed upper"
        endif
      endelse

      i++
    endwhile


    IF (m EQ 0) then begin
      m=i
      r2=0
    endif
    ;To avoid infinite loops
    Notm=m
    ;;;;;;;;;;;;;;

    if (uppermiss EQ 0) then begin
      m=i

      while (lowerint EQ 0) do begin


        ;      print,m
        IF (round(lower[m]*res_int) EQ round(ray[m]*res_int)) THEN BEGIN

          ;print,m

          lowerint=1

          r2=x[m]



        endif else begin
          if (m GT (Notm+round((1/reso))*10))Then begin

            lowerint=1

            m=Notm

            r1=rmiss

            ;print,"I missed lower"
          endif

        endelse

        m++
      endwhile

    endif


    ;print,'hello'

    ;Here we define the limits of the Haze
    haze_amp = haze_amp   ;this is how many times the amplitud the haze extends to in the +z direction (it is normally 1)

    h_haze = 2. * haze_amp * A
    ;h haze is not the altitud of the haze but the vertical distance the light transcurs in it.
    ;By not substracting the thickness of the ring what I'm doing is overlapping the haze and the ring when calculating the optical depth
    ;That's not right

    r_haze1 = (-haze_amp * A  + y_int)/s + p

    r_haze2=  ( haze_amp * A  + y_int)/s + p

    ; print,r_haze1
    ; print,r_haze2

    ;this is the fraction of the particles fragmented
    ;f=1./100.   ; for dps I used

    ;this is the ratio of the radius. It is the whole particle over fragmented particle so it is always greater than 1. This is also the cuberoot of the number
    ;of fragments in which each particle shatters. We use 1000. because a min in the outer a ring in 1 cm and A max is 20 cm. We do 10 m / 10 cm (there's probably (pun)
    ;a good stadistical way of doing this)

    ratiorad =  ratiorad    ; for dsp I used

    denominator = (ratiorad - 1.)

    ;This are the properties of the surface mass density

    peak=peak
    dev  = dev
    Amp  = amp


    optd[uint(round(p/step))] = Hamp * (ratiorad/denominator) * (1/h_haze) * abs(1/(Cos(phi)*Cos(B)))  * sqrt(!DPi/2) * dev * Amp * (erf((r_haze2 - peak)/(sqrt(2) * dev)) - erf((r_haze1 - peak)/(sqrt(2) * dev))) $
      +Hamp * abs(1/(d_nor*Cos(B)*Cos(phi))) *(r2-r1)  - Hamp * sqrt(!Dpi/2) * abs(1/(Cos(B)*Cos(phi)*d_nor*denominator)) * dev * Amp * (erf((r2  -   peak)/(sqrt(2) * dev)) - erf((r1    - peak)/(sqrt(2) * dev))) ;$
    ;+(ratiorad/denominator) * (1/h_haze) * (1/Cos(Beff))  * sqrt(!DPi/2) * dev * Amp * (erf((r_haze2  -   peak)/(sqrt(2) * dev)) - erf((r2-peak)/(sqrt(2) * dev)))
;print, Hamp * (1/(d_nor*Cos(Beff))) *(r2-r1)

  endfor
  ;Now that we've ended the for loop, we can calculate chi square or plot the results
  
  tao[-1] = tao[-2] ;To do this we have to smooth out the endpoints of out data. The endpoints are problematic because of the binning process
  tao[0]  = tao[1]  ;Because these points are not to be compared with the model, I can just interpolate their values.
  
  optd[0] = optd[1] ;Since the simulation is an iterative process the first value is generally off and I interpolate it (altenatively I could also remove it. Probablu what I should do)


  if keyword_set(plot) then begin

    graphic3   = plot(r, optd, xrange=[131650.,132000.], yrange=[0.,10.], thick=2, min_value=0.1)
    graphic4   = oplotmimas(radius,data,star='+ Haze ('+stars_format+')',tao=tao, title=1, angles = [B, phi], phase = phase_print)

  endif
 

  taosigma[-1] = taosigma[-2]
  taosigma[0]  = taosigma[1]
 ; print, taosigma


  chi = goodfit2(radius, tao, r, optd, points_in_model, taosigma=taosigma)
  ;print, 'chi_gampeg032 =', chi_gampeg



end
