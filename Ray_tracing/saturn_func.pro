;Outline
;A version of saturn_star that is written as a function to be minimized with amoeba it takes the constant beta (beta), a viscosity (v), and a surface mass density (sigma).
; 
;Inputs
; v     -> kinematic viscosity of the rings in wave region.
; sigma -> surface mass density of the A-ring in te wave region
; d     -> thickness of the rings
; beta  -> proportionality constant between the torque on the self-gravity wakes and the optical depth due to the haze 

;Output
;chi^2 -> a sum of chi^2 for all the stars in the data set we want to fit.

;Example
;IDL> saturn_func([1.48, 0.0623., 370],plot=1) epsilon1 phase theory
;IDL> saturn_func([1.39, 0.0576., 367],plot=1) epsilon1 phase obs   (1.3905730     0.057039060       366.76825) if magnus2 is used (the other dataset that does substract the background from Irr) chi¨2 4.5189223, if the sqrt(slope) is used in haze_creation, then 4.5171900
;IDL> saturn_func([1.37, 0.057., 366],plot=1) epsilon0.5 phase obs

;IDL> IDL> amoeba(0.001 ,function_name = 'saturn_func', scale=[.1,0.005, 1.], P0 = [3.2, 0.054, 364], function_value = f)
;IDL> IDL> amoeba(0.001 ,function_name = 'saturn_func', scale=[.1,0.005, 1.], P0 = [2.98, 0.056., 365], function_value = f)
;IDL> IDL> amoeba(0.001 ,function_name = 'saturn_func', scale=[.1,0.005, 1.], P0 = [1.53, 0.037., 365], function_value = f)
;IDL> amoeba(0.001 ,function_name = 'saturn_func', scale=[.1,0.005, 1.], P0 = [1.53, 0.037., 365], function_value = f)

;Authors
;Daniel Sega

;Change Log:

;9 / 26 / 2020 -> Created
;6 /  8 / 2022 -> added create haze keyword where the haze layer is created from first principles. 

function SATURN_FUNC, P0, name=name, plots = plots, isotro = isotro, highb = HighB, fit_rv = fit_rv, create_haze = create_haze
 ;tic
 ;fit_rv = 0
 fit_phase = 1
 create_haze = 1
 isotro=1
 epsilon=1
 
beta =  p0[0]
v    =  P0[1]
sigma=  P0[2]
W    = 58.5

restore, 'stars_robust.sav' ;this is the list of occultations
if keyword_set(highb) then restore, 'stars_robust.sav'

;stars= stars[-20:-1]
;B    = B[-20:-1]
;phi= phi[-20:-1]
;  lon_mimas = lon_mimas[-20:-1]
;  up_down = up_down[-20:-1]
 ; z_mimas = z_mimas[-20:-1]

;stars= stars[3:-1]
;B    = B[3:-1]
;phi= phi[3:-1]
;  lon_mimas = lon_mimas[3:-1]
;  up_down = up_down[3:-1]
;  z_mimas = z_mimas[3:-1]


noccs = n_elements(stars)
chis  = fltarr(noccs) ; this is an array to store the chi squared for each occultation  

if keyword_set(name) then begin
;  stars = stars(where(strmatch(stars, name) EQ 1))
  suba = where(strmatch(stars, name) EQ 1)
   stars = stars[suba]
  noccs = n_elements(stars)
  chis  = fltarr(noccs) ; this is an array to store the chi squared for each occultation
  
  B = B[suba]
  phi = phi[suba]
 ; stars = stars[suba]
  lon_mimas = lon_mimas[suba]
  up_down = up_down[suba]
  z_mimas = z_mimas[suba]
  
  B_copy = B
  phi_copy=phi

endif

B_copy = B
phi_copy=phi
phases = make_array(noccs)
theory = make_array(noccs)


for star_i=0, noccs-1 do begin
  
  beta =  p0[0]
   y_add = 0
B = b_copy
phi = phi_copy

  
    ;Definitions of some constants
  
    G  = 6.67408e-11      ;SI
    cD = 1.50217913e-7    ;this is the curly D in SCL, (1/s^2)
    mu = 1.2957087e-4     ;(1/s) vertical frequency of particles at resonance
    muM= .771634e-4       ;(1/s) vertical frequency of mimas
    aM = 185539.          ;m  mimas semi-major axis
    iM = 1.574 * !Dpi/180 ;rads, mimas inclination
    A_mimas = (sin(im)*am) +1.5 ;The 1.5 is to match the maximum measured elevation from the SPICE data and prevent a crash for asin()
    Omega= 1.28897e-4       ;keplerian frequency of ring particles at resonance
    kappa= 1.2822e0-4        ; Linmbad frequency of ring particles at resonance
    ;orbital paramenters from Mimas taken from Nasa webpage
  
  
    j=complex(0,1)
    ;initialize your variable array. Yes we have to do it again because each occ is optimized for a particular t-vector
  
    reso  =   0.001d ;km ;This is the resolution of the x axis, not of the ray-tracing
    km    = 300
    start = -21
    x     = dindgen(round(km*(1/reso)), Increment=reso, start=start) ;radial axis, in Km. t=0 at resonance
  
  
    d=0.015
   ; d is the vertical thinckness of the ring in Km
  
    stars_format = stars[star_i]
    ;stars_format = Rename_fn(stars[star_i])
    
    star_file = file_search(filepath(stars_format+ '*', root_dir = ['magnus2']))
  
    ;print,star_file
  
    ;Open a file with the data already converted into optical depth (tao), and a radius (radius) vector
    Restore, star_file
   
    i_geo  = star_i
    i_mimas= star_i
    
 

  
     ;print, n_elements()
    ;print,rev[i_geo]
  vuelta = 0
    flip1 = 0
    flip2 = 0
    fase = 0
    
  ;  B[i_geo] = abs(b[i_geo])
    if B[i_geo[0]] LT 0 then begin     
     if (phi[i_geo[0]] LE 270 and phi[i_geo[0]] GE 90) then begin
       flip1 = !dpi
       endif         
    endif
    
     if  B[i_geo[0]] gT 0 then begin
                  if (phi[i_geo[0]] gE 270 or phi[i_geo[0]] lE 90) then begin 
            flip1 = !dpi
       endif
    
    ; if (phi[i_geo[0]] lt 180) then lon_mimas[i_mimas] = -lon_mimas[i_mimas]

    
          ;phi[i_geo] = phi[i_geo] ;+ !dpi 
    endif


    B   = ( B[i_geo] /  180d)   *(!Dpi)
    phi = (phi[i_geo]/  180d) * (!Dpi)
  
    B   = B[0]
    phi = phi[0]
  
  
  
    ;For some occs on this catalog of geometries, B is given in a negative number, which means the ring is fliped upsidedown
    ;for the values used to recreate that occultation. Because of how the code is written, we rather not change the B angle but
    ;add a phase of Pi to the wave (it produces the same result)
    fase = 0
    sign = 1
    if (UP_DOWN[i_mimas[0]] EQ -1) then begin
      sign = -1
      fase = !dpi
      
    endif                                  ;all of this is because the sign of the velocity of mimas comes into play
  
 ;  if z_mimas[i_mimas[0]] lt 0 and UP_DOWN[i_mimas[0]] EQ 1) then fase2 = 2*!dpi
  
    phase = 4*(360 -lon_mimas[i_mimas[0]])*!dpi/180d + sign*asin(z_mimas[i_mimas[0]]/A_mimas) + fase - !dpi/4 + flip1;eqn 47a in SCL\
    phase = phase[0]
    phases[star_i] = phase
   
    while (phase LT 0) do begin
      phase = phase + 2*!dpi
    endwhile
  
    B = abs(B)
    Beff=abs(atan(tan(B)/cos(phi))) ;this is the effective angle (alpha in Jerousek 2012)
  
    d_nor = d
    res_int = round(1/reso)
    
  
    if Beff GT 1. then begin ;0.872          ; for occultations of phi GT 47 degrees, the thickness of the ring changes the optical depth very little
      d_nor = 10d    ;100             ; (smaller than the noise in the data of AlpVir008E). After this the effect of the thickness in the model is that of 
      Res_int = 10d    ;2 '50         ; increasing the resolution (see Gresh 1986 et al., their model has this feature as well),
      ;Resolution of intercept            ; so we take a computational shortcut and use lower resolution (0.001 Km is too much resoution for normal occultations)
   
    endif
  
  
    if ~keyword_set(sigma) then sigma = 363.                        ;sigma is the surface mass density in  Kg/m^2 ;Also, this is a best fit for gampeg32
    
    if keyword_set(fit_rv) then begin
       COMMON SHAREname, nombre
       nombre = stars[star_i]
       
      rv = amoeba(0.01 ,function_name = 'saturn_rv', scale=[.5], P0 = [131902.0], function_value = f)      
      rv = rv[0]
    endif else begin
    rv    = 131902.0 ;resonant radius (km)
  endelse
  
  if keyword_set(fit_phase) then begin
    COMMON SHAREname, nombre
    restore, 'phases1.sav'
    
    if keyword_set(name) then begin
      stars = stars[suba]
      theory = theory[suba]
      phases = phases[suba]
    endif
    nombre = stars[star_i]
   theory[star_i] = phase
  ;  phase = amoeba(0.01 ,function_name = 'saturn_phase', scale=[.5], P0 = [phase], function_value = f)
    ;phase = amoeba(0.01 ,function_name = 'saturn_phase', scale=[.03], P0 = [phase], function_value = f)

    ;phases[star_i] = phase[0]
    phase = phases[star_i]
  ; print,'fitted =',phase
  ; print,(theory[star_i]-phase)*180/!dpi
  endif
  
    v=v ;viscosity (m^2/s)
    ; alpha = G*sigma*rv*.4
    ;SCL theory damping
  ;print,'rv =',rv
    dampingfactor = (1./3.)*(cD^2)*mu/(2*!Dpi*G*sigma)^3             ;this is in SI (meters. remember t and rv are in km)
    damp          =  Exp(-((x^3/rv^2)*1000.*(dampingfactor)*(v)))    ;the 1000 is to transform x into meters since the viscosity is in meters
  
    ;This is the slope of the ring, and some of the processes of the simulation depend on it.
  
  
    ;The slope of the ray is related to the opening angle by
    s = Tan(beff) ; Oh no. Reply: Beff is okay
  
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
  
    e=(2d*!DPI*G*sigma/(rv*1000d*cD)) ;the 1000 is to transform rv into m
 
  
    ;Calculate the damping due to the torque on the self_gravity wakes
  
    if keyword_set(create_haze) then begin
          res = .4d
  
      upper1= A * damp * (1./(sqrt(!DPI))) * Exp(j*(phase + !DPI/4.)) * (Exp(j*x^2/(2.*e*rv*rv))*(0.5*sqrt(!DPI/2.)-0.5*j*sqrt(!DPI/2) + fresnel_complex(x/(sqrt(2.*e)*rv))))
      yhaze = haze_func(upper1, damp, res, sigma, phase, x, epsilon)
      
      if stars[star_i] eq 'BetCen077I' then begin
      ; if Beff GT 0.85 then begin
        yhaze[*,0]+= d_nor
        yhaze[*,1]-= d_nor
      ;  beta = beta/(2*d_n)
        
      endif
     
   ;  if stars[star_i] eq 'BetCen077I' then begin
   ;    yadd = 10.
   ;    yhaze[*,0] = yhaze[*,0] + yadd
   ;    yhaze[*,1] = yhaze[*,1] - yadd
   ;  endif

      ;save, yhaze, filename = 'yhaze.sav'
    endif
    
    ;Here comes the wave profile at maximum amplitude at resonance (this sets the phase to -Pi/4). The sqrt of pi in the amplitude comess from the fresnel integral, don't panick.
    ;it is not part of the amplitude
    upper= A * damp * (1/(sqrt(!DPI))) * Exp(j*(phase + !DPI/4)) * (Exp(j*x^2/(2*e*rv*rv))*(0.5*sqrt(!DPI/2) - 0.5*j*sqrt(!DPI/2) + fresnel_complex(x/(sqrt(2*e)*rv)))) +d_nor/2.
    lower= A * damp * (1/(sqrt(!DPI))) * Exp(j*(phase + !DPI/4)) * (Exp(j*x^2/(2*e*rv*rv))*(0.5*sqrt(!DPI/2) - 0.5*j*sqrt(!DPI/2) + fresnel_complex(x/(sqrt(2*e)*rv)))) -d_nor/2.
 
  
  
    slope     = -(shift(upper, 1) - upper)/reso
    slope[0]  = 0
    slope=abs(slope)
  
   if keyword_set(isotro) then Hamp  =  meantao(radius,tao,tao=1)*Sin(B) else $
    Hamp  =  meantao(radius,tao,tao=1)*Sin(Beff)
  
    
    ;Initialize the arrays in which you're going to put your results for the optical depth.
    ;How long is the array depends on the geometry, assuming we always simulate light rays covertin 200 km from resonance  
  
    n       = 200
    step    = 0.4
    points_in_model = round(n/step) ;number of rays in the model
    optd=FLTARR(round(n * (1/step)) + 1)
    ;initialize the array where you are going to dump your results
    ;and the x value for each result (we take the x(radious - rv) where the light intersects the upper ring
    r= -dindgen(points_in_model, increment=step, start=start) +rv
    ;initialize all the indices
  
    i= 0d ;this is so that the programs find where the ring intercepts the light ray quickly
  
    m = 0
    m = long(m)
    p = 0
    p = long(p)
  
    ;This are for security
    Noti=0l
    Notm=0l
    r1=0
    tao_slope0= 0l
  
  
    y_int = real_part(start * s)
  
    FOR p=0.0, n, step do begin
  
  
      ;next, plot the ray
      ray= -s*(x -p) + y_int
  ;       oplot,x,ray
  
  
      ;set all the while loops to false
      upperint=0
      lowerint=0
      uppermiss=0
      ;this is to exit the loop
      Noti=i
  
      tao_slope = make_array(n_elements(x), value=0.)
      ;tao_slope= 0
  
      while (upperint EQ 0) do begin
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
            tao_slope =  tao_slope0
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
  
          tao_slope[m] = slope[m] * beta * reso
          ;print, m
  
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
  
              tao_slope =  tao_slope0
  
            endif
  
  
          endelse
  
  
  
          m++
        endwhile
      endif
  
      tao_slope0 = tao_slope
  
      ;print,'hello',m-i
  
      ;Here we define the limits of the Haze
      haze_amp = 0.8 ;+ y_add /A  ;this is how many times the amplitud the haze extends to in the +z direction (it is normally 1)
        
      h_haze = 2. * haze_amp * A ;+ y_add
      ;h haze is not the altitude of the haze 
      r_haze1 = (-haze_amp * A  + y_int)/s + p
  
      r_haze2=  ( haze_amp * A  + y_int)/s + p
      
        r_haze1_i = floor(mean(where(round(x*1000)/1000. EQ round(r_haze1*1000)/1000.)))
        r_haze2_i = floor(mean(where(round(x*1000)/1000. EQ round(r_haze2*1000)/1000.)))
 if ~keyword_set(create_haze) then begin
  yhaze = make_array(n_elements(x),2)
  yhaze[*,1]=haze_amp*A
  yhaze[*,0]=-haze_amp*A
 endif


      
  if keyword_set(create_haze) then begin
              
      junk    = min(abs(-s*(x-p) + y_int - yhaze[*,0]), index)
      r_haze1_i =  index

      junk = min(abs(-s*(x-p) + y_int  - yhaze[*,1]), index)
      if index lt r_haze1_i then begin
        r_haze2_i = r_haze1_i 
        r_haze1_i = index
        endif else r_haze2_i = index
endif   

  
  
      if (r_haze1 LT x[0]) or (r_haze1_i eq r_haze2_i) then begin
        optd_haze   = 0
        optd_minus  = 0
  
      endif else begin 

  
  if keyword_set(isotro) then   optd_haze  = INT_TABULATED(x[r_haze1_i : r_haze2_i], abs((1/((yhaze[r_haze1_i:r_haze2_i,1]-yhaze[r_haze1_i:r_haze2_i,0])*Cos(B)*Cos(phi))) * beta * abs(slope[r_haze1_i : r_haze2_i]))) else $
          optd_haze  = int_tabulated(x[r_haze1_i : r_haze2_i],  abs((hamp/((yhaze[r_haze1_i:r_haze2_i,1]-yhaze[r_haze1_i:r_haze2_i,0])*Cos(Beff))) * beta* abs(slope[r_haze1_i : r_haze2_i])))
        
         
        if (r1 EQ r2) OR (r2 EQ d_nor/tan(beff)) then begin
          optd_minus = 0
          r2=d_nor/tan(beff)
          r1=0
        endif
      endelse
      optd_minus=0
         if keyword_set(isotro) then optd[uint(round(p/step))] = optd_haze  + Hamp * abs(1./(d_nor*Cos(B)*Cos(phi))) * (r2-r1) else optd[uint(round(p/step))] = optd_haze + Hamp * abs(1./(d_nor*Cos(Beff))) * (r2-r1)

     if stars[star_i] eq 'BetCen077I' then optd[uint(round(p/step))] = optd_haze + Hamp 
    endfor
    ;Now that we've ended the for loop, we can calculate chi square or plot the results
    
    if stars[star_i] eq 'BetCen077I' then begin
      optd[where(optd eq 0)] =  Hamp 
      optd[where(optd gt 5)] =  Hamp
      endif 
  
    tao[-3] = tao[-4] ;To do this we have to smooth out the endpoints of out data. The endpoints are problematic because of the binning process
    tao[-2]  = tao[-3]  ;Because these points are not to be compared with the model, I can just interpolate their values.
    tao[-1] = tao[-2] ;To do this we have to smooth out the endpoints of out data. The endpoints are problematic because of the binning process
    tao[0]  = tao[1]  ;Because these points are not to be compared with the model, I can just interpolate their values.
  
    Irr[0] = Irr[1]
    Irr[-1] = Irr[-2]
  
    optd[7] = optd[8] ;Since the simulation is an iterative process the first values are generally off and I interpolate it (altenatively I could also remove it) Note the points at the endges of the model are not used in the fit but I still want to smooth the edges
    optd[6] = optd[7] 
    optd[5] = optd[6] 
    optd[4] = optd[5]
    optd[3] = optd[4]
    optd[2] = optd[3]   
    optd[1] = optd[2] 
    optd[0] = optd[1] 
    
    optd[-5] = optd[-2]
    optd[-4] = optd[-5]
    optd[-3] = optd[-4]
    optd[-2] = optd[-3]
    optd[-1] = optd[-2]
 ; joshplot=1
    if keyword_set(plots) then begin

      modelplot   = plot(r-rv, optd, xrange=[-200.,50.], yrange=[0.,5.], thick=2, name="This Work",linestyle=0, dimensions = [1366,768], buffer= 1)
          
      graphic4   = oplotmimas(radius,data,star='('+stars_format+')',tao=tao, title=1, uperror=taoplus,downerror=taominus,modelplot=modelplot);,scl=1);, error=taosigma) ;, angles = [B, phi], phase = phase_print)
  
      graphic4   = plot(r-rv, optd, xrange=[-200.,50.], yrange=[0.,5.], thick=2, name="New Model",overplot=1,linestyle=0, dimensions = [1366,768])

      graphic4.save, stars_format+'.png'
  
    endif
    
   if keyword_set(joshplot) then begin
      optional2   = plot(-x[0:where(x eq 180)], slope*2 + 0.6, xrange=[-200.,50.], yrange=[0.,2.], thick=4, name="Slope",$
        linestyle=4, dimensions = [1820,720],AXIS_STYLE=1,color='g',layout=[2,1,2],margin=[0, .15, .25, .05],font_size=28, xtitle='Distance from resonance [Km]',ytickformat="(A1)",xtickvalues=[-150,-100,-50,0])
        
      optional1   = plot(-x[0:where(x eq 180)], Abs(upper1[0:where(x eq 180)]) + .6, xrange=[-200.,50.], yrange=[0.,2.], thick=4, name="Amplitude",overplot=1,linestyle=5,color='b')
      a=axis('Y', location=[95,00], tickdir=1,textpos=1,tickvalues=[0,1], target=nop,COORD_TRANSFORM=[0, .6/1], TICKFONT_SIZE=20 ,AXIS_RANGE=[0,3],title='Amplitude [Km]',gridstyle=5,color='b',subgridstyle=5)
      a=axis('Y', location=[40,0], tickdir=1,textpos=1,tickvalues=[0,0.3], target=nop,COORD_TRANSFORM=[0, .2/1], TICKFONT_SIZE=20 ,AXIS_RANGE=[0,3],title='Slope',gridstyle=3,color='g',subgridstyle=3)
       graphic  = oplotmimas(radius,data,star='('+stars_format+')',tao=tao, title=1,uperror=taoplus,downerror=taominus,layout=1,modelplot=optional2, optional=optional1);, error=taosigma) ;, angles = [B, phi], phase = phase_print)
     graphic4   = plot(r-rv, optd, xrange=[-200.,50.], yrange=[0.,2.], thick=4, name="This Work",overplot=1,linestyle=0, dimensions = [1366,768])
     
     restore, "../scl.sav"
     newm= plot(rscl-rv, optd,thick=3, color='g', name='SCL Theory',linestyle=5,overplot=1)


      leg  = legend(target = [graphic4,newm, graphic],position=[.23,.93], font_size = 20)
     ; a=axis('Y', location=[75,0], tickdir=1,textpos=1,tickvalues=[0,1], target=nop,COORD_TRANSFORM=[0, .6/1], TICKFONT_SIZE=18 ,AXIS_RANGE=[0,3],title='Amplitude [Km]',gridstyle=5,color='b',subgridstyle=5)
     ; a=axis('Y', location=[50,0], tickdir=1,textpos=1,tickvalues=[0,0.3], target=nop,COORD_TRANSFORM=[0, .2/1], TICKFONT_SIZE=18 ,AXIS_RANGE=[0,3],title='Slope',gridstyle=3,color='g',subgridstyle=3)
      graphic4.save, stars_format+'.png'
    endif
      chi = goodfit3(radius, Irr/I0, r, exp(-optd), points_in_model, stars_format, taosigma = Sqrt(Irr)/I0)

chis[star_i] = chi
;print, stars[star_i],chi
endfor

Chi_Total = total(chis)
print,beta, chi_total/(498.*59. + 497), keyword_set(fit_phase), keyword_set(isotro),epsilon

save,theory, phases,stars, filename = 'phases2.sav'

;toc
Return, Chi_total/(498.*59. + 497) ;these numbers are the degrees of freedom. There is one occultation EpsCas104E where we had to remove a point and hence it has one less degree of freedom

end