;Outline
;A version of saturn_star that is written to plot 6 other occultations of different geometries. The geometries correspond with the subsets of the UVIS dataser derfined in Sega 2024. The raytracing for all the geometries is done in the code. Calls on stars_robust.sav and 
;phases1.sav (which is a table for the oberved phases of the wave for each occulation) and the "magnus" directory.
;
;;Requirements
;It calls on func_shugampeg.pro,
;stars_robust.sav (set to look for it in the root directory for your idl instance) and
;phases1.sav (which is a table for the oberved phases of the wave for each occulation) (set to look for it in the root directory for your idl instance)
;goodfit3,pro   --> in this directory
;oplotmimas_fig13.pro --> in this directory
;haze_func.pro --> in the Create_haze directory
;linspace.pro  --> in the root of the main branch it github
;tsum.pro       --> ''
;and the "magnus" directory.(set to look for it in the root directory for your idl instance)
;
;Log.
;Jan 2024 -> Created. Daniel Sega

pro Figure13;, P0, name=name, plots = plots, isotro = isotro, highb = HighB, fit_rv = fit_rv, create_haze = create_haze


  junk = plot([1,1],/nodata, dimensions = [1820,1020],axis_style=0,buffer=1)
  ;tic
  ;fit_rv = 0
  fit_phase = 1
  create_haze = 1
  epsilon=1
  isotro = 1
  name=1

  W    = 58.5

  restore, 'stars_robust.sav' ;this is the list of occultations
  if keyword_set(highb) then restore, 'stars_robust.sav'

  suba = [where(strmatch(stars, 'EpsCas104E*') eq 1),where(strmatch(stars, 'EtaLup034E*') eq 1),where(strmatch(stars, 'KapCMa168I*') eq 1),where(strmatch(stars, 'BetCen077I*') eq 1),where(strmatch(stars, 'KapCen035E*') eq 1),where(strmatch(stars, 'AlpVir232E*') eq 1)];,where(strmatch(stars, 'GamPeg032*') eq 1)]

  ;  stars = stars(where(strmatch(stars, name) EQ 1))    noccs = n_elements(stars)
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



  B_copy = B
  phi_copy=phi
  phases = make_array(noccs)
  theory = make_array(noccs)

  graphs = make_array(12)


  for star_i=0, noccs-1 do begin


    P0=[1.39899, 0.058184460, 365.33575]

    beta =  p0[0]
    v    =  P0[1]
    sigma=  P0[2]

    beta =  p0[0]
    y_add = 0
    B = b_copy
    phi = phi_copy
    isotro=1


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


    d=0.03
    ; d is the vertical thinckness of the ring in Km

    stars_format = stars[star_i]
    ;stars_format = Rename_fn(stars[star_i])

    star_file = file_search(filepath(stars_format+ '*', root_dir = ['magnus']))

    ;print,star_file

    ;Open a file with the data already converted into optical depth (tao), and a radius (radius) vector
    Restore, star_file
    i_geo  = star_i
    i_mimas= star_i

    vuelta = 0
    flip1 = 0
    flip2 = 0
    fase = 0

    ;  B[i_geo] = abs(b[i_geo])
    if B[i_geo[0]] LT 0 then begin
      ;    up_down[i_mimas] = -1*up_down[i_mimas]
      ;   z_mimas[i_mimas] = -z_mimas[i_mimas]
      ; lon_mimas[i_mimas] = -lon_mimas[i_mimas]
      ;   flip1 = !dpi
      ;   B[i_geo[0]] = abs(b[i_geo[0]])
      ; lon_mimas[i_mimas] = lon_mimas[i_mimas] * (-1)
      ; flip1 = !dpi
      ;      vuelta = 1
      ;if ((phi[i_geo[0]] LE 90) OR ((phi[i_geo[0]] LE 360 AND phi[i_geo[0]] GE 270)))  then begin
      ;             strings = stars[star_i]
      ;      if strmatch(strings.Substring(-1),'E') then flip2 = !dpi
      if (phi[i_geo[0]] LE 270 and phi[i_geo[0]] GE 90) then begin
        ; print,  stars[star_i]
        flip1 = !dpi
      endif
      ;  endif
      ;  if (phi[i_geo[0]] GT 180) then up_down[i_mimas] = -1*up_down[i_mimas]
      ;   if (phi[i_geo[0]] LE 270 AND phi[i_geo[0]] GE 90) then up_down[i_mimas] = -1*up_down[i_mimas]

    endif

    if  B[i_geo[0]] gT 0 then begin

      if (phi[i_geo[0]] gE 270 or phi[i_geo[0]] lE 90) then begin
        flip1 = !dpi

      endif


    endif




    B   = ( B[i_geo] /  180d)   *(!Dpi)
    phi = (phi[i_geo]/  180d) * (!Dpi)

    B   = B[0]
    phi = phi[0]
    Beff=abs(atan(tan(abs(B))/cos(phi))) ;this is the effective angle (alpha in Jerousek 2012)


    angles =[B,phi,Beff]


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
    ; print,stars[star_i],UP_DOWN[i_mimas[0]]

    ;print,'hrllo'     ;phase = phase[star_i]
    ;  phase_print = phases
    while (phase LT 0) do begin
      phase = phase + 2*!dpi
    endwhile
    ;  print,phase_print
    ;  print,phase
    ;  phase_print = phase_print*180/!dpi

    ;  print,'phase =',phase_print, '   star = ', stars[star_i], '  phi = ', phi*180/!dpi


    B = abs(B)
    Beff=abs(atan(tan(B)/cos(phi))) ;this is the effective angle (alpha in Jerousek 2012)

    d_nor = d
    res_int = round(1/reso)


    if Beff GT 1. then begin ;0.872          ; for occultations of phi GT 47 degrees, the thickness of the ring changes the optical depth very little
      d_nor = 10d    ;100             ; (smaller than the noise in the data of AlpVir008E). After this the effect of the thickness in the model is that of
      Res_int = 10d    ;2 '50         ; increasing the resolution (see Gresh 1986 et al., their model has this feature as well),
      ;Resolution of intercept            ; so we take a computational shortcut and use lower resolution (0.001 Km is too much resoution for normal occultations)




    endif
    ;Hamp is the normal optical depth outsiede of the wave for this occultation. We get it with the granola bar model


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

    upper= A * damp * (1/(sqrt(!DPI))) * Exp(j*(phase + !DPI/4)) * (Exp(j*x^2/(2*e*rv*rv))*(0.5*sqrt(!DPI/2) - 0.5*j*sqrt(!DPI/2) + fresnel_complex(x/(sqrt(2*e)*rv)))) +d_nor/2.
    lower= A * damp * (1/(sqrt(!DPI))) * Exp(j*(phase + !DPI/4)) * (Exp(j*x^2/(2*e*rv*rv))*(0.5*sqrt(!DPI/2) - 0.5*j*sqrt(!DPI/2) + fresnel_complex(x/(sqrt(2*e)*rv)))) -d_nor/2.
    ;upper0=(A/(sqrt(!DPI)))*Exp(j*( !DPI/4.)) *damp* (Exp(j*x^2./(2.*e*rv*rv))*(0.5*sqrt(!DPI/2.)-0.5*j*sqrt(!DPI/2.) + fresnel_complex(x/(sqrt(2.*e)*rv))))



    slope     = -(shift(upper, 1) - upper)/reso
    slope[0]  = 0
    ;slope1 = -(shift(slope, 1) - slope)/reso
    slope=abs(slope)

    if keyword_set(isotro) then Hamp  =  meantao(radius,tao,tao=1)*Sin(B) else $
      Hamp  =  meantao(radius,tao,tao=1)*Sin(Beff)

    n       = 200
    step    = 0.1
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

      ;   beta = beta/2*h



      h_haze = 2. * haze_amp * A ;+ y_add


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

    optd[7] = optd[8]
    optd[6] = optd[7] ;Since the simulation is an iterative process the first values are generally off and I interpolate it (altenatively I could also remove it) Note the points at the endges of the model are not used in the fit but I still want to smooth the edges
    optd[5] = optd[6]
    optd[4] = optd[5]
    optd[3] = optd[4]
    optd[2] = optd[3]
    optd[1] = optd[2]
    optd[0] = optd[1]

    ;optd[-3] = optd[-4]
    ; optd[-5] = optd[-6]
    optd[-5] = optd[-2]
    optd[-4] = optd[-5]
    optd[-3] = optd[-4]
    optd[-2] = optd[-3]
    optd[-1] = optd[-2]
    ; joshplot=1

    riso = r
    optdiso = optd

    sm=5
    if star_i eq 5 then sm=70
    ; if star_i eq 6 then sm=70

    if star_i le 2 then begin
      if star_i eq 0 then iso  = plot(r-rv, smooth(optd,sm), xrange=[-200.,50.], yrange=[0.,5.], thick=2, name="This Work",linestyle=0, dimensions = [1366,768], buffer=1,layout=[3,2,star_i+4],/current,AXIS_STYLE=1, xtitle='Distance from resonance [km]', ytitle='Optical Depth', font_size = 14, XTICKFONT_SIZE=16, YTICKFONT_SIZE=16,margin=[.12, .15, .19, .05],/nodata) else $
        iso  = plot(r-rv, smooth(optd,sm), xrange=[-200.,50.], yrange=[0.,5.], thick=2, name="This Work",linestyle=0, dimensions = [1366,768], buffer= 0,layout=[3,2,star_i+4],/current,AXIS_STYLE=1,margin=[0, .15, .25, .05], xtitle='Distance from resonance [km]', font_size = 14, YTICKFONT_SIZE=16, XTICKFONT_SIZE=16,/nodata)


    endif else begin

      if star_i eq 3 then   iso  = plot(r-rv, smooth(optd,sm), xrange=[-200.,50.], yrange=[0.,5.], thick=2, name="This Work",linestyle=0, dimensions = [1366,768], buffer=1,layout=[3,2,star_i-2],/current,AXIS_STYLE=1,xtickformat="(A1)", ytitle='Optical Depth', YTICKFONT_SIZE=16,XTICKFONT_SIZE=16,margin=[.12, 0, .19, .15],/nodata) else $
        iso  = plot(r-rv, smooth(optd,sm), xrange=[-200.,50.], yrange=[0.,5.], thick=2, name="This Work",linestyle=0, dimensions = [1366,768], buffer= 0,layout=[3,2,star_i-2],/current,AXIS_STYLE=1,margin=[0, 0, .25, .2],xtickformat="(A1)", YTICKFONT_SIZE=16,XTICKFONT_SIZE=16)

    endelse
    ;  graphic4.save, stars_format+'.png'


    ;;;;\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    ;;;;\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    ;;;;\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    ;;;;\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


    fit_phase = 1
    create_haze = 1
    isotro=1
    epsilon=0.5

    ;    P0=[3.003670, 0.05656125, 365.51056]

    beta = 0
    v    = .028
    sigma=  P0[2]
    W    = 58.5


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


    d=0.03
    ; d is the vertical thinckness of the ring in Km
    ;I'm reading off the parameter of Mimia's orbit
    ;Restore, Filepath('Mimas_lon_z_danielOX.sav', subdir = ['Mimas_data'])
    ; Restore, 'Mimas_lon_z_danielOX.sav'
    ;  Restore, 'mimascoor.sav'

    ;I'm reading off the geometry of the occultation
    ; ;Restore, 'geometries.sav'
    ;Restore, 'geometry.sav'


    ;Finally, I'm reading off the value for the normal optical depth of the hump according to the granola bar model (Jerousek 2016)
    ;  Restore, 'granola_haze.sav'

    ;  i_geo  = where(stars[star_i] EQ rev)
    ; i_mimas= i_geo + 1

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
    ; print,stars[star_i],UP_DOWN[i_mimas[0]]

    ;print,'hrllo'     ;phase = phase[star_i]
    ;  phase_print = phases
    while (phase LT 0) do begin
      phase = phase + 2*!dpi
    endwhile
    ;  print,phase_print
    ;  print,phase
    ;  phase_print = phase_print*180/!dpi

    ;  print,'phase =',phase_print, '   star = ', stars[star_i], '  phi = ', phi*180/!dpi


    B = abs(B)
    Beff=abs(atan(tan(B)/cos(phi))) ;this is the effective angle (alpha in Jerousek 2012)

    d_nor = d
    res_int = round(1/reso)


    if Beff GT 1. then begin ;0.872          ; for occultations of phi GT 47 degrees, the thickness of the ring changes the optical depth very little
      d_nor = 10d    ;100             ; (smaller than the noise in the data of AlpVir008E). After this the effect of the thickness in the model is that of
      Res_int = 10d    ;2 '50         ; increasing the resolution (see Gresh 1986 et al., their model has this feature as well),
      ;Resolution of intercept            ; so we take a computational shortcut and use lower resolution (0.001 Km is too much resoution for normal occultations)




    endif
    ;Hamp is the normal optical depth outsiede of the wave for this occultation. We get it with the granola bar model


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
      phase = phases[star_i]

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

    e=(2d*!DPI*G*sigma/(rv*1000d*cD)) ;the 1000 is to transform rv into m

    upper= A * damp * (1/(sqrt(!DPI))) * Exp(j*(phase + !DPI/4)) * (Exp(j*x^2/(2*e*rv*rv))*(0.5*sqrt(!DPI/2) - 0.5*j*sqrt(!DPI/2) + fresnel_complex(x/(sqrt(2*e)*rv)))) +d_nor/2.
    lower= A * damp * (1/(sqrt(!DPI))) * Exp(j*(phase + !DPI/4)) * (Exp(j*x^2/(2*e*rv*rv))*(0.5*sqrt(!DPI/2) - 0.5*j*sqrt(!DPI/2) + fresnel_complex(x/(sqrt(2*e)*rv)))) -d_nor/2.
    ;upper0=(A/(sqrt(!DPI)))*Exp(j*( !DPI/4.)) *damp* (Exp(j*x^2./(2.*e*rv*rv))*(0.5*sqrt(!DPI/2.)-0.5*j*sqrt(!DPI/2.) + fresnel_complex(x/(sqrt(2.*e)*rv))))



    slope     = -(shift(upper, 1) - upper)/reso
    slope[0]  = 0
    ;slope1 = -(shift(slope, 1) - slope)/reso
    slope=abs(slope)

    if keyword_set(isotro) then Hamp  =  meantao(radius,tao,tao=1)*Sin(B) else $
      Hamp  =  meantao(radius,tao,tao=1)*Sin(Beff)

    n       = 200
    step    = 0.1
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

      ;   beta = beta/2*h



      h_haze = 2. * haze_amp * A ;+ y_add

      r_haze1 = (-haze_amp * A  + y_int)/s + p
      r_haze2=  ( haze_amp * A  + y_int)/s + p

      r_haze1_i = floor(mean(where(round(x*1000)/1000. EQ round(r_haze1*1000)/1000.)))
      r_haze2_i = floor(mean(where(round(x*1000)/1000. EQ round(r_haze2*1000)/1000.)))
      if ~keyword_set(create_haze) then begin
        yhaze = make_array(n_elements(x),2)
        yhaze[*,1]=haze_amp*A
        yhaze[*,0]=-haze_amp*A
      endif

      ;junk    = min(abs(-s*(x-p) + y_int - A* damp), index)
      ;r_haze1 = x[index]

      ;junk = min(abs(-s*(x-p) + y_int + A* damp), index)
      ;r_haze2 = x[index]

      ; if stars[star_i] eq 'BetCen077I' then begin
      ;     yadd = 10.
      ;yhaze[*,0] = yhaze[*,0] + yadd
      ;yhaze[*,1] = yhaze[*,1] - yadd
      ;      endif


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
        ;  x_inter     = findgen(r_haze2-r_haze1,increment = reso, start=r_haze1)
        ; slope_inter = slope[0] + (slope[1] - slope[0])/(x[1] - x[0]) * (r_haze2-r_haze1)
        ; INT_TABULATED(x_inter[0 : -1],(ratiorad/denominator) * (1/h_haze) * abs((1/(Cos(B)*Cos(phi)))) * beta * slope_inter[0: -1])
        optd_haze   = 0
        optd_minus  = 0

      endif else begin


        if keyword_set(isotro) then   optd_haze  = INT_TABULATED(x[r_haze1_i : r_haze2_i], abs((1/((yhaze[r_haze1_i:r_haze2_i,1]-yhaze[r_haze1_i:r_haze2_i,0])*Cos(B)*Cos(phi))) * beta * abs(slope[r_haze1_i : r_haze2_i]))) else $
          optd_haze  = int_tabulated(x[r_haze1_i : r_haze2_i],  abs((hamp/((yhaze[r_haze1_i:r_haze2_i,1]-yhaze[r_haze1_i:r_haze2_i,0])*Cos(Beff))) * beta* abs(slope[r_haze1_i : r_haze2_i])))
        ;optd_haze  = INT_TABULATED(x[r_haze1_i : r_haze2_i], abs((1/(Cos(B)*Cos(phi))) * beta * abs(slope[r_haze1_i : r_haze2_i])))
        ;  print,optd_haze


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

    optd[7] = optd[8]
    optd[6] = optd[7] ;Since the simulation is an iterative process the first values are generally off and I interpolate it (altenatively I could also remove it) Note the points at the endges of the model are not used in the fit but I still want to smooth the edges
    optd[5] = optd[6]
    optd[4] = optd[5]
    optd[3] = optd[4]
    optd[2] = optd[3]
    optd[1] = optd[2]
    optd[0] = optd[1]

    ;optd[-3] = optd[-4]
    ; optd[-5] = optd[-6]
    optd[-5] = optd[-2]
    optd[-4] = optd[-5]
    optd[-3] = optd[-4]
    optd[-2] = optd[-3]
    optd[-1] = optd[-2]
    ; joshplot=1

    sm=5
    if star_i eq 5 then sm=50

    ;  if star_i eq 6 then sm=20

    if star_i le 2 then begin

      ;fonts=
      aniso   = plot(r-rv, smooth(optd,sm), xrange=[-200.,50.], yrange=[0.,max([max(optd),max(tao), max(optdiso)])+.6], thick=2, name="SCL",linestyle=2, dimensions = [1366,768], buffer=1,overplot=1, color='green',layout=[3,2,star_i+4],/current)
      ; if star_i eq 6 then sm=5

      graphic4   = oplotmimas_fig13(radius,data,star='('+stars_format+')',tao=tao, title=1, uperror=taoplus,downerror=taominus,angles=angles);,scl=1);, error=taosigma) ;, angles = [B, phi], phase = phase_print)

      graphic5   = plot(r-rv, smooth(optd,sm), xrange=[-200.,50.], yrange=[0.,max([max(optd),max(tao), max(optdiso)])+.6], thick=2, name="New Model",overplot=1,linestyle=2, dimensions = [1366,768], color='green',layout=[3,2,star_i+4])

      graphic5   = plot(riso-rv, optdiso, xrange=[-200.,50.], yrange=[0.,max([max(optd),max(tao), max(optdiso)])+.6], thick=2, name="New Model",overplot=1,linestyle=0, dimensions = [1366,768],layout=[3,2,star_i+4])


      if star_i le 1 then begin
        if star_i eq 0 then         junk =  text(0.05+star_i*.29, 0.43, '(d)',font_size=12)
        if star_i eq 1 then         junk =  text(0.05+star_i*.29, 0.43, '(e)',font_size=12)



        junk =  text(0.05+star_i*.29, 0.40, 'B  =' + string(angles[0]*180./!dpi,FORMAT='(F6.2)')+'°',font_size=14)
        junk =  text(0.05+star_i*.29, 0.37, '$\phi$  =' + string(angles[1]*180./!dpi,FORMAT='(F6.2)')+'°',font_size=14)
        junk =  text(0.05+star_i*.29, 0.34, '$B_{eff}$ =' + string(angles[2]*180./!dpi,FORMAT='(F6.2)')+'°',font_size=14)
      endif else begin

        junk =  text(0.1+star_i*.29, 0.43, '(f)',font_size=12)


        junk =  text(0.1+star_i*.29, 0.40, 'B  =' + string(angles[0]*180./!dpi,FORMAT='(F6.2)')+'°',font_size=14)
        junk =  text(0.1+star_i*.29, 0.37, '$\phi$  =' + string(angles[1]*180./!dpi,FORMAT='(F6.2)')+'°',font_size=14)
        junk =  text(0.1+star_i*.29, 0.34, '$B_{eff}$ =' + string(angles[2]*180./!dpi,FORMAT='(F6.2)')+'°',font_size=14)

      endelse


      ; if star_i eq 6 then sm=20

    endif else begin

      aniso   = plot(r-rv, smooth(optd,sm), xrange=[-200.,50.], yrange=[0.,max([max(optd),max(tao), max(optdiso)])+.6], thick=2, name="SCL",linestyle=2, dimensions = [1366,768], buffer=1,overplot=1, color='green',layout=[3,2,star_i-2],/current)

      graphic4   = oplotmimas_fig13(radius,data,star='('+stars_format+')',tao=tao, title=1, uperror=taoplus,downerror=taominus,angles=angles);,scl=1);, error=taosigma) ;, angles = [B, phi], phase = phase_print)
      ; if star_i eq 6 then sm=5
      graphic5   = plot(r-rv, smooth(optd,sm), xrange=[-200.,50.], yrange=[0.,max([max(optd),max(tao), max(optdiso)])+.6], thick=2, name="New Model",overplot=1,linestyle=2, dimensions = [1366,768], color='green',/current)

      graphic5   = plot(riso-rv, smooth(optdiso,sm), xrange=[-200.,50.], yrange=[0.,max([max(optd),max(tao), max(optdiso)])+.6], thick=2, name="New Model",overplot=1,linestyle=0, dimensions = [1366,768],/current)

      if star_i le 4 then begin
        if star_i eq 3 then        junk =  text(0.05+(star_i-3)*.29, 0.83, '(a)',font_size=12)
        if star_i eq 4 then        junk =  text(0.05+(star_i-3)*.29, 0.83, '(b)',font_size=12)


        junk =  text(0.05+(star_i-3)*.29, 0.80, 'B  =' + string(angles[0]*180./!dpi,FORMAT='(F6.2)')+'°',font_size=14)
        junk =  text(0.05+(star_i-3)*.29, 0.77, '$\phi$  =' + string(angles[1]*180./!dpi,FORMAT='(F6.2)')+'°',font_size=14)
        junk =  text(0.05+(star_i-3)*.29, 0.74, '$B_{eff}$ =' + string(angles[2]*180./!dpi,FORMAT='(F6.2)')+'°',font_size=14)
      endif else begin


        junk =  text(0.095+(star_i-3)*.29, 0.83, '(c)',font_size=12)
        junk =  text(0.095+(star_i-3)*.29, 0.80, 'B  =' + string(angles[0]*180./!dpi,FORMAT='(F6.2)')+'°',font_size=14)
        junk =  text(0.095+(star_i-3)*.29, 0.77, '$\phi$  =' + string(angles[1]*180./!dpi,FORMAT='(F6.2)')+'°',font_size=14)
        junk =  text(0.095+(star_i-3)*.29, 0.74, '$B_{eff}$ =' + string(angles[2]*180./!dpi,FORMAT='(F6.2)')+'°',font_size=14)

      endelse

    endelse


    print,phase






    ;graphic4   = oplotmimas_fig13(radius,data,star='('+stars_format+')',tao=tao, title=1, uperror=taoplus,downerror=taominus,modelplot=modelplot);,scl=1);, error=taosigma) ;, angles = [B, phi], phase = phase_print)

    ;graphic4   = plot(r-rv, optd, xrange=[-200.,50.], yrange=[0.,5.], thick=2, name="New Model",overplot=1,linestyle=0, dimensions = [1366,768])


    chi = goodfit3(radius, Irr/I0, r, exp(-optd), points_in_model, stars_format, taosigma = Sqrt(Irr)/I0)

    chis[star_i] = chi
    ;print, stars[star_i],chi
  endfor
  leg  = legend(position=[.65,.48], font_size = 17,target=[iso,aniso,graphic4],linestyle=6)
  graphic5.save, stars_format+'matrix2.png'


  Chi_Total = total(chis)
  print,beta, chi_total/(498.*59. + 497), keyword_set(fit_phase), keyword_set(isotro),epsilon

  save,theory, phases,stars, filename = 'phases2.sav'

  ;toc

end