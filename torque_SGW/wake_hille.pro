
Pro Wake_hillE

  ;Description
  ;Version of wake_rot24.pro with an ellipsoidal shape for the wake. It compues the equilibrium angle for the ellipsoidal shape

tic
  flag  = 1 ;to change initial conditions,1 is such that wakes start flat at a trhough. 0 is suh that the wake starts with no omega but inclined as the amplitude of the slope of the wave
  flag1 = 1 ;to turn on or off the wake potential. 1 is off

  rho = 492
  sigma = rho * 0.5 * sqrt(!dpi)
  lambda_c = 58.
  k = 2*!dpi/lambda_c
  h   = 4d
  l   = 232d
  W   = 18d

  n0    = 30
  ;n0=3
  a0    = .1
  amax  = 18
  amin  = 13e-3
  ;amin=3.2e-3
  q     = -3.1
  
  amin = 13e-3
  amax = 10
  q = -3.15
  n0 = 25
  a0 = .1

  ;rho_p = 0.4
  rho_p = 400.

  as = linspace(amin,amax,100000)

  n      = n0*(as/a0)^q
  print,mean(n)
  a_mean = total(n*as)/total(as)
  
 print,300*2/(H*sqrt(!dpi))

  rho_dist = rho_p * !dpi* (4./3.) * as^3 *  n0*(as/a0)^q
  print,mean(rho_dist)

  rho_mean =   495d*3/5
   rho_mean =   385d
 ; rho_mean = 85.
  ;rho_mean = mean(rho_dist)*3
  ;rho_mean =   360d
  ;rho_mean =   250d

  if flag1 then begin
    ;rho_mean = 1000d
    ;rho_mean = 40d
    rho_mean = 425*3./5.*!dpi/4.
    ;rho_mean = 150.
  endif

  j=dcomplex(0,1)

  ;taux = 5.33333*10^7 Sin[(30*!dpi/180 + 0.0980198 thetax*thetay]
  ;tauy = 533333. Cos[(5 \[Pi])/36 + thetax*thetay]

  rv  = 131902000.
  G   = 6.67e-11
  Ms  = 5.6832e26
  Omega = 0.000128897
  mu = 1.2957087e-4

  s = 1.50524
  ;s = 1.
  ; s = .5
  phi = 40./180. *!dpi

  ;vkdr = -0.0000642821d*3
  vkdr = -s*Omega


  muminus = 2*G*Ms/(rv^3) - 2*(!dpi)*G*(sigma * k )* H*(cos(phi)^2)
  muplus  =   G*Ms/(rv^3) + 2*(!dpi)*G*(sigma * k )* H*(sin(phi)^2)

  if flag1 then begin

    muminus = 2*G*Ms/(rv^3) ;- 2*(!dpi)*G*(sigma * k )* H*(cos(phi)^2)
    muplus  =   G*Ms/(rv^3) ;+ 2*(!dpi)*G*(sigma * k )* H*(sin(phi)^2)

    sigma = 0

  endif

  ;4*Pi*G*k*(140)*H*(500)*l^3*H*W/24 (-Cos(Phi)^2 * Cos(x) * Sin(x) + Sin(Phi)^2 *Cos(x) * Sin(x) + Sin(Phi) * Cos (Phi) (Cos[x]^2 - Sin[x]^2))

  omega = mu

  period = 48493d
  ; dt= .01d
  dt= 50d

  Ts = 2d
  time = period/dt*Ts


  db = .002d
  startb = -13
  b = dindgen(abs(startb)/db ,start = startb, increment = db)
  intw = linspace(-w/2,w/2,80)
  inth = linspace(-.7*W,.7*W,60)
  intl = linspace(-l/2.,l/2,80)


  alphax = make_array(time + 1d, /DOUBLE)
  alphay = make_array(time + 1d, /DOUBLE)
  alphaz = make_array(time + 1d, /DOUBLE)

  omegax = make_array(time + 1d, /DOUBLE)
  omegay = make_array(time + 1d, /DOUBLE)
  omegaz = make_array(time + 1d, /DOUBLE)

  thetax = make_array(time + 1d, /DOUBLE)
  thetay = make_array(time + 1d, /DOUBLE)
  thetaz = make_array(time + 1d, /DOUBLE)



  slope = 0d
  m = rho*H*l*W

  Ixx = 1/20.*m*(H^2 + l^2) ; 2.66773e7
  Iyy = 1/20.*m*(H^2 + W^2) ;  277333.
  Izz = 1/20.*m*(W^2 + l^2) ;2.69333e7 for l=200


  seconds = findgen(time,increment = dt)

  ; dhdot0 =  omega*(slope*cos(0))
  dhdot0 =  -omega*(slope*sin(0))

  omegay[1] = 0.
  thetay[1] =slope*cos(phi) ;omegay[1]*dt                                                 rho_mean * W * Abs(-omegax[i-1] - dhdothalf*cos(thetax[i-1])^2*sin(thetaz[i-1] + phi)) * (-omegax[i-1] - dhdothalf*cos(thetax[i-1])^2*sin(thetaz[i-1] + phi)) * (l/2.)^4 *dt/(2*Ixx)                  ;

  omegax[1] = 0.
  thetax[1] = -slope*sin(phi) ; omegax[1]*dt
  thetay[1] =  -slope*cos(phi)*cos(thetax[1]) ;omegay[1]*dt                                                 rho_mean * W * Abs(-omegax[i-1] - dhdothalf*cos(thetax[i-1])^2*sin(thetaz[i-1] + phi)) * (-omegax[i-1] - dhdothalf*cos(thetax[i-1])^2*sin(thetaz[i-1] + phi)) * (l/2.)^4 *dt/(2*Ixx)                  ;

  omegaz[1] = ((Ixx  -Iyy)/Izz) * omegax[1]*omegay[1]/omega ; (Ixx-Iyy)/Izz * omegay[1]*omegax[1]*dtdeltav = (-omegax - dhdot*cos(thetax)^2*sin(thetaz + phi))*l/2. + (omegay  + dhdot*cos(thetay)^2*cos(thetaz + phi))*W/2.
  ;omegay[i] =  omegay[i-1] -2*cos(omega*i*dt)*slope*(omega^2)*(rho*h*l*w*w*w/24)*(cos(phi)- sin(phi)*thetaz[i-1])/(Iyy)*dt
  ;omegax[i] =  omegax[i-1] +2*cos(omega*i*dt)*slope*(omega^2)*(rho*h*w*l*l*l/24)*(sin(phi) +cos(phi)*thetaz[i-1])/(Ixx)*dt + (Iyy - Izz)/Ixx*omegay[i-1]*omegaz[i-1]
  ;plot, seconds*(2*!dpi/omega)^(-1), thetay,xrange=[0,25],yrange=[-3,3]

  if flag then begin


    dhdot0 =  omega*(slope*cos(0))

    omegay[1] =  -dhdot0*cos(phi)
    thetay[1] =  0;-slope*cos(phi) ;omegay[1]*dt                                                 rho_mean * W * Abs(-omegax[i-1] - dhdothalf*cos(thetax[i-1])^2*sin(thetaz[i-1] + phi)) * (-omegax[i-1] - dhdothalf*cos(thetax[i-1])^2*sin(thetaz[i-1] + phi)) * (l/2.)^4 *dt/(2*Ixx)                  ;

    omegax[1] = -dhdot0*sin(phi)
    thetax[1] =  0;-slope*sin(phi) ; omegax[1]*dt

    omegaz[1] = ((Ixx  -Iyy)/Izz) * omegax[1]*omegay[1]/omega ; (Ix

  endif

  thetaz[1]  = 0. ; omegaz[1]*dt

  ; The first rotation of the major axis of the wake is a z rotation by phi angles.

  Rotation = [[Cos(phi), Sin(phi), 0], [-Sin(phi), Cos(phi), 0], [0, 0, 1]]

  x = Rotation # [1d, 0, 0]
  y = Rotation # [0, 1d, 0]
  z = [0, 0, 1d]

  x0s     = make_array(time + 1d)
  x1s     = make_array(time + 1d)
  x2s     = make_array(time + 1d)
  y0s     = make_array(time + 1d)
  y1s     = make_array(time + 1d)
  y2s     = make_array(time + 1d)
  z0s     = make_array(time + 1d)
  z1s     = make_array(time + 1d)
  z2s     = make_array(time + 1d)
  phis    = make_array(time + 1d)
  thetaxs = make_array(time + 1d)
  thetays = make_array(time + 1d)
  vzrecord  = make_array(time + 1d)

  phis[0:2] = phi

  for i=2d, time do begin


    if y[2] gt .97 then y[2] = .97d
    if y[2] lt -.97 then y[2] = -.97d
    if x[2] gt .97 then x[2] = .97d
    if x[2] lt -.97 then x[2] = -.97d
    if z[2] gt .97 then z[2] = .97d
    if z[2] lt -.97 then z[2] = -.97d



    if x[1] gt .97 then x[1] = .97d
    if x[1] lt -.97 then x[1] = -.97d
    if y[1] gt .97 then y[1] = .97d
    if y[1] lt -.97 then y[1] = -.97d
    if z[1] gt .97 then z[1] = .97d
    if z[1] lt -.97 then z[1] = -.97d


    if x[0] gt .97 then x[0] = .97d
    if x[0] lt -.97 then x[0] = -.97d
    if z[0] gt .97 then z[0] = .97d
    if z[0] lt -.97 then z[0] = -.97d
    if y[0] gt .97 then y[0] = .97d
    if y[0] lt -.97 then y[0] = -.97d

    x0s[i-2] = x[0]
    x1s[i-2] = x[1]
    x2s[i-2] = x[2]
    y0s[i-2] = y[0]
    y1s[i-2] = y[1]
    y2s[i-2] = y[2]
    z0s[i-2] = z[0]
    z1s[i-2] = z[1]
    z2s[i-2] = z[2]


    dhdothalf =  -omega*(slope*sin(omega*dt*((i-1) + 1/2.)))
    dhdot     =  -omega*(slope*sin(omega*(i-1)*dt))

    dhdotdothalf = -(omega^2)*slope*cos(omega*dt*(i-1./2.))
    dhdotdot     = -(omega^2)*slope*cos(omega*dt*(i-1))
    if flag then begin


      dhdothalf    =  omega  *  slope*cos(omega*dt*(i-1./2.))
      dhdot        =  omega  *  slope*cos(omega*dt*(i-1))

      dhdotdothalf = -(omega^2)*slope*sin(omega*dt*(i-1./2.))
      dhdotdot     = -(omega^2)*slope*sin(omega*dt*(i-1))

    endif



    torqueyhalf  = -2*dhdotdot*(rho*h*l*w^3/24.) * z[2]*x[0] /(Iyy) * dt/2. - !dpi/120d * muminus*(rho*h*l*w^3) * x[0]*z[0] / (Iyy) * dt/2. + !dpi/120d * muplus*(rho*h*l*w^3) * x[1]*z[1] / (Iyy) * dt/2. + 4*(!dpi)*G*(sigma)*(k)* H * (rho*h*l*w*w*w/24.) * Sin(Phi) * Cos(Phi) * (z[1]*x[0] + z[0]*x[1]) / (Iyy) * dt/2.
    torquexhalf  =  2*dhdotdot*(rho*h*w*l^3/24.) * z[2]*y[0] /(Ixx) * dt/2. + !dpi/120d * muminus*(rho*h*w*l^3) * y[0]*z[0] / (Ixx) * dt/2. - !dpi/120d * muplus*(rho*h*w*l^3) * y[1]*z[1] / (Ixx) * dt/2. - 4*(!dpi)*G*(sigma)*(k)* H * (rho*h*w*l*l*l/24.) * Sin(Phi) * Cos(Phi) * (z[1]*y[0] + z[0]*y[1]) / (Ixx) * dt/2.
    torquezhalf  = -2*dhdotdot*(rho*h*w*l^3/24.) * x[2]*y[0] /(Izz) * dt/2. - !dpi/120d * muminus*(rho*h*w*l^3) * y[0]*x[0] / (Izz) * dt/2. + !dpi/120d * muplus*(rho*h*w*l^3) * x[1]*y[1] / (Izz) * dt/2. + 4*(!dpi)*G*(sigma)*(k)* H * (rho*h*w*l*l*l/24.) * Sin(Phi) * Cos(Phi) * (x[1]*y[0] + x[0]*y[1]) / (Izz) * dt/2.


    inertiayhalf = ((Izz - Ixx)/Iyy)*omegaz[i-1]*omegax[i-1]*dt/2.
    inertiaxhalf = ((Iyy - Izz)/Ixx)*omegay[i-1]*omegaz[i-1]*dt/2.
    inertiazhalf = ((Ixx - Iyy)/Izz)*omegay[i-1]*omegax[i-1]*dt/2.

    stub = hillcol4dz(b, l,w,h,s,rho,rho_mean, x[0], x[1], x[2], y[0], y[1], y[2], z[0], z[1], z[2], omegax[i-1],omegay[i-1],omegaz[i-1], ind=i, phis=phis, dhdot = dhdot, thetays = thetays, thetaxs=thetaxs, intl = IntL,intw=intw,inth=inth,vzrecord=vzrecord)

    dragyhalf    =    stub[0]*dt/(2*Iyy) ;- torquecol(l,w,h,s,rho_mean, y[0], y[1], y[2], omegarot = omegaz[i-1], z2 = z[2], z1 = z[1],        ind=i, phis=phis, oldphi = phis[i-1])*z[2]*dt/(2*Izz)
    dragxhalf    =    stub[1]*dt/(2*Ixx)          ;note that delta v is positive if particles are "going up" (positive z coordinate) wrt the particles in the wake positive y coordinate
    dragzhalf    =    stub[2]*dt/(2*Izz)   ;* sin(thetaz[i-1] + phi);*cos(thetay[i-1])*cos(thetax[i-1])*(sin(thetay[i-1]) + cos(thetax[i-1]))

    print, torquezhalf/(dt)*Izz
    print,'col', dragzhalf/(dt)*Izz
    print,phis[i]*180/!dpi


    omegayhalf =  omegay[i-1] + torqueyhalf + inertiayhalf + dragyhalf
    omegaxhalf =  omegax[i-1] + torquexhalf + inertiaxhalf + dragxhalf
    omegazhalf =  omegaz[i-1] + torquezhalf + inertiazhalf + dragzhalf


    thetayhalf  =  omegay[i-1]*dt/2.
    thetaxhalf  =  omegax[i-1]*dt/2.
    thetazhalf  =  omegaz[i-1]*dt/2.

    Rotation  = [[1, thetazhalf, -thetayhalf], [-thetazhalf, 1, thetaxhalf], [thetayhalf, -thetaxhalf, 1]]

    xhalf = Matrix_multiply(Rotation, x)
    yhalf = Matrix_multiply(Rotation, y)
    zhalf = Matrix_multiply(Rotation, z)




    ;  print,zhalf[2]
    ;  print, xhalf[1]


    ;print,torquezhalf
    ;       print,omegay[i-1]

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ; dhdot =  omega*(slope*cos(omega*i*dt))




    torquey = - 2*dhdotdothalf*(rho*h*l*w*w*w/24.) * zhalf[2] * xhalf[0]/(Iyy) * dt - !dpi/120d * muminus*(rho*h*l*w^3) * xhalf[0]*zhalf[0] / (Iyy) * dt + !dpi/120d * muplus*(rho*h*l*w^3) * xhalf[1]*zhalf[1] / (Iyy) * dt + 4*(!dpi)*G*(sigma)*(k)* H * (rho*h*l*w*w*w/24.) * Sin(Phi) * Cos(Phi) * (zhalf[1]*xhalf[0] + zhalf[0]*xhalf[1]) / (Iyy) * dt
    torquex =   2*dhdotdothalf*(rho*h*w*l*l*l/24.) * zhalf[2] * yhalf[0]/(Ixx) * dt + !dpi/120d * muminus*(rho*h*w*l^3) * yhalf[0]*zhalf[0] / (Ixx) * dt - !dpi/120d * muplus*(rho*h*w*l^3) * yhalf[1]*zhalf[1] / (Ixx) * dt - 4*(!dpi)*G*(sigma)*(k)* H * (rho*h*w*l*l*l/24.) * Sin(Phi) * Cos(Phi) * (zhalf[1]*yhalf[0] + zhalf[0]*yhalf[1]) / (Ixx) * dt
    torquez = - 2*dhdotdothalf*(rho*h*w*l*l*l/24.) * xhalf[2] * yhalf[0]/(Izz) * dt - !dpi/120d * muminus*(rho*h*w*l^3) * yhalf[0]*xhalf[0] / (Izz) * dt + !dpi/120d * muplus*(rho*h*w*l^3) * xhalf[1]*yhalf[1] / (Izz) * dt + 4*(!dpi)*G*(sigma)*(k)* H * (rho*h*w*l*l*l/24.) * Sin(Phi) * Cos(Phi) * (xhalf[1]*yhalf[0] + xhalf[0]*yhalf[1]) / (Izz) * dt

    inertiay = ((Izz - Ixx)/Iyy) * omegazhalf*omegaxhalf*dt
    inertiax = ((Iyy - Izz)/Ixx) * omegayhalf*omegazhalf*dt
    inertiaz = ((Ixx  -Iyy)/Izz) * omegayhalf*omegaxhalf*dt

    stub = hillcol4dz(b, l,w,h,s,rho,rho_mean, xhalf[0], xhalf[1], xhalf[2], yhalf[0], yhalf[1], yhalf[2], zhalf[0], zhalf[1], zhalf[2], omegaxhalf, omegayhalf, omegazhalf, ind=i, phis=phis, dhdot = dhdot, thetays = thetays, thetaxs=thetaxs, intl = IntL,intw=intw,inth=inth,vzrecord=vzrecord)

    dragy    =  stub[0]*dt/(Iyy)
    dragx    =  stub[1]*dt/(Ixx)                  ;note that delta v is positive if particles are "going up" (positive z coordinate) wrt the particles in the wake positive y coordinate
    dragz    =  stub[2]*dt/(Izz)

    print, torquez/(dt)*Izz
    print,'col', dragz/(dt)*Izz
    print,phis[i]*180/!dpi

    omegay[i] =  omegay[i-1] + torquey + inertiay + dragy
    omegax[i] =  omegax[i-1] + torquex + inertiax + dragx
    omegaz[i] =  omegaz[i-1] + torquez + inertiaz + dragz

    ; print,'torquez  =',torquez * Izz/dt


    thetay[i] =  thetay[i-1] + omegayhalf*dt
    thetax[i] =  thetax[i-1] + omegaxhalf*dt
    thetaz[i]  = thetaz[i-1] + omegazhalf*dt


    alphay[i] = (torquey + inertiay + dragy)/dt
    alphax[i] = (torquex + inertiax + dragx)/dt
    alphaz[i] = (torquez + inertiaz + dragz)/dt


    thetayhalf  =  omegayhalf*dt
    thetaxhalf  =  omegaxhalf*dt
    thetazhalf  =  omegazhalf*dt

    Rotation  = [[1, thetazhalf, -thetayhalf], [-thetazhalf, 1, thetaxhalf], [thetayhalf, -thetaxhalf, 1]]

    x = Matrix_multiply(Rotation, x)
    y = Matrix_multiply(Rotation, y)
    z = Matrix_multiply(Rotation, z)

    if alphax[i] gt 3e-7 then begin
      print,'hello'
    endif

    ;print,'dragz', dragz
    ;print, 'tidal', - 3 * muminus*(rho*h*w*l^3/24.) * yhalf[0]*xhalf[0] / (Izz) * dt
    ;print,'secods', seconds[i-1]*(2*!dpi/omega)^(-1)

    ;if abs(x[1]) gt 1 then x[1]= .98
    ;if abs(z[1]) gt 1 then z[1]=zcopy[1]
    ;if abs(z[2]) gt 1 then z[2]=zcopy[2]




  endfor

  ;  j = s*omega*b^2

  ;define energy

  ; E = 1/2.*(s*omega*b)^2

  ;define eccentricity

  ;ec = Sqrt(1 + 2*E*j^2/(GM)^2)

  ; define argument of periapsis

  ;po = acos(b*(GM/j^2)*(ec^2 - 1)/ec)
  ;cos1 = x1s*cos(po) + xs*sin(po)
  ;vtheta = GM*(1 + ec*Cos1)/j



  dhdot  = -omega*(slope*sin(omega*seconds))
  if flag then dhdot  = omega*(slope*cos(omega*seconds))

  deltav = (-omegax + dhdot*y0s*z2s + vkdr *y0s*z1s) * (l/2.)  + (omegay + dhdot*x0s*z2s + vkdr * x0s*z1s) * (W/2.)

  deltavx =( omegaz + dhdot*y0s*x2s + vkdr *y0s*x1s) * (l/2.)

  ;omegay[i] =  omegay[i-1] -2*cos(omega*i*dt)*slope*(omega^2)*(rho*h*l*w*w*w/24)*(cos(phi)- sin(phi)*thetaz[i-1])/(Iyy)*dt
  ;omegax[i] =  omegax[i-1] +2*cos(omega*i*dt)*slope*(omega^2)*(rho*h*w*l*l*l/24)*(sin(phi) +cos(phi)*thetaz[i-1])/(Ixx)*dt + (Iyy - Izz)/Ixx*omegay[i-1]*omegaz[i-1]
  ;plot, seconds*(2*!dpi/omega)^(-1), thetay,xrange=[0,25],yrange=[-3,3]

  deltavk = - omegax * z2s * l/2 + omegay * W/2* z2s + omegaz * l/2. * x2s +  dhdot * y0s * l/2 + dhdot * x0s * W/2 ;+ vzrecord

  deltavj = - omegax * z1s * l/2 + omegay * W/2* z1s + omegaz * l/2. * x1s +  vkdr *y0s* l/2 + vkdr  * x0s * W/2

  ; deltavi = (omegaz*cos(thetay)*sin(thetaz + phi)  - vkdr *cos(thetax)*sin(phi + thetaz])) * (l/2.) + (omegay*cos(thetay)*sin(thetax)  - vkdr* cos(thetay)*cos(thetaz + phi))*W/2 + (-omegax*cos(thetay)*sin(thetax)   - vkdr*sin(thetaz + phi)*cos(thetay))*l/2



  ; plot,  seconds*(2*!dpi/omega)^(-1), deltav;,xrange=[0,2]

  ; plot,  seconds*(2*!dpi/omega)^(-1), sqrt(deltavk^2 + deltavj^2);,xrange=[0,2]
  ; oplot,  seconds*(2*!dpi/omega)^(-1),  sqrt(deltav^2 + deltavx^2), color=1000
  ; oplot, seconds*(2*!dpi/omega)^(-1) ,1e-2*sin(omega*seconds),color=50000

  wi,1,wsize = [1880,980]
  !p.multi=[1,1,1]

  ;plot,[1,2]
  plot,[1,2]
  deltavk= deltavk[0:-2]

  ; i = where(abs(shift(phis,1) - phi) gt 3)

!p.multi[0,0]
  cgplot,seconds*(2*!dpi/omega)^(-1),   phis*180/!dpi ,color=cgcolor('blue'),charsize=2.8, title='Time evolution of self-gravity pitch-angle', xtitle='Time (Orbital Periods)', ytitle='Degrees',yrange=[0,50],xrange=[0,i/time];,yrange=[-2,2]
  cgoplot,seconds*(2*!dpi/omega)^(-1),   y2s, color = cgcolor('red')
  cgoplot,seconds*(2*!dpi/omega)^(-1),   z2s,color= cgcolor('green')
  ;  cgplot,seconds*(2*!dpi/omega)^(-1),   alphax ,color=cgcolor('blue'),charsize=2.8, xtitle='Time (Orbital Periods)', ytitle='Angular acceleraion [1/s^2]',xrange=[0,2];,yrange=[-2,2]
  ;  cgoplot,seconds*(2*!dpi/omega)^(-1),  alphay, color = cgcolor('red')
  ;  cgoplot,seconds*(2*!dpi/omega)^(-1),  alphaz,color= cgcolor('green')
  ;    cgplot, seconds*(2*!dpi/omega)^(-1),    omegax, color=cgcolor('blue'),charsize=2.8, xtitle='Time (Orbital Periods)', ytitle='Angular velocity [1/s]',xrange=[0,30],yrange=[-6e-4,15e-4]
  ;    cgoplot, seconds*(2*!dpi/omega)^(-1),   omegay ,color= cgcolor('red')
  ;    cgoplot, seconds*(2*!dpi/omega)^(-1),   omegaz,color= cgcolor('green')
  ;  cgAxis, YAxis=1, YRange=[0, 0.01],charsize=2.8, /Save, ytitle='Relative speed [m/s]'
  cgplot,seconds*(2*!dpi/omega)^(-1),  deltav*100, nsum=1, xtitle='Time (Orbital Periods)',  ytitle='Relative speed [cm/s]',charsize=2.8,xrange=[0,Ts],yrange=[-5,5];
  cgLegend, Color=['Blue','Red', 'Green'], Location=[0.12, 0.92],  Titles=['!4h!Xx', '!4h!Xy','!4h!Xz'], Length=0.075, /Box, VSpace=2.75,/Background, bg_color='white',charsize=2.2
  cgLegend, Color=['Blue','Red', 'Green'], Location=[0.2, 0.48],  Titles=['!4a!Xx', '!4a!Xy','!4a!Xz'], Length=0.075, /Box, VSpace=2.75,/Background, bg_color='white',charsize=2.2
  cgLegend, Color=['Black'], Location=[0.85, 0.41],  Titles=['!4D!Xv'], Length=0.075, /Box, VSpace=2.75,/Background, bg_color='white',charsize=2.2

  print,mean(phis*180/!dpi)
  ;  deltavsk  = max(abs(deltavk[period/dt *2:-1]),i)
  deltavsk1 = max(abs(deltavk),i)
  timesk =  sqrt(mean((deltavk)^2 ))

  ;   print,deltavsk,deltavsk1,timesk

  func   = (dhdot[0:-2]*y0s[0:-2] - omegax[0:-2]*z2s[0:-2])*l/2 + omegaz * l/2. * x2s
  ;ts_hants(func, amplitudes=amp,phases=phases,num_period=period/dt,seno=1)
  write_png,'ellipseequilibrio.png',TVRD(/TRUE)
  SAVE, /VARIABLES, FILENAME = 'ellipseequilibrio.sav'
  ;E, /VARIABLES, FILENAME = 'seethefolderdt5.sav'
  ;save, deltav, z2s,y2s,x2s,thetaxs,thetays,phis, filename ='flag1longvariables.sav'

  ;cgLegend, Color=['Blue','Red', 'Green','Black'], Location=[0.1, 0.1],  Titles=['!4h!Xx', '!4h!Xy','!4h!Xz', '|!4D!Xv|'], Length=0.075, /Box, VSpace=2.75,/Background, bg_color='white',charsize=2.2

  ;plot, seconds*(2*!dpi/omega)^(-1),    omegay, xrange=[0,3];,yrange=[-2,2]
  ;oplot, seconds*(2*!dpi/omega)^(-1),   omegax ,color=1000.
  ;oplot, seconds*(2*!dpi/omega)^(-1),   omegaz,color=50000.

  ;plot, seconds*(2*!dpi/omega)^(-1),   x0s ,yrange=[-2,2]
  ;oplot, seconds*(2*!dpi/omega)^(-1),  x1s ,color=1000.
  ;oplot, seconds*(2*!dpi/omega)^(-1),  x2s ,color=50000.

  ; print, sin(omega*seconds[where(deltavk eq max(deltavk))])


  ;plot,
  ;  oplot, seconds*(2*!dpi/omega)^(-1) ,1e-3*cos(omega*seconds),color=1000
  z
  zeroes = [12083*2, 12083*4, 12083*6, 12083*8, 12083*10, 12083*12, 12083*14]

end
