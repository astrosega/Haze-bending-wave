
;Program that inegrates Hill equations of motion. The goal is to inegrate the trajectories to study the rotational evolution of self-gravity wakes.
;
;Log:
;
;3/30/2023 -> created
;
;
;
;FUNCTION differential, X, Y

;Ax = -(-2 Vy[t] + X[t] (-3 + 3/(X[t]^2 + Y[t]^2)^(3./2.)))
;Ay = -( 2 Vx[t] +      (3*Y)/(X^2 + Y^2)^(3./2.))

;deltaX/h = Vx
;deltay/h = Vy

;  RETURN, [-0.5 * Y[0], 4.0 - 0.3 * Y[1] - 0.1 * Y[0]]
;END

function Hillcol4Deq2, bs, lc, w, h, s, rho, rho_mean, x0, x1, x2, y0, y1, y2, z0, z1, z2, omegax,omegay,omegaz, ind=ind, phis=phis, oldphi = oldphi, dhdot = dhdot, thetaxs = thetaxs, thetays = thetays, intw=intW, IntH=inth, intl=intl, vzrecord=vzrecord
  omegax1=omegax
  omegay1=omegay
  omegaz1=omegaz

  std = 10d
  rv  = 131902000.
  G   = 6.67e-11
  M   = rho*h*w*lc          ;for l=240, delta a 58, 54 for l=200
  omega = 1.28897e-4
  GM = G*M
  vkdr = -s*Omega

  phi  =  atan(-y0,y1)
  phis[ind] = phi


  tetax = asin(y2)
  thetaxs[ind] = tetax ; acos(x0)
  tetax = abs(tetax)

  if x2 gt 1 then x2 = 1
  if x2 lt -1 then x2 = -1
  phi1   = atan(-x0,x1)
  tetay = !dpi/2 - acos(x2)
  thetays[ind] = tetay; acos(y1)

  if ~finite(y0) then print, 'it has happened'
  ;  print,'y0', y0
  ;  print,'y1', y1
  ;print,'y2', y2
  ;  print,'z0', z0
  ;  print,'z1', z1
  ;  print,'z2', z2
  ;  print,'x0', x0
  ;  print,'x1', x1
  ;  print,'x2', x2
  ;  print, y0^2 + y1^2 + y2^2
  ;  print, z0^2 + z1^2 + z2^2
  ;  print, x0^2 + x1^2 + x2^2
  ;endif
  fudge1x = tsum(intw,Exp(-(z2*(H/2)/(std/2)  + intw*x2/(std/2))^2))
  fudge1y = tsum(intl, Exp(-(z2*(H/2)/(std/2) + intl*y2/(std/2))^2))
  fudge1z = tsum(inth,Exp(-(x2*(W/2)/(std/2) + inth*z2/(std/2))^2))


  if phi  lt 0 then phi  = !dpi + phi
  ; if phi1 lt 0 then phi1 = !dpi + phi


  ;  tic


  Ms = 5.6832e26

  rh = (rho*20.*H*W/(3*Ms))^(1/3.)*131902000.
  ;rh = 8
  ;m = 492*200*20*5/(Ms +  492*200*20*5)

  if ~finite(1/cos(tetax)^4) then tetax =0.01d

  Lcp = Lc;*abs(cos(tetax))
  Wcp = (W*abs(cos(tetay)))
  L = Lc/(rh)
  m=3d

  !p.multi=[0,0,0]

  dt = .05d
  ;dt = 5d
  time = 2*!dpi/dt


  vx = make_array(time+1,n_elements(bs),9, /double)
  vy = make_array(time+1,n_elements(bs),9,/double)
  vz = make_array(time+1,n_elements(bs),9,/double)
  x =  make_array(time+1,n_elements(bs),9,/double)
  y =  make_array(time+1,n_elements(bs),9,/double)
  z =  make_array(time+1,n_elements(bs),9,/double)
  q0 = make_array(time+1,n_elements(bs),9,/double)
  q1 = make_array(time+1,n_elements(bs),9,/double)
  v0 = make_array(time+1, n_elements(bs),9,/double)
  v1 = make_array(time+1, n_elements(bs),9,/double)

  t = findgen(time+2)
  torque = make_array(n_elements(bs))

  gforcex = make_array(n_elements(bs),9)
  gforcey = make_array(n_elements(bs),9)
  gforcez = make_array(n_elements(bs),9)
  q = make_array(n_elements(bs),9)
  vp = make_array(n_elements(bs),9)

  i= 0d
  hit = 0
  x[0,*,*]  = [bs,bs,bs,bs,bs,bs,bs,bs,bs]
  y[0,*,*]  = -30d
  z[0,*,0]  = .21d
  z[0,*,1]  = .24d
  z[0,*,2]  = .27d
  z[0,*,3]  =  .3d
  z[0,*,4]  =   .335d
  z[0,*,5]  =   .36d
  z[0,*,6]  =   .39d
  z[0,*,7]  =   .42d
  z[0,*,8]  =   .45d
  
;  z[0,*,0]  = .23d
;  z[0,*,1]  = .26d
;  z[0,*,2]  = .29d
;  z[0,*,3]  =  .32d
;  z[0,*,4]  =   .35d
;  z[0,*,5]  =   .38d
;  z[0,*,6]  =   .4d
;  z[0,*,7]  =   .43d
;  z[0,*,8]  =   .47d

  vx[0,*,*] = 0d
  vy[0,*,*] = -s*x[0,*,*];(3/2)*x[0];+ 3/r
  vz[0,*,*]  = 0d



  for i=1d, time do begin

    if i gt 53 and i lt 63 then  begin
      vx[i,where(bs lt -5.5),*]  = 0
      vy[i,where(bs lt -5.5),*]  = 0
      vz[i,where(bs lt -5.5),*]  = 0
      ; x[i, where(bs lt -5),*] = 30
      ; y[i, where(bs lt -5),*] = 30
      ; z[i, where(bs lt -5),*]  = 30
    endif


    gforcex[*,*] =3*x[i-1,*,*] + (m*(-(L/2) + x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*] + Sqrt((x[i-1,*,*]*x0 + x1*y[i-1,*,*] + x2*z[i-1,*,*])^2 + (L/2 - x[i-1,*,*]*y0 - y[i-1,*,*]*y1 - y2*z[i-1,*,*])^2 + (x[i-1,*,*]*z0 + y[i-1,*,*]*z1 + z[i-1,*,*]*z2)^2))*((y0 + (2*x0*(x[i-1,*,*]*x0 + x1*y[i-1,*,*] + x2*z[i-1,*,*]) + 2*y0*(L/2 + x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*]) + 2*z0*(x[i-1,*,*]*z0 + y[i-1,*,*]*z1 + z[i-1,*,*]*z2))$
      /(2*Sqrt((x[i-1,*,*]*x0 + x1*y[i-1,*,*] + x2*z[i-1,*,*])^2 + (L/2 + x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*])^2 + (x[i-1,*,*]*z0 + y[i-1,*,*]*z1 + z[i-1,*,*]*z2)^2)))/(-(L/2) + x[i-1,*,*]*y0 + y[i-1,*,*]*y1 +y2*z[i-1,*,*] + Sqrt((x[i-1,*,*]*x0 + x1*y[i-1,*,*] + x2*z[i-1,*,*])^2 + (L/2 - x[i-1,*,*]*y0 - y[i-1,*,*]*y1 - y2*z[i-1,*,*])^2 + (x[i-1,*,*]*z0 + y[i-1,*,*]*z1 + z[i-1,*,*]*z2)^2))$
      - ((y0 + (2*x0*(x[i-1,*,*]*x0 + x1*y[i-1,*,*] + x2*z[i-1,*,*]) - 2*y0*(L/2 - x[i-1,*,*]*y0 - y[i-1,*,*]*y1 - y2*z[i-1,*,*]) + 2*z0*(x[i-1,*,*]*z0 + y[i-1,*,*]*z1 + z[i-1,*,*]*z2))/(2*Sqrt((x[i-1,*,*]*x0 + x1*y[i-1,*,*] + x2*z[i-1,*,*])^2 + (L/2 - x[i-1,*,*]*y0 - y[i-1,*,*]*y1 - y2*z[i-1,*,*])^2 + (x[i-1,*,*]*z0 + y[i-1,*,*]*z1 + z[i-1,*,*]*z2)^2)))$
      *(L/2 + x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*] + Sqrt((x[i-1,*,*]*x0 + x1*y[i-1,*,*] + x2*z[i-1,*,*])^2 + (L/2 + x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*])^2 + (x[i-1,*,*]*z0 + y[i-1,*,*]*z1 + z[i-1,*,*]*z2)^2)))$
      /(-(L/2) + x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*] + Sqrt((x[i-1,*,*]*x0 + x1*y[i-1,*,*] + x2*z[i-1,*,*])^2 + (L/2 - x[i-1,*,*]*y0 - y[i-1,*,*]*y1 - y2*z[i-1,*,*])^2 + (x[i-1,*,*]*z0 + y[i-1,*,*]*z1 + z[i-1,*,*]*z2)^2))^2))/(L/2 + x[i-1,*,*]*y0 + y[i-1,*,*]*y1 +y2*z[i-1,*,*] + Sqrt((x[i-1,*,*]*x0 + x1*y[i-1,*,*] + x2*z[i-1,*,*])^2 + (L/2 + x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*])^2 + (x[i-1,*,*]*z0 + y[i-1,*,*]*z1 + z[i-1,*,*]*z2)^2))

    gforcey[*,*] = (m*(-(L/2) + x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*] + Sqrt((x[i-1,*,*]*x0 + x1*y[i-1,*,*] + x2*z[i-1,*,*])^2 + (L/2 - x[i-1,*,*]*y0 - y[i-1,*,*]*y1 - y2*z[i-1,*,*])^2 + (x[i-1,*,*]*z0 + y[i-1,*,*]*z1 + z[i-1,*,*]*z2)^2)) $
      *((y1 + (2*x1*(x[i-1,*,*]*x0 + x1*y[i-1,*,*] + x2*z[i-1,*,*]) + 2*y1*(L/2 + x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*]) + 2*z1*(x[i-1,*,*]*z0 + y[i-1,*,*]*z1 + z[i-1,*,*]*z2))/(2*Sqrt((x[i-1,*,*]*x0 + x1*y[i-1,*,*] + x2*z[i-1,*,*])^2 + (L/2 + x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*])^2 + (x[i-1,*,*]*z0 + y[i-1,*,*]*z1 + z[i-1,*,*]*z2)^2)))$
      /(-(L/2) + x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*] + Sqrt((x[i-1,*,*]*x0 + x1*y[i-1,*,*] + x2*z[i-1,*,*])^2 + (L/2 - x[i-1,*,*]*y0 - y[i-1,*,*]*y1 - y2*z[i-1,*,*])^2 + (x[i-1,*,*]*z0 + y[i-1,*,*]*z1 + z[i-1,*,*]*z2)^2)) $
      - ((y1 + (2*x1*(x[i-1,*,*]*x0 + x1*y[i-1,*,*] + x2*z[i-1,*,*]) - 2*y1*(L/2 - x[i-1,*,*]*y0 - y[i-1,*,*]*y1 - y2*z[i-1,*,*]) + 2*z1*(x[i-1,*,*]*z0 + y[i-1,*,*]*z1 + z[i-1,*,*]*z2))/(2*Sqrt((x[i-1,*,*]*x0 + x1*y[i-1,*,*] + x2*z[i-1,*,*])^2 + (L/2 - x[i-1,*,*]*y0 - y[i-1,*,*]*y1 - y2*z[i-1,*,*])^2 + (x[i-1,*,*]*z0 + y[i-1,*,*]*z1 + z[i-1,*,*]*z2)^2)))$
      *(L/2 + x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*] + Sqrt((x[i-1,*,*]*x0 + x1*y[i-1,*,*] + x2*z[i-1,*,*])^2 + (L/2 + x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*])^2 + (x[i-1,*,*]*z0 + y[i-1,*,*]*z1 + z[i-1,*,*]*z2)^2)))/(-(L/2) + x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*] + Sqrt((x[i-1,*,*]*x0 + x1*y[i-1,*,*] + x2*z[i-1,*,*])^2 + (L/2 - x[i-1,*,*]*y0 - y[i-1,*,*]*y1 - y2*z[i-1,*,*])^2 + (x[i-1,*,*]*z0 + y[i-1,*,*]*z1 + z[i-1,*,*]*z2)^2))^2))$
      /(L/2 + x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*] + Sqrt((x[i-1,*,*]*x0 + x1*y[i-1,*,*] + x2*z[i-1,*,*])^2 + (L/2 + x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*])^2 + (x[i-1,*,*]*z0 + y[i-1,*,*]*z1 + z[i-1,*,*]*z2)^2))


    gforcez[*,*] = (m*(-(L/2) + x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*] + Sqrt((x[i-1,*,*]*x0 + x1*y[i-1,*,*] + x2*z[i-1,*,*])^2 + 1./4.*(L - 2*(x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*]))^2 + (x[i-1,*,*]*z0 + y[i-1,*,*]*z1 + z[i-1,*,*]*z2)^2))$
      *(-(((L/2 + x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*] + Sqrt((x[i-1,*,*]*x0 + x1*y[i-1,*,*] + x2*z[i-1,*,*])^2 + (L/2 + x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*])^2 + (x[i-1,*,*]*z0 + y[i-1,*,*]*z1 + z[i-1,*,*]*z2)^2))*(y2 + (-L*y2 + 2*x2*(x[i-1,*,*]*x0 + x1*y[i-1,*,*] + x2*z[i-1,*,*]) + 2*y2*(x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*]) + 2*z2*(x[i-1,*,*]*z0 + y[i-1,*,*]*z1 + z[i-1,*,*]*z2))$
      /(2*Sqrt((x[i-1,*,*]*x0 + x1*y[i-1,*,*] + x2*z[i-1,*,*])^2 + 1/4.*(L - 2*(x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*]))^2 + (x[i-1,*,*]*z0 + y[i-1,*,*]*z1 + z[i-1,*,*]*z2)^2))))$
      /(-(L/2) + x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*] + Sqrt((x[i-1,*,*]*x0 + x1*y[i-1,*,*] + x2*z[i-1,*,*])^2 + 1/4.*(L - 2*(x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*]))^2 + (x[i-1,*,*]*z0 + y[i-1,*,*]*z1 + z[i-1,*,*]*z2)^2))^2) + (y2 + (x2*(x[i-1,*,*]*x0 + x1*y[i-1,*,*] + x2*z[i-1,*,*]) + y2*(L/2 + x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*]) + z2*(x[i-1,*,*]*z0 + y[i-1,*,*]*z1 + z[i-1,*,*]*z2))$
      /Sqrt((x[i-1,*,*]*x0 + x1*y[i-1,*,*] + x2*z[i-1,*,*])^2 + (L/2 + x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*])^2 + (x[i-1,*,*]*z0 + y[i-1,*,*]*z1 + z[i-1,*,*]*z2)^2))$
      /(-(L/2) + x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*] + Sqrt((x[i-1,*,*]*x0 + x1*y[i-1,*,*] + x2*z[i-1,*,*])^2 + 1./4.*(L - 2*(x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*]))^2 + (x[i-1,*,*]*z0 + y[i-1,*,*]*z1 + z[i-1,*,*]*z2)^2))))$
      /(L/2 + x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*] + Sqrt((x[i-1,*,*]*x0 + x1*y[i-1,*,*] + x2*z[i-1,*,*])^2 + (L/2 + x[i-1,*,*]*y0 + y[i-1,*,*]*y1 + y2*z[i-1,*,*])^2 + (x[i-1,*,*]*z0 + y[i-1,*,*]*z1 + z[i-1,*,*]*z2)^2))


    vxhalf = vx[i-1, *,*] + gforcex*dt/2. +2.*vy[i-1, *,*]*dt/2.
    vyhalf = vy[i-1, *,*] + gforcey*dt/2. -2.*vx[i-1, *,*]*dt/2.
    vzhalf = vz[i-1, *,*] + gforcez*dt/2.


    xhalf  =  x[i-1, *,*] + vx[i-1, *,*]*dt/2
    yhalf  =  y[i-1, *,*] + vy[i-1, *,*]*dt/2
    zhalf  =  z[i-1, *,*] + vz[i-1, *,*]*dt/2

    gforcehalfx =3*xhalf + (m*(-(L/2) + xhalf*y0 + yhalf*y1 + y2*zhalf + Sqrt((xhalf*x0 + x1*yhalf + x2*zhalf)^2 + (L/2 - xhalf*y0 - yhalf*y1 - y2*zhalf)^2 + (xhalf*z0 + yhalf*z1 + zhalf*z2)^2))*((y0 + (2*x0*(xhalf*x0 + x1*yhalf + x2*zhalf) + 2*y0*(L/2 + xhalf*y0 + yhalf*y1 + y2*zhalf) + 2*z0*(xhalf*z0 + yhalf*z1 + zhalf*z2))$
      /(2*Sqrt((xhalf*x0 + x1*yhalf + x2*zhalf)^2 + (L/2 + xhalf*y0 + yhalf*y1 + y2*zhalf)^2 + (xhalf*z0 + yhalf*z1 + zhalf*z2)^2)))/(-(L/2) + xhalf*y0 + yhalf*y1 +y2*zhalf + Sqrt((xhalf*x0 + x1*yhalf + x2*zhalf)^2 + (L/2 - xhalf*y0 - yhalf*y1 - y2*zhalf)^2 + (xhalf*z0 + yhalf*z1 + zhalf*z2)^2))$
      - ((y0 + (2*x0*(xhalf*x0 + x1*yhalf + x2*zhalf) - 2*y0*(L/2 - xhalf*y0 - yhalf*y1 - y2*zhalf) + 2*z0*(xhalf*z0 + yhalf*z1 + zhalf*z2))/(2*Sqrt((xhalf*x0 + x1*yhalf + x2*zhalf)^2 + (L/2 - xhalf*y0 - yhalf*y1 - y2*zhalf)^2 + (xhalf*z0 + yhalf*z1 + zhalf*z2)^2)))$
      *(L/2 + xhalf*y0 + yhalf*y1 + y2*zhalf + Sqrt((xhalf*x0 + x1*yhalf + x2*zhalf)^2 + (L/2 + xhalf*y0 + yhalf*y1 + y2*zhalf)^2 + (xhalf*z0 + yhalf*z1 + zhalf*z2)^2)))$
      /(-(L/2) + xhalf*y0 + yhalf*y1 + y2*zhalf + Sqrt((xhalf*x0 + x1*yhalf + x2*zhalf)^2 + (L/2 - xhalf*y0 - yhalf*y1 - y2*zhalf)^2 + (xhalf*z0 + yhalf*z1 + zhalf*z2)^2))^2))/(L/2 + xhalf*y0 + yhalf*y1 +y2*zhalf + Sqrt((xhalf*x0 + x1*yhalf + x2*zhalf)^2 + (L/2 + xhalf*y0 + yhalf*y1 + y2*zhalf)^2 + (xhalf*z0 + yhalf*z1 + zhalf*z2)^2))

    gforcehalfy = (m*(-(L/2) + xhalf*y0 + yhalf*y1 + y2*zhalf + Sqrt((xhalf*x0 + x1*yhalf + x2*zhalf)^2 + (L/2 - xhalf*y0 - yhalf*y1 - y2*zhalf)^2 + (xhalf*z0 + yhalf*z1 + zhalf*z2)^2)) $
      *((y1 + (2*x1*(xhalf*x0 + x1*yhalf + x2*zhalf) + 2*y1*(L/2 + xhalf*y0 + yhalf*y1 + y2*zhalf) + 2*z1*(xhalf*z0 + yhalf*z1 + zhalf*z2))/(2*Sqrt((xhalf*x0 + x1*yhalf + x2*zhalf)^2 + (L/2 + xhalf*y0 + yhalf*y1 + y2*zhalf)^2 + (xhalf*z0 + yhalf*z1 + zhalf*z2)^2)))$
      /(-(L/2) + xhalf*y0 + yhalf*y1 + y2*zhalf + Sqrt((xhalf*x0 + x1*yhalf + x2*zhalf)^2 + (L/2 - xhalf*y0 - yhalf*y1 - y2*zhalf)^2 + (xhalf*z0 + yhalf*z1 + zhalf*z2)^2)) $
      - ((y1 + (2*x1*(xhalf*x0 + x1*yhalf + x2*zhalf) - 2*y1*(L/2 - xhalf*y0 - yhalf*y1 - y2*zhalf) + 2*z1*(xhalf*z0 + yhalf*z1 + zhalf*z2))/(2*Sqrt((xhalf*x0 + x1*yhalf + x2*zhalf)^2 + (L/2 - xhalf*y0 - yhalf*y1 - y2*zhalf)^2 + (xhalf*z0 + yhalf*z1 + zhalf*z2)^2)))$
      *(L/2 + xhalf*y0 + yhalf*y1 + y2*zhalf + Sqrt((xhalf*x0 + x1*yhalf + x2*zhalf)^2 + (L/2 + xhalf*y0 + yhalf*y1 + y2*zhalf)^2 + (xhalf*z0 + yhalf*z1 + zhalf*z2)^2)))/(-(L/2) + xhalf*y0 + yhalf*y1 + y2*zhalf + Sqrt((xhalf*x0 + x1*yhalf + x2*zhalf)^2 + (L/2 - xhalf*y0 - yhalf*y1 - y2*zhalf)^2 + (xhalf*z0 + yhalf*z1 + zhalf*z2)^2))^2))$
      /(L/2 + xhalf*y0 + yhalf*y1 + y2*zhalf + Sqrt((xhalf*x0 + x1*yhalf + x2*zhalf)^2 + (L/2 + xhalf*y0 + yhalf*y1 + y2*zhalf)^2 + (xhalf*z0 + yhalf*z1 + zhalf*z2)^2))

    gforcehalfz = (m*(-(L/2) + xhalf*y0 + yhalf*y1 + y2*zhalf + Sqrt((xhalf*x0 + x1*yhalf + x2*zhalf)^2 + 1./4.*(L - 2*(xhalf*y0 + yhalf*y1 + y2*zhalf))^2 + (xhalf*z0 + yhalf*z1 + zhalf*z2)^2))$
      *(-(((L/2 + xhalf*y0 + yhalf*y1 + y2*zhalf + Sqrt((xhalf*x0 + x1*yhalf + x2*zhalf)^2 + (L/2 + xhalf*y0 + yhalf*y1 + y2*zhalf)^2 + (xhalf*z0 + yhalf*z1 + zhalf*z2)^2))*(y2 + (-L*y2 + 2*x2*(xhalf*x0 + x1*yhalf + x2*zhalf) + 2*y2*(xhalf*y0 + yhalf*y1 + y2*zhalf) + 2*z2*(xhalf*z0 + yhalf*z1 + zhalf*z2))$
      /(2*Sqrt((xhalf*x0 + x1*yhalf + x2*zhalf)^2 + 1/4.*(L - 2*(xhalf*y0 + yhalf*y1 + y2*zhalf))^2 + (xhalf*z0 + yhalf*z1 + zhalf*z2)^2))))$
      /(-(L/2) + xhalf*y0 + yhalf*y1 + y2*zhalf + Sqrt((xhalf*x0 + x1*yhalf + x2*zhalf)^2 + 1/4.*(L - 2*(xhalf*y0 + yhalf*y1 + y2*zhalf))^2 + (xhalf*z0 + yhalf*z1 + zhalf*z2)^2))^2) + (y2 + (x2*(xhalf*x0 + x1*yhalf + x2*zhalf) + y2*(L/2 + xhalf*y0 + yhalf*y1 + y2*zhalf) + z2*(xhalf*z0 + yhalf*z1 + zhalf*z2))$
      /Sqrt((xhalf*x0 + x1*yhalf + x2*zhalf)^2 + (L/2 + xhalf*y0 + yhalf*y1 + y2*zhalf)^2 + (xhalf*z0 + yhalf*z1 + zhalf*z2)^2))$
      /(-(L/2) + xhalf*y0 + yhalf*y1 + y2*zhalf + Sqrt((xhalf*x0 + x1*yhalf + x2*zhalf)^2 + 1./4.*(L - 2*(xhalf*y0 + yhalf*y1 + y2*zhalf))^2 + (xhalf*z0 + yhalf*z1 + zhalf*z2)^2))))$
      /(L/2 + xhalf*y0 + yhalf*y1 + y2*zhalf + Sqrt((xhalf*x0 + x1*yhalf + x2*zhalf)^2 + (L/2 + xhalf*y0 + yhalf*y1 + y2*zhalf)^2 + (xhalf*z0 + yhalf*z1 + zhalf*z2)^2))





    vx[i,*,*] = vx[i-1,*,*] + gforcehalfx * dt + 2*vyhalf * dt
    x[i, *,*] = x[i-1,*,*] + vxhalf * dt

    vy[i, *,*] = vy[i-1, *,*] + gforcehalfy * dt - 2*vxhalf * dt
    y[i, *,*] =   y[i-1, *,*] + vyhalf * dt

    vz[i, *,*] = vz[i-1, *,*] + gforcehalfz * dt
    z[i, *,*] =  z[i-1, *,*] + vzhalf * dt



  endfor

  if phi lt 0 then phi = !dpi + phi
  q0 = (x*x0 + y*x1 + z*x2)*rh
  q1 = (x*y0 + y*y1 + z*y2)*rh
  v0 = (vx*x0 + vy*x1 + vz*x2)*rh*omega
  ;v0 = (vx*x0 + vy*x1+ vz*x2)*rh*omega
  v1 = (vx*y0 + vy*y1 + vz*y2)*rh*omega


  ;phi1


  junk =  where_xyz(abs(q1) lt lc/2 and abs(q0) lt w/2., xind=xind, yind=yind,zind=zind, n)
  xind = xind[where(-shift(yind,1) + yind gt 0)]
  yind = yind[where(-shift(yind,1) + yind gt 0)]
  ;zind = zind[where(-shift(yind,1) + yind gt 0)]





  ; fudge1z = shift(fudge1z,1)


  ;index = where(v0 gt 0)

  ;bs1 = bs1[index]
  ;v0 = v0[index]
  ;q1 = q1[index]

  ; z[where(z lt H/2*z2)
  ; fudge1z = [2*tsum(intz,Exp(-(x2*(W/2)/(std/2) + intz*z2/(std/2))^2)),

  ;toc

  if n_elements(yind) lt 2 then begin

    if ~finite(1/tan(tetax)^4) then tetax =0.001d
    torquez = -4* rho_mean * Abs( omegaz1 + dhdot*y0*x2) * ( omegaz1 + dhdot*y0*x2) *fudge1z  * std^2 * (std^2 - Exp(-((l * tan(tetax))/std)^2) * (std^2 + l^2 * tan(tetax)^2))/(32*tan(tetax)^4)

  endif else begin

    bs1 = bs[yind,*]*rh
    findex = where(shift(bs1,-1)-bs1 lt -1,n)
    if n eq 0 then findex=n_elements(yind) -1
    bs2 = bs1[0:findex[0]]
    ;v0 = v0[*,*,0]
    ;q1 = q1[*,*,0]
    v0 = v0[xind[0:findex[0]],yind[0:findex[0]],zind[0:findex[0]]]
    q1 = q1[xind[0:findex[0]],yind[0:findex[0]],zind[0:findex[0]]]
    vzrecord[ind] = mean(abs(vz[xind[0:findex[0]],yind[0],*]))*rh*omega

    zin   =  z[*,yind[0:findex[0]],*]*rh
    zin   = reform(zin[0,*,*])
    ;  z    = z[xind[0:findex[0]],0:-1,*]
    ;z2   = for j = 0, n_elements(xind[0:findex[0]]) = z[xind[0:findex[0]],yind[0:findex[0]],*]*rh
    ;
    z00   = z[*,*,0]*rh
    z11   = z[*,*,1]*rh
    z22   = z[*,*,2]*rh
    z3   = z[*,*,3]*rh
    z4   = z[*,*,4]*rh
    z5   = z[*,*,5]*rh
    z6   = z[*,*,6]*rh
    z7   = z[*,*,7]*rh
    z8   = z[*,*,8]*rh

    n = n_elements(bs2)
    zf = make_array(n_elements(bs2),9)


    zf[*,0]   = z00[xind[0:n-1],yind[0:n-1]]
    zf[*,1]   = z11[xind[0:n-1],yind[0:n-1]]
    zf[*,2]   = z22[xind[0:n-1],yind[0:n-1]]
    zf[*,3]   = z3[xind[0:n-1],yind[0:n-1]]
    zf[*,4]   = z4[xind[0:n-1],yind[0:n-1]]
    zf[*,5]   = z5[xind[0:n-1],yind[0:n-1]]
    zf[*,6]   = z6[xind[0:n-1],yind[0:n-1]]
    zf[*,7]   = z7[xind[0:n-1],yind[0:n-1]]
    zf[*,8]   = z8[xind[0:n-1],yind[0:n-1]]


    heff = make_array(n_elements(bs2))
    for j=0, n_elements(bs2)-1 do begin

      junk = min(abs(zf[j,*] - h/2.),indi)
      heff[j] = zin[j,indi]
    endfor
    fudge1z = abs([std/2*sqrt(!dpi)/2*erf(x2*W/std + heff*z2/(std/2.)) - std/2*sqrt(!dpi)/2*erf(x2*W/std)  + abs(std/2*sqrt(!dpi)/2*erf(x2*W/std  -heff*z2/(std/2.))) - std/2*sqrt(!dpi)/2*erf(x2*W/std)])
    ;fudge1z = std*sqrt(!dpi)/2*erf(x2*W/std + heff*z2/(std/2.))
    fudge = Exp(-(q1*tan(tetax)/(std/2))^2)

    if ~finite(1/tan(tetax)^4) then tetax =0.001d
    ; torquex =  4*fudge1x*tsum(bs[yind], abs(v0[xind,yind]*(-signum(x2))*(abs(cos(phi)*z0) + abs(sin(phi)*z1)) + (dhdot*y0*z2 - omegax1) * q1[xind,yind]/abs(cos(tetax))) * (v0[xind,yind]*(-signum(x2))*(abs(cos(phi)*z0) + abs(sin(phi)*z1))     + (dhdot*y0*z2 - omegax1) * q1[xind,yind]/abs(cos(tetax))) * rho_mean * fudge * q1[xind,yind]/abs(cos(tetax)))
    ; torquez = -4*tsum(bs2, fudge1z*abs(v0 * signum(z2) + (dhdot*y0*x2 + omegaz1) * q1/abs(cos(tetax))) * (v0* signum(z2)     + (dhdot*y0*x2 + omegaz1) * q1/abs(cos(tetax))) * rho_mean * fudge * q1/abs(cos(tetax)))

    torquez = -4*tsum(bs2, fudge1z*abs(v0*signum(z2)*signum(x1) + (dhdot*y0*x2 + omegaz1) * q1) * (v0* signum(z2)*signum(x1) + (dhdot*y0*x2 + omegaz1) * q1) * rho_mean * fudge * q1)

    ;print,'b',max(bs2[where(q1 gt 0,num)])
    ;print,'n',num
  endelse

  q0 = (x*z0 + y*z1 + z*z2)*rh
  q1 = (x*y0 + y*y1 + z*y2)*rh
  v0 = (vx*z0 + vy*z1 + vz*z2)*rh*omega
  ;v0 = (vx*x0 + vy*x1 +vz*x2)*rh*omega
  v1 = (vx*y0 + vy*y1 + vz*y2)*rh*omega

  junk =  where_xyz(abs(q1) lt lc/2 and abs(q0) lt h/2, xind=xind, yind=yind,zind=zind, n)
  xind = xind[where(-shift(yind,1) + yind gt 0)]
  yind = yind[where(-shift(yind,1) + yind gt 0)]

  if n_elements(yind) lt 2 then begin
    if ~finite(1/tan(tetax)^4) then tetax =0.001d
    torquex =  4* rho_mean * Abs(-omegax1 + dhdot*y0*z2) * (-omegax1 + dhdot*y0*z2) * fudge1x  * std^2 * (std^2 - Exp(-((l * tan(tetax))/std)^2) * (std^2 + l^2 * tan(tetax)^2))/(32*tan(tetax)^4) ;(l/2)^2

  endif else begin
    if ~finite(1/tan(tetax)^4) then tetax =0.001d

    bs1    = bs[yind,*]*rh
    findex = where(shift(bs1,-1)-bs1 lt -1,n)
    if n eq 0 then findex=n_elements(yind) -1

    bs2    = bs1[0:findex[0]]
   ; v0     = v0[*,*,0]
   ; q1     = q1[*,*,0]
    v0     = v0[xind[0:findex[0]],yind[0:findex[0]],zind[0:findex[0]]] ;Zind wast here on the one that worked first
    q1     = q1[xind[0:findex[0]],yind[0:findex[0]],zind[0:findex[0]]] ;Zind wast here on the one that worked first

    fudge = Exp(-(q1*tan(tetax)/(std/2))^2)



    torquex =  4*fudge1x*tsum(bs2, abs(v0*signum(x2)*(-signum(z1)) + (dhdot*y0*z2 - omegax1) * q1) * (v0*signum(x2)*(-signum(z1)) + (dhdot*y0*z2 - omegax1) * q1) * rho_mean * fudge * q1)



    ;print, torquez
  endelse

  q0 = (x*z0 + y*z1 + z*z2)*rh
  q1 = (x*x0 + y*x1 + z*x2)*rh

  ;q0 = (x*z0 + y*z1)*rh
  ;q1 = (x*x0 + y*x1)*rh

  v0 = (vx*z0 + vy*z1 + vz*z2)*rh*omega
  v1 = (vx*x0 + vy*x1 + vz*x2)*rh*omega


  junk =  where_xyz(abs(q1) lt w/2 and abs(q0) lt h/2,xind=xind, yind=yind, zind=zind, n)
  xind = xind[where(-shift(yind,1) + yind gt 0)]
  yind = yind[where(-shift(yind,1) + yind gt 0)]
  ; zind = zind[where(-shift(yind,1) + yind gt 0)]



  if n_elements(yind) lt 2 then begin
    if ~finite(1/tan(tetay)^4) then tetay =0.001d
    torquey =  -4 * rho_mean * Abs( omegay1 + dhdot*x0*z2) *  (omegay1 + dhdot*x0*z2) * fudge1y * std^2 * (std^2 - Exp(-(w*tan(tetay)/std)^2) * (std^2 + w^2 * tan(tetay)^2))/(32*tan(tetay)^4)
  endif else begin
    if ~finite(1/tan(tetay)^4) then tetay =0.001d

    bs1 = bs[yind,*]*rh
    findex = where(shift(bs1,-1)-bs1 lt -1,n)
    if n eq 0 then findex=n_elements(yind) -1
    bs2 = bs1[0:findex[0]]
 
    v0 = v0[xind[0:findex[0]],yind[0:findex[0]],yind[0:findex[0]],zind[0:findex[0]]]
    q1 = q1[xind[0:findex[0]],yind[0:findex[0]],yind[0:findex[0]],zind[0:findex[0]]]

    fudge = Exp(-(q1*tan(tetax)/(std/2))^2)

    ;v0 = v0[xind[0:findex[0]],yind[0:findex[0]],*]
    ;q1 = q1[xind[0:findex[0]],yind[0:findex[0]],*]




    fudge = Exp(-(q1*tan(tetay)/(std/2))^2)
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;bs[yind]*rh WASNT THE CASE FOR THE SUBMITTED PAPER. RATHER, just bs[yind]
    torquey =  -4*fudge1y*tsum(bs2, abs(v0*signum(z1)*signum(y2) + (dhdot*x0*z2 + omegay1) * q1) * (v0*signum(z1)*signum(y2) + (dhdot*x0*z2 + omegay1) * q1) * rho_mean * fudge * q1)

    ;torquey =  -4*fudge1y*tsum(bs2, abs(v0*signum(z1)*signum(y2) + (dhdot*x0*z2 + omegay1) * q1/abs(cos(tetay))) * (v0*signum(z1)*signum(y2) + (dhdot*x0*z2 + omegay1) * q1/abs(cos(tetay))) * rho_mean * fudge * q1/abs(cos(tetay)))


    ; torquey =  -4*fudge1y*tsum(bs2, abs(v0*signum(y2)*(abs(sin(phi)*z0) + abs(cos(phi)*z1)) + (dhdot*x0*z2 + omegay1) * q1/abs(cos(tetay))) * (v0*signum(y2)*(abs(sin(phi)*z0) + abs(cos(phi)*z1)) + (dhdot*x0*z2 + omegay1) * q1/abs(cos(tetay))) * rho_mean * fudge * q1/abs(cos(tetay)))
    ;torquey =  -4*fudge1y*tsum(bs2, abs(v0*signum(y2) + (dhdot*x0*z2 + omegay1) * q1/abs(cos(tetay))) * (v0*signum(y2) + (dhdot*x0*z2 + omegay1) * q1/abs(cos(tetay))) * rho_mean * fudge * q1/abs(cos(tetay)))

  endelse


  ;  toc
  ;plothill,x,y,bs;,norange=1

  ;if torquez gt 1e7 then begin
  ; print,'stop'
  ;endif
  ;  print,torquez
  ;  print,'phi',phi*180/!dpi
  return,[torquey,torquex,torquez]
end