;Function that runs the dynamical torque program (wake_rot27.pro) for 0.6 periods, in a flat ring (no bending wave) in order to see at what angles does the wake equilibrate.

;Inputs
;s ---> shear rate
;rho_mean ---> mean spatial density of particles in the ring.
;l,h,w ---> dimensions of the self-gravity wake (length, highth, width)
;omega ---> orbital frequency
;phi (optional)---> initial orientation
;
;Output 
;Phis --> mean of the phi value in the last 6 timesteps of the simulation.
;
;Change Log
;10/22/2024 -> Created by Daniel Sega

function equilibrium_H, s, rho_mean,l,h,w,omega, phi=phi

tic
flag  = 1 ;to change initial conditions,1 is such that wakes start flat at a trhough. 0 is suh that the wake starts with no omega but inclined as the amplitude of the slope of the wave
flag1 = 1 ;to turn on or off the wake potential. 1 is off
;phi=25*!dpi/180
sigma = 492. ;* 0.5 * sqrt(!dpi)
lambda_c = 58.
k = 2*!dpi/lambda_c
;h   = 4d;4
;l   = 232d
;W   = 18d;18

n0    = 9.5
a0    = 10e-2
amax  = 9.5
amin  = 11e-3
q     = -3

;rho_p = 0.4
rho_p = 400.

as = dindgen(10000, start = amin, increment = amax/10000)

n      = n0*(as/a0)^q
a_mean = total(n*as)/total(as)

rho_dist = rho_p * !dpi* (4./3.) * as^3 *  n0*(as/a0)^q

j=dcomplex(0,1)

;taux = 5.33333*10^7 Sin[(30*!dpi/180 + 0.0980198 thetax*thetay]
;tauy = 533333. Cos[(5 \[Pi])/36 + thetax*thetay]

rv  = 131902000.
G   = 6.67e-11
Ms  = 5.6832e26
;Omega = 0.000128897
rho = 492.
mu = 1.2957087e-4

;s = 1.
; s = .5
if ~keyword_set(phi) then phi = 25./180. *!dpi

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

Ixx = 2.66773e7
Iyy = 277333.
Izz = 2.69333e7

omega = mu

period = 48492d
; dt= .01d
;dt= 48.49d
 dt = 3250d

Ts = 0.6
time = period/dt*Ts


db = .004d
startb = -13    ;normally 13, 5 for l=77
b = dindgen(abs(startb)/db ,start = startb, increment = db)
intw = linspace(-w/2,w/2,50)
inth = linspace(-.7*W,.7*W,30)
intl = linspace(-l/2.,l/2,60)


alphax = make_array(time + 1d, /DOUBLE)
alphay = make_array(time + 1d, /DOUBLE)
alphaz = make_array(time + 1d, /DOUBLE)

omegax = make_array(time + 1d, /DOUBLE)
omegay = make_array(time + 1d, /DOUBLE)
omegaz = make_array(time + 1d, /DOUBLE)

thetax = make_array(time + 1d, /DOUBLE)
thetay = make_array(time + 1d, /DOUBLE)
thetaz = make_array(time + 1d, /DOUBLE)

lo = l ;lo is for the purpose of distringuishing between the L of the emmited potential and the L of the wake that feels the potential. The theoretical right thing is that al the L's are the SAME


slope = 0d
m = rho*H*l*W

Ixx = 1/12.*m*(H^2 + l^2) ; 2.66773e7
Iyy = 1/12.*m*(H^2 + W^2) ;  277333.
Izz = 1/12.*m*(W^2 + l^2) ;2.69333e7 for l=200


seconds = findgen(time,increment = dt)

; dhdot0 =  omega*(slope*cos(0))
dhdot0 =  -omega*(slope*sin(0))

omegay[1] = 0.
thetay[1] =slope*cos(phi) ;omegay[1]*dt                                                                  ;

omegax[1] = 0.
thetax[1] = -slope*sin(phi) ; omegax[1]*dt
thetay[1] =  -slope*cos(phi)*cos(thetax[1]) ;omegay[1]*dt                                                          ;

omegaz[1] = ((Ixx  -Iyy)/Izz) * omegax[1]*omegay[1]/omega ; 

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


  if y[2] gt .99 then y[2] = .99d
  if y[2] lt -.99 then y[2] = -.99d
  if x[2] gt .99 then x[2] = .99d
  if x[2] lt -.99 then x[2] = -.99d
  if z[2] gt .99 then z[2] = .99d
  if z[2] lt -.99 then z[2] = -.99d



  if x[1] gt .99 then x[1] = .99d
  if x[1] lt -.99 then x[1] = -.99d
  if y[1] gt .99 then y[1] = .99d
  if y[1] lt -.99 then y[1] = -.99d
  if z[1] gt .99 then z[1] = .99d
  if z[1] lt -.99 then z[1] = -.99d


  if x[0] gt .99 then x[0] = .99d
  if x[0] lt -.99 then x[0] = -.99d
  if z[0] gt .99 then z[0] = .99d
  if z[0] lt -.99 then z[0] = -.99d
  if y[0] gt .99 then y[0] = .99d
  if y[0] lt -.99 then y[0] = -.99d

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



  torqueyhalf  = -2*dhdotdot*(rho*h*l*w^3/24.) * z[2]*x[0] /(Iyy) * dt/2. - 2 * muminus*(rho*h*l*w^3/24.) * x[0]*z[0] / (Iyy) * dt/2. + 2 * muplus*(rho*h*l*w^3/24.) * x[1]*z[1] / (Iyy) * dt/2. + 4*(!dpi)*G*(sigma)*(k)* H * (rho*h*l*w*w*w/24.) * Sin(Phi) * Cos(Phi) * (z[1]*x[0] + z[0]*x[1]) / (Iyy) * dt/2.
  torquexhalf  =  2*dhdotdot*(rho*h*w*l^3/24.) * z[2]*y[0] /(Ixx) * dt/2. + 2 * muminus*(rho*h*w*l^3/24.) * y[0]*z[0] / (Ixx) * dt/2. - 2 * muplus*(rho*h*w*l^3/24.) * y[1]*z[1] / (Ixx) * dt/2. - 4*(!dpi)*G*(sigma)*(k)* H * (rho*h*w*l*l*lo/24.) * Sin(Phi) * Cos(Phi) * (z[1]*y[0] + z[0]*y[1]) / (Ixx) * dt/2.
  torquezhalf  = -2*dhdotdot*(rho*h*w*l^3/24.) * x[2]*y[0] /(Izz) * dt/2. - 2 * muminus*(rho*h*w*l^3/24.) * y[0]*x[0] / (Izz) * dt/2. + 2 * muplus*(rho*h*w*l^3/24.) * x[1]*y[1] / (Izz) * dt/2. + 4*(!dpi)*G*(sigma)*(k)* H * (rho*h*w*l*l*lo/24.) * Sin(Phi) * Cos(Phi) * (x[1]*y[0] + x[0]*y[1]) / (Izz) * dt/2.


  inertiayhalf = ((Izz - Ixx)/Iyy)*omegaz[i-1]*omegax[i-1]*dt/2.
  inertiaxhalf = ((Iyy - Izz)/Ixx)*omegay[i-1]*omegaz[i-1]*dt/2.
  inertiazhalf = ((Ixx - Iyy)/Izz)*omegay[i-1]*omegax[i-1]*dt/2.

  stub = hillcol4deq2(b, l,w,h,s,rho,rho_mean, x[0], x[1], x[2], y[0], y[1], y[2], z[0], z[1], z[2], omegax[i-1],omegay[i-1],omegaz[i-1], ind=i, phis=phis, dhdot = dhdot, thetays = thetays, thetaxs=thetaxs, intl = IntL,intw=intw,inth=inth,vzrecord=vzrecord)

  dragyhalf    =    stub[0]*dt/(2*Iyy) ;- torquecol(l,w,h,s,rho_mean, y[0], y[1], y[2], omegarot = omegaz[i-1], z2 = z[2], z1 = z[1],        ind=i, phis=phis, oldphi = phis[i-1])*z[2]*dt/(2*Izz)
  dragxhalf    =    stub[1]*dt/(2*Ixx)          ;note that delta v is positive if particles are "going up" (positive z coordinate) wrt the particles in the wake positive y coordinate
  dragzhalf    =    stub[2]*dt/(2*Izz)   ;* sin(thetaz[i-1] + phi);*cos(thetay[i-1])*cos(thetax[i-1])*(sin(thetay[i-1]) + cos(thetax[i-1]))

  ; print, torquezhalf/dt*Izz
  ; print,'col', dragzhalf/dt*Izz
  ;print,phis[i]*180/!dpi


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




  torquey = - 2*dhdotdothalf*(rho*h*l*w*w*w/24.) * zhalf[2] * xhalf[0]/(Iyy) * dt - 2 * muminus*(rho*h*l*w^3/24.) * xhalf[0]*zhalf[0] / (Iyy) * dt + 2 * muplus*(rho*h*l*w^3/24.) * xhalf[1]*zhalf[1] / (Iyy) * dt + 4*(!dpi)*G*(sigma)*(k)* H * (rho*h*lo*w*w*w/24.) * Sin(Phi) * Cos(Phi) * (zhalf[1]*xhalf[0] + zhalf[0]*xhalf[1]) / (Iyy) * dt
  torquex =   2*dhdotdothalf*(rho*h*w*l*l*l/24.) * zhalf[2] * yhalf[0]/(Ixx) * dt + 2 * muminus*(rho*h*w*l^3/24.) * yhalf[0]*zhalf[0] / (Ixx) * dt - 2 * muplus*(rho*h*w*l^3/24.) * yhalf[1]*zhalf[1] / (Ixx) * dt - 4*(!dpi)*G*(sigma)*(k)* H * (rho*h*w*lo*l*l/24.) * Sin(Phi) * Cos(Phi) * (zhalf[1]*yhalf[0] + zhalf[0]*yhalf[1]) / (Ixx) * dt
  torquez = - 2*dhdotdothalf*(rho*h*w*l*l*l/24.) * xhalf[2] * yhalf[0]/(Izz) * dt - 2 * muminus*(rho*h*w*l^3/24.) * yhalf[0]*xhalf[0] / (Izz) * dt + 2 * muplus*(rho*h*w*l^3/24.) * xhalf[1]*yhalf[1] / (Izz) * dt + 4*(!dpi)*G*(sigma)*(k)* H * (rho*h*w*lo*l*l/24.) * Sin(Phi) * Cos(Phi) * (xhalf[1]*yhalf[0] + xhalf[0]*yhalf[1]) / (Izz) * dt

  inertiay = ((Izz - Ixx)/Iyy) * omegazhalf*omegaxhalf*dt
  inertiax = ((Iyy - Izz)/Ixx) * omegayhalf*omegazhalf*dt
  inertiaz = ((Ixx  -Iyy)/Izz) * omegayhalf*omegaxhalf*dt

  stub = hillcol4deq2(b, l,w,h,s,rho,rho_mean, xhalf[0], xhalf[1], xhalf[2], yhalf[0], yhalf[1], yhalf[2], zhalf[0], zhalf[1], zhalf[2], omegaxhalf, omegayhalf, omegazhalf, ind=i, phis=phis, dhdot = dhdot, thetays = thetays, thetaxs=thetaxs, intl = IntL,intw=intw,inth=inth,vzrecord=vzrecord)

  dragy    =  stub[0]*dt/(Iyy)
  dragx    =  stub[1]*dt/(Ixx)                  ;note that delta v is positive if particles are "going up" (positive z coordinate) wrt the particles in the wake positive y coordinate
  dragz    =  stub[2]*dt/(Izz)


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


endfor
toc
;cgplot,seconds*(2*!dpi/omega)^(-1),   phis*180/!dpi,color= cgcolor('green')

;save,phis,filename='wakeflatring.sav'
return, mean(phis[-6:-1])
end