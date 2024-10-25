Pro figure10
;Description
;Program that runs wpl232dt200 for 1 period runs for a number of slope values. It produces figure 10 in Sega et al 2024. The wake dimensions are L=232m, W=18m, H=4m and the timesstep is 200 seconds (twice as much as in figure 8 and 9's wpl232dt97.pro andn nwpl232dt97.pro) 


  h   = 4.
  l   = 232.
  W   = 18.
  by =  linspace(38d, l/2 +10, 10000)
  b = linspace(38d, l/2 +10, 10000)
  intw = linspace(-w/2,w/2,40)
  inth = linspace(-h/2.,h/2,20)
  intl = linspace(-l/2.,l/2,60)


  slopes = [0.001d,0.0025d, 0.01d, 0.02d,.03d,.04d ,0.05d,.06d ,0.07d,.08d,.09d, .1d,.11d,.12d,.13d,.14d,.15d,.16d,.17d, .18d,.19d, .2d,.21d, .22d,.23d,.24d, .25d, .26d,.27d, .28d,.29d, .3d]
  ;  slopes = [0,0.01, 0.03, 0.04,.05,.07,.09, .1,.13,.15,.17,.19,.21,.23, .25,.27,.28]
  thetas = linspace(0,!dpi,15)
  ;thetas = thetas[-1]
  flag = 1 ;to change initial conditions,1 is such that wakes start flat at a trhough. 0 is such that the wake starts with no omega but inclined as the amplitude of the slope
  if ~keyword_set(thetas) then thetas = [0]

  flag1 = 0 ;to turn on or off the wake potential. 1 is off

  Cstimes  = make_array(n_elements(thetas))
  Csdeltavs = make_array(n_elements(thetas))

  Cstimesk2  = make_array(n_elements(thetas))
  Csdeltavsk2 = make_array(n_elements(thetas))

  Cstimesk  = make_array(n_elements(thetas))
  Csdeltavsk = make_array(n_elements(thetas))

  Cstimesk1  = make_array(n_elements(thetas))
  Csdeltavsk1 = make_array(n_elements(thetas))

  result      = make_array(n_elements(thetas))

  rms      = make_array(n_elements(thetas), n_elements(slopes))
  modes      = make_array(n_elements(thetas), n_elements(slopes))
  rms2      = make_array(n_elements(thetas), n_elements(slopes))
  modes2     = make_array(n_elements(thetas), n_elements(slopes))

  ind2 = 0

  Foreach theta, thetas do begin

    deltavs  = make_array(n_elements(slopes))
    times    = make_array(n_elements(slopes))
    deltavsk2 = make_array(n_elements(slopes))
    timesk2   = make_array(n_elements(slopes))
    deltavsk = make_array(n_elements(slopes))
    deltavsk1 = make_array(n_elements(slopes))
    timesk   = make_array(n_elements(slopes))
    timesk1   = make_array(n_elements(slopes))
    ind = 0

    Foreach slope, slopes do begin

      ;ORIENTATION AND KEPLERIAM SHEAR AND WAKEPOTENTIAL IN CARTISIAN
      ;This is the ones with the torque bieng a sin and EVERYTHING including a proper treatement of the orientation of the wake with iterative infenitesimal rotations (which are commutativ and make the problem much easier)
      ;We take a mean particle size in the gap Using Richard Jerousek 2018 (thesis)

      rho = 492.
      sigma = 492. * 0.5 * sqrt(!dpi)
      lambda_c = 58.
      k = 2*!dpi/lambda_c
      h   = 4.
      l   = 232.
      W   = 18.

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

      rho_mean =   385*1.3
      ;rho_mean = mean(rho_dist)

      if flag1 then begin
        rho_mean = 385.
      endif

      j=dcomplex(0,1)

      ;taux = 5.33333*10^7 Sin[(30*!dpi/180 + 0.0980198 thetax*thetay]
      ;tauy = 533333. Cos[(5 \[Pi])/36 + thetax*thetay]

      rv  = 131902000.
      G   = 6.67e-11
      Ms  = 5.6834e26
      Omega = 0.000128897
      rho = 492.
      mu = 1.2957087e-4

      s = 1.50524
      ; s = 1.
      ; s = .5
      phi = 25./180. *!dpi

      ;vkdr = -0.0000642821d*3
      vkdr = -s*Omega


      muminus = 2*G*Ms/(rv^3) - 2*(!dpi)*G*(sigma * k )* H*(cos(phi)^2)
      muplus  =   G*Ms/(rv^3) + 2*(!dpi)*G*(sigma * k )* H*(sin(phi)^2)

      if flag1 then begin

        muminus = 2*G*Ms/(rv^3) ;- 2*(!dpi)*G*(sigma * k )* H*(cos(phi)^2)
        muplus  =   G*Ms/(rv^3) ;+ 2*(!dpi)*G*(sigma * k )* H*(sin(phi)^2)

        sigma = 0

      endif


      omega = mu

      period = 48493d
      dt= 200d

      time = period/dt*1.01



      alphax = make_array(time + 1d, /DOUBLE)
      alphay = make_array(time + 1d, /DOUBLE)
      alphaz = make_array(time + 1d, /DOUBLE)

      omegax = make_array(time + 1d, /DOUBLE)
      omegay = make_array(time + 1d, /DOUBLE)
      omegaz = make_array(time + 1d, /DOUBLE)

      thetax = make_array(time + 1d, /DOUBLE)
      thetay = make_array(time + 1d, /DOUBLE)
      thetaz = make_array(time + 1d, /DOUBLE)



      ; slope = 0.03
      m = rho*H*l*W

      Ixx = 1/12.*m*(H^2 + l^2) ; 2.66773e7
      Iyy = 1/12.*m*(H^2 + W^2) ;  277333.
      Izz = 1/12.*m*(W^2 + l^2) ;2.69333e7 for l=200


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


        dhdot0 =  omega*(slope)*cos(theta)

        omegay[1] =  -dhdot0*cos(phi)
        thetay[1] =  -atan(sin(slope*sin(theta))*cos(phi),sqrt(cos(slope*sin(theta))^2*cos(phi)^2 + sin(phi)^2))  ;omegay[1]*dt                                                 rho_mean * W * Abs(-omegax[i-1] - dhdothalf*cos(thetax[i-1])^2*sin(thetaz[i-1] + phi)) * (-omegax[i-1] - dhdothalf*cos(thetax[i-1])^2*sin(thetaz[i-1] + phi)) * (l/2.)^4 *dt/(2*Ixx)                  ;

        omegax[1] =  -dhdot0*sin(phi)
        thetax[1] =  -atan(sin(slope*sin(theta))*sin(phi),sqrt(cos(slope*sin(theta))^2*sin(phi)^2 + cos(phi)^2)) ; omegax[1]*dt

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

      for i=2d, time do begin

        num = .95

        if y[2] gt num then y[2] = num
        if y[2] lt -num then y[2] = -num
        if x[2] gt num then x[2] = num
        if x[2] lt -num then x[2] = -num
        if z[2] gt num then z[2] = num
        if z[2] lt -num then z[2] = -num



        if x[1] gt num then x[1] = num
        if x[1] lt -num then x[1] = -num
        if y[1] gt num then y[1] = num
        if y[1] lt -num then y[1] = -num
        if z[1] gt num then z[1] = num
        if z[1] lt -num then z[1] = -num


        if x[0] gt num then x[0] = num
        if x[0] lt -num then x[0] = -num
        if z[0] gt num then z[0] = num
        if z[0] lt -num then z[0] = -num
        if y[0] gt num then y[0] = num
        if y[0] lt -num then y[0] = -num

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

          ;theta=0

          dhdothalf    =  omega  *  slope*cos(omega*dt*(i-1./2.) + theta)
          dhdot        =  omega  *  slope*cos(omega*dt*(i-1) + theta)

          dhdotdothalf = -(omega^2)*slope*sin(omega*dt*(i-1./2.) + theta)
          dhdotdot     = -(omega^2)*slope*sin(omega*dt*(i-1) + theta)

        endif



        torqueyhalf  = -2*dhdotdot*(rho*h*l*w^3/24.) * z[2]*x[0] /(Iyy) * dt/2. - 2 * muminus*(rho*h*l*w^3/24.) * x[0]*z[0] / (Iyy) * dt/2. + 2 * muplus*(rho*h*l*w^3/24.) * x[1]*z[1] / (Iyy) * dt/2. + 4*(!dpi)*G*(sigma)*(k)* H * (rho*h*l*w*w*w/24.) * Sin(Phi) * Cos(Phi) * (z[1]*x[0] + z[0]*x[1]) / (Iyy) * dt/2.
        torquexhalf  =  2*dhdotdot*(rho*h*w*l^3/24.) * z[2]*y[0] /(Ixx) * dt/2. + 2 * muminus*(rho*h*w*l^3/24.) * y[0]*z[0] / (Ixx) * dt/2. - 2 * muplus*(rho*h*w*l^3/24.) * y[1]*z[1] / (Ixx) * dt/2. - 4*(!dpi)*G*(sigma)*(k)* H * (rho*h*w*l*l*l/24.) * Sin(Phi) * Cos(Phi) * (z[1]*y[0] + z[0]*y[1]) / (Ixx) * dt/2.
        torquezhalf  = -2*dhdotdot*(rho*h*w*l^3/24.) * x[2]*y[0] /(Izz) * dt/2. - 2 * muminus*(rho*h*w*l^3/24.) * y[0]*x[0] / (Izz) * dt/2. + 2 * muplus*(rho*h*w*l^3/24.) * x[1]*y[1] / (Izz) * dt/2. + 4*(!dpi)*G*(sigma)*(k)* H * (rho*h*w*l*l*l/24.) * Sin(Phi) * Cos(Phi) * (x[1]*y[0] + x[0]*y[1]) / (Izz) * dt/2.

        ;print,x[1]


        inertiayhalf = ((Izz - Ixx)/Iyy)*omegaz[i-1]*omegax[i-1]*dt/2.
        inertiaxhalf = ((Iyy - Izz)/Ixx)*omegay[i-1]*omegaz[i-1]*dt/2.
        inertiazhalf = ((Ixx - Iyy)/Izz)*omegay[i-1]*omegax[i-1]*dt/2.

        dragyhalf    =     torquecol2(by,l,w,h,s,rho,rho_mean, x[0], x[1], x[2], y[0], y[1], y[2], z[0], z[1], z[2], omegarot = omegay[i-1],        ind=i, phis=phis, dhdot = dhdot, thetays = thetays, intl=intl  )*dt/(2*Iyy) ;- torquecol(l,w,h,s,rho_mean, y[0], y[1], y[2], omegarot = omegaz[i-1], z2 = z[2], z1 = z[1],        ind=i, phis=phis, oldphi = phis[i-1])*z[2]*dt/(2*Izz)
        dragxhalf    =     torquecol1(b, l,w,h,s,rho,rho_mean, x[0], x[1], x[2], y[0], y[1], y[2], z[0], z[1], z[2], omegarot = omegax[i-1],  xi=1, ind=i, phis=phis, dhdot = dhdot, thetaxs = thetaxs, intw = intw)*dt/(2*Ixx)          ;note that delta v is positive if particles are "going up" (positive z coordinate) wrt the particles in the wake positive y coordinate
        dragzhalf    =     torquecol1(b, l,w,h,s,rho,rho_mean, x[0], x[1], x[2], y[0], y[1], y[2], z[0], z[1], z[2], omegarot = omegaz[i-1],        ind=i, phis=phis, dhdot = dhdot, thetaxs = thetaxs, inth = inth)*dt/(2*Izz)   ;* sin(thetaz[i-1] + phi);*cos(thetay[i-1])*cos(thetax[i-1])*(sin(thetay[i-1]) + cos(thetax[i-1]))


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



        torquey = - 2*dhdotdothalf*(rho*h*l*w*w*w/24.) * zhalf[2] * xhalf[0]/(Iyy) * dt - 2 * muminus*(rho*h*l*w^3/24.) * xhalf[0]*zhalf[0] / (Iyy) * dt + 2 * muplus*(rho*h*l*w^3/24.) * xhalf[1]*zhalf[1] / (Iyy) * dt + 4*(!dpi)*G*(sigma)*(k)* H * (rho*h*l*w*w*w/24.) * Sin(Phi) * Cos(Phi) * (zhalf[1]*xhalf[0] + zhalf[0]*xhalf[1]) / (Iyy) * dt
        torquex =   2*dhdotdothalf*(rho*h*w*l*l*l/24.) * zhalf[2] * yhalf[0]/(Ixx) * dt + 2 * muminus*(rho*h*w*l^3/24.) * yhalf[0]*zhalf[0] / (Ixx) * dt - 2 * muplus*(rho*h*w*l^3/24.) * yhalf[1]*zhalf[1] / (Ixx) * dt - 4*(!dpi)*G*(sigma)*(k)* H * (rho*h*w*l*l*l/24.) * Sin(Phi) * Cos(Phi) * (zhalf[1]*yhalf[0] + zhalf[0]*yhalf[1]) / (Ixx) * dt
        torquez = - 2*dhdotdothalf*(rho*h*w*l*l*l/24.) * xhalf[2] * yhalf[0]/(Izz) * dt - 2 * muminus*(rho*h*w*l^3/24.) * yhalf[0]*xhalf[0] / (Izz) * dt + 2 * muplus*(rho*h*w*l^3/24.) * xhalf[1]*yhalf[1] / (Izz) * dt + 4*(!dpi)*G*(sigma)*(k)* H * (rho*h*w*l*l*l/24.) * Sin(Phi) * Cos(Phi) * (xhalf[1]*yhalf[0] + xhalf[0]*yhalf[1]) / (Izz) * dt

        inertiay = ((Izz - Ixx)/Iyy) * omegazhalf*omegaxhalf*dt
        inertiax = ((Iyy - Izz)/Ixx) * omegayhalf*omegazhalf*dt
        inertiaz = ((Ixx  -Iyy)/Izz) * omegayhalf*omegaxhalf*dt

        dragy    =  torquecol2(by, l,w,h,s,rho,rho_mean, xhalf[0], xhalf[1], xhalf[2], yhalf[0], yhalf[1], yhalf[2], zhalf[0], zhalf[1], zhalf[2], omegarot = omegayhalf,        ind=i, phis=phis, dhdot = dhdothalf, thetays = thetays, intl = intl)*dt/(Iyy)
        dragx    =  torquecol1(b,  l,w,h,s,rho,rho_mean, xhalf[0], xhalf[1], xhalf[2], yhalf[0], yhalf[1], yhalf[2], zhalf[0], zhalf[1], zhalf[2], omegarot = omegaxhalf, xi=1,  ind=i, phis=phis, dhdot = dhdothalf, thetaxs = thetaxs ,intw = intw)*dt/(Ixx)                  ;note that delta v is positive if particles are "going up" (positive z coordinate) wrt the particles in the wake positive y coordinate
        dragz    =  torquecol1(b,  l,w,h,s,rho,rho_mean, xhalf[0], xhalf[1], xhalf[2], yhalf[0], yhalf[1], yhalf[2], zhalf[0], zhalf[1], zhalf[2], omegarot = omegazhalf,        ind=i, phis=phis, dhdot = dhdothalf, thetaxs = thetaxs, inth = inth)*dt/(Izz)


        omegay[i] =  omegay[i-1] + torquey + inertiay + dragy
        omegax[i] =  omegax[i-1] + torquex + inertiax + dragx
        omegaz[i] =  omegaz[i-1] + torquez + inertiaz + dragz


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
      dhdot  = -omega*(slope*sin(omega*seconds))
      if flag then dhdot  = omega*(slope*cos(omega*seconds + theta))

      deltav = (-omegax + dhdot*y0s*z2s + vkdr *y0s*z1s) * (l/2.)  + (omegay + dhdot*x0s*z2s + vkdr * x0s*z1s) * (W/2.)

      deltavx =( omegaz + dhdot*y0s*x2s + vkdr *y0s*x1s) * (l/2.)

      deltavk = - omegax * z2s * l/2 + omegay * W/2* z2s + omegaz * l/2. * x2s +  dhdot * y0s * l/2 + dhdot * x0s * W/2

      deltavj = - omegax * z1s * l/2 + omegay * W/2* z1s + omegaz * l/2. * x1s +  vkdr *y0s*z1s* l/2 + vkdr  * x0s * W/2

      !p.multi=[1,1,3]
      deltavk= deltavk[0:-2]


      func   = (dhdot[0:-2]*y0s[0:-2] - omegax[0:-2]*z2s[0:-2])*l/2 + omegaz * l/2. * x2s


      zeroes = [12083*2, 12083*4, 12083*6, 12083*8, 12083*10, 12083*12, 12083*14]

      deltavsk2[ind] = max(abs(deltavk[period/dt *1: -1]),i)
      deltavs[ind] = max(abs(deltav),i)
      ; deltavs2[ind] = max(abs(deltav[0:period/dt *3]),i)
      ; deltavs3[ind] = max(abs(deltav[period/dt *3:period/dt*4]),i)
      ;deltavs[ind] = mean(sqrt(deltavk^2))
      timesk2[ind] =  sqrt(mean((deltavk[period/dt *1: -1])^2 ));seconds[i]*(2*!dpi/omega)^(-1)
      times[ind] =  sqrt(mean((deltav)^2 ))


      deltavsk1[ind] = max(abs(deltavk[0: period/dt *1]),i)
      deltavsk[ind] = max(abs(deltavk),i)
      timesk1[ind] =  sqrt(mean((deltavk[0: period/dt *1])^2 ));seconds[i]*(2*!dpi/omega)^(-1)
      timesk[ind] =  sqrt(mean((deltavk)^2 ))
      rms[ind2,ind] =    timesk1[ind]
      modes[ind2,ind] =    deltavsk1[ind]
      rms2[ind2,ind] =    timesk[ind]
      modes2[ind2,ind] =    deltavsk[ind]
      ind+=1
    endforeach
    Csdeltavs[ind2] = correlate(slopes,deltavs)
    Cstimes[ind2]  = correlate(slopes, times)
    Csdeltavsk2[ind2] = correlate(slopes,deltavsk2)
    Cstimesk2[ind2]  = correlate(slopes, timesk2)

    Csdeltavsk[ind2] = correlate(slopes,deltavsk)
    Cstimesk[ind2]  = correlate(slopes, timesk)
    Csdeltavsk1[ind2] = correlate(slopes,deltavsk1)
    Cstimesk1[ind2]  = correlate(slopes, timesk1)

    print, rms

    result[ind2] = Regress(slopes, times, CONST=const)



    ind2+=1

  endforeach


  print, deltavs
  print, deltavsk2
  ;print, deltavs2
  ;print, deltavs3

  print, 'Csdeltavs' , Csdeltavs
  print, 'Csdeltavsk2', Csdeltavsk2
  print, 'Cstimes'   , Cstimes
  print, 'Cstimesk2'  , Cstimesk2

  print, 'Csdeltavsk' , Csdeltavsk
  print, 'Csdeltavsk1', Csdeltavsk1
  print, 'Cstimesk'   , Cstimesk
  print, 'Cstimesk1'  , Cstimesk1

  print, 'slopes', result


  print, times
  print, timesk2
  save, rms, modes,rms2,modes2, filename = 'linearfigure10.sav'

  wi,2
  !p.multi=[0,0,0]
  cgplot,slopes,rms[0,*]
  for i=1,n_elements(rms[*,0])-1 do begin
    cgoplot,slopes,rms[i,*]
  endfor
 rms = modes
  slopes = [0.001d,.0025d, 0.01d, 0.02d,.03d,.04d ,0.05d,.06d ,0.07d,.08d,.09d, .1d,.11d,.12d,.13d,.14d,.15d,.16d,.17d, .18d,.19d, .2d,.21d, .22d,.23d,.24d, .25d, .26d,.27d, .28d,.29d, .3d]
  !p.multi = [0,0,0]
  cgplot,slopes,rms[0,*],psym=16,symsize=2,yrange=[0,.02]
  for i=1,n_elements(rms[*,0])-1 do begin
    cgoplot,slopes,rms[i,*],psym=16,symsize=2,yrange=[0,.1]
  endfor

  ensembleavg = make_Array(n_elements(rms[0,*]))
  sigma = make_Array(2,n_elements(rms[0,*]))
  errors = make_Array(n_elements(rms[0,*]))

  for i=0,n_elements(rms[0,*])-1 do begin
    ensembleavg[i] = mean(rms[*,i])
    sigma[*,i] = percentiles(rms[*,i], conflimit=0.6)
    errors[i] = stddev(rms[*,i])
  endfor


  wi,2,wsize = [1720,1100]

  ;cgplot,slopes,ensembleavg,psym=16,symsize=2,yrange=[-.001,.022],xrange=[-.01, .31],err_yhigh =  ,err_ylow = errors

  cgplot,slopes,ensembleavg,psym=16,symsize=2,yrange=[-.001,.022],xrange=[-.01, .31],err_yhigh = abs(sigma[1,*]-ensembleavg) ,err_ylow = abs(sigma[0,*]-ensembleavg),xtitle='!4Slope',ytitle='!4Maximum speed of the wake relative to the ring',font=1,charsize=4.2
  cgplot,slopes,ensembleavg,psym=16,symsize=2,yrange=[-.001,.022],xrange=[-.01, .31],err_yhigh = errors ,err_ylow =errors,xtitle='!4Slope',ytitle='!4Maximum speed of the wake relative to the ring',font=1,charsize=4.2

  fig = errorplot(slopes,ensembleavg*100,errors*100,yrange=[-.1,2.2],xrange=[-.01, .31],xtitle='!4Slope [-]',ytitle='!4Maximum collisional speed [cm/s]',symbol='o',linestyle=' ',ERRORBAR_COLOR="k",sym_color = 'k',sym_filled = 1)
  fig.thick = 3
  fig.sym_thick = 2
  fig.font_size = 32


  ;leg  = legend(position=[.8,.8], font_size = 16)


  b = Regress(slopes, ensembleavg, CONST=constb);,measure_errors = errors)

  a = Regress(alog10(slopes), alog10(ensembleavg), CONST=const);,measure_errors = alog10(errors))

  a1 = Regress(sqrt(slopes), ensembleavg, CONST=const1,measure_errors = errors);abs(sigma[1,*]-ensembleavg))
  b1 = Regress((slopes)^(1./4), ensembleavg,measure_errors = errors);abs(sigma[1,*]-ensembleavg)) the (slopes)^(1./4) fit at first may seem better than the sqrt fit, but it assumes particles can get lifted at close to cero relative speeds. The sqrt fit requires a threshold which is more consistent with microgravity experiments. See Sega et al 2024

  slopes = [0,slopes]
  errors = [.001,errors]
  ensembleavg = [0,ensembleavg]
  ;c = Regress(slopes, exp(ensembleavg), CONST=constc,measure_errors = exp(errors))

  ;slopes[0]=0.001

  cgoplot, slopes,sqrt(slopes)*(a1[0])+const1,color=cgcolor('red'),thick=2

  ;fig2 = plot(slopes,(sqrt(slopes)*(a1[0])+const1)*100,color='r',thick=2, name = 's$\sqrt{\middh\slashdr\mid  }+ b$',overplot=1)
    fig2 = plot(slopes,(sqrt(slopes)*(a1[0])+const1)*100,color='r',thick=2, name = '$2.42\sqrt{slope} + 0.47$',overplot=1)
 ;  fig2 = plot(slopes,((slopes)^(1./4)*(b1[0]))*100,color='r',thick=2, name = '$2.42(slope)^{0.25}$',overplot=1)

  

  leg  = legend(position=[.85,.5], font_size = 25, target = [fig2],linestyle=6)

  le_string = string(83B)


  cgLegend, Color=['Red'], Location=[0.2, 0.88],  Titles=['!9r'+'!X'+'!3slope + b!X'], Length=0.075, /Box, VSpace=3,/Background, bg_color = 'white', charsize=2.5

  print,const1[0]
  print,a1[0]

end