; Function that takes the dimensions and density of the wake (l, w, h, rho), and the shear rate (s), the angular velocity (omegaz) and a vector with the ortientation of the z' Eulerian axis (z2),
; and the angle phi (that can be a number or a 1D array) and returns the torque on the wake.

;If at the upper lever omegax is and z1 are introduce instead, and x=1, then we get the torque wrt to the x axis.

;Chage Log
;11/01/21 Created - Daniel Sega

;Example call
;l=300, w=20, h=4, s=1.5, rho=.5, phi = dindgen(15, Increment=.1)
;tor = torquecol(l,w,h,s,rho,phi)
;
;Example use in saturn_wake program:
;torquecol(l,w,h,s,rho,phi, omegaz = omegaz[i-1], z2 = z[2])
;torquecol(l,w,h,s,rho,phi, omegaz = omegax[i-1], z2 = z[1])

function torquecol3, l, w, h, s, rho_mean, x0, x1, x2, y0, y1, y2, z0, z1, z2, omegarot=omegarot, Xi=xi, ind=ind, phis=phis, oldphi = oldphi

  rho = 500.
  rv  = 131902000.
  G   = 6.67e-11
  M   = rho*h*w*l
  omega = 1.28897e-4
  GM = G*M
  ;copy = y1

  ;  if (-y0 lt 0) and (y1 lt 0) then begin
  ;phi = 2*!dpi + atan(-y0,y1)
  phi = atan(-y0,y1)


  ;  if oldphi - phi gt 2 then phi = -2*!dpi + atan(-y0,y1)
  ;  if oldphi - phi lt 2 then phi = 2*!dpi + atan(-y0,y1)
  ;   endif else begin
  ;phi = atan(-y0,y1)

  if abs(y0) lt 1e-3 then phi = atan(-y0/y1)
  ;if abs(phis[ind-1] - phi) gt 1 then phi = -abs(atan(-y0,y1)) -!dpi/2

  ;  if (oldphi - phi gt 3) and (ind gt 10) then begin
  ;   phis[ind] = 2*!dpi + atan(-y0,y1)
  ;  endif
  ;  if oldphi - phi lt -3 and (ind gt 10) then begin
  ;    phis[ind] = -2*!dpi + atan(-y0,y1)
  ;  endif else begin

  ; if abs(phi*180/!dpi) lt .1 then begin
  ;   print, 'why'
  ; endif


  ; if (abs(phis[ind-1] - phi) gt 1) and abs(phi) lt .1 then begin
  ; phis[ind] = -!dpi
  ;endif else begin
  phis[ind] = phi
  ;endelse
  ;if   (phi gt 2) and (abs(phis[ind-1] - phi) gt 3)  then phis[ind] = phi - 2*!dpi

  ;phis[ind] = phi ;this is a one call line, it just works for certain parameters in the system
  ;  endelse
  ;
  ; print, phis[ind]*180/!dpi


  phi = abs(phi)
  ;if copy gt 1 then phi = asin( 2- copy )
  ;if copy lt-1 then phi = asin( -2+ copy )

  tetax = atan(y2,z2)
  ; tetax = abs(tetax)

  if ~keyword_Set(z1) then z1=0

  ;print,phi
  ;print,tetax

  sign = -1
  if z2 gt 0 then sign = 1

  if keyword_set(xi) and (z1 gt 0) then sign = 1
  ;if keyword_set(yi) and (z1 lt 0) then sign = 1

  ;set of impact parameters

  ;b = dindgen(2000, increment = .1)
  b = linspace(10,l,300)

  ;phi = abs(phi)
  ;define angular momemtum

  j = s*omega*b^2

  ;define energy

  E = 1/2.*(s*omega*b)^2

  ;define eccentricity

  ec = Sqrt(1 + 2*E*j^2/(GM)^2)

  ; define argument of periapsis

  po = acos(b*(GM/j^2)*(ec^2 - 1)/ec)

  ;define perpendicular (to the wake) velocity.

  vtheta = GM*(1 + ec*Cos(!dpi/2. - phi - po))/j

  ;find the integration limits

  ;r = j^2/GM * 1/(1 + ec*Cos(!dpi/2. -phi - po))

  ;cos1 = y0*cos(po) + y1*sin(po)

  ;print,cos[30]
  ;print,Cos(!dpi/2. - phi - po[30])

  if n_elements(phi) eq 1 then begin

    i = where(j^2/GM * 1/(1 + ec*Cos(!dpi/2. -phi - po)) - l/2.*cos(tetax) lt .1)
    k = where(j^2/GM * 1/(1 + ec*Cos(!dpi/2.  -phi - po)) - W/2*cos(phi) lt .1)

    if k[-1] eq -1 then k[-1] = 0
    if k[-1] eq i[-1] then k[-1] = 10

    if k[-1] ge i[-1] then return, 0


    ;vtheta = GM*(1 + ec*Cos1)/j

    ;torque = int_tabulated(vtheta[0:i[0]]^2 * rho * H * j[0:i[0]]^2/GM * 1/(1 + ec[0:i[0]]*cos(phi - s[0:i[0]])), b[0:i[0]])
    if n_elements(omegarot) gt 1 then omegarot = omegarot[-1]
    
    return, j[i[-1]]^2/GM *1/(1 + ec[i[-1]]*cos(- po[i[-1]]))


    if keyword_set(xi) then begin

      torque = 4*tsumcumul((vtheta[k[-1]:i[-1]] + sign * omegarot * j[k[-1]:i[-1]]^2/GM *1/(1 + ec[k[-1]:i[-1]]*cos(!dpi/2. - phi - po[k[-1]:i[-1]])))^2 * rho_mean * H * j[k[-1]:i[-1]]^2/GM * 1/(1 + ec[k[-1]:i[-1]]*cos(!dpi/2. - phi - po[k[-1]:i[-1]])), b[k[-1]:i[-1]])*(cos(phi)*z0 + sin(phi)*z1)
      torque= abs(torque[-1])

    endif else begin

      torque = 4*tsumcumul((vtheta[k[-1]:i[-1]] + sign * omegarot * j[k[-1]:i[-1]]^2/GM *1/(1 + ec[k[-1]:i[-1]]*cos(!dpi/2. - phi - po[k[-1]:i[-1]])))^2 * rho_mean * W * j[k[-1]:i[-1]]^2/GM * 1/(1 + ec[k[-1]:i[-1]]*cos(!dpi/2. - phi - po[k[-1]:i[-1]])), b[k[-1]:i[-1]])*(-cos(phi)*x0 - sin(phi)*x1)
      torque= abs(torque[-1])

    endelse



    ; torque = 4*tsumcumul((vtheta[k[-1]:i[-1]])^2 * rho_mean * H * j[k[-1]:i[-1]]^2/GM * 1/(1 + ec[k[-1]:i[-1]]*cos(!dpi/2. - phi - po[k[-1]:i[-1]])), b[k[-1]:i[-1]])
    ;torque= torque[-1]


    ;print,torque1

    ;torque = 4*tsumcumul((vtheta[k[-1]:i[-1]] + sign * omegarot * j[k[-1]:i[-1]]^2/GM *1/(1 + ec[k[-1]:i[-1]]*cos1[k[-1]:i[-1]]))^2 * rho_mean * H * j[k[-1]:i[-1]]^2/GM * 1/(1 + ec[k[-1]:i[-1]]*cos1[k[-1]:i[-1]]), b[k[-1]:i[-1]])
    ;torque= torque[-1]
    ;print,torque

  endif else begin

    torque       = make_array(n_elements(phi))
    ;torquesemiap = make_array(n_elements(phi))

    for ind=0, n_elements(phi) -1 do begin

      vtheta = GM*(1 + ec*Cos(!dpi/2. - phi[ind] - po))/j
      i = where(j^2/GM * 1/(1 + ec*Cos(!dpi/2. -phi[ind] - po)) - l/2.*cos(tetax) lt .1)
      k = where(j^2/GM * 1/(1 + ec*Cos(!dpi/2. -phi[ind] - po)) - W/2.*cos(phi) lt .1)

      torque[ind] = 2*int_tabulated(vtheta[k[-1]:i[-1]]^2 * rho_mean * H * j[k[-1]:i[-1]]^2/GM * 1/(1 + ec[k[-1]:i[-1]]*cos(!dpi/2. - phi[ind] - po[k[-1]:i[-1]])), b[k[-1]:i[-1]])

      ;torquesemiap[ind] = 2 * int_tabulated(rho * H * GM * (1 + ec[k[-1]:i[-1]]*Cos(!dpi/2. -phi[ind] - po[k[-1]:i[-1]])), b[k[-1]:i[-1]])
    endfor

  endelse

  torqueapprox = 4* rho_mean * H * GM * ((Sin(phi) + Cos(phi))*(b[i[-1]] - 30) + s^2*omega^2*(b[i[-1]]^4 - 30^4)*Sin(phi)/(4*GM)); + GM/(10*s^4*Omega^4)*(1/30^5 - 1/b[i[-1]]^5)*Sin(phi) + (Sin(phi) + Cos(phi))* GM/(4*s^4*Omega^4)*(1/30^2 - 1/b[i[-1]]^2))

  ;print,torque
  ;print,torquesemiap

  ;print,torque/torqueapprox

  return, torque
end