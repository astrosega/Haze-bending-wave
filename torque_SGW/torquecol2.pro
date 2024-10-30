;torqyecol but for the torque abou the y axis.
;The trajectories for the collisions are obtained using a gravitational focusing approach outlines in section 2 of Sega et al 2024, ICARUS

;If at the upper level omegax and z1 are introduced instead, and the keyworkd x=1, then we get the torque wrt to the x axis.

;Chage Log
;11/01/21 Created - Daniel Sega

;Example call
;l=300, w=20, h=4, s=1.5, rho=.5, phi = dindgen(15, Increment=.1)
;tor = torquecol(l,w,h,s,rho,phi)
;
;Example use in saturn_wake program:
;torquecol2(l,w,h,s,rho,rho_mean, x[0], x[1], x[2], y[0], y[1], y[2], z[0], z[1], z[2], omegarot = omegay[i-1],        ind=i, phis=phis, dhdot = dhdot)

function torquecol2, b, l, w, h, s, rho, rho_mean, x0, x1, x2, y0, y1, y2, z0, z1, z2, omegarot=omegarot, ind=ind, phis=phis, dhdot = dhdot, thetays = thetays, intL=intl

std = 10d
  rv  = 131902000.
  G   = 6.67e-11
  M   = rho*h*w*l
  omega = 1.28897e-4
  GM = G*M
  vkdr = -s*Omega


  phi   = atan(x1,-x0)
  tetay = !dpi/2 - acos(x2)
  ;tetay = abs(tetay)
  
  if keyword_set(thetays) then thetays[ind] = tetay; acos(y1)
 ; tetaz = atan(-y0,y1)
 
 if ~keyword_set(intl) then intl = linspace(-l/2,l/2,60)  
 
    if phi lt 0 then phi = !dpi + phi


 ;b0 = linspace(.01d,0.99d,250)
 ;b1 = linspace(1d,4.99d,250)
 ;b2 = linspace(5d, W, 500)

 ;b=[b0,b1,b2]
 
 ;b = linspace(w/2*sin(phi)*abs(cos(tetay)) + 10, l/2 + 50,2000)
; b = linspace(w/2*sin(phi) + 10, w/2*sin(phi) + 2.5*10,2000)
  
 b = linspace(39d, 70 + 50*sin(phi), 10000)


 ;b = linspace(5d, W, 5000)

  ;define angular momemtum

  j = s*omega*b^2

  ;define energy

  E = 1/2.*(s*omega*b)^2

  ;define eccentricity

  ec = Sqrt(1 + 2*E*j^2/(GM)^2)

  ; define argument of periapsis

  po = acos(b*(GM/j^2)*(ec^2 - 1)/ec)

  ;define perpendicular (to the face of the wake's projection in the plane being hit) velocity.

  vtheta = GM*(1 + ec*Cos(!dpi/2 - phi - po))/j

  ;find the integration limitsoldphi

  r = j^2/GM * 1/(1 + ec*Cos(!dpi/2. -phi - po))


    i = where(abs(j^2/GM * 1/(1 + ec*Cos(!dpi/2 - phi - po)) - W/2.*abs(cos(tetay))) lt 2)
   ; k = where(j^2/GM * 1/(1 + ec*Cos(!dpi/2 - phi - po)) - H/2*abs(cos(phi)) lt .01)
    k = where(abs(j^2/GM * 1/(1 + ec*Cos(!dpi/2 - phi - po)) - 2) lt 2)
    
  k[-1] = 0 ;this commented out for the plots in the draft
  
    
 
if k[-1] eq -1 then k[-1] = 0


;index = where( i - shift(i,1) gt 2)
;i = i[0:index[0]-1]
  
 ;if n_elements(k) gt 1 then begin
;index = where( k - shift(k,1) gt 2)
;k = k[0:index[0]-1]
;endif

fudge1 = tsum(intl, Exp(-(z2*(H/2)/(std/2) + intl*y2/(std/2))^2))

;fudge1 = tsumcumul(intl, Exp(-(z2*(H/2)/(std/2) + intl*y2/(std/2))^2))
;fudge1 = fudge1[-1]
;print,fudge1

    if k[-1] ge i[-1] then begin
               if ~finite(1/tan(tetay)^4) then tetay =0.001d
      ;print,  -4 * rho_mean * Abs( omegarot + dhdot*x0*z2) *  (omegarot + dhdot*x0*z2) * fudge1 * std^2 * (std^2 - Exp(-(w*tan(tetay)/std)^2) * (std^2 + w^2 * tan(tetay)^2))/(32*tan(tetay)^4)
      return, -4 * rho_mean * Abs( omegarot + dhdot*x0*z2) *  (omegarot + dhdot*x0*z2) * fudge1 * std^2 * (std^2 - Exp(-(w*tan(tetay)/std)^2) * (std^2 + w^2 * tan(tetay)^2))/(32*tan(tetay)^4)
    endif

    fudge = Exp(-(r[k[-1]:i[-1]]*tan(tetay)/(std/2))^2)

    


;print,r[i[-1]]/abs(cos(tetay))*2 - w

    
 ; print,  b[i[-1]]/(W*.5*x0*z1) 
;return, b[i[-1]]/(W*.5**z1*x0)
;print,  (j[i[-1]]^2/GM *1/(1 + ec[i[-1]]*cos( - po[i[-1]])))

     ;torque =  -4*int_tabulated(b[k[-1]:i[-1]], abs(vtheta[k[-1]:i[-1]]*(sin(phi)*z0 + cos(phi)*z1) + (dhdot*x0*z2 + omegarot) * r[k[-1]:i[-1]]/abs(cos(tetay))) * (vtheta[k[-1]:i[-1]]*(sin(phi)*z0 + cos(phi)*z1)     + (dhdot*x0*z2 + omegarot) * r[k[-1]:i[-1]]/abs(cos(tetay))) * rho_mean * fudge * L * r[k[-1]:i[-1]]/abs(cos(tetay)))
     
       torque =  -4*int_tabulated(b[k[-1]:i[-1]], abs(vtheta[k[-1]:i[-1]]*signum(y2)*(abs(sin(phi)*z0) + abs(cos(phi)*z1)) + (dhdot*x0*z2 + omegarot) * r[k[-1]:i[-1]]/abs(cos(tetay))) * (vtheta[k[-1]:i[-1]]*signum(y2)*(abs(sin(phi)*z0) + abs(cos(phi)*z1)) + (dhdot*x0*z2 + omegarot) * r[k[-1]:i[-1]]/abs(cos(tetay))) * rho_mean * fudge * fudge1 * r[k[-1]:i[-1]]/abs(cos(tetay)))
              ; torque =  -4*tsum(b[k[-1]:i[-1]], abs(vtheta[k[-1]:i[-1]]*(sin(phi)*z0 + cos(phi)*z1) + (dhdot*x0*z2 + omegarot) * r[k[-1]:i[-1]]/abs(cos(tetay))) * (vtheta[k[-1]:i[-1]]*(sin(phi)*z0 + cos(phi)*z1) + (dhdot*x0*z2 + omegarot) * r[k[-1]:i[-1]]/abs(cos(tetay))) * rho_mean * fudge * fudge1 * r[k[-1]:i[-1]]/abs(cos(tetay)))

;print, torque

  ;return,  -rho_mean * Abs( omegarot + dhdot*x0*z2) *  (omegarot + dhdot*x0*z2) * (W/2.)^4 * fudge1
  return,  torque
end
