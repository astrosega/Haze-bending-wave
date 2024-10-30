; Function that takes the dimensions and density of the wake (l, w, h, rho), and the shear rate (s), the angular velocity (omegaz) and a vector with the ortientation of the z' Eulerian axis (z2),
; and the angle phi (that can be a number or a 1D array) and returns the torque on the wake.
;The trajectories for the collisions are obtained using a gravitational focusing approach outlines in section 2 of Sega et al 2024, ICARUS

;If at the upper lever omegax is and z1 are introduce instead, and x=1, then we get the torque wrt to the x axis.

;Chage Log
;11/01/21 Created - Daniel Sega
;01/24/22 Added the rho paramter to make the Hill density a knob of the simulation

;Example call
;l=300, w=20, h=4, s=1.5, rho=.5, phi = dindgen(15, Increment=.1)
;tor = torquecol(l,w,h,s,rho,phi)
;
;Example use in saturn_wake program:
;torquecol1(l,w,h,s,rho,rho_mean, x[0], x[1], x[2], y[0], y[1], y[2], z[0], z[1], z[2], omegarot = omegax[i-1],  xi=1, ind=i, phis=phis, dhdot = dhdot)

function torquecol1, b, l, w, h, s, rho, rho_mean, x0, x1, x2, y0, y1, y2, z0, z1, z2, omegarot=omegarot, Xi=xi, ind=ind, phis=phis, oldphi = oldphi, dhdot = dhdot, thetaxs = thetaxs, intw=intW, IntH=inth, tsuma = tsuma

  if ~keyword_Set(xi) then xi=0
  

;tic

  std = 10d
  rv  = 131902000.
  G   = 6.67e-11
  M   = rho*h*w*l          ;for l=240, delta a 58, 54 for l=200
  ;M   = 900*7.5^3*!dpi*4/3 ;delta A is 29,65, or 20 if radius=5m
  ;M   = 900*10^3*!dpi*4/3 
  omega = 1.28897e-4
  GM = G*M
  vkdr = -s*Omega


;if ~keyword_set(intw) and xi  then intw = linspace(-w/2,w/2,60)
;if ~keyword_set(inth) and ~xi then inth = linspace(-h/2,h/2,40)  
  ;copy = y1

  ;  if (-y0 lt 0) and (y1 lt 0) then begin
  ;phi = 2*!dpi + atan(-y0,y1)
  phi  =  atan(-y0,y1)


  ;  if oldphi - phi gt 2 then phi = -2*!dpi + atan(-y0,y1)
  ;  if oldphi - phi lt 2 then phi = 2*!dpi + atan(-y0,y1)
  ;   endif else begin
  ;phi = atan(-y0,y1)

  ; if abs(y0) lt 1e-3 then phi = atan(-y0,y1)
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
;  if phi gt 3 then begin
;    print, 'it has happened'
;    print, phi
;  endif
;if phi lt 0 then begin
;  phi = atan(-y0/y1)
;  if phi gt .5 then phi = atan(-y0,y1)
;endif


;print,phi

  ; if (abs(phis[ind-1] - phi) gt 1) and abs(phi) lt .1 then begin
  ; phis[ind] = -!dpi
  ;endif else begin
  phis[ind] = phi
 
 ; b = linspace(38d, 70 + 50*sin(phi), 10000)
  ;endelse
  ;if   (phi gt 2) and (abs(phis[ind-1] - phi) gt 3)  then phis[ind] = phi - 2*!dpi

  ;phis[ind] = phi ;this is a one call line, it just works for certain parameters in the system
  ;  endelse
  ;
  ; print, phis[ind]*180/!dpi


 ; phi = abs(phi)
  ;if (phi lt 0) and (phi gt -1) then phi=!dpi + phi
  ;if copy gt 1 then phi = asin( 2- copy )
  ;if copy lt-1 then phi = asin( -2+ copy )
  tetax = asin(y2)
    if keyword_set(thetaxs) then thetaxs[ind] = tetax ; acos(x0)
  tetax = abs(tetax)
  
  if ~finite(y0) then begin
    print, 'it has happened'
    print,'y0', y0
    print,'y1', y1
    ;print,'y2', y2
    print,'z0', z0
    print,'z1', z1
    print,'z2', z2
    print,'x0', x0
    print,'x1', x1
    print,'x2', x2
    print, y0^2 + y1^2 + y2^2
    print, z0^2 + z1^2 + z2^2
    print, x0^2 + x1^2 + x2^2
  endif
  
  ;print,'x0', x0
  ;print,'x1', x1
  ;print,'x2', x2
  
  
  
;if z1 lt 0 then begin
;  print, 'check'
;endif



  ;print,phi
  ;print,tetax
  ;if keyword_set(yi) and (z1 lt 0) then sign = 1

  ;set of impact parameters
  
  start = .001d
  ;start = h/2d
  ;start = w
  
  ;if (abs(tan(phi)) lt w/l) then start = w/2 * abs(cos(phi))
   ;if (abs(tan(phi)) lt h/l) then start = h/2 * abs(cos(phi))
    if phi lt 0 then phi = !dpi + phi
 ;  b = linspace(l/2*sin(phi), l, 5000)
 ; b0 = linspace(start,2.999d,10000)
 ; b1 = linspace(3d,19.999d,5000)
 ;   if sin(phi)*l/2 lt 10*2.5  then b = linspace(10d*1.7d,10*2.5 , 5000) else   b = linspace(10d*1.7d, l/2*sin(phi)*l/2, 5000)
     
  b = linspace(39d, 70 + 50*sin(phi), 10000)

  
 ; b=[b0,b1,b2]
 ; b = linspace(55d*.5, l/2 +10, 2000)
   ; b = linspace(20, l/2 +10, 5000)

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

  r = j^2/GM * 1/(1 + ec*Cos(!dpi/2. - phi - po))
  
  ;vtheta1 = j/r

  ;cos1 = y0*cos(po) + y1*sin(po)
;if ind gt 48000.*2/5. and x1 then begin
;  print,'skopy'
;endif
  ;print,cos[30]
  ;print,Cos(!dpi/2. - phi - po[30])
  

                        i = where(abs(j^2/GM * 1/(1 + ec*Cos(!dpi/2. - phi - po)) - l/2.*abs(cos(tetax))) lt 3) ;sqrt((l/2.*abs(cos(tetax)))^2 + w^2/4)
                        k = where(abs(j^2/GM * 1/(1 + ec*Cos(!dpi/2. - phi - po)) - h/2 )   lt 2) ;*abs(cos(phi)) or w/2*sqrt(2) and lt 3
if xi then k = where(abs(j^2/GM * 1/(1 + ec*Cos(!dpi/2. -phi - po)) - 2)  lt 1) ;*abs(cos(phi)

k[-1] = 0

;if k[0] eq -1 then begin
;print,'k',k
;endif

;if i[0] eq -1 then begin
 ; print,'i',i
;endif


if k[-1] eq -1 then k[-1] = 0

index = where( i - shift(i,1) gt 10, n)
if n gt 1 then begin
        index = where( i - shift(i,1) gt 10)
        i = i[0:index[0]-1]
endif

 index = where( k - shift(k,1) gt 10, n)
if n gt 1 then begin
      index = where( k - shift(k,1) gt 10)
      k = k[0:index[0]-1]
endif


; k[-1] = 0
    ;if k[-1] eq -1 then k[-1] = 0
    if (k[-1] ge i[-1]) and xi  then begin
     ; print, 'oh noo',xi
      fudge1 = Int_tabulated(intw,Exp(-(z2*(H/2)/(std/2)  + intw*x2/(std/2))^2))
      
   ;   if keyword_set(tsuma) then fudge1 = tsum(intw,Exp(-(z2*(H/2)/(std/2)  + intw*x2/(std/2))^2)) else fudge1 = Int_tabulated(linspace(-w/2,w/2,40),Exp(-(z2*(H/2)/(std/2)  + linspace(-w/2,w/2,40)*x2/(std/2))^2))


     
     if ~finite(1/tan(tetax)^4) then tetax =0.001d
      
      return,  4* rho_mean * Abs(-omegarot + dhdot*y0*z2) * (-omegarot + dhdot*y0*z2) * fudge1  * std^2 * (std^2 - Exp(-((l * tan(tetax))/std)^2) * (std^2 + l^2 * tan(tetax)^2))/(32*tan(tetax)^4) ;(l/2)^2
      
    endif
    if (k[-1] ge i[-1]) and ~xi then begin
       fudge1 = Int_tabulated(inth,Exp(-(x2*W/H + inth*z2/(std/2))^2))
      ;fudge1 = tsum(inth,Exp(-(x2*(W/2)/(std/2) + inth*z2/(std/2))^2))
    ; if keyword_set(tsuma) then fudge1 = tsum(inth,Exp(-(x2*(W/2)/(std/2) + inth*z2/(std/2))^2)) else fudge1 = Int_tabulated(linspace(-H/2,H/2,40),Exp(-(x2*(W/2)/(std/2) + linspace(-H/2,H/2,40)*z2/(std/2))^2))

 
      
          if ~finite(1/tan(tetax)^4) then tetax =0.001d
            return, -4* rho_mean * Abs( omegarot + dhdot*y0*x2) * ( omegarot + dhdot*y0*x2) *fudge1  * std^2 * (std^2 - Exp(-((l * tan(tetax))/std)^2) * (std^2 + l^2 * tan(tetax)^2))/(32*tan(tetax)^4)
      ;return, -4* rho_mean * Abs( omegarot + dhdot*y0*x2) * ( omegarot + dhdot*y0*x2) *fudge1  * int_tabulated(intl2,exp(-(intl2*tan(tetax)/std)^2) *intl2^3)
      
    endif
  
  fudge = Exp(-(r[k[-1]:i[-1]]*tan(tetax)/(std/2))^2)
  ;fudge = Exp(-(r[k[0]:i[-1]]*tan(tetax)/(std/2))^2)
    
    if xi then begin

     ;   fudge1 = Int_tabulated(intw,Exp(-(z2*(H/2)/(std/2)  + intw*x2/(std/2))^2))
       fudge1 = tsum(intw,Exp(-(z2*(H/2)/(std/2)  + intw*x2/(std/2))^2))
  ;     if keyword_set(tsuma) then fudge1 = tsum(intw,Exp(-(z2*(H/2)/(std/2)  + intw*x2/(std/2))^2)) else fudge1 = Int_tabulated(linspace(-w/2,w/2,40),Exp(-(z2*(H/2)/(std/2)  + linspace(-w/2,w/2,40)*x2/(std/2))^2))

       ;fudge1 = Int_tabulated(linspace(-w/2,w/2,40),Exp(-(z2*(H/2)/(std/2)  + linspace(-w/2,w/2,40)*x2/(std/2))^2))

 
        torque =  4*fudge1*int_tabulated(b[k[-1]:i[-1]], abs(vtheta[k[-1]:i[-1]]*(-signum(x2))*(abs(cos(phi)*z0) + abs(sin(phi)*z1)) + (dhdot*y0*z2 - omegarot) * r[k[-1]:i[-1]]/abs(cos(tetax))) * (vtheta[k[-1]:i[-1]]*(-signum(x2))*(abs(cos(phi)*z0) + abs(sin(phi)*z1))     + (dhdot*y0*z2 - omegarot) * r[k[-1]:i[-1]]/abs(cos(tetax))) * rho_mean * fudge * r[k[-1]:i[-1]]/abs(cos(tetax)))

      ;  torque =  4*fudge1*int_tabulated(b[k[-1]:i[-1]], abs(vtheta[k[-1]:i[-1]]*(cos(phi)*z0 + sin(phi)*z1) + (dhdot*y0*z2 - omegarot) * r[k[-1]:i[-1]]/abs(cos(tetax))) * (vtheta[k[-1]:i[-1]]*(cos(phi)*z0 + sin(phi)*z1)     + (dhdot*y0*z2 - omegarot) * r[k[-1]:i[-1]]/abs(cos(tetax))) * rho_mean * fudge * r[k[-1]:i[-1]]/abs(cos(tetax)))
      ; torque =  4*fudge1*int_tabulated(b[k[-1]:i[-1]], abs(vtheta[k[-1]:i[-1]]*(cos(phi)*z0 + sin(phi)*z1) + (dhdot*y0*z2 - omegarot) * r[k[-1]:i[-1]]/abs(cos(tetax))) * (vtheta[k[-1]:i[-1]]*(cos(phi)*z0 + sin(phi)*z1)     + (dhdot*y0*z2 - omegarot) * r[k[-1]:i[-1]]/abs(cos(tetax))) * rho_mean * fudge * r[k[-1]:i[-1]]/abs(cos(tetax)))
      
             ;torque =  4*int_tabulated(b[k[0]:i[-1]], abs(vtheta[k[0]:i[-1]]*(cos(phi)*z0 + sin(phi)*z1) + (dhdot*y0*z2 - omegarot) * r[k[0]:i[-1]]/abs(cos(tetax))) * (vtheta[k[0]:i[-1]]*(cos(phi)*z0 + sin(phi)*z1)     + (dhdot*y0*z2 - omegarot) * r[k[0]:i[-1]]/abs(cos(tetax))) * rho_mean * fudge * fudge1 * r[k[0]:i[-1]]/abs(cos(tetax)))


    
     endif else begin

      ;fudge1 = int_tabulated(inth,Exp(-(x2*(W/2)/(std/2) + inth*z2/(std/2))^2))
      fudge1 =          tsum(inth,Exp(-(x2*(W/2)/(std/2) + inth*z2/(std/2))^2))
      ; fudge1 = Int_tabulated(linspace(-H/2,H/2,40),Exp(-(x2*(W/2)/(std/2) + linspace(-H/2,H/2,40)*z2/(std/2))^2))
   ;  if keyword_set(tsuma) then fudge1 = tsum(inth,Exp(-(x2*(W/2)/(std/2) + inth*z2/(std/2))^2)) else fudge1 = Int_tabulated(linspace(-H/2,H/2,40),Exp(-(x2*(W/2)/(std/2) + linspace(-H/2,H/2,40)*z2/(std/2))^2))


     torque = -4*fudge1*int_tabulated(b[k[-1]:i[-1]], abs(vtheta[k[-1]:i[-1]]*signum(z2)*(abs(cos(phi)*x0) + abs(sin(phi)*x1)) + (dhdot*y0*x2 + omegarot) * r[k[-1]:i[-1]]/abs(cos(tetax))) * (vtheta[k[-1]:i[-1]]*signum(z2)*(abs(cos(phi)*x0) + abs(sin(phi)*x1))     + (dhdot*y0*x2 + omegarot) * r[k[-1]:i[-1]]/abs(cos(tetax))) * rho_mean * fudge * r[k[-1]:i[-1]]/abs(cos(tetax)))
      
     ; torque = -4*int_tabulated(b[k[0]:i[-1]], abs(vtheta[k[0]:i[-1]]*(cos(phi)*x0 + sin(phi)*x1) + (dhdot*y0*x2 + omegarot) * r[k[0]:i[-1]]/abs(cos(tetax))) * (vtheta[k[0]:i[-1]]*(cos(phi)*x0 + sin(phi)*x1)     + (dhdot*y0*x2 + omegarot) * r[k[0]:i[-1]]/abs(cos(tetax))) * rho_mean * fudge * fudge1 * r[k[0]:i[-1]]/abs(cos(tetax)))
    ;   torque = -4*fudge1*int_tabulated(b[k[-1]:i[-1]], abs(vtheta[k[-1]:i[-1]]*(cos(phi)*x0 + abs(sin(phi)*x1)) + (dhdot*y0*x2 + omegarot) * r[k[-1]:i[-1]]/abs(cos(tetax))) * (vtheta[k[-1]:i[-1]]*(cos(phi)*x0 + sin(phi)*x1)     + (dhdot*y0*x2 + omegarot) * r[k[-1]:i[-1]]/abs(cos(tetax))) * rho_mean * fudge * r[k[-1]:i[-1]]/abs(cos(tetax)))
               ;torque = -4*fudge1*tsum(b[k[-1]:i[-1]], abs(vtheta[k[-1]:i[-1]]*(cos(phi)*x0 + sin(phi)*x1) + (dhdot*y0*x2 + omegarot) * r[k[-1]:i[-1]]/abs(cos(tetax))) * (vtheta[k[-1]:i[-1]]*(cos(phi)*x0 + sin(phi)*x1)     + (dhdot*y0*x2 + omegarot) * r[k[-1]:i[-1]]/abs(cos(tetax))) * rho_mean * fudge * r[k[-1]:i[-1]]/abs(cos(tetax)))



     ; vth = vtheta[i[-1]]
     ; print, 'vth',vth
     ; print, 'vkdr', vkdr*l/2*y0
      
    endelse
    
;
;if ~xi and torque gt 0 then begin
 ; print,'torque', torque
 ;print, 'omegarot', omegarot
;   print,ind
; endif
 ; print,'torque', torque, xi
 ;print,r[i[-1]]/abs(cos(tetax))*2 - l, 'l'
 
; toc
  return, torque
end
