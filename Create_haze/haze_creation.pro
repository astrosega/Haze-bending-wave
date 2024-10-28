;Program that uses a SCL model being kicked at different points in its
;trayectory to create a Haze
;it takes an angle phi, which is the phase at which the SHO is
;initially kicked
pro haze_creation, phi, theta, omega, i, sigma, nhits, epsilon, x, upper, slope, damp, k, A = A, koeffs = koeffs,s=s, plotf = plotf, fase = fase, landr = landr,c=c

  if ~keyword_set(A) then begin
     Hsigma=(sigma/1000.)
     A=(476.92d/(Sqrt(Hsigma)))/1000d ;amplitud of the wave
  endif
  
  if ~keyword_set(fase) then fase = 0
  
  if keyword_set(landr) then nhits = 2

                                ;First I draw the wave

  G  = 6.67408e-11              ;SI
  cD = 1.50217913e-7            ;this is the curly D in SCL, (1/s^2)
  mu = 1.2957087e-4             ;(1/s) vertical frequency of particles
  muM= .771634e-4               ;(1/s) vertical frequency of mimas
  aM = 185539.                  ;Km  mimas semi-major axis
  iM = 1.574 * !Dpi/180         ;rads, mimas inclination
  rv    = 131902.0


  J4 = -935.83e-6               ;

  J2 = 16290.71e-6              ;
  M = 5.6834e26                 ;
  r = 185539000                 ;
  R = 60330000                  ;

  omegaM = 0.0000771634
  muM    = 0.000077365
                                ;orbital paramenters from Mimas taken from Nasa webpage

  C = make_array(2,nhits)
  
  mu1 = mu - epsilon*k[i] * sqrt(slope[i]) * s*sqrt(2.)/2. + 0.00000475

  B = [[cos(phi), sin(phi)],[-mu1*sin(phi), mu1*cos(phi)]]
  R = [A*damp[i]*cos(phi), -omega[0]*A*damp[i]*sin(phi) + epsilon*sqrt(slope[i])*s*sqrt(2.)/2.+ 0.00000475]

  C[* , 0] = Invert(B) ## R     ; I verified this is the right operator for this order of opertations
  
  if keyword_set(landr) and nhits gt 1 then begin
    mu2 = mu + epsilon*k[i] * sqrt(slope[i]) * s*sqrt(2.)/2.+ 0.00000475
    B = [[cos(phi), sin(phi)],[-mu2*sin(phi), mu2*cos(phi)]]
    R = [A*damp[i]*cos(phi), -omega[0]*A*damp[i]*sin(phi) + epsilon*sqrt(slope[i])*s*sqrt(2.)/2.+ 0.00000475]

    C[* , 1] = Invert(B) ## R     ; I verified this is the right operator for this order of opertations
    nhits = 1
  endif

for j = 1, nhits-1 do begin ;This considers the scenario where a ring particles get hit multiple times. Does not change the results in Sega et al. 2024.

;print,  k[i]*!dpi/mu * epsilon * slope[i]*s
;print, C[0, j-1]

phi = phi + !dpi ;+ k[i]*2*!dpi/mu *

alt = 0

B = [[cos(phi), sin(phi)],[-mu*sin(phi), mu*cos(phi)]]
 R = [A*damp[i]*cos(phi),-omega*A*damp[i]*sin(phi)*(1 +epsilon)*(2./3.) + (1./3.)*(1 -2*epsilon) * mu*(-C[0,j-1]*sin(phi) +C[1,j-1]*cos(phi)) ]

C[*,j] = Invert(B) ## R

;phase = atan(-C[1],C[0])

;print, A, 'a'
;print,c[0,j],c[1,j]
;print,r
endfor
q = c
 
 if keyword_set(landr) then nhits = 2

  if keyword_set(plotf) then begin

     cgoplot,theta/(2*!dpi),A*damp[i]*cos(theta),color=cgcolor('red'),yrange = [-1.5,1.5], background = cgcolor('white')
     cgoplot, theta/(2*!dpi), c[0,0]*cos(theta) + c[1,0]*sin(theta),color=cgcolor('blue')
                             
  endif
  if abs(i) lt 1e-3 then i = 1                        
  phase  = tsum(x[0:i],k[0:i])
  koeffs = [c[0,0]*cos(phase- !dpi/4 + fase) + c[1,0]*sin(phase- !dpi/4 + fase)]
  
  if keyword_set(landr) then  koeffs = [c[0,0]*cos(phase- !dpi/4 + fase) + c[1,0]*sin(phase- !dpi/4 + fase),c[0,1]*cos(phase - !dpi/4 + fase) + c[1,1]*sin(phase- !dpi/4 + fase)]
  
  
                            
  if nhits gt 1  then begin
   if keyword_set(plotf) then cgoplot, theta/(2*!dpi),c[0,1]*cos(theta) + c[1,1]*sin(theta),color=cgcolor('blue')
      koeffs = [c[0,0]*cos(phase - !dpi/4 + fase) + c[1,0]*sin(phase - !dpi/4 + fase),c[0,1]*cos(phase - !dpi/4 + fase) + c[1,1]*sin(phase- !dpi/4 + fase)]
   endif
  if nhits gt 2 then begin
    if keyword_set(plotf) then cgoplot, theta/(2*!dpi), c[0,2]*cos(theta) + c[1,2]*sin(theta),color=cgcolor('cyan')
  koeffs = [c[0,0]*cos(phase- !dpi/4 + fase) + c[1,0]*sin(phase- !dpi/4 + fase),c[0,1]*cos(phase- !dpi/4 + fase) + c[1,1]*sin(phase- !dpi/4 + fase),c[0,2]*cos(phase- !dpi/4 + fase) + c[1,2]*sin(phase- !dpi/4 + fase)]
  endif
  if nhits gt 3 then cgoplot, theta/(2*!dpi), c[0,3]*cos(theta) + c[1,3]*sin(theta),color=cgcolor('green')
  if nhits gt 4 then cgoplot, theta/(2*!dpi), c[0,4]*cos(theta) + c[1,4]*sin(theta), color = cgcolor('red')
  if nhits gt 5 then cgoplot, theta/(2*!dpi), c[0,5]*cos(theta) + c[1,5]*sin(theta)
  if nhits gt 6 then cgoplot, theta/(2*!dpi), c[0,6]*cos(theta) + c[1,6]*sin(theta)
  if nhits gt 7 then cgoplot, theta/(2*!dpi), c[0,7]*cos(theta) + c[1,7]*sin(theta)
  if nhits gt 8 then cgoplot, theta/(2*!dpi), c[0,8]*cos(theta) + c[1,8]*sin(theta)
  if nhits gt 9 then cgoplot, theta/(2*!dpi), c[0,9]*cos(theta) + c[1,9]*sin(theta)

end



