;It calculates the optical depth between the radial distances r1 and r2 given Intensity (Irradiance) and radius.
;To do so it rebins the original data by using occdatabin.pro

;Parameters
;--------------
;Radius: array

;data: array

;reso: radial resolution of output optical depth, data and radius array in Km. float or int

;taomax:  set to one in order to also get the maximum allowed value for tao given the selected resolution

;taoplus: set to one in order to also get the upper error bar for tao

;i: this is a dummy to be able to pass the index value if program is run on a loop that is going trough a list of occultations

;Output
;-----------------
;binned :(3 by (size of binned data array)) containing the binned radius (binned[*,0]), the optical depth (binned[*,1]) and the irradiance (binned[*,2]) of the star.

;if taomax = 1 then it will also contain the maximum allowed value for tao, given the selected resolution, in binned[*, 3]

;If taoplus = 1 then it will also contain the upper error bar for the optical depth in binned[*,3], and the difference between the optical depth and the upper errobar in
;binned[*, 4] 

;Updates
;07/12/2016 - Created by DDSN
;6/26/2017  - Added taoplus
;6/27/2017  - Added the possibility of passing a b value.

Function odepth_occdatabin, rad, data, reso, r1, r2, Taomax=taomax, Taoplus=taoplus, b=b, i=i, af=af, laplace = laplace


  ;First find the resolution of the desired region

  r = findradius(rad, r1,r2)
  points_in_bin = floor(abs((r[0] - r[1]) / ((r2-r1))) * reso) ;the region has to be such that the variation in the rad can be taken to be constant
  print,points_in_bin
  bins     = n_elements(data)/points_in_bin
  reradius = rebin(rad[0:bins*points_in_bin-1], bins)
  ;plot,reradius
  re = findradius(reradius, r1,r2)
  reradius = reradius[re[0]:re[1]]
  ;print,reradius

  ;we only want to rebin a particular region, so we create subvectors

  radreg  = rad[r[0]:r[1]]
  datreg  =   data[r[0]:r[1]]
  
  if r[0] eq r[1] then begin
    print,'no data at this redius'
    return,[[1,2,3],[1,2,3],[1,2,3]]
  endif

  occdatabin, reradius, radreg, datreg, dbin
  
   redata = dbin

  ;;;;;;; finding the background and maximum irradiance
  ;Section 3.3 of JOsh et all 2010 AstJnr 140
  
  

  e    = findradius(rad, 137000,137500)
if keyword_set(af) then e=findradius(rad, 133450,133700)
if keyword_set(laplace) then begin
  e=findradius(rad, 119860,119950)
endif

  Imax = mean(data[e[0]:e[1]])
  print,'imax',imax

  ;this is a section of the B ring that is opaque

  rb = findradius(rad,107870,107900) ;Range from Colwell 2010
  
  if (rb[0] EQ rb[1]) AND ~keyword_set(b) then begin
    print, 'no b ring data in this occultation. Using b default value for AlpVir occultations'
    b = 0.260163
    
    if keyword_set(i) then print, 'no b in', i ;this is so that I can see which occs have no b
  
  endif
  
  if ~keyword_set(b) then b  = mean(data[rb[0]:rb[1]])
  
print,'b',b
  ;;;;;;;;
  
  redata = dbin
  ;The stars aparent brightness will be given by the unoculted data minus the background.

  I0  = Imax-b
  
  if keyword_set(i) AND I0 LE 10 then print, 'Low I0 in ', i 
  
  Irr = redata-b

  ;plot,redata
  
  ;this is needed for very low angle optical depths, since the 'data' may equal 0: the ring gets completely opaque and tao diverges

  Irr[where(Irr LE b, /NULL)] = b

  If KEYWORD_SET(taomax) then taomax = -alog(sqrt(Irr)/I0)
  

  ;r=mimas53(reradius)
  ;plot,reradius[r[0]:r[1]],redata[r[0]:r[1]]

  tao     = -alog(Irr/I0)
  taoplus =  alog(I0) - alog(I0*Exp(-tao) - sqrt(I0*Exp(-tao) + b)) ;Colwell et al. 1990. upper limit of error bar
  taominus = alog(I0) - alog(I0*Exp(-tao) + sqrt(I0*Exp(-tao) + b)) 


  binned = [[reradius],[tao],[redata]]

  If KEYWORD_SET(taomax) then begin

    binned = [[reradius],[tao],[redata],[taomax]]

  Endif
  
  If KEYWORD_SET(taoplus) then begin
    
    sigma = (I0*Exp(-tao)-Sqrt(Irr))
    sigma[where(sigma LE 0)] = b
    taosigma = alog(Irr/sigma)
    taoerror = 1./Sqrt(irr)
  
  I0 = Make_array(n_elements(reradius), value = I0)
    
    binned = [[reradius],[tao],[Irr],[taominus],[taoplus],[sigma],[taoerror], [I0]]
    
    ;Note that for the paper I passed redata instead of Irr which means I didn't substract the background from the signal that I  compare to the modeled I0*exp(-optd)
    

  Endif


  return, binned

end