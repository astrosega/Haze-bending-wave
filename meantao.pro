;Outline
;It computes the average optical depth outside the wave region. From radial values 131940,131910km. If optical depth is not provided it attempts to compute it from the irradiance vector provided

;Inputs
;radius -> radiaul vector of the data
;data   -> data vector (irradiance)
;tao    -> keyword. set if you pass the optical depth instead of irradiance.

;Output
;g -> mean optical depth

function meantao, radius, data, tao=tao, long = long;tao is a flag

if keyword_set(tao) then begin
  r=findradius(radius,131940,131910)
  
  if keyword_set(long) then r=findradius(radius,131600,131700) 
  ;print,r
  g=mean(data[r[0]:r[1]], /Nan)
  ;print,data[r[0]:r[1]]
  
  return, g



endif
 
tao=odepth(radius,data)
r=findradius(radius,131650,131700)
 
 
  g=mean(tao[r[0]:r[1]])
  
  return, g
  
  end
