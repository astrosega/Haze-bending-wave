;It calculates the optical depth given Intensity (Irradiance) and radius.

;Parameters
;--------------
;Radius: array

;Data: array
;- Irradiance

;taomax: boolean
;if set to 1, then instead of defining a minimum irradiance, we define a maximum optical depth using the function taomax

;Output
;-----------------
;Optical depth

;Updates
;07/12/2015 - Created by DDSN
;11/12/2016 - Log creared 

FUNCTION odepth, radius,data, taomax=taomax,af=af

;encke gap is between 133450 and 133700


if keyword_set(af) then e=findradius(radius, 133450,133700) else e=findradius(radius, 137000,139000)

;This is between the A ring and F ring

;e=findradius(radius, 137000,139000)

Imax=mean(data[e[0]:e[1]])
print,'imax',imax

;this is a section of the b ring that is opaque

  rb = findradius(radius,107870,107900) ;Range from Colwell 2010

if rb[0] EQ rb[1] then print,'Check if there is B-ring data in this occultation'

 b=mean(data[rb[0]:rb[1]])

print,'b',b

;The stars aparent brightness will be given by the unoculted data minus the background.

I0=Imax-b

;print,I0
;I0=612
;data[where(data LE 0, /NULL)]=0.000000000000001

Irr = data - b

;Irr[where(Irr LE 1, /NULL)] = 1
;this is needed for very low angle optical depths, since the 'data' may equal 0: the ring gets completely opaque and tao diverges

Irr[where(Irr LE b, /NULL)] = b


;data[where(data LE 0, /NULL)]=1

;tao=-alog((data-b)/I0)

tao=alog(I0/Irr)

if KEYWORD_Set(taomax) then begin

taomaxi = taomax(radius, data)

tao[where(tao GE taomaxi)] = taomaxi
endif

return, tao

end