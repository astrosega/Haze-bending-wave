;It compares the results from the model and the occultation and calculates the square of the diference of both. So the program tells you how good
; the fits is, i.e how much the simulated optical depth differs from the observed.

;Inputs
;radius -> (vector in the cassini data files giving the radial distance)
;data -> (vector in the cassini data files giving the irradiance)
;r -> radius vector in the simulation
;optd - simulated optical depth

;Output
;totales -> the sum of all the elements of the array made by substracting observed with simulated data

;Authors
;Daniel Sega

;Change Log:

;10/ 8 / 2016 -> changed differencia from the absolute value of the difference to the square of the difference,
;then I changed it back

function goodfit3,  radius, tao, r, optd, points_in_model, star, Taosigma=taosigma

  taofit=FLTARR(points_in_model)

  ;this is to match both x axis for the occultation and the data. E.I it makes sure that the theoretical data's indices matches the indices of the
  ;r(simulation) vector. Remember that data and radius are huge and span all the ring, so their indices don't match the indices of the vectors created for the sim
  ;
  ;The simulation has more points than the data so some points are removed. What the loop does is that it asks taosum -> "what's your value at
  ;radius=(so and so) kilometers, which corresponds to r[k] and to radius[srad]. Then it says, "well, put that value in taofit[k].
  ;The consequence is that I can now use r(simulation) as an x axis for the observed data, and compare.


  for k=0,points_in_model-1 do begin
    dump=min(abs(radius-r[k]),srad)
    taofit[k]=tao[srad]

    if star eq 'EpsCas104E' and srad eq 471 then begin
      taofit[k] = optd[k]
      taosigma[471] = 1
      
    endif
  endfor


  ;graphic6=plot(taofit)
  ;graphic7=plot(taosum)
  ;;
  if keyword_Set(taosigma) then diferencia=((taofit-optd)/taosigma)^2 $
  else begin
    diferencia=(taofit-optd)^2
  endelse
  ;plot,diferencia
  totales=total(diferencia)

  return, totales

end
