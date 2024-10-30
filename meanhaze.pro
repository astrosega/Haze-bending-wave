;Outline
;It computes the average optical depth within the wave region. From radial values 131812km ,131852km. If optical depth is not provided it attempts to compute it from the irradiance vector provided
;Used for Table 1 in Sega et al 2024. Icarus

;Inputs
;radius -> radiaul vector of the data
;data   -> data vector (irradiance)
;tao    -> keyword. set if you pass the optical depth instead of irradiance.

;Output
;g -> mean optical depth

function meanhaze, radius, data, tao=tao ;tao is a flag


  if keyword_set(tao) then begin
    r=findradius(radius,131902 - 90,131902-50)
    ;print,r
    g=mean(data[r[0]:r[1]], /Nan)
    ;print,data[r[0]:r[1]]

    return, g



  endif

  tao=odepth(radius,data)
  r=findradius(radius,131902 - 90,131902-50)


  g=mean(tao[r[0]:tao[1]])

  return, g

end
