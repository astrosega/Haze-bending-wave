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