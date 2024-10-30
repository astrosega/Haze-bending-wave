FUNCTION findradius, radius,r1,r2


  dump=min(abs(radius-r1),a)

  dump=min(abs(radius-r2),b)


  if (a GT b) then begin
    
    return, [b,a]

  endif else begin

   return, [a,b]

  endelse

  return, 0

end