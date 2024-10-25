FUNCTION MIMAS53, radius, long = long


dump=min(abs(radius-131600),a)
 
dump=min(abs(radius-132100),b)


;dump=min(abs(radius-131000),a)

;dump=min(abs(radius-133000),b)


if (a GT b) then begin

return, [b,a]

endif else begin
  
return,[a,b]

endelse

return, 0

end