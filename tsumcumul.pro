;+
; PURPOSE:
;  This function uses the trapezoidal rule to perform a running
;  integration of tabulated data. It is adapted from the tsum routine
;  in the IDL astronomy user's library, but returns the integral at
;  each point of x, and not just over one large region.
;
; CATEGORY:
;  Numerical recipes
;
; INPUTS:
;  x: Input x values
;  y: Input y values
;
; OUTPUTS:
;  Integral(y dx) from x[0] to x[i]. A vector of the same length as x
;  and y. (not really. The lenght is n_elements(x)-1)
;
; MODIFICATION HISTORY:
;  September 2009: Lifted from IDL astronomy user's library by Chris Beaumont
;- April 2020; It added the option to linearly extrapolate the last element of the sum to make the output actually then same lenght of x and y ny Daniel Sega
function tsumcumul, x, y, Extrapolate=extrapolate
  compile_opt idl2
  On_error,2

  npar = N_params()
  if npar ne 2 then begin
    print, 'tsumcumul calling sequence:'
    print, 'result = tsumcumul(x,y)'
    return, !values.f_nan
  endif

  if n_elements(x) ne n_elements(y) then $
    message, 'X and Y must have the same number of elements'

  yy = y
  xx = x
  npts = n_elements(x)

  ; Compute areas of trapezoids and sum result
  xdif = xx[1:*] - xx
  yavg =  ( yy[0:npts-2] + yy[1:npts-1] ) / 2.
  sum = total( xdif*yavg, /cumul )
  
  ;;
  if keyword_set(extrapolate) then sum = append(sum, [2*sum[-1] - sum[-2]])

  return, sum
end