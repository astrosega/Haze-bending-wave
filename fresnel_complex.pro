;This function returns the value for the fresnel complex integral

;Input
;t (variable, this can be just x, as the integral is usually defined on text books)

;Output
;F(t) where F is the complex fresnel integral

;Authors
;Daniel Sega

FUNCTION FRESNEL_COMPLEX, t



  z=(dCOMPLEX((t/(Sqrt(2))),(t/(sqrt(2)))))

  z1=dCOMPLEX(1,1)

  F=(erf(z)*(sqrt(!dpi/2.)*(1./z1)))

  return,F

end