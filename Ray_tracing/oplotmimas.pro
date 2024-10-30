function oplotmimas, radius, data, nsum = nsum, rebin=rebin, star=star, title=title, Tao=taosum, ANGLES=angles, PHASE=phase, uperror=uperror,downerror=downerror, scl=SCL,layout=layout,modelplot=modelplot,optional=optional,aniso=aniso, joshplot=joshplot

;A plot warpper made to overplot the data to the model's prediction
;Takes data(irradiance) and radius and plots optical depth. It is made to overplot over the result of the simulation. 
;If you specify Angles and phase it will print the angles of the occultation to the plot. If error is provided, it will plot error bars.
;If star is provided it places the name of the star on the title of the plot.

rv = 131902.0

  IF ~KEYWORD_SET(nsum) THEN nsum = 50
 
  IF ~KEYWORD_SET(star) THEN star = ''
  
  IF ~KEYWORD_SET(taosum) THEN begin
  
   IF KEYWORD_SET(rebin) THEN begin
  
     taosum = odepth_occdatabin(radius,data, rebin)
  
    bins = n_elements(data)/rebin
  
    radius=rebin(radius[0:bins*rebin-1], bins)
  
    ENDIF ELSE begin
    
     tao = odepth(radius,data)
     taosum=smooth(tao,nsum)
     endelse
  
  endif

  r   = mimas53(radius)


 
  ;graphic=plot(radius[r[0]:r[1]],taosum[r[0]:r[1]],xrange=[131700,132000],yrange=[0,5],xtitle='Distance from Saturn´s center [Km]',ytitle='Optical Depth',color='red',overplot=1)
  
  if KEYWORD_SET(title) then begin
    
    if Keyword_SET(SCL) then begin

     ; restore, "scl.sav"

     if ~keyword_set(layout) then  newm= plot(rscl-rv, optd, overplot=1,thick=1, color='g', name='SCL Theory',linestyle=5,/current)
      
      endif
      
      if keyword_set(layout) then begin
        
       if keyword_set(joshplot) then begin
        
        graphic = plot(radius[r[0]+30:r[1]] -rv, taosum[r[0]+30:r[1]], yrange=[0,2.], color='red', font_size = 22, $
          XTICKFONT_SIZE=24, YTICKFONT_SIZE=24, Name = 'Data' ,dim=[1820,980],/current,overplot=1,transparency=30)
    

          
        graphic = plot(radius[r[0]:r[1]] -rv, taosum[r[0]:r[1]], yrange=[0,2.], xtitle='Distance from resonance [km]', color='red', font_size = 22, $
        XTICKFONT_SIZE=24, YTICKFONT_SIZE=24, Name = 'Data' ,dim=[1820,980],/current,layout=[2,1,1],margin=[.13, .15, .05, .05],ytitle='Optical Depth')
        
        
                    if Keyword_set(uperror) then begin
              
                      downerror[0]  = downerror[1]   ;this smooths the ends that diverge due to the binning process.
                      downerror[-1] = downerror[-2]
              
                      uperror[0]  = uperror[1]   ;this smooths the ends that diverge due to the binning process.
                      uperror[-1] = uperror[-2]
              
              
                      graf2 =fillplot(radius[r[0]:r[1]]-rv,[[uperror[r[0]:r[1]]],[downerror[r[0]:r[1]]]], overplot = 1,thick=1,fill_color="red",fill_transparency=80.,color ='red')
              
              
              
                    endif
                    
                    graphic = plot(radius[r[0]:r[1]] -rv, taosum[r[0]:r[1]], yrange=[0,5], xtitle='Distance from resonance [km]', ytitle='Optical Depth', color='red', font_size = 15, $
                      XTICKFONT_SIZE=24, YTICKFONT_SIZE=24, Name = 'Data',overplot=1,thick=1.5,layout=layout)

                      
 endif
        
              endif else begin
                
                if Keyword_set(uperror) then begin

                  downerror[0]  = downerror[1]   ;this smooths the ends that diverge due to the binning process.
                  downerror[-1] = downerror[-2]

                  uperror[0]  = uperror[1]   ;this smooths the ends that diverge due to the binning process.
                  uperror[-1] = uperror[-2]


                  graf2 =fillplot(radius[r[0]:r[1]]-rv,[[uperror[r[0]:r[1]]],[downerror[r[0]:r[1]]]], overplot = 1,thick=1,fill_color="red",fill_transparency=80.,linestyle=6)



                endif
  
  graphic = plot(radius[r[0]:r[1]] -rv, taosum[r[0]:r[1]], yrange=[0,5], color='red', Name = 'Data',/overplot,thick=.9 )
  graphic = plot(radius[r[0]:r[1]] -rv, taosum[r[0]:r[1]], yrange=[0,5], xtitle='Distance from resonance [km]', ytitle='Optical Depth', color='red', font_size = 14, XTICKFONT_SIZE=24, YTICKFONT_SIZE=24, Name = 'Data',overplot=1,thick=1)
  
if keyword_set(layout) then graphic = plot(radius[r[0]:r[1]] -rv, taosum[r[0]:r[1]], yrange=[0,5], xtitle='Distance from resonance [km]', ytitle='Optical Depth', color='red', font_size = 14, XTICKFONT_SIZE=24, YTICKFONT_SIZE=24, Name = 'Data',overplot=1,thick=2)

  
; graphic.title.font_size = 24
 endelse

    if keyword_set(layout) and keyword_set(joshplot) then leg  = legend(position=[.66,.93], font_size = 20,target=[modelplot,optional,graphic]) else if keyword_set(aniso) then leg  = legend(position=[.8,.8], font_size = 20,target=[modelplot,modelplot2,graphic]) else if keyword_set(modelplot) then leg  = legend(position=[.8,.8], font_size = 20,target=[modelplot,graphic])

  if KEYWORD_SET(angles) and (keyword_set(modelplot)  or keyword_set(aniso)) then begin
    
   junk =  text(0.2, 0.80, 'B  =' + string(angles[0]*180./!dpi,FORMAT='(F6.2)')+'°',font_size=18)
   junk =  text(0.2, 0.75, '$\phi$  =' + string(angles[1]*180./!dpi,FORMAT='(F6.2)')+'°',font_size=18)
   junk =  text(0.2, 0.7, '$B_{eff}$ =' + string(angles[2]*180./!dpi,FORMAT='(F6.2)')+'°',font_size=18)

   
   endif

   
  if KEYWORD_SET(phase) then junk =  text(0.2, 0.70, 'Phase   =' + string(phase))
  
;  if Keyword_set(uperror) then begin
    
 ;  downerror[0]  = downerror[1]   ;this smooths the ends that diverge due to the binning process.
 ;  downerror[-1] = downerror[-2]
    
  ;  uperror[0]  = uperror[1]   ;this smooths the ends that diverge due to the binning process.
   ; uperror[-1] = uperror[-2]
     
     
    ;  graf2 =fillplot(radius[r[0]:r[1]]-rv,[[uperror[r[0]:r[1]]],[downerror[r[0]:r[1]]]], overplot = 1,thick=1,fill_color="red",fill_transparency=70.,color ='cornflower',name='Haze')

           

  ;endif

  if keyword_set(scl) and ~keyword_set(layout) then newm = plot(rscl-rv, optd, overplot=1,thick=3, color='g', name='SCL Theory',linestyle=5)

  
  endif else begin
    
    graphic=plot(radius[r[0]:r[1]-rv],taosum[r[0]:r[1]],yrange=[0,10],xtitle='Distance from Saturn´s center [km]',ytitle='Optical Depth',color='red',overplot=1,$
      XTICKFONT_SIZE=24, YTICKFONT_SIZE=24 , Font_Size = 24)
  
  endelse





  ;backColor = cgColor("White")


;  oplot,radius[r[0]:r[1]],tao[r[0]:r[1]],nsum = nsum,color=1000
  ;background=backcolor, color=0
   ; ytitle='Optical Depth',$
    ;font=50,title='Observed optical depth of Mimas 5:3 Bending Wave'
    
    return, graphic

end