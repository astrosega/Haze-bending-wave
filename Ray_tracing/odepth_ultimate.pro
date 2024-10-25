;An iteration over odepth_occdatabin.pro that runs the program and saves the arrays for the occultations used in saturn_ultimate
;
;Log:
;7/10/2017 -> Created by DDSN
;

PRO odepth_ultimate

term  = '_data_rad_lon_nobin_cpck_'
ext  = '.sav'
output = '_Mimas53BW_400m_res'

restore, 'stars.sav'
           
Noccs = n_elements(stars)


star_rev = STRARR(Noccs)

;stars = STRMID(stars,0,9) + strmid(stars, 10, 1)
;stars = stars.insert(" ", 3)
;stars = stars.insert(" ", 7)
;stars = stars.insert(" ", 11)
;stars = stars.tolower()
;stars = stars.capwords()
;stars = stars.compress()

;m = where(strmatch(stars, 'AlpSco*'))
;stars[m]= stars[m].insert('B0', 6)

;rename, m, stars

;stars[m] = stars[m].remove(9,10)

;print, n_elements(stars)

;m = where(strmatch(stars, '*Cma*'))
;stars[m]= stars[m].replace('m', 'M')





for i=0, (noccs-1) do begin
  print,stars[i]
  stars[i] = rename_fn(stars[i])
  star_rev[i] = file_search(filepath(stars[i]+ '*', subdir = ['occs']))
  print, star_rev[i]
endfor

;print, star_rev




           
r1   =   131690.
r2   =   131950.
reso =        0.4 ;Km

For i=0, noccs-1 do begin
  Restore, star_rev[i]
  
  ;some of the values for the background have to be estimated because these data files lack data on the b_ring. The b data is estimates by looking
  ;for occultations of the same star with a closerby rev
  
  if     strmatch(stars[i], 'TheAra40*')   then binned = odepth_occdatabin(radius, data, reso, r1, r2, taoplus = 1, i=i, b=0.075) $
  else if strmatch(stars[i], 'TheCar190*') then binned = odepth_occdatabin(radius, data, reso, r1, r2, taoplus = 1, i=i, b=0.146) $
  else if strmatch(stars[i], 'AlpCru100*') then binned = odepth_occdatabin(radius, data, reso, r1, r2, taoplus = 1, i=i, b=0.7)$
  else if strmatch(stars[i], 'DelSco236*') then binned = odepth_occdatabin(radius, data, reso, r1, r2, taoplus = 1, i=i, b=0.27)$
  else if strmatch(stars[i], 'GamPeg36*')  then binned = odepth_occdatabin(radius, data, reso, r1, r2, taoplus = 1, i=i, b=0.44)$
  else if strmatch(stars[i], 'GamAra37*')  then binned = odepth_occdatabin(radius, data, reso, r1, r2, taoplus = 1, i=i, b=0.10)$
  else if strmatch(stars[i], 'EpsCas104*') then binned = odepth_occdatabin(radius, data, reso, r1, r2, taoplus = 1, i=i, b=0.10)$
  else if strmatch(stars[i], 'GamPeg172*') then binned = odepth_occdatabin(radius, data, reso, r1, r2, taoplus = 1, i=i, b=0.44)$
  else if strmatch(stars[i], 'AlpVir008I*')then binned = odepth_occdatabin(radius, data, reso, r1, r2, taoplus = 1, i=i, b=0.44)$
  else if strmatch(stars[i], 'KapCMa168*') then binned = odepth_occdatabin(radius, data, reso, r1, r2, taoplus = 1, i=i, b=0.44)$
  else if strmatch(stars[i], 'KapCen042I*')then binned = odepth_occdatabin(radius, data, reso, r1, r2, taoplus = 1, i=i, b=0.44)$

  else binned = odepth_occdatabin(radius, data, reso, r1, r2, taoplus = 1, i=i)
  
  radius   = binned[*, 0]
  tao      = binned[*, 1]
  irr      = binned[*, 2]
  taominus = binned[*, 3]
  taoplus  = binned[*, 4]
  Isigma   = binned[*, 5]
  I0       = binned[0, 7]
  
  save, radius, tao, taosigma, Irr, I0, filename=stars[i] + output + ext

Endfor

end
