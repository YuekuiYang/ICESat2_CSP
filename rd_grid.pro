pro rd_grid
;
; This program reads in the output from precip_grid.pro calculates the precip amount and displays the result
;
domec_lat = -75.0 - 6.0 / 60.0
domec_lon = 123.0 + 24.0 / 60.0


input_file = 'Antarctic_IWC_Apr-Oct2021'

nbs = 0L
nr = 0L
imax = 0
jmax = 0
lat1 = 0.0
lat2 = 0.0
lon1 = 0.0
lon2 = 0.0

domec_obs = 0L
bsnow_cnt = 0L
bs_freq = 0.0

cap_bscat = fltarr(700,100000)

get_lun,lun
openr,lun,input_file
;readu,lun,imax,jmax,nr,nbs
readu,lun,imax,jmax,nr
grid_obs = fltarr(imax,jmax)
grid_cap = fltarr(imax,jmax)
grid_IWC = fltarr(imax,jmax)
grid_od = fltarr(imax,jmax)
readu,lun,lat1,lat2,lon1,lon2
readu,lun,grid_obs
readu,lun,grid_IWC
readu,lun,grid_cap
readu,lun,domec_obs,bsnow_cnt,bs_freq
close,lun

print,'nr, nbs = ', nr, nbs
bsnow_freq = float(bsnow_cnt) / float(domec_obs)
print,'Dome C blowing snow frequency = ', bsnow_freq, bsnow_cnt, domec_obs

grid_lat = fltarr(imax,jmax)
grid_lon = fltarr(imax,jmax)
latwid = 2.0
lonwid = 4.0

mindiff_lat = 999.0
mindiff_lon = 999.0

for j=0,jmax-1 do begin
  lat = lat1 + float(j) * latwid
    lat_diff = abs(lat - domec_lat)
    if (lat_diff lt mindiff_lat) then begin
       mindiff_lat = lat_diff
       minj = j
    endif
  for i=0,imax-1 do begin
    grid_lat(i,j) = lat
    lon = lon1 + float(i) * lonwid
    grid_lon(i,j) = lon
    lon_diff = abs(lon - domec_lon)
    if (lon_diff lt mindiff_lon) then begin
       mindiff_lon = lon_diff
       mini = i
    endif
  endfor
endfor


fall_speed = 0.10   ; m/s
hours = (4.0 * 31.0 + 3.0 * 30.0) * 24.0 ; hours for April 1 - October 31

ptr = where(grid_obs eq 0.0)
grid_obs(ptr) = 1.0
grid_IWC = grid_IWC / grid_obs

rate = grid_IWC * fall_speed  ; g/m^2/s
kg_per_hr = rate/1000.0*3600.0 ; kg/m^2/hr
kg_total = kg_per_hr * hours

roe = 1000.0 ; kg/m^3
meters = kg_total / roe
precip_mm = meters * 1000.0

help,precip_mm
print,'Dome C i,j = ', mini,minj
print,'Dome C precip (mm) = ', precip_mm(mini,minj)


device,decomposed=0
loadct,39,file='/home/spalm/colors.palm.tbl
loadct,38,file='/home/spalm/colors.palm.tbl
!P.Background=255

dat_lat_min = -87.0
dat_lat_max = -70.0
dat_lon_min = 0.0
dat_lon_max = 180.0
dat_lat_min = lat1
dat_lat_max = lat2
dat_lon_min = lon1
dat_lon_max = lon2

map_lat_min = -90.0
map_lat_max = -65.0
map_lon_min = -180.0
map_lon_max = 180.0
rot = 0
lonavg = (map_lon_min + map_lon_max) / 2.0

;goto,win6

ptr = where (grid_IWC gt 5000.0 or grid_IWC lt 0.)
demimg = grid_IWC 
demimg(ptr) = 0.0


 window,1,xsize=800,ysize=750,retain=2,title='Precip mm, scale: 2 - 16'
  map_set,map_lat_min,lonavg, rot, latdel = 5, londel = 20, glinestyle=1, /lambert,$
     glinethick=2, $
     limit=[map_lat_min,map_lon_min,map_lat_max,map_lon_max], label=2,charsize=2.0,  $
        /continents,mlinethick=2,/hires, color=0

  img = precip_mm
  img = bytscl(img,min=2.0,max=16.0)
  
  img = demimg
  img = bytscl(img,min=2000.0,max=4000.0)


  mapped_image = map_image (img,stx,sty,xs,ys, latmin=dat_lat_min, $
    latmax=dat_lat_max, lonmin=dat_lon_min, lonmax=dat_lon_max, compress=1, /bilinear)

ptr = where (mapped_image gt 120 and mapped_image lt 126)
mapped_image(ptr) = 0
  tv,mapped_image,stx,sty,xsize=xs,ysize=ys

  map_continents,/hires,mlinethick=3, /coasts
  map_grid, glinethick=5, latdel=5, latlab=0.0, londel=22.5, lonlab=-67.0, charsize=3, label=2, charthick=3

  xyouts,domec_lon,domec_lat,'X',charsize=3,color=255, charthick=5

stop

 window,2,xsize=800,ysize=750,retain=2, title='Precip Frequency, scale: 0.4 - 0.8'
  map_set,map_lat_min,lonavg, rot, latdel = 5, londel = 20, glinestyle=1, /lambert,$
     glinethick=2, $
     limit=[map_lat_min,map_lon_min,map_lat_max,map_lon_max], label=2,charsize=2.0,  $
        /continents,mlinethick=2,/hires, color=0


  precip_freq = fltarr(imax,jmax)
  ptr = where (grid_obs gt 0.0)
  precip_freq(ptr) = grid_cap(ptr) / grid_obs(ptr)

  img = precip_freq
  img = bytscl(img,min=0.4,max=0.8)

  mapped_image = map_image (img,stx,sty,xs,ys, latmin=dat_lat_min, $
    latmax=dat_lat_max, lonmin=dat_lon_min, lonmax=dat_lon_max, compress=1, /bilinear)

  tv,mapped_image,stx,sty,xsize=xs,ysize=ys

  map_continents,/hires,mlinethick=3, /coasts
  map_grid, glinethick=5, latdel=5, latlab=0.0, londel=22.5, lonlab=-67.0, charsize=3, label=2, charthick=3

window,3,xsize=800,ysize=750,retain=2, title='Number Observations, scale: 0 - 200'
  map_set,map_lat_min,lonavg, rot, latdel = 5, londel = 20, glinestyle=1, /lambert,$
     glinethick=2, $
     limit=[map_lat_min,map_lon_min,map_lat_max,map_lon_max], label=2,charsize=2.0,  $
        /continents,mlinethick=2,/hires, color=0

  img = grid_obs
  img = bytscl(img,min=0.0,max=200.0)

  mapped_image = map_image (img,stx,sty,xs,ys, latmin=dat_lat_min, $
    latmax=dat_lat_max, lonmin=dat_lon_min, lonmax=dat_lon_max, compress=1, /bilinear)

  tv,mapped_image,stx,sty,xsize=xs,ysize=ys

  map_continents,/hires,mlinethick=3, /coasts
  map_grid, glinethick=5, latdel=5, latlab=0.0, londel=22.5, lonlab=-67.0, charsize=3, label=2, charthick=3

window,4,xsize=800,ysize=750,retain=2, title='CAP Count, scale: 10 - 100'
  map_set,map_lat_min,lonavg, rot, latdel = 5, londel = 20, glinestyle=1, /lambert,$
     glinethick=2, $
     limit=[map_lat_min,map_lon_min,map_lat_max,map_lon_max], label=2,charsize=2.0,  $
        /continents,mlinethick=2,/hires, color=0

  img = grid_cap
  img = bytscl(img,min=10.0,max=100.0)

  mapped_image = map_image (img,stx,sty,xs,ys, latmin=dat_lat_min, $
    latmax=dat_lat_max, lonmin=dat_lon_min, lonmax=dat_lon_max, compress=1, /bilinear)

  tv,mapped_image,stx,sty,xsize=xs,ysize=ys

  map_continents,/hires,mlinethick=3, /coasts
  map_grid, glinethick=5, latdel=5, latlab=0.0, londel=22.5, lonlab=-67.0, charsize=3, label=2, charthick=3

win5:

window,5,xsize=800,ysize=750,retain=2, title='Column Optical Depth, scale 0 - 0.50'
  map_set,map_lat_min,lonavg, rot, latdel = 5, londel = 20, glinestyle=1, /lambert,$
     glinethick=2, $
     limit=[map_lat_min,map_lon_min,map_lat_max,map_lon_max], label=2,charsize=2.0,  $
        /continents,mlinethick=2,/hires, color=0


  cod = fltarr(imax,jmax)
  ptr = where (grid_cap gt 0.0)
  cod(ptr) = grid_od(ptr) / grid_cap(ptr)

  img = cod
  img = bytscl(img,min=0.0,max=0.50)

  mapped_image = map_image (img,stx,sty,xs,ys, latmin=dat_lat_min, $
    latmax=dat_lat_max, lonmin=dat_lon_min, lonmax=dat_lon_max, compress=1, /bilinear)

  tv,mapped_image,stx,sty,xsize=xs,ysize=ys

  map_continents,/hires,mlinethick=3, /coasts
  map_grid, glinethick=5, latdel=5, latlab=0.0, londel=22.5, lonlab=-67.0, charsize=3, label=2, charthick=3

win6:

navg = 40
image = fltarr(700,10000)

j = 0
for i=0,nbs-navg,navg do begin
   image(*,j) = total(cap_bscat(*,i:i+navg-1),2) / float(navg)
   j++
endfor

window,6,xsize=j,ysize=700
img = bytscl(image(*,0:j-1),min=-1.0e-7, max=1.50e-5)
tv,rotate(img,3)


end
