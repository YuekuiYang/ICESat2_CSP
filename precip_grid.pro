pro precip_grid, restart

; This program reads in ATL09 files for a given time period and produces a gridded file of Ice Water Content (IWC) from which 
; precipitation is computed as described in Palm and Yang 2025
;restart = 0 is the normal way to run
;restart = 1 is used to resume the program after it had failed (usually do to an internet connection drop)

threshold = 2.00e-6 ; value above which CAP is likely. This is compared with the average signal of the 1 km above the surface.

output_file = 'Antarctic_IWC_Apr-Oct2021'  ; this is the output file that contains the gridded IWC
get_lun,output_lun

domec_lat = -75.0 - 6.0 / 60.0
domec_lon = 123.0 + 24.0 / 60.0
domec_i = 76
domec_j = 6

domec_obs = 0L
bsnow_cnt = 0L
bs_freq = 0.0

;set up the grid

nr = 0L

;define area of grid
lat1 = -88.0
lon1 = -180.0
lat2 = -65.0
lon2 = 180.0
latwid = 2.0 ;grid latitude width
lonwid = 4.0 ;grid longitude width

imax = fix((lon2 - lon1) / lonwid) + 1
jmax = fix((lat2 - lat1) / latwid) + 1
grid_dat = fltarr(imax,jmax)
grid_obs = fltarr(imax,jmax)
grid_cap = fltarr(imax,jmax)
grid_lat = fltarr(imax,jmax)
grid_lon = fltarr(imax,jmax)

;define the grid lat/lons
for j=0,jmax-1 do begin
  lat = lat1 + float(j) * latwid
  for i=0,imax-1 do begin
    grid_lat(i,j) = lat
    grid_lon(i,j) = lon1 + float(i) * lonwid
  endfor
endfor


; This is used to resume the program in case it fails before finishing.
; normally this is not used
nrecs_restart = 0
if (restart) then begin
  openr,output_lun,output_file
  readu,output_lun,imax,jmax,nr
  readu,output_lun,lat1,lat2,lon1,lon2
  readu,output_lun,grid_obs
  readu,output_lun,grid_dat
  readu,output_lun,grid_cap
  close,output_lun

  print,'nr = ', nr
  nrecs_restart = nr
endif

profn = 1  ;ATL09 contains 3 backscatter profiles. 1 and 3 are the best (higher S/N)

prof = 'profile_' + strcompress(string(profn,format='(i1)'),/remove_all)
prof0 = prof
prof = prof + '/high_rate'

;define ATL09 parameters needed
var_name = strarr(12)
var_name(0) = prof + '/cab_prof'
var_name(1) = prof + '/latitude'
var_name(2) = prof + '/longitude'
var_name(3) = prof + '/delta_time'
var_name(4) = 'ancillary_data/atlas_sdp_gps_epoch'
var_name(8) = prof + '/surface_bin'
var_name(9) = prof + '/solar_elevation'
var_name(10) = prof + '/dem_h'
var_name(11) = prof + '/bsnow_h'
prof = prof0 + '/low_rate'
var_name(5) = prof + '/mol_backscatter'


nr = 0L
filename = ''
openr,1,'flist_Apr-Oct2021'  ;this file contains the list of ATL09 files to process

if (restart) then begin
  for i=1,nrecs_restart do begin
    readf,1,filename
  endfor
  nr = long(nrecs_restart)
endif


while (NOT EOF(1)) do begin

; /css/icesat-2/ATLAS/ATL09.006/2021.04.01/ATL09_20210401001541_01101101_006_02.h5
  readf,1,filename
  nr++

  print,'Processing file ', filename
  file_id = H5F_OPEN(filename)

; get calibrated attenuated backscatter (CAB)
  dataset_id = H5D_OPEN(file_id, var_name(0))
  cab = H5D_READ(dataset_id)
  H5D_close,dataset_id

;get latitude
  dataset_id = H5D_OPEN(file_id, var_name(1))
  latitude = H5D_READ(dataset_id)
  H5D_close,dataset_id

;get longitude
  dataset_id = H5D_OPEN(file_id, var_name(2))
  longitude = H5D_READ(dataset_id)
  H5D_close,dataset_id

; get surface bin
  dataset_id = H5D_OPEN(file_id, var_name(8))
  surf_bin = H5D_READ(dataset_id)
  H5D_close,dataset_id

;get solar elevation
  dataset_id = H5D_OPEN(file_id, var_name(9))
  solar_elev = H5D_READ(dataset_id)
  H5D_close,dataset_id

;get DEM
  dataset_id = H5D_OPEN(file_id, var_name(10))
  dem = H5D_READ(dataset_id)
  H5D_close,dataset_id

;get blowing snow data
  dataset_id = H5D_OPEN(file_id, var_name(11))
  bsnow_h = H5D_READ(dataset_id)
  H5D_close,dataset_id

;get molecular backscatter
  dataset_id = H5D_OPEN(file_id, var_name(5))
  mol_back = H5D_READ(dataset_id)
  H5D_close,dataset_id


  H5F_close,file_id

;The molecular backscatter is at 1 second resolution, while the backscatter profiles are at 25 Hz
;The below interpolates the molecular backscatter to 25 Hz
  r = size(latitude)
  molec_backs = congrid(mol_back,700,r(1),/interp)
  help, molec_backs

;cycle through the grid 

  for ig=0,imax-1 do begin
    for jg=0,jmax-1 do begin


      latmin = grid_lat(ig,jg) - (latwid/2.0)
      latmax = grid_lat(ig,jg) + (latwid/2.0)
      lonmin = grid_lon(ig,jg) - (lonwid/2.0)
      lonmax = grid_lon(ig,jg) + (lonwid/2.0)

      ptr = where(latitude gt latmin and latitude lt latmax and longitude gt lonmin and longitude lt lonmax and dem gt 2500.0)
      bscat = cab(*,ptr)
      lat = latitude(ptr)
      lon = longitude(ptr)
      sbin = surf_bin(ptr)
      selev = solar_elev(ptr)
      bsnow = bsnow_h(ptr)
      molecbacks = molec_backs(*,ptr)

      r = size(lat)
      nrecs = r(1)
      if (nrecs gt 100) then begin ;there must be at least 100 profiles in a given grid box to proceed
        grid_obs(ig,jg)++

        np = 0
        sum = 0.0
        avg_sig = fltarr(2000)

        bs_cnt = 0L
        bsh_sum = 0.0
        for i=0,nrecs-1 do begin
          if (sbin(i) lt 700 and sbin(i) gt 500) then begin ; must have surface detected
             b2 = sbin(i) - 2 ; two bins above surface bin to be safe (not to include ground signal)
             b1 = sbin(i) -  33 ; this is 1 km above the ground
             if (bsnow(i) le 500.0) then begin ; if true, there is a blowing snow layer present. Set b2 to just above its top
               bs_cnt++
               bsh_sum = bsh_sum + bsnow(i)
               b2 = b2 - fix(bsnow(i)/30.0) - 1
             endif
             sig = 0.0
             scnt = 0.0
             for b=b1,b2 do begin
               if (bscat(b,i) lt 2.0e-4) then begin  ; make sure we are not including any ground signal or cloud
                 sig = sig + bscat(b,i)
                 scnt++
               endif
             endfor
             avg_sig(np) = 0.0
             if (scnt gt 0.0) then $
                avg_sig(np) = sig / scnt ; average backscatter

             if (avg_sig(np) gt threshold) then begin ; test against threshold
               sum = sum + avg_sig(np)
               np++
               nbins = b2 - 200 
             endif
          endif
        endfor

        if (ig eq domec_i and jg eq domec_j) then begin  ;this was added to compute the blowing snow frequency at Dome C
           domec_obs++
           bs_height = 0.0
           if (bs_cnt gt 0) then $
             bs_height = bsh_sum / float(bs_cnt)
           if (bs_cnt gt 500 and bs_height gt 30.0) then bsnow_cnt++
           bs_freq = float(bsnow_cnt) / float(domec_obs)
        endif

; np is now the number of backscatter profiles identified as CAP for this grid box

        max_sig = 0.0
        avg_sig_gbox = 0.0
        med_sig_gbox = 0.0
        if (np gt 25) then begin   ; must have more than 25 CAP profiles for a given grid box to proceed
          avg_sig_gbox = sum / float(np)
          max_sig = max(avg_sig(0:np-1),maxi)
          if (maxi gt 0 and maxi lt np-1) then $
            max_sig = (max_sig + avg_sig(maxi-1) + avg_sig(maxi+1)) / 3.0
            med_sig_gbox = median(avg_sig(0:np-1))
        endif

        avg_selev = total(selev) / float(nrecs)

        print,'avg_sig_gbox = ', avg_sig_gbox, max_sig, avg_selev

        IWC = 0.0
        S = 30.0
        if (avg_sig_gbox gt threshold) then begin
         ; likely CAP in this grid box. Compute IWC
           grid_cap(ig,jg)++
           Bscatt = avg_sig_gbox
           extinc = S * Bscatt
           IWC = 186.0 * extinc^1.15   ;g/m^3  ;This eqn from Heymsfield et al., 2014. See that paper or Palm and Yang 2025
        endif

        print,'IWC, bs_freq, ig, jg = ', IWC, bs_freq, ig, jg
        grid_dat(ig,jg) = grid_dat(ig,jg) + IWC

      endif ; if (nrecs gt 100) then begin

    endfor
  endfor

  if ((nr mod 50) eq 0) then begin
     openw,output_lun,output_file
     writeu,output_lun,imax,jmax,nr
     writeu,output_lun,lat1,lat2,lon1,lon2
     writeu,output_lun,grid_obs
     writeu,output_lun,grid_dat
     writeu,output_lun,grid_cap
     writeu,output_lun,domec_obs,bsnow_cnt,bs_freq
     close,output_lun
  endif

endwhile

openw,output_lun,output_file
writeu,output_lun,imax,jmax,nr
writeu,output_lun,lat1,lat2,lon1,lon2
writeu,output_lun,grid_obs
writeu,output_lun,grid_dat
writeu,output_lun,grid_cap
writeu,output_lun,domec_obs,bsnow_cnt,bs_freq
close,output_lun


end

