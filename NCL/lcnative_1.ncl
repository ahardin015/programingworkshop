;*******************************************
; lcnative_1.ncl
;*******************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
begin
;*******************************************
; open file and read in data
;*******************************************
  f    = addfile ("pre.8912.mon.nc", "r")
  P    = f->pre
  LAT2D= f->lat
  LON2D= f->lon
;********************************************
;reorder arrays
;********************************************
  p = P(time|:,ycoord|:,xcoord|:)
  lon2d=LON2D(ycoord|:,xcoord|:)
  lat2d=LAT2D(ycoord|:,xcoord|:)
;********************************************
; get some parameters
;********************************************
  nlat = dimsizes(lat2d&ycoord)
  nlon = dimsizes(lon2d&xcoord)
;********************************************
; create plot
;********************************************
  wks  = gsn_open_wks ("ps", "lcnative")          ; open workstation
  gsn_define_colormap (wks,"gui_default")         ; choose color map
  
  res                        = True               ; plot mods desired
  res@cnFillOn               = True               ; color fill  
  res@cnLinesOn              = False              ; no contour lines
  res@cnLineLabelsOn         = False              ; no contour labels
  res@gsnSpreadColors        = True               ; use total colormap
  res@gsnSpreadColorStart    = 4 
  res@gsnSpreadColorEnd      = -1
  res@cnInfoLabelOn          = False              ; no contour info label
  
  res@mpGridLineDashPattern  = 2                  ; lat/lon lines as dashed
  res@pmTickMarkDisplayMode  = "Always"           ; turn on tickmarks
  
  res@tiMainString           = "Native Lambert Conformal Grid"
  res@tiMainFontHeightF      = 0.020              ; smaller title
  
  res@gsnAddCyclic           = False              ; regional data
; !!!!! any plot of data that is on a native grid, must use the "corners"
; method of zooming in on map.
  
  res@mpLimitMode       = "Corners"            ; choose range of map
  res@mpLeftCornerLatF  = lat2d(0,0)
  res@mpLeftCornerLonF  = lon2d(0,0)
  res@mpRightCornerLatF = lat2d(nlat-1,nlon-1)
  res@mpRightCornerLonF = lon2d(nlat-1,nlon-1)
  
; The following 4 pieces of information are REQUIRED to properly display
; data on a native lambert conformal grid.  This data should be specified
; somewhere in the model itself.
  
  res@mpProjection        = "LambertConformal"
  res@mpLambertParallel1F = 30.
  res@mpLambertParallel2F = 55.
  res@mpLambertMeridianF  = 45.
  
; usually, when data is placed onto a map, it is TRANSFORMED to the specified
; projection. Since this model is already on a native lambert conformal grid,
; we want to turn OFF the tranformation.
  
  res@tfDoNDCOverlay = True
  
  
  plot = gsn_csm_contour_map(wks,p(0,:,:),res)     ; Draw contours over a map.
  
end
