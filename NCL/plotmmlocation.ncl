load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
do j=1,6
if j .eq. 6 then
    j=j+1
end if
fname = "p"+j+"_100613.qcd"
col = numAsciiCol(fname)
row = numAsciiRow(fname)
print(col)
print(row)
data = asciiread(fname,(/row,col/),"float")

id = data(0,0)
time = data(:,1)
lat = data(:,2)
lon = data(:,3)*-1.0
tfast = data(:,4)














;;;;;;;;;;;;;;;;
;Plot Map
;;;;;;;;;;;;;;;

wks = gsn_open_wks("pdf","MM"+j)

res = True
res@mpFillOn = False
res@gsnFrame = False
res@gsnMaximize = True

res@mpOutlineBoundarySets = "AllBoundaries"
  res@mpDataBaseVersion     = "Ncarg4_1"
  res@mpDataSetName         = "Earth..4"
  res@mpMinLatF = min(lat)-0.25  
  res@mpMaxLatF = max(lat) 
  res@mpMinLonF = min(lon)-0.25 
  res@mpMaxLonF = max(lon) 
  mkres               = True
  mkres@gsMarkerIndex = 16 
  mkres@gsMarkerColor = "Green"
map = gsn_csm_map(wks,res)
do i=0,row-1,100
gsn_polymarker(wks,map,lon(i),lat(i),mkres)
end do
  mkres@gsMarkerIndex = 2 
  mkres@gsMarkerColor = "Red"
gsn_polymarker(wks,map,-101.709,35.233,mkres)
frame(wks)

plt = gsn_csm_xy(wks,time(32000::),tfast(32000::),res)
frame(wks)
delete(data)
delete(time)
delete(lon)
delete(lat)
delete(tfast)
delete(wks)
end do
end

