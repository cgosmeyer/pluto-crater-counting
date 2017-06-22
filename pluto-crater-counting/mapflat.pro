;+
;Transforms New Horizons JPEGS of Pluto and Charon into 'flat' transformed
;maps.
;
;Author:
;
;    C.M. Gosmeyer, Oct. 2015
;    
;Use:
;    
;    idl> mapfloat
;    
;Outputs: 
;
;    
;References:
;    GRIDDATA: http://www.exelisvis.com/docs/GRIDDATA.html
;    TRIANGULATE: http://www.exelisvis.com/docs/TRIANGULATE.html
;    MAP_SET: http://www.exelisvis.com/docs/MAP_SET_Procedure.html
;
;Notes:
;    Griddata should provide the output for triangulate.
;    Do I do this on the individual squares of images, or on my mosaics?
;    
;-

; How do I specify the coordinates using pixels??? 
; Make use of F to make a transformation?
; Knowing angular size of Pluto.
; (To break into lat/long need know tilt, angular size of disk, and rotation.)
; Place a grid on top of Pluto, so have a lat/long of each pixel.
; Then use a transformation to project it flat.
;  ->Use MAP_SET !! 
; Or looking at Sean's thing.... get the lat/long, throw through
; map_set, then griddata to grid the points.
; 
; I'll never be able to do a 'perfect' transformation. I'll just add a
; weight for the craters counted at the edges. Make it a function of 
; longitude from the center.




;Sphere From Cartesian Coordinates
;Result = GRIDDATA( X, Y, Z,  F,  /SPHERE )

;Sphere From Spherical Coordinates
;Result = GRIDDATA( Lon,  Lat,  F,  /SPHERE, METHOD='NearestNeighbor', TRIANGLES=array, /SPHERE) 

; In order to display the flattened image:
; See Example 2 in http://www.exelisvis.com/docs/trigrid.html#T_809226861_1095045
;f = SIN(lon * !DTOR)^2 * COS(lat * !DTOR)
;TRIANGULATE, X, Y, Triangles [, B] [, CONNECTIVITY=variable] [, SPHERE=variable [, /DEGREES]] [, FVALUE=variable] [, REPEATS=variable] [, TOLERANCE=value]


;; Figure out 'long and lat' of each pixel.
;; Overplot a grid.

;; IDEA:
;; 1. find x,y center of planet and radius in pixels
;; 2. Loop through each pixel, until hit radius, assigning a lat, long to each pixel.
;; (or is 2. essentially what griddata does?)

;; PLUTO: 4330 pixel diam; 330 l, 310 r, 310 t, 368 b
;;     Crater center at 2110, 3290
;; CHaron: 1620 pixel diam; 178 l, 200 r, 167 t, 200 b    
;;     to center use crater left rim at 868, 1070

;; This did something. Not sure what, but it did something

read_jpeg, '../data/pluto.jpg', planet, ctable, COLORS=!D.N_COLORS-1, /GRAYSCALE
length = (size(planet))[1]  ;long
height = (size(planet))[2]  ;lat
;; change these to get different 'clips' of planet in order to do
;; a higher-res transformed image.
;; Pluto should be divided 4 x 4, so in longitude and latitude space, 
;;     0-45, 45-90, 90-135, 135-180 degrees.
;;     might be able ot get away with not transforming middle four regions?
length_clip_l = 0.
length_clip_r = long(float(length) / 4.) * 1.
height_clip_b = 0.
height_clip_t = long(float(height) / 4.) * 1.
planet_clip = planet[length_clip_l:length_clip_r, height_clip_b:height_clip_t]



Window, 1, xsize=length, ysize=height ;, xsize=(length_clip_r - length_clip_l), ysize=(height_clip_t - height_clip_b)
;; miller_cylindrical may be the best.

map_set, 0,0, /aitoff, /isotropic, /grid, clip=0, limit=[-90,-90,90,90] ;, scale=100E6  ;/robinson, ;; this for setting up map display?                                             
result = map_image(planet, Startx, Starty, compress=1, missing=0, $
	               LATMIN=-90, LATMAX=90, LONMIN=-90, LONMAX=90)  ;; does the long, lat conversion?
	               ;LATMIN=-90, LATMAX=-45, LONMIN=0, LONMAX=45)

;tv, result, Startx, Starty     
display, result       
write_jpeg, '../outputs/charon_flat.jpg', result


;; Or try:
; http://www.exelisvis.com/docs/map_proj_init.html
; http://www.exelisvis.com/docs/MAP_PROJ_IMAGE.html
Window, 1, xsize=length, ysize=height 
;; resize image data
;planet_resize = REBIN(planet, 180, 180 )
sMap = MAP_PROJ_INIT('azimuthal')
mapStruct = MAP_PROJ_INIT( 'gnomonic', /GCTP )
;coordtrans = MAP_PROJ_INVERSE (planet, MAP_STRUCTURE=mapStruct )
result = MAP_PROJ_IMAGE(planet, [-90, -90, 90, 90], IMAGE_STRUCTURE=sMap, $  ;[-90, -90, 90, 90]
	MAP_STRUCTURE=mapStruct, /bilinear) ;, missing=0 ) 
display, result
write_jpeg, '../outputs/pluto_flat.jpg', result



;; or again try:

data = planet
cgMap_Set, /Cylindrical, 0, 180, Position=[0.1, 0.1, 0.9, 0.8], Limit=[-90, 90, 90, 270]
Triangulate, lons, lats, triangles
gridded = Trigrid(lons, lats, data, triangles, [0.5,0.5], [0,-90,360,90], XGrid=x, YGrid=y)
cgContour, gridded, x, y, Levels=levels, /CELL_FILL, C_COLORS=c_colors, /OVERPLOT
cgMap_Grid, LatDel=20, LonDel=40, /Box_Axes



;;; maybe use example to build a sphere, whose lat and lon can be projected onto image
;;; http://www.exelisvis.com/docs/MAP_PATCH.html



R_x = length / 2.
R_y = height / 2.
longitude=[]
latitude=[]
for x=0, length do begin
   for y=0, height do begin
        longitude[x] = x / (R_x)
        latitude[y] = 2. * np.atan(np.exp(y/R_y) - np.pi / 2.
    endfor
endfor

; 3-D coords:
S = 1. # mapped radius
x_trans = S * np.cos(latitude) * np.cos(longitude)
y_trans = S * np.cos(latitude) * np.sin(longitude)
z_trans = S * np.sin(latitude)

s = griddata(longitude, latitude, planet[*], /sphere, /degrees)
Window, 1, xsize=length, ysize=height ;, xsize=(length_clip_r - length_clip_l), ysize=(height_clip_t - height_clip_b)
;; miller_cylindrical may be the best.

map_set, 0,0, /aitoff, /isotropic, /grid, clip=0, limit=[-90,-90,90,90] ;, scale=100E6  ;/robinson, ;; this for setting up map display?                                             
result = map_image(s, Startx, Starty, compress=1, missing=0, $
	               LATMIN=-90, LATMAX=90, LONMIN=-90, LONMAX=90)  ;; does the long, lat conversion?
display, result       
write_jpeg, '../outputs/charon_flat.jpg', result




; need transform pluto[0], pluto[1] to long, lat
longitude = x / R_x
latitude = 2. * np.atan(np.exp(y/R_y)) - np.pi / 2.

s = griddata(size(pluto[0]), size(pluto[1]), pluto[*], ...)


;; Using the solution above, grid pixels to sphere.
s = griddata(lon[hvlow_a], lat[hvlow_a], rfuva_sec[hvlow_a], $
	/sphere, /degrees, xout=findgen(90)/90.*360.-180., yout=findgen(90)/90.*180.-90., /grid)
display, s, findgen(90)/90.*360.-180., findgen(90)/90.*180.-90.
contour, s, findgen(90)/90.*360.-180., findgen(90)/90.*180.-90., col=cgcolor('red'), nlevels=5, /over

;; Plot map
map_set,0,0,/aitoff,/grid,/isotropic,glines=0,/noborder,chars=1.5,/cont

