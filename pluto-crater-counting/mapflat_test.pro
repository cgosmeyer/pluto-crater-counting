PRO mapflat_test


read_jpeg, '../data/charon.jpg', planet, ctable, COLORS=!D.N_COLORS-1, /GRAYSCALE
length = (size(planet))[1]  ;long
height = (size(planet))[2]  ;lat


R_x = length / 2.
R_y = height / 2.
longitude=fltarr(length)
latitude=fltarr(height)
for x=0, R_x do begin
   for y=0, R_y do begin
        longitude[x] = x / (R_x)
        latitude[y] = 2. * atan(exp(y/R_y) - (3.14159)/ 2.)
        ;print, '___'
        ;print, x, y
        ;print, longitude[x], latitude[y]
    endfor
endfor

; 3-D coords:
S = 100. ; mapped radius
x_trans = S * cos(latitude) * cos(longitude)
y_trans = S * cos(latitude) * sin(longitude)
z_trans = S * sin(latitude)
help, x_trans
print, x_trans
print, y_trans

sg = griddata(x_trans, y_trans, planet[*]) ;, /sphere, /degrees, /grid, xout=findgen(90)/90.*360.-180., yout=findgen(90)/90.*180.-90.)
;Window, 1, xsize=length, ysize=height ;, xsize=(length_clip_r - length_clip_l), ysize=(height_clip_t - height_clip_b)
;; miller_cylindrical may be the best.
help, sg

map_set, 0,0, /aitoff, /isotropic, /grid, clip=0, limit=[-90,-90,90,90] ;, scale=100E6  ;/robinson, ;; this for setting up map display?                                             
result = map_image(sg, Startx, Starty, compress=1, missing=0, $
	               LATMIN=-90, LATMAX=90, LONMIN=-90, LONMAX=90)  ;; does the long, lat conversion?
display, result       
write_jpeg, '../outputs/charon_flat.jpg', result





END