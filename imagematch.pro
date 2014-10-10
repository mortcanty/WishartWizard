; docformat = 'rst'
; imagematch.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.

;+
; :Description:
; Function to match two polSAR images, that is,
; warp image (fid2) will exactly overlap base image (fid1)
; If PTS keyword is set, the provided ground reference
; points wil be used, otherwise they will be dertived from
; the four corners of the two images
; :Params:
;      fid1:  in, required
;         reference image
;      fid2:  in, out, required
;         target image
;      dims2: out, required
;         matched spatial subset of target 
;      map_info: in, out, required
;         geocoding information        
; :KEYWORDS:
;     pts: in, optional
;         ground reference points for matching
; :USES:
;    ENVI
; :RETURNS:
;        1 on success, else 0
; :Author:
;    Mort Canty (2014)
;-
function imagematch, fid1, fid2, dims2, map_info, pts=pts

  Compile_Opt idl2
;  Standard error handling.
  Catch, theError
  IF theError NE 0 THEN BEGIN
    Catch, /CANCEL
    void = Error_Message()
    RETURN, 0
  ENDIF
  
  envi_file_query, fid1, dims=dims1
  envi_file_query, fid2, dims=dims2
  cols1 = dims1[2]-dims1[1]+1
  rows1 = dims1[4]-dims1[3]+1
  map_info = envi_get_map_info(fid=fid1)
  envi_convert_file_coordinates, fid1, $
    dims1[1], dims1[3], e, n, /to_map
  map_info.mc = [0D,0D,e,n]
  if n_elements(pts) eq 0 then begin
; get tie points for geographic matching
    pts = fltarr(4,4)
    x1 = dims1[1]
    y1 = dims1[3]    
    iProj = envi_get_projection(fid=fid1)
    oProj = envi_get_projection(fid=fid2)     
    envi_convert_file_coordinates, fid1, x1, y1, e, n, /to_map
    envi_convert_projection_coordinates, e,n,iProj,e,n,oProj
    envi_convert_file_coordinates, fid2, x2, y2, e, n
    pts[*,0] = [x1,y1,x2,y2]   
    x1 = dims1[2]
    y1 = dims1[3]
    envi_convert_file_coordinates, fid1, x1, y1, e, n, /to_map
    envi_convert_projection_coordinates, e,n,iProj,e,n,oProj
    envi_convert_file_coordinates, fid2, x2, y2, e, n    
    pts[*,1] = [x1,y1,x2,y2]    
    x1 = dims1[1]
    y1 = dims1[4]
    envi_convert_file_coordinates, fid1, x1, y1, e, n, /to_map
    envi_convert_projection_coordinates, e,n,iProj,e,n,oProj
    envi_convert_file_coordinates, fid2, x2, y2, e, n
    pts[*,2] = [x1,y1,x2,y2]    
    x1 = dims1[2]
    y1 = dims1[4]
    envi_convert_file_coordinates, fid1, x1, y1, e, n, /to_map
    envi_convert_projection_coordinates, e,n,iProj,e,n,oProj
    envi_convert_file_coordinates, fid2, x2, y2, e, n
    pts[*,3] = [x1,y1,x2,y2]
  endif
; co-register
; -----workaround------------------
  x0 = map_info.mc[2]
  y0 = map_info.mc[3]
  xsize = map_info.ps[0]*cols1
  ysize = map_info.ps[1]*rows1
; ---------------------------------
  ENVI_DOIT, 'ENVI_REGISTER_DOIT', b_fid=fid1, w_fid=fid2, $
    w_dims=dims2, $
    pts=pts, r_fid=r_fid, $
    method=4, $
    x0 = x0, y0 = y0, $
    xsize = xsize, ysize = ysize, $
    /in_memory
  fid2 = r_fid 
  dims2 = [-1,0,cols1-1,0,rows1-1]
  
  return, 1
  
end