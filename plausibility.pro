; docformat = 'rst'
; plausibility.pro
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
; Procedure to check the plausibilzty of reference/warp gcps
;          Li, H., Manjunath, B. S., and Mitra, S. K. (1995).
;          A contour-based approach ;to multisensor image
;          registration. IEEE Transactions on Image Processing,
;          4(3), 320–334.
; :Params:
;      pts:  in, required
;         gcps in ENVI format
;      outpts:  in, out, required
;         gcps with outliers removed
; :returns:
;      1 on success else 0        
; :USES:
;    ENVI
; :Author:
;    Mort Canty (2014)
;-
function plausibility, pts, outpts

  Compile_Opt idl2
  ;  Standard error handling.
  Catch, theError
  IF theError NE 0 THEN BEGIN
    Catch, /CANCEL
    void = Error_Message()
  ENDIF
  
  base_gcps = pts[0:1,*]
  warp_gcps = pts[2:3,*]
  n_gcps = n_elements(pts[0,*])
  ratios = fltarr(n_gcps*(n_gcps-1)/2,3)
  k=0L
  for i=0,n_gcps-1 do for j=i+1,n_gcps-1 do begin
    den = norm(warp_gcps[*,i]-warp_gcps[*,j])
    if den gt 0 then ratios[k,*] = [norm(base_gcps[*,i]-base_gcps[*,j])/den,i,j] $
    else ratios[k,*]=[10,i,j]
    k=k+1
  endfor
  ratios[*,0]=alog(ratios[*,0])
  hist = histogram(ratios[*,0],nbins=100,min=-1.0,max=1.0,reverse_indices=R)
; just take the ratios in the histogram maximum bin
  i = (where(hist eq max(hist)))[0]
  max_indices = R[R[i] : R[i+1]-1]
  ratios = ratios[max_indices,*] 
; extract the GCPs from the array of ratios
  k = 0
  for i=0,n_gcps-1 do begin
    indices = where(round(ratios[*,1]) eq i,count1)
    indices = where(round(ratios[*,2]) eq i,count2)
    if (count1 ne 0) and (count2 ne 0) then begin
      outpts[0:1,k]=base_gcps[*,i]
      outpts[2:3,k]=warp_gcps[*,i]
      k=k+1
    endif
  endfor
  if k ge 4 then begin
    outpts = outpts[*,0:k-1]
    return, 1
  end else return, 0  
  
end  