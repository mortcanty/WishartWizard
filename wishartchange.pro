; docformat = 'rst'
; wishartchange.pro
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
;    Polarimetric change detection for polarised SAR images
;    (multi-look covariance representation)  
;    Perform test for equality of two complex 1x1, 2x2 or 3x3
;    variance-covariance matrices representing
;      full polarimetry (9 bands)
;      dual polarimetry (4 bands
;      single polarimetry (1 band)
;    Based on a Matlab script by Allan Nielsen.      
; :Params:
;    fid1, fid2:
;       file ids, image1, image2
;    dims1, dims2:
;       spatial dimensions, image1, image2   
;    n1, n2:
;       ENLs for image1, image2           
;    bands:  
;       number of polarimetric bands        
; :KEYWORDS:
;    changestat: out
;       change statistic
;    change prob: out
;       change probability   
;    span1: out
;       span of image1        
; :Uses:
;    ENVI 
; :RETURNS:
;    1 on success else 0 
; :Author:
;    Mort Canty (2014)        
;-
FUNCTION wishartchange, fid1, fid2, dims1, dims2, n1, n2, bands, changestat=Z, changeprob=CP, span1=span1
; calculate determinants and change statistic distribution parameters
  case bands of
    9: begin
      k1 = n1*envi_get_data(fid=fid1,dims=dims1,pos=0)  ;c11
      a1 = envi_get_data(fid=fid1,dims=dims1,pos=1)  ;c12
      im = envi_get_data(fid=fid1,dims=dims1,pos=2)
      a1 = n1*complex(a1,im)
      rho1 = envi_get_data(fid=fid1,dims=dims1,pos=3) ;c13
      im = envi_get_data(fid=fid1,dims=dims1,pos=4)
      rho1 = n1*complex(rho1,im)
      xsi1 = n1*envi_get_data(fid=fid1,dims=dims1,pos=5) ;c22
      b1 = envi_get_data(fid=fid1,dims=dims1,pos=6)   ;c23
      im = envi_get_data(fid=fid1,dims=dims1,pos=7)
      b1 = n1*complex(b1,im)
      zeta1 = n1*envi_get_data(fid=fid1,dims=dims1,pos=8);c33
      det1 = k1*xsi1*zeta1 + 2*real_part(a1*b1*conj(rho1)) - xsi1*(abs(rho1)^2) - k1*(abs(b1)^2) - zeta1*(abs(a1)^2)
      span1 = k1 + 2*xsi1 + zeta1
      k2 = n2*envi_get_data(fid=fid2,dims=dims2,pos=0)  ;c11
      a2 = envi_get_data(fid=fid2,dims=dims2,pos=1)  ;c12
      im = envi_get_data(fid=fid2,dims=dims2,pos=2)
      a2 = n2*complex(a2,im)
      rho2 = envi_get_data(fid=fid2,dims=dims2,pos=3) ;c13
      im = envi_get_data(fid=fid2,dims=dims2,pos=4)
      rho2 = n2*complex(rho2,im)
      xsi2 = n2*envi_get_data(fid=fid2,dims=dims2,pos=5) ;c22
      b2 = envi_get_data(fid=fid2,dims=dims2,pos=6)   ;c23
      im = envi_get_data(fid=fid2,dims=dims2,pos=7)
      b2 = n2*complex(b2,im)
      zeta2 = n2*envi_get_data(fid=fid2,dims=dims2,pos=8);c33
      det2 = k2*xsi2*zeta2 + 2*real_part(a2*b2*conj(rho2)) - xsi2*(abs(rho2)^2) - k2*(abs(b2)^2) - zeta2*(abs(a2)^2)
      k3    = k1 + k2
      a3    = a1 + a2
      rho3  = rho1 +  rho2
      xsi3  = xsi1 +  xsi2
      b3    = b1 +    b2
      zeta3 = zeta1 + zeta2
      det3 = k3*xsi3*zeta3 + 2*real_part(a3*b3*conj(rho3)) - xsi3*(abs(rho3)^2) - k3*(abs(b3)^2) - zeta3*(abs(a3)^2)
      p = 3
      f = p^2
      cst = p*((n1+n2)*alog(n1+n2)-n1*alog(n1)-n2*alog(n2))
      rho = 1. - (2.*p^2-1.)*(1./n1 + 1./n2 - 1./(n1+n2))/(6.*p)
      omega2 = -(p*p/4.)*(1. - 1./rho)^2 + p^2*(p^2-1.)*(1./n1^2 + 1./n2^2 - 1./(n1+n2)^2)/(24.*rho^2)
    end
    4: begin
      k1 = n1*envi_get_data(fid=fid1,dims=dims1,pos=0) ;c11
      a1 = envi_get_data(fid=fid1,dims=dims1,pos=1) ;c12
      im = envi_get_data(fid=fid1,dims=dims1,pos=2)
      a1 = n1*complex(a1,im)
      xsi1 = n1*envi_get_data(fid=fid1,dims=dims1,pos=3);c22
      det1 = k1*xsi1 - abs(a1)^2
      span1 = k1 + xsi1
      k2 = n2*envi_get_data(fid=fid2,dims=dims2,pos=0) ;c11
      a2 = envi_get_data(fid=fid2,dims=dims2,pos=1) ;c12
      im = envi_get_data(fid=fid2,dims=dims2,pos=2)
      a2 = n2*complex(a2,im)
      xsi2 = n2*envi_get_data(fid=fid2,dims=dims2,pos=3);c22
      det2 = k2*xsi2 - abs(a2)^2
      k3    = k1 + k2
      a3    = a1 + a2
      xsi3  = xsi1 +  xsi2
      det3 = k3*xsi3 - abs(a3)^2
      p = 2
      cst = p*((n1+n2)*alog(n1+n2)-n1*alog(n1)-n2*alog(n2))
      f = p^2;
      rho = 1-(2*f-1)*(1./n1+1./n2-1./(n1+n2))/(6.*p);
      omega2 = -f/4.*(1-1./rho)^2 + f*(f-1)*(1./n1^2+1./n2^2-1./(n1+n2)^2)/(24.*rho^2);
    end
    1: begin
      k1 = n1*envi_get_data(fid=fid1,dims=dims1,pos=0) ;c11
      span1 = k1
      det1 = k1
      k2 = n2*envi_get_data(fid=fid2,dims=dims2,pos=0) ;c11
      det2 = k2
      k3 = k1 + k2
      det3 = k3
      p = 1
      cst = p*((n1+n2)*alog(n1+n2)-n1*alog(n1)-n2*alog(n2))
      f = p^2;
      rho = 1-(2.*f-1)*(1./n1+1./n2-1./(n1+n2))/(6.*p);
      omega2 = -f/4.*(1-1./rho)^2+f*(f-1)*(1./n1^2+1./n2^2-1./(n1+n2)^2)/(24.*rho^2);
    end
    else: begin
      message, 'Incorrect input data, aborting...'
      return, 0
    end
  endcase
  envi_file_mng, id = fid2, /remove
  ; check fo bad data
  idx = where(det1 le 0,count)
  if count gt 0 then begin
    print, 'Warning: det(image1) has non-positive values'
    det1[idx] = 0.0001;(machar()).eps
  endif
  idx = where(det2 le 0,count)
  if count gt 0 then begin
    print, 'Warning: det(image2) has non-positive values'
    det2[idx] = 0.0001;(machar()).eps
  endif
  idx = where(det3 le 0,count)
  if count gt 0 then begin
    print, 'Warning: det(image1+image2) has non-positive values'
    det3[idx] = 0.0001;(machar()).eps
  endif
  ; change statistic
  lnQ = cst+n1*alog(det1)+n2*alog(det2)-(n1+n2)*alog(det3)
  Z = -2*rho*lnQ
  ; change probability
  CP =  (1.-omega2)*chisqr_pdf(Z,f)+omega2*chisqr_pdf(Z,f+4)
  
  return, 1

END