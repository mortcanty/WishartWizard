; docformat = 'rst'
; wishart_wizard_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.

PRO wishart_wizard_run_define_buttons, buttonstate
   ENVI_DEFINE_MENU_BUTTON, buttonstate, $
      VALUE = 'Wishart Wizard', $
      REF_VALUE = 'Polarimetric Tools', $
      EVENT_PRO = 'wishart_wizard_run', $
      UVALUE = 'WISHART_WIZARD',$
      POSITION = 'after'
END

pro reset_sens, status
  widget_control, status.ingest2ID, sensitive = 0 
  widget_control, status.saveviewID, sensitive = 0
  widget_control, status.save1ID, sensitive = 0
  widget_control, status.save2ID, sensitive = 0
  widget_control, status.coregisterID, sensitive = 0
  widget_control, status.estenlID, sensitive = 0
  widget_control, status.viewspan1ID, sensitive = 0
  widget_control, status.viewspan2ID, sensitive = 0
  widget_control, status.viewenl1ID, sensitive = 0
  widget_control, status.viewenl2ID, sensitive = 0
  widget_control, status.toggleID, sensitive = 0
  widget_control, status.changemapID, sensitive = 0
end  

function checkforsarscape
  common osb, osb,workdir,dem_file,rg_looks
  COMPILE_OPT IDL2
; Standard error handling.
  Catch, theError
  IF theError NE 0 THEN BEGIN
    Catch, /CANCEL
    RETURN, 0
  ENDIF
  osb = SARscapeBatch()
  rg_looks=3L
  workdir = 'none'
  return, 1
end

; ------------------
; GUI event handlers
; ------------------

pro wishart_cleanup, tlb
  widget_control, tlb, get_Uvalue=status,/no_copy
  if n_elements(status) eq 0 then return
  ptr_free,status.map_info
  ptr_free,status.dims1
  ptr_free,status.bnames
  ptr_free,status.image1
  ptr_free,status.image2
  ptr_free,status.enlim1
  ptr_free,status.enlim2
  print, 'Done ---------------------------'
end

pro wishart_wizard_event, event
   COMPILE_OPT IDL2
   widget_control, event.top, get_Uvalue=status
   if (TAG_NAMES(event, /STRUCTURE_NAME) eq 'WIDGET_DRAW') and (status.current_view gt 0) then begin
      X = event.x
      X = X>0
      X = X<status.cols1
      Y = status.rows1-event.y-1
      Y = Y>0
      Y = Y<status.rows1
      case status.current_view of
        1: begin
             data = (*status.image1)[X,Y]
             widget_control, status.xytextID, set_value='('+strtrim(X,2)+', '+strtrim(Y,2)+'): '+strtrim(data,2)
           end  
        2: begin
             data = (*status.image2)[X,Y]
             widget_control, status.xytextID, set_value='('+strtrim(X,2)+', '+strtrim(Y,2)+'): '+strtrim(data,2)
           end 
        3: begin            
             data = (*status.enlim1)[X,Y]
             widget_control, status.xytextID, set_value='('+strtrim(X,2)+', '+strtrim(Y,2)+'): <'+strtrim(data,2)+'>
           end   
        4: begin
             data = (*status.enlim2)[X,Y]
             widget_control, status.xytextID, set_value='('+strtrim(X,2)+', '+strtrim(Y,2)+'): <'+strtrim(data,2)+'>
           end
        5: begin
             data = (*status.Z)[X,Y]
             widget_control, status.xytextID, set_value='('+strtrim(X,2)+', '+strtrim(Y,2)+'): '+strtrim(data,2)
           end          
        6: begin
             data = (*status.CP)[X,Y]
             widget_control, status.xytextID, set_value='('+strtrim(X,2)+', '+strtrim(Y,2)+'): '+strtrim(data,2)
           end  
        7: begin
             data = (*status.CP)[X,Y]
             widget_control, status.xytextID, set_value='('+strtrim(X,2)+', '+strtrim(Y,2)+'): '+strtrim(data,2)
           end              
        else: return  
      endcase
      
   endif   
end

; File menu handlers SARscape

pro Onsetwork, event
  common osb, osb,workdir,dem_file,rg_looks
  widget_control, event.top, get_Uvalue=status, /no_copy  
  COMPILE_OPT IDL2  
; Standard error handling.
  Catch, theError
  IF theError NE 0 THEN BEGIN
    Catch, /CANCEL
    void = Error_Message()
    widget_control,event.top,set_Uvalue=status,/no_copy
    RETURN
  ENDIF
  workdir = dialog_pickfile(/directory)
  if workdir ne '' then begin
    print, 'Working directory: '+workdir
    widget_control, status.statetextID,set_value=workdir
    ok = SARscape_set_working_in_actual_default(workdir)
    widget_control, status.slcID,sensitive=1
    widget_control, status.clearworkID,sensitive=1
    cd, workdir
  endif  
  widget_control,event.top,set_Uvalue=status,/no_copy
end

pro Onclearwork, event
  common osb, osb,workdir,dem_file,rg_looks 
  widget_control, event.top, get_Uvalue=status, /no_copy  
  COMPILE_OPT IDL2  
; Standard error handling.
  Catch, theError
  IF theError NE 0 THEN BEGIN
    Catch, /CANCEL
    void = Error_Message()
    widget_control,event.top,set_Uvalue=status,/no_copy
    RETURN
  ENDIF
  cd, workdir
  files = file_search()
  if n_elements(files) gt 1 then foreach file, files do file_delete, file 
  print, workdir+' was cleared'
   widget_control, status.statetextID,set_value='Cleared: '+workdir
  widget_control,event.top,set_Uvalue=status,/no_copy
end

pro Onshowwork, event
  common osb, osb,workdir,dem_file,rg_looks 
  widget_control, event.top, get_Uvalue=status, /no_copy
  COMPILE_OPT IDL2  
; Standard error handling.
  Catch, theError
  IF theError NE 0 THEN BEGIN
    Catch, /CANCEL
    void = Error_Message()
    widget_control,event.top,set_Uvalue=status,/no_copy
    RETURN
  ENDIF
  widget_control, status.statetextID,set_value=workdir
  widget_control,event.top,set_Uvalue=status,/no_copy
end

pro ONradarsat2, event
  common osb, osb,workdir,dem_file,rg_looks
  COMPILE_OPT IDL2  
; Standard error handling.
  Catch, theError
  IF theError NE 0 THEN BEGIN
    Catch, /CANCEL
    void = Error_Message()
    RETURN
  ENDIF  
  ok = osb.setupmodule(module='IMPORTRADARSAT2FORMAT')
  ok = osb.setparam('parameter_extraction_data_type','radarsat2_slc')
;  ok = osb.setparam('parameter_extraction_data_version','default') 
  infile = dialog_pickfile(/read,filter='product.xml')
  if infile eq '' then return
  ok = osb.setparam('input_file_list',infile) 
  ok = osb.setparam('output_file_list',workdir+'\product') 
  print, 'Importing RadarSat-2 quadpol imagery ...'
  ok = osb.executeprogress(ErrMsg=err)
  if ok then print, 'Import radarsat-2 output files in '+workdir+'product'+'_...' $
  else begin
    print, 'Import Radarsat-2 failed:'+err
    return
  endelse    
  ok = osb.setupmodule(module='POLPOLARIMETRICFEATURES')
  infiles = workdir+'\product_'+['HH_slc','HV_slc','VH_slc','VV_slc'] 
  ok = osb.setparam('in_hh_file_name',infiles[0])
  ok = osb.setparam('in_hv_file_name',infiles[1])
  ok = osb.setparam('in_vh_file_name',infiles[2])
  ok = osb.setparam('in_vv_file_name',infiles[3])
  out_root_name = workdir+'product'
  ok = osb.setparam('out_root_name',out_root_name)
  ok = osb.setparam('span','NotOK')
  ok = osb.setparam('polrat','NotOK')
  ok = osb.setparam('lindepolrat','NotOK')
  ok = osb.setparam('ppd_hh_vv','NotOK')
  ok = osb.setparam('coherence_hh_vv','NotOK')
  ok = osb.setparam('ppd_hh_hv','NotOK')
  ok = osb.setparam('coherence_hh_hv','NotOK')
  ok = osb.setparam('ppd_hv_vv','NotOK')
  ok = osb.setparam('coherence_hv_vv','NotOK')
  ok = osb.setparam('norm_flag','NotOK')
; get the multilooking factors
  envi_open_file, infiles[0], r_fid=fid, /invisible
  envi_file_query, fid, pixel_size=ps
  envi_file_mng,id=fid,/remove 
  az_looks = round(rg_looks*ps[0]/ps[1])
  grid = [rg_looks*ps[0], az_looks*ps[1]]
  ok = osb.setparam('az_looks',az_looks)
  ok = osb.setparam('rg_looks',rg_looks)
  print, 'Calculating polarimetric features (covariance matrix elements)...'
  print, 'az looks: '+strtrim(az_looks,2)
  print, 'rg looks: '+strtrim(rg_looks,2)
  print, 'grid size: '+strtrim(grid[0],2)+', '+strtrim(grid[1],2)
  ok = osb.executeprogress(ErrMsg=err)  
  if ok then print, 'Covariance matrix element files written to '+out_root_name+'_...' $
  else print, 'Polarimetric features Radarsat-2 failed:'+err
end  

pro ONtsx, event
  common osb, osb,workdir,dem_file,rg_looks
  COMPILE_OPT IDL2
  ; Standard error handling.
  Catch, theError
  IF theError NE 0 THEN BEGIN
    Catch, /CANCEL
    void = Error_Message()
    RETURN
  ENDIF
  ok = osb.setupmodule(module='IMPORTTSXFORMAT')
    infile = dialog_pickfile(/read,filter='*.xml')
  if infile eq '' then return
  outfile = workdir+'\product'
  ok = osb.setparam('input_file_list',infile)
  ok = osb.setparam('output_file_list',outfile)
  ok = osb.setparam('parameter_extraction_data_type','TSX1_SSC_SM_QUAD')
;  ok = osb.setparam('parameter_extraction_data_version','default')
  ok = osb.executeprogress(ErrMsg=err)
  if ok then print, 'Import TSX output file: '+outfile $
    else begin
      print, 'Import TSX failed '+err 
      return
    endelse 
  ok = osb.setupmodule(module='POLPOLARIMETRICFEATURES')
  infiles = workdir+'\product_'+['HH_slc','HV_slc','VH_slc','VV_slc']
  ok = osb.setparam('in_hh_file_name',infiles[0])
  ok = osb.setparam('in_hv_file_name',infiles[1])
  ok = osb.setparam('in_vh_file_name',infiles[2])
  ok = osb.setparam('in_vv_file_name',infiles[3])
  out_root_name = workdir+'product'
  ok = osb.setparam('out_root_name',out_root_name)
  ok = osb.setparam('span','NotOK')
  ok = osb.setparam('polrat','NotOK')
  ok = osb.setparam('lindepolrat','NotOK')
  ok = osb.setparam('ppd_hh_vv','NotOK')
  ok = osb.setparam('coherence_hh_vv','NotOK')
  ok = osb.setparam('ppd_hh_hv','NotOK')
  ok = osb.setparam('coherence_hh_hv','NotOK')
  ok = osb.setparam('ppd_hv_vv','NotOK')
  ok = osb.setparam('coherence_hv_vv','NotOK')
  ok = osb.setparam('norm_flag','NotOK')
  ; get the multilooking factors
  envi_open_file, infiles[0], r_fid=fid, /invisible
  envi_file_query, fid, pixel_size=ps
  envi_file_mng,id=fid,/remove
  az_looks = round(rg_looks*ps[0]/ps[1])
  grid = [rg_looks*ps[0], az_looks*ps[1]]
  ok = osb.setparam('az_looks',az_looks)
  ok = osb.setparam('rg_looks',rg_looks)
  print, 'Calculating polarimetric features (covariance matrix elements)...'
  print, 'az looks: '+strtrim(az_looks,2)
  print, 'rg looks: '+strtrim(rg_looks,2)
  print, 'grid size: '+strtrim(grid[0],2)+', '+strtrim(grid[1],2)
  ok = osb.executeprogress(ErrMsg=err)
  if ok then print, 'Covariance matrix element files written to '+out_root_name+'_...' $
  else print, 'Polarimetric features TSX failed:'+err
end

pro ONtsxsl, event
  common osb, osb,workdir,dem_file,rg_looks
  ; Standard error handling.
  Catch, theError
  IF theError NE 0 THEN BEGIN
    Catch, /CANCEL
    void = Error_Message()
    RETURN
  ENDIF
  ok = osb.setupmodule(module='IMPORTTSXFORMAT')
  infile = dialog_pickfile(/read,filter='*.xml')
  if infile eq '' then return
  ok = osb.setparam('input_file_list',infile)
  outfile = workdir+file_basename(infile,'.xml')
  ok = osb.setparam('output_file_list',outfile)
  ok = osb.setparam('parameter_extraction_data_type','TSX1_SSC_SL_SINGLE')
  ok = osb.executeprogress(ErrMsg=err)
  if ok then print, 'Import TSX output file: '+outfile+'_HH_slc' $
    else begin
      print, 'Import TSX failed '+err 
      return
    endelse     
  ok = osb.setupmodule(module='BASEMULTILOOKING')   
  infile = outfile+'_HH_slc'
  outfile = outfile+'_HH_pwr'
; get the multilooking factors
  envi_open_file, infile, r_fid=fid, /invisible
  envi_file_query, fid, pixel_size=ps
  envi_file_mng,id=fid,/remove
  az_looks = round(rg_looks*ps[0]/ps[1])
  grid = [rg_looks*ps[0], az_looks*ps[1]]
  print, 'az looks: '+strtrim(az_looks,2)
  print, 'rg looks: '+strtrim(rg_looks,2)
  print, 'grid size: '+strtrim(grid[0],2)+', '+strtrim(grid[1],2)
  ok = osb.setparam('input_file_list',infile)
  ok = osb.setparam('output_file_list',outfile)
  ok = osb.setparam('azimuth_multilook',az_looks)
  ok = osb.setparam('range_multilook',rg_looks) 
  ok = osb.executeprogress(ErrMsg=err)
  if ok then print, 'Single polarization file written to '+outfile $
  else print, 'Multilooking TSX failed: '+err
end

pro ONcsk, event
  common osb, osb,workdir,dem_file,rg_looks
  ; Standard error handling.
  Catch, theError
  IF theError NE 0 THEN BEGIN
    Catch, /CANCEL
    void = Error_Message()
    RETURN
  ENDIF
  ok = osb.setupmodule(module='IMPORTCSKFORMAT')
  infile = dialog_pickfile(/read,filter='*.h5')
  if infile eq '' then return
  ok = osb.setparam('input_file_list',infile)
  outfile = workdir+file_basename(infile,'.h5')
  ok = osb.setparam('output_file_list',outfile)
  ok = osb.executeprogress(ErrMsg=err)
  infile = outfile+'_vv_slc'
  outfile = outfile+'_vv_pwr'
  if ok then print, 'Import CMS output file: '+infile $
  else begin
    print, 'Import CSK failed'
    return
  endelse
  ok = osb.setupmodule(module='BASEMULTILOOKING') 
  envi_open_file, infile, r_fid=fid, /invisible
  envi_file_query, fid, pixel_size=ps
  envi_file_mng,id=fid,/remove
  az_looks = round(rg_looks*ps[0]/ps[1])
  grid = [rg_looks*ps[0], az_looks*ps[1]]
  print, 'az looks: '+strtrim(az_looks,2)
  print, 'rg looks: '+strtrim(rg_looks,2)
  print, 'grid size: '+strtrim(grid[0],2)+', '+strtrim(grid[1],2)
  ok = osb.setparam('input_file_list',infile)
  ok = osb.setparam('output_file_list',outfile)
  ok = osb.setparam('azimuth_multilook',az_looks)
  ok = osb.setparam('range_multilook',rg_looks)
  ok = osb.executeprogress(ErrMsg=err)
  if ok then print, 'Single polarization file written to '+outfile $
  else print, 'Multilooking CSK failed: '+err
end  

pro ONdem, event
  common osb, osb,workdir,dem_file,rg_looks
  COMPILE_OPT IDL2
  ; Standard error handling.
  Catch, theError
  IF theError NE 0 THEN BEGIN
    Catch, /CANCEL
    void = Error_Message()
    RETURN
  ENDIF
  dem_file = dialog_pickfile(/read,title='Enter a digital elevation model file',filter='*_dem')
  if dem_file eq '' then begin
    print, 'A DEM file must be selected'
    message, 'A DEM file must be selected'
  endif
  widget_control, event.top, get_Uvalue=status
  print, 'DEM file: '+dem_file
  widget_control, status.statetextID,set_value=dem_file
  widget_control, status.geocodeID, sensitive=1
end

pro ONgeocode, event
  common osb, osb,workdir,dem_file,rg_looks
  COMPILE_OPT IDL2
; Standard error handling.
  Catch, theError
  IF theError NE 0 THEN BEGIN
    Catch, /CANCEL
    void = Error_Message()
    RETURN
  ENDIF
  cd, workdir
  ok = osb.setupmodule(module='BASICGEOCODING')
  input_file_list = dialog_pickfile(title='Select matrix element bands for geocoding', $
                    /read, $
                    /multiple_files, $
                    filter=['*_hh','*_hv','*_vh','*_vv','*_HH_pwr','*_hh_pwr','*_VV_pwr','*_vv_pwr'])               
  output_file_list = input_file_list + '_geo' 
  envi_open_file, input_file_list[0], r_fid=fid, /invisible
  envi_file_query, fid, pixel_size=ps
  envi_file_mng,id=fid,/remove
  geocode_grid_size = round(total(ps)*10.0/2.0)/10.0
  print, 'Geocoding and terrain matching...'
  print, 'geocode grid size: '+strtrim(geocode_grid_size,2)
  ok = osb.setparam('input_file_list',input_file_list)
  ok = osb.setparam('output_file_list',output_file_list)
  ok = osb.setparam('dem_file_name',dem_file)
  ok = osb.setparam('geocode_grid_size_x',geocode_grid_size)
  ok = osb.setparam('geocode_grid_size_y',geocode_grid_size)
  ok = osb.setparam('geocode_resampling_type','bilinear_interpolation')
;  ok = osb.setparam('geocode_resampling_type','3th_order_cc')
  ok = osb.executeprogress(ErrMsg=err)
  if ok then print, 'Geocoded file(s) written to '+workdir $
  else print, 'Geocoding failed:'+err
end  

; File menu handlers Wishart

pro ONingest1, event
   widget_control, event.top, get_Uvalue=status, /no_copy   
   COMPILE_OPT IDL2
  ; Standard error handling.
   Catch, theError
   IF theError NE 0 THEN BEGIN
     Catch, /CANCEL
     void = Error_Message()
     widget_control,event.top,set_Uvalue=status,/no_copy
     RETURN
   ENDIF      
   envi_select, /band_only, dims=dims1, fid=fid1, pos=pos, title='Choose (spatial subset of) any covariance band for 1st image'
   envi_file_query, fid1, fname=fname1
;   fname1 = dialog_pickfile( title='Choose any covariance band for 1st image', $
;     filter=['*_geo','*.tif','*.bin'])
;   envi_open_file, fname1, r_fid=fid1
   if (fid1 eq -1) then begin
     widget_control,event.top,Set_Uvalue=status,/no_copy
     return
   endif  
   dir = file_dirname(fname1) 
   if dir eq '.' then message, 'Cannot ingest from ENVI memory'
   print, 'Ingesting from directory '+dir+' ...'
   widget_control,status.statetextID,set_value='Ingesting ...'
   if polsaringest(dir,dims=dims1,outim=outim,bnames=bnames,map_info=map_info) then begin
     status.bnames = ptr_new(bnames)
     status.dims1 = ptr_new(dims1)
     status.cols1 = dims1[2]-dims1[1]+1
     status.rows1 = dims1[4]-dims1[3]+1
     envi_convert_file_coordinates, fid1, dims1[1], dims1[3], e, n, /to_map
     map_info.mc = [0D,0D,e,n]
     status.map_info = ptr_New(map_info)
     status.image1 = ptr_new(outim)
     envi_file_mng,/remove,id=fid1
     print, 'First image (subset) ingested'
;   reset menu sensitivities
     reset_sens, status
     widget_control, status.ingest2ID, sensitive = 1
     widget_control, status.viewspan1ID, sensitive = 1
     widget_control, status.saveviewID, sensitive = 1
     widget_control, status.save1ID, sensitive = 1
     widget_control, status.load2ID, sensitive = 1
   end else begin
     reset_sens, status
     envi_file_mng,/remove,id=fid1
     widget_control,status.statetextID,set_value='Image ingest failed' 
     print, 'Image ingest failed' 
     widget_control,event.top,set_Uvalue=status,/no_copy 
     return
   endelse                   
   widget_control,event.top,set_Uvalue=status,/no_copy 
   onViewspan1, event   
end

pro ONingest2, event
   widget_control, event.top, get_Uvalue=status, /no_copy
   
   COMPILE_OPT IDL2
   ; Standard error handling.
   Catch, theError
   IF theError NE 0 THEN BEGIN
     Catch, /CANCEL
     void = Error_Message()
     widget_control,event.top,set_Uvalue=status,/no_copy
     RETURN
   ENDIF  
   fname2 = dialog_pickfile( title='Choose any covariance band for 2nd image', $
     filter=['*_geo','*.tif','*.bin'])
   envi_open_file, fname2, r_fid=fid2
   if (fid2 eq -1) then begin
      widget_control,event.top,Set_Uvalue=status,/no_copy
      return
   endif
   dir = file_dirname(fname2) 
   if dir eq '.' then message, 'Cannot ingest from ENVI memory' 
   envi_file_mng,/remove,id=fid2
   print, 'Ingesting from directory '+dir+' ...'
   widget_control,status.statetextID,set_value='Ingesting ...'
   if polsaringest(dir,dims=dims2,outim=outim,bnames=bnames,map_info=map_info2) then begin
     bands = n_elements(*status.bnames)
     if bands ne n_elements(bnames) then begin
       print, 'Band mismatch'
       widget_control,status.statetextID,set_value='Image ingest failed'
       return
     endif
;   setup for matching using four corners only
     envi_enter_data, (*status.image1)[*,*,0], map_info=*status.map_info, r_fid=fid1     
     envi_enter_data, outim, r_fid=fid2, map_info=map_info2
     fid2a = fid2 
     if imagematch(fid1, fid2, dims2, map_info) then begin
       img = fltarr(status.cols1,status.rows1,bands)
       for k=0,bands-1 do img[*,*,k] = envi_get_data(dims=dims2,fid=fid2,pos=k)
       status.image2 = ptr_new(img)
       status.dims1 = ptr_new(dims2)
     end else begin
        widget_control,status.statetextID,set_value='Image match failed' 
        print, 'Image match failed' 
        widget_control,event.top,set_Uvalue=status,/no_copy 
        return
     endelse   
     print, 'Second image ingested'
     envi_file_mng,/remove,id=fid1
     envi_file_mng,/remove,id=fid2
     envi_file_mng,/remove,id=fid2a 
     widget_control, status.viewspan2ID, sensitive = 1
     widget_control, status.toggleID, sensitive = 1
     widget_control, status.save2ID, sensitive = 1
     widget_control, status.estenlID, sensitive = 1
     widget_control, status.coregisterID, sensitive= 1
     widget_control, status.changemapID, sensitive = 1
     widget_control,event.top,set_Uvalue=status,/no_copy 
     onViewspan2, event
     return
   end else begin
     widget_control,status.statetextID,set_value='Image ingest failed' 
     print, 'Image ingest failed' 
   endelse                   
   widget_control,event.top,set_Uvalue=status,/no_copy  
end

pro ONload1, event
  COMPILE_OPT IDL2
  widget_control, event.top, get_Uvalue=status, /no_copy
  envi_select, title='Choose first polSAR image', $
    fid=fid1, dims=dims1, pos=pos, /no_spec
  if (fid1 eq -1) then begin
    print, 'cancelled'
    widget_control,event.top,set_Uvalue=status,/no_copy 
    return
  endif
  envi_file_query, fid1, fname=fname, bnames=bnames
  map_info = envi_get_map_info(fid=fid1)
  status.cols1 = dims1[2]-dims1[1]+1
  status.rows1 = dims1[4]-dims1[3]+1
  status.bnames = ptr_new(bnames)
  envi_convert_file_coordinates, fid1, dims1[1], dims1[3], e, n, /to_map
  map_info.mc = [0D,0D,e,n]
  status.map_info = ptr_New(map_info) 
  bands = n_elements(pos)
  image = fltarr(status.cols1,status.rows1,bands)
  for i=0,bands-1 do $
    image[*,*,i] = envi_get_data(fid=fid1,dims=dims1,pos=pos[i])
  idx = where(finite(image,/NAN),count)
  if count gt 0 then image[idx]=0.0  
  status.image1 = ptr_new(image)
  status.dims1 = ptr_new([-1,0,status.cols1-1,0,status.rows1-1])
  print, 'input file '+fname
  widget_control, status.ingest2ID, sensitive = 1
  widget_control, status.viewspan1ID, sensitive = 1
  widget_control, status.saveviewID, sensitive = 1
  widget_control, status.save1ID, sensitive = 1
  widget_control, status.load2ID, sensitive = 1
  widget_control,event.top,set_Uvalue=status,/no_copy 
  onViewspan1, event
end  
  
pro ONload2, event
  COMPILE_OPT IDL2
  ; Standard error handling.
  Catch, theError
  IF theError NE 0 THEN BEGIN
    Catch, /CANCEL
    void = Error_Message()
    widget_control,event.top,set_Uvalue=status,/no_copy
    RETURN
  ENDIF  
  widget_control, event.top, get_Uvalue=status, /no_copy
  envi_select, title='Choose second polSAR image', $
    fid=fid, dims=dims,pos=pos, /no_spec
  if (fid eq -1) then begin
    print, 'cancelled'
    widget_control,event.top,set_Uvalue=status,/no_copy 
    return
  endif
  envi_file_query, fid, fname=fname, bnames=bnames, map_info=map_info
  cols = dims[2]-dims[1]+1
  rows = dims[4]-dims[3]+1
  bands = n_elements(pos)
  if (cols ne status.cols1) or (rows ne status.rows1) then begin
    void = dialog_message('size mismatch',/error)
    message, 'size mismatch',/continue
    widget_control,event.top,set_Uvalue=status,/no_copy 
    return
  endif
  image = fltarr(status.cols1,status.rows1,bands)
  for i=0,bands-1 do $
    image[*,*,i] = envi_get_data(fid=fid,dims=dims,pos=pos[i])
  idx = where(finite(image,/NAN),count)
  if count gt 0 then image[idx]=0.0  
  status.image2 = ptr_new(image)
  print, 'input file '+fname
  widget_control, status.viewspan2ID, sensitive = 1
  widget_control, status.toggleID, sensitive = 1
  widget_control, status.save2ID, sensitive = 1
  widget_control, status.estenlID, sensitive = 1
  widget_control, status.coregisterID, sensitive= 1
  widget_control, status.changemapID, sensitive = 1
  widget_control,event.top,set_Uvalue=status,/no_copy 
  onViewspan2, event
end   

pro ONsave1, event
  COMPILE_OPT IDL2
  widget_control,event.top,get_Uvalue=status,/no_copy 
; output destination
  base = widget_auto_base(title='PolSAR Output, first image')
  sb = widget_base(base, /row, /frame)
  wp = widget_outfm(sb, uvalue='outf', /auto)
  result1 = auto_wid_mng(base)
  if (result1.accept eq 0) then begin
    print, 'Output cancelled'
    widget_control,event.top,set_Uvalue=status,/no_copy
    return
  end
  if (result1.outf.in_memory eq 1) then begin
    envi_enter_data, *status.image1+0, descrip = 'First image', $
    bnames = *status.bnames, $
    map_info= *status.map_info
  end else begin
    openw, unit, result1.outf.name, /get_lun
    writeu, unit, *status.image1
    envi_setup_head,fname=result1.outf.name, ns=status.cols1, nl=status.rows1, nb=n_elements(*status.bnames), $
      data_type=4, $
      interleave=0, $
      file_type=0, $
      bnames=*status.bnames, $
      map_info=*status.map_info, $
      descrip='polSAR image1: '+file_basename(result1.outf.name), /write,/open
    print, 'File created ', result1.outf.name
    free_lun, unit
  endelse
  widget_control,event.top,set_Uvalue=status,/no_copy
end

pro ONsave2, event
  COMPILE_OPT IDL2
  widget_control,event.top,get_Uvalue=status,/no_copy 
; output destination
  base = widget_auto_base(title='PolSAR Output, second image')
  sb = widget_base(base, /row, /frame)
  wp = widget_outfm(sb, uvalue='outf', /auto)
  result1 = auto_wid_mng(base)
  if (result1.accept eq 0) then begin
    print, 'Output cancelled'
    widget_control,event.top,set_Uvalue=status,/no_copy
    return
  end
  if (result1.outf.in_memory eq 1) then begin
    envi_enter_data, *status.image2+0, descrip = 'Second image', $
    bnames = *status.bnames, $
    map_info= *status.map_info
  end else begin
    openw, unit, result1.outf.name, /get_lun
    writeu, unit, *status.image2
    envi_setup_head,fname=result1.outf.name, ns=status.cols1, nl=status.rows1, nb=n_elements(*status.bnames), $
      data_type=4, $
      interleave=0, $
      file_type=0, $
      bnames=*status.bnames,$
      map_info=*status.map_info, $
      descrip='polSAR image2: '+file_basename(result1.outf.name), /write,/open
    print, 'File created ', result1.outf.name
    free_lun, unit
  endelse
  widget_control,event.top,set_Uvalue=status,/no_copy
end

pro ONsaveview, event
  COMPILE_OPT IDL2
  widget_control,event.top,get_Uvalue=status,/no_copy
  case status.current_view of
    1: begin
      span = total( (*status.image1)[*,*,[0,5,8]], 3 )
      envi_enter_data, span, descrip = 'span of first image', $
        bnames = ['span of first image'], $
        map_info= *status.map_info
    end
    2: begin
      span = total( (*status.image2)[*,*,[0,5,8]], 3 )
      envi_enter_data, span, descrip = 'span of second image', $
        bnames = ['span of second image'], $
        map_info= *status.map_info
    end
    3: begin
      enl = (*status.enlim1)
      envi_enter_data, enl, descrip = 'ENL of first image', $
        bnames = ['ENL of first image'], $
        map_info= *status.map_info
    end
    4: begin
      enl = (*status.enlim2)
      envi_enter_data, enl, descrip = 'ENL of second image', $
        bnames = ['ENL of second image'], $
        map_info= *status.map_info
    end
    5: begin
      Z = (*status.Z)
      envi_enter_data, Z, descrip = 'Change statistic', $
        bnames = ['Change statistic'], $
        map_info= *status.map_info
    end
    6: begin
      CP = (*status.CP)
      envi_enter_data, CP, descrip = 'Change probability', $
        bnames = ['Change probability'], $
        map_info= *status.map_info
    end
    7: begin
      CP = (*status.CP)
      idx = where(CP gt (1.0-status.signif),count)
      img = bytarr(status.cols1,status.rows1)
      if count gt 0 then img[idx] = 1
      envi_enter_data, img, descrip = 'Change map', $
        file_type = 3, $
        num_classes = 2, $
        lookup = [[0,0,0],[255,0,0]], $
        class_names=['no change','change'], $
        bnames = ['change map'], $
        map_info= *status.map_info
    end
  endcase 
  widget_control,event.top,set_Uvalue=status,/no_copy
end

pro ONquit, event
  widget_control,event.top,get_Uvalue=status,/no_copy
  ptr_free,status.map_info
  ptr_free,status.dims1
  ptr_free,status.bnames
  ptr_free,status.image1
  ptr_free,status.image2
  ptr_free,status.enlim1
  ptr_free,status.enlim2
  widget_control, event.top, /destroy
  print, 'Done ---------------------------'
end

; Looks menu handlers

pro ONrglooks, event
  common osb, osb,workdir,dem_file,rg_looks
  COMPILE_OPT IDL2
  base = widget_auto_base(title='Number of range looks')
  we = widget_param(base, dt=1, field=5, floor=1, $
    default=rg_looks, uvalue='param', xsize=32, /auto)
  result = auto_wid_mng(base)
  if not (result.accept eq 0) then rg_looks=result.param
end

pro ONsetenl1, event
  COMPILE_OPT IDL2
  widget_control,event.top,get_Uvalue=status,/no_copy
  
  base = widget_auto_base(title='Equivalent number of looks, first image')
  we = widget_param(base, dt=4, field=5, floor=1, $
    default=status.enl1, uvalue='param', xsize=32, /auto)
  result = auto_wid_mng(base)
  if (result.accept eq 0) then status.enl1=float(8) else status.enl1=float(result.param)
  widget_control, status.enl1textID, set_value='ENL1:  '+strtrim(status.enl1,2)
  
  widget_control,event.top,set_Uvalue=status,/no_copy  
end

pro ONsetenl2, event
  COMPILE_OPT IDL2
  widget_control,event.top,get_Uvalue=status,/no_copy
  
  base = widget_auto_base(title='Equivalent number of looks, second image')
  we = widget_param(base, dt=4, field=5, floor=1, $
    default=status.enl2, uvalue='param', xsize=32, /auto)
  result = auto_wid_mng(base)
  if (result.accept eq 0) then status.enl2=float(8) else status.enl2=float(result.param)
  widget_control, status.enl2textID, set_value='ENL2:  '+strtrim(status.enl2,2)
  
  widget_control,event.top,set_Uvalue=status,/no_copy  
end
  
; View menu handlers  

pro ONviewspan1, event
  COMPILE_OPT IDL2
  widget_control,event.top,get_Uvalue=status,/no_copy
  
  widget_control, status.drawID, draw_xsize=status.cols1, draw_ysize=status.rows1
  wset,status.wID
  if n_elements(*status.bnames) eq 9 then begin
    span = total( (*status.image1)[*,*,[0,5,8]], 3 )
  end else if n_elements(*status.bnames) eq 4 then begin
    span = total( (*status.image1)[*,*,[0,3]], 3 )
  end else begin
    span = (*status.image1)[*,*,0]
  endelse
  tvscl, alog(span+0.001)
  widget_control,status.statetextID,set_value='Log span of Image1'
  widget_control, status.saveviewID, sensitive = 1
  status.current_view = 1
  
  widget_control,event.top,set_Uvalue=status,/no_copy
end

pro ONviewspan2, event
  COMPILE_OPT IDL2
  widget_control,event.top,get_Uvalue=status,/no_copy
  
  widget_control, status.drawID, draw_xsize=status.cols1, draw_ysize=status.rows1
  wset,status.wID
  if n_elements(*status.bnames) eq 9 then begin
    img = (*status.image2)[*,*,0] + (*status.image2)[*,*,5] + (*status.image2)[*,*,8]
  end else if n_elements(*status.bnames) eq 4 then begin
    img = (*status.image2)[*,*,0] + (*status.image2)[*,*,3]
  end else begin
    img = (*status.image2)[*,*,0]
  endelse
  tvscl, alog(img+0.001)
  widget_control,status.statetextID,set_value='Log span of Image2'
  widget_control, status.saveviewID, sensitive = 1
  status.current_view = 2
  
  widget_control,event.top,set_Uvalue=status,/no_copy
end

pro ONviewenl1, event
  COMPILE_OPT IDL2
  widget_control,event.top,get_Uvalue=status,/no_copy
  
  widget_control, status.drawID, draw_xsize=status.cols1, draw_ysize=status.rows1
  wset,status.wID
  img = (*status.enlim1)[*,*,0]
  tvscl, img
  widget_control,status.statetextID,set_value='ENL of Image1'
  widget_control, status.saveviewID, sensitive = 1
  status.current_view = 3
  
  widget_control,event.top,set_Uvalue=status,/no_copy
end

pro ONviewenl2, event
  COMPILE_OPT IDL2
  widget_control,event.top,get_Uvalue=status,/no_copy
  
  widget_control, status.drawID, draw_xsize=status.cols1, draw_ysize=status.rows1
  wset,status.wID
  img = (*status.enlim2)[*,*,0]
  tvscl, img
  widget_control,status.statetextID,set_value='ENL of Image2'
  widget_control, status.saveviewID, sensitive = 1
  status.current_view = 4
  
  widget_control,event.top,set_Uvalue=status,/no_copy
end

pro ONviewstat, event
  COMPILE_OPT IDL2
  widget_control,event.top,get_Uvalue=status,/no_copy
  
  widget_control, status.drawID, draw_xsize=status.cols1, draw_ysize=status.rows1
  wset,status.wID
  img = (*status.Z)[*,*,0]
  tvscl, img>0.0
  widget_control,status.statetextID,set_value='Change statistic'
  widget_control, status.saveviewID, sensitive = 1
  status.current_view = 5
  
  widget_control,event.top,set_Uvalue=status,/no_copy
end

pro ONviewprob, event
  COMPILE_OPT IDL2
  widget_control,event.top,get_Uvalue=status,/no_copy
  
  widget_control, status.drawID, draw_xsize=status.cols1, draw_ysize=status.rows1
  wset,status.wID
  img = (*status.CP)[*,*,0]
  tvscl, img>0.0
  widget_control,status.statetextID,set_value='Change statistic'
  widget_control, status.saveviewID, sensitive = 1
  status.current_view = 6
  
  widget_control,event.top,set_Uvalue=status,/no_copy
end

pro ONviewchange, event
  COMPILE_OPT IDL2
  widget_control,event.top,get_Uvalue=status,/no_copy
  
  widget_control, status.drawID, draw_xsize=status.cols1, draw_ysize=status.rows1
  wset,status.wID
  img = (*status.CP)[*,*,0]
  idx = where(img gt (1.0-status.signif), count)
  img = bytarr(status.cols1*status.rows1,3)
  img1 = bytscl(alog((*status.image1)[*,*,0]+0.001))
  img[*,0] = img1
  img[*,1] = img1
  img[*,2] = img1
  if count gt 0 then begin
    img[idx,0] = 255
    img[idx,1] = 0
    img[idx,1] = 0
  endif
  tvscl, reform(img,status.cols1,status.rows1,3), true=3
  widget_control,status.statetextID,set_value='Change map'
  widget_control, status.saveviewID, sensitive = 1
  status.current_view = 7
  
  widget_control,event.top,set_Uvalue=status,/no_copy
end

; Run menu handlers

pro ONcoregister, event
  COMPILE_OPT IDL2
  widget_control,event.top,get_Uvalue=status,/no_copy
  print, 'running ...'
  widget_control,status.statetextID,set_value='Running ...'
; setup to get tie-points
  envi_enter_data, alog((*status.image1)[*,*,0]+0.001), map_info=*status.map_info, r_fid=fid1
  envi_enter_data, alog((*status.image2)[*,*,0]+0.001), r_fid=fid2
  ENVI_DOIT, 'ENVI_AUTO_TIE_POINTS_DOIT', base_fid=fid1, warp_fid=fid2, out_tie_points_array=pts
  if not plausibility(pts,pts) then begin
    widget_control,status.statetextID,set_value='Insufficient control points, match failed'
    print, 'Image co-register failed'
    envi_file_mng,/remove,id=fid1
    envi_file_mng,/remove,id=fid2 
    widget_control,event.top,set_Uvalue=status,/no_copy
    return
  endif
  envi_file_mng,/remove,id=fid1
  envi_file_mng,/remove,id=fid2 
; setup for matching 
  envi_enter_data, (*status.image1)[*,*,0], map_info=*status.map_info, r_fid=fid1
  envi_enter_data, (*status.image2), r_fid=fid2  
  fid2a = fid2
  if imagematch(fid1, fid2, dims2, status.map_info, pts=pts) then begin
     bands = n_elements(*status.bnames)
     img = fltarr(status.cols1,status.rows1,bands)
     for k=0,bands-1 do img[*,*,k] = envi_get_data(dims=dims2,fid=fid2,pos=k)
     status.image2 = ptr_new(img)
  end else begin
     widget_control,status.statetextID,set_value='Image match failed' 
     print, 'Image co-register failed' 
  endelse 
  envi_file_mng,/remove,id=fid1
  envi_file_mng,/remove,id=fid2
  envi_file_mng,/remove,id=fid2a  
  widget_control,status.statetextID,set_value='Co-register completed'
  print, 'Image co-register completed'
  widget_control,event.top,set_Uvalue=status,/no_copy
end

pro ONestenl1, event
  COMPILE_OPT IDL2
  widget_control,event.top,get_Uvalue=status,/no_copy
  base = widget_auto_base(title='Window size')
  wg = widget_slist(base, list = ['3x3','7x7','11x11'], $
    default = 1, uvalue='list', /auto)
  result = auto_wid_mng(base)
  if (result.accept eq 0) then begin
    widget_control,event.top,set_Uvalue=status,/no_copy
    return
  end else case result.list of
       0: wsize = 3
       1: wsize = 7
       2: wsize = 11
  endcase     
  print, 'running ...'
  widget_control,status.statetextID,set_value='Running ENL'  
; setup  
  envi_enter_data, (*status.image1)+0, r_fid=fid
  if enlml(fid, *status.dims1, n_elements(*status.bnames), enl=enl, wsize=wsize) then begin
    envi_file_mng,/remove,id=fid
    status.enlim1=ptr_new(enl)
    widget_control, status.viewenl1ID, sensitive = 1
    print, 'done'
    widget_control,status.statetextID,set_value='ENL done'
    hist = histogram(enl,min=0.0,max=20.0,nbins=200,OMIN=mn)
    hist[0]=0.0
    void = max(hist, mxpos)
    mode = float(mn + mxpos)/10.0
    print, 'mode at: '+strtrim(mode,2)
    envi_plot_data,findgen(200)/10.0,hist[0:*],$
      title='ENL1',$
      base=base
    widget_control, base, group_leader=event.top
  end else begin
    envi_file_mng,/remove,id=fid
    print, 'ENL calculation failed'
    widget_control,status.statetextID,set_value='ENL calculation failed'
  endelse
  widget_control,event.top,set_Uvalue=status,/no_copy
  ONviewenl1, event
end

pro ONestenl2, event
  COMPILE_OPT IDL2
  widget_control,event.top,get_Uvalue=status,/no_copy
  base = widget_auto_base(title='Window size')
  wg = widget_slist(base, list = ['3x3','7x7','11x11'], $
    default = 1, uvalue='list', /auto)
  result = auto_wid_mng(base)
  if (result.accept eq 0) then begin
    widget_control,event.top,set_Uvalue=status,/no_copy
    return
  end else case result.list of
  0: wsize = 3
  1: wsize = 7
  2: wsize = 11
endcase
  print, 'running ...'
  widget_control,status.statetextID,set_value='Running ENL'
; setup  
  envi_enter_data, *status.image2+0, r_fid=fid
  if enlml(fid, *status.dims1, n_elements(*status.bnames), enl=enl, wsize=wsize) then begin
    envi_file_mng,/remove,id=fid
    status.enlim2=ptr_new(enl)
    widget_control, status.viewenl2ID, sensitive = 1
    print, 'done'
    widget_control,status.statetextID,set_value='ENL done'
    hist = histogram(enl,min=0.0,max=20.0,nbins=200,OMIN=mn)
    hist[0]=0.0
    void = max(hist, mxpos)
    mode = float(mn + mxpos)/10.0
    print, 'mode at: '+strtrim(mode,2)
    widget_control, status.enl2textID, set_value='ENL2:  '+strtrim(status.enl2,2)
    envi_plot_data,findgen(200)/10.0,hist[0:*],$
      title='ENL2',$
      base=base
    widget_control, base, group_leader=event.top    
  end else begin
    envi_file_mng,/remove,id=fid
    print, 'ENL calculation failed'
    widget_control,status.statetextID,set_value='ENL calculation failed'
  endelse
  widget_control,event.top,set_Uvalue=status,/no_copy
  ONviewenl2,event
end

pro ONchangemap, event
  COMPILE_OPT IDL2
  widget_control,event.top,get_Uvalue=status,/no_copy
  print, 'running ...'
  widget_control,status.statetextID,set_value='Running ...'
; setup
  envi_enter_data, *status.image1+0, r_fid=fid1
  envi_enter_data, *status.image2+0, r_fid=fid2
  bands = n_elements(*status.bnames)
  if wishartchange(fid1, fid2, *status.dims1, *status.dims1, status.enl1, status.enl2, bands, changestat=Z, changeprob=CP) then begin
     status.Z = ptr_new(Z)
     CP = median(CP,3)
     status.CP = ptr_new(CP) 
     print, 'Change map done'
     widget_control, status.statetextID,set_value='Change map done' 
     widget_control, status.viewstatID, sensitive = 1
     widget_control, status.viewprobID, sensitive = 1
     widget_control, status.viewchangeID, sensitive = 1
     widget_control, status.slider_sigID, sensitive = 1
  end else begin
     print, 'Change map failed'
     widget_control,status.statetextID,set_value='Change map failed' 
  endelse
  envi_file_mng,/remove,id=fid1
  envi_file_mng,/remove,id=fid2
  widget_control,event.top,set_Uvalue=status,/no_copy
  ONviewchange, event
end

; Help menu handlers

pro onHelp, event
  COMPILE_OPT IDL2
  help = file_which('wisharthelp.pdf',/include_current_dir)
  ONLINE_HELP, 'MadView', BOOK=help
end

PRO onAbout, event
  COMPILE_OPT IDL2
  dummy=dialog_message(['Complex Wishart Change Detection Wizard, Version 1.0 ',$
    'M. Canty, 2014'],/state)
END

; Status bar handlers

pro ONtoggle, event
  COMPILE_OPT IDL2
  widget_control,event.top,get_Uvalue=status
  if status.toggle eq 1 then begin
    status.toggle = 0
    status.current_view = 2
    Onviewspan2, event
  end else begin
    status.toggle = 1
    status.current_view = 1
    Onviewspan1, event
  endelse
  widget_control,event.top,set_Uvalue=status,/no_copy
end

pro ONslider_sig, event
  COMPILE_OPT IDL2
  widget_control,event.top,get_Uvalue=status,/no_copy
  status.signif = float(event.value)/100.0
  widget_control,event.top,set_Uvalue=status,/no_copy
  ONviewchange, event
end

;+
; :Description:
;       GUI for calculating change maps for polSAR imagery
;       with the complex Wishart change detection algorithm  
; :Params:
;       event:  in, optional 
;          required if called from ENVI              
; :Uses:
;       polsaringest::
;       enlml::
;       imagematch::
;       wishartchange:: 
;       ENVI::            
;       COYOTE
; :Author:
;       Mort Canty (2014)        
;-
PRO wishart_wizard_run, event

  COMPILE_OPT IDL2
; Standard error handling.
  Catch, theError
  IF theError NE 0 THEN BEGIN
    Catch, /CANCEL
    void = Error_Message()
    RETURN
  ENDIF

  print, '--------------------------------'
  print, 'Complex Wishart Change Detection'
  print, systime(0)
  print, '--------------------------------'
  
; widget creation and registration

  tlb           = widget_base(/column,title=" The Wishart Wizard: Complex Wishart polSAR Change Detection",mbar=menubarID,xoffset=10,yoffset=80)
  
  fileID        = widget_button(menubarID,value='File',/menu)
    ingestID      = widget_button(fileID,value='Ingest',/menu)
      ingest1ID     = widget_button(ingestID, value='First covariance matrix elements', event_pro='ONingest1')
      ingest2ID     = widget_button(ingestID, value='Second covariance matrix elements', event_pro='ONingest2',Sensitive=0)
    saveID        = widget_button(fileID,value='Save',/menu)   
      save1ID       = widget_button(SaveID, value='First covariance image', event_pro='ONsave1', Sensitive=0)
      save2ID       = widget_button(SaveID, value='Second covariance image', event_pro='ONsave2',Sensitive=0)      
    loadID        = widget_button(fileID,value='Reload',/menu)  
      load1ID        = widget_button(loadID, value='First covariance image', event_pro='ONload1')
      load2ID        = widget_button(loadID, value='Second covariance image', event_pro='ONload2',Sensitive=0)
    saveviewID    = widget_button(fileID, value='Save view to ENVI', event_pro='ONsaveview',Sensitive=0)
    quitID        = widget_button(fileID,value='Quit', event_pro='ONquit', /separator)
    
  sarscapeID    = widget_button(menubarID,value='(SARscape)',/menu,Sensitive = checkForSarscape())
    workdirID     = widget_button(sarscapeID,value='Working directory',/menu)
      setworkID     = widget_button(workdirID, value='Set working directory',event_pro='ONsetwork')
      clearworkID   = widget_button(workdirID, value='Clear working directory',event_pro='ONclearwork',sensitive=0)
      showworkID    = widget_button(workdirID, value='Show working directory',event_pro='ONshowwork')
    rglooksID     = widget_button(sarscapeID, value='Set range looks', event_pro='Onrglooks')
    setdemID      = widget_button(sarscapeID, value='Set DEM', event_pro='ONdem',sensitive=1)
    slcID         = widget_button(sarscapeID, value='Import SLC', /menu, sensitive=0)
      radarsat2ID   = widget_button(slcID, value='RADARSAT-2 QuadPol', event_pro='ONradarsat2')
      tsxID         = widget_button(slcID, value='TerraSAR-X QuadPol', event_pro='ONtsx')
      tsxslID       = widget_button(slcID, value='TerraSAR-X Spotlight SinglePol', event_pro='ONtsxsl')
      cskID         = widget_button(slcID, value='COSMO-SkyMed Singlepol', event_pro='ONcsk')
    geocodeID     = widget_button(sarscapeID, value='Geocode', event_pro='ONgeocode',sensitive=0)    
  
    
  
  viewID        = widget_button(menubarID, value = 'View', /menu)
    viewspanID    = widget_button(viewID, value = 'Span', /menu)
      viewspan1ID   = widget_button(viewspanID, value = 'First covariance image', event_pro = 'ONviewspan1', sensitive=0)
      viewspan2ID   = widget_button(viewspanID, value = 'Second covariance image', event_pro = 'ONviewspan2', sensitive=0)
    viewenlID   = widget_button(viewID, value = 'ENL', /menu)
      viewenl1ID    = widget_button(viewenlID, value = 'First covariance image', event_pro = 'ONviewenl1', sensitive=0)
      viewenl2ID    = widget_button(viewenlID, value = 'Second covariance image', event_pro = 'ONviewenl2', sensitive=0)
    viewstatID    = widget_button(viewID, value = 'Change statistic', event_pro = 'ONviewstat', sensitive=0)
    viewprobID    = widget_button(viewID, value = 'Change probability', event_pro = 'ONviewprob', sensitive=0)
    viewchangeID  = widget_button(viewID, value = 'Change map', event_pro = 'ONviewchange', sensitive=0)
  
  runID         = widget_button(menubarID, value='Run', /menu)
    coregisterID  = widget_button(runID, value='Co-register', event_pro = 'ONcoregister', sensitive=0)
    estenlID        = widget_button(runID, value='Calculate ENL', /menu,sensitive=0)
      estenl1ID     = widget_button(estenlID, value='First covariance image', event_pro = 'ONestenl1')
      estenl2ID     = widget_button(estenlID, value='Second covariance image', event_pro = 'ONestenl2')
    setenl1ID     = widget_button(runID, value='Set ENL1', event_pro='ONsetenl1')
    setenl2ID     = widget_button(runID, value='Set ENL2', event_pro='ONsetenl2')  
    changemapID   = widget_button(runID, value='Change detection', event_pro = 'ONchangemap', sensitive=0) 
  
  helpID        = widget_button(menubarID,value='Help',/menu)
    mvHelpID      = widget_button(helpID,value='Open help file',event_pro='ONhelp')
    mvAboutID     = widget_button(helpID,value='About',event_pro='ONabout')
  
  
  widget_control,tlb,get_Uvalue=status 
  statetextID   = widget_label(tlb,value='The Wishart Wizard',/align_left,/dynamic_resize)
  drawID        = widget_draw(tlb,xsize=850,ysize=850, X_Scroll_size=800, Y_Scroll_size=800, /scroll, /motion_events)
  statusbaseID  = widget_base(tlb,/row,title='status') 
  space1ID      = widget_label(statusbaseID,value='  ')
  enl1textID    = widget_label(statusbaseID,value='  ENL1: 8.0         ', xsize = 80, /frame)
  enl2textID    = widget_label(statusbaseID,value='  ENL2: 8.0         ', xsize= 80, /frame)
  space3ID      = widget_label(statusbaseID,value='     ')
  toggleID      = widget_button(statusbaseID, value = 'Toggle Span', event_pro = 'ONtoggle', sensitive = 0, /align_center)
  space4ID      = widget_label(statusbaseID,value='     ')
  sliderbaseID  = widget_base(statusbaseID,xsize=100)
  slider_sigID  = widget_slider(sliderbaseID,value=1,min=1,max=5,event_pro='ONslider_sig', frame=2, title = ' Significance (%)', Sensitive=0)
  space2ID      = widget_label(statusbaseID,value='  ')
  xytextID      = widget_label(statusbaseID,value=' (0,0) 0            ',xsize=150)
  
  widget_control, tlb, /realize
  widget_control, drawID, get_value=wID
  
; widget communication structure
  
  status =  { $ 
             clearworkID:clearworkID, $
             slcID: slcID, $
             setdemID: setdemID, $
             geocodeID: geocodeID, $
             saveviewID: saveviewID, $
             ingest1ID: ingest1ID, $
             ingest2ID: ingest2ID, $
             load1ID: load1ID, $
             load2ID: load2ID, $
             save1ID: save1ID, $
             save2ID: save2ID, $
             coregisterID: coregisterID, $
             estenlID: estenlID, $
             viewspan1ID: viewspan1ID, $
             viewspan2ID: viewspan2ID, $
             viewenl1ID: viewenl1ID, $
             viewenl2ID: viewenl2ID, $
             viewstatID: viewstatID, $
             viewprobID: viewprobID, $
             viewchangeID: viewchangeID, $
             toggleID: toggleID, $
             changemapID: changemapID, $
             drawID: drawID, $
             wID:      wID, $
             xytextID: xytextID, $
             statetextID: statetextID, $
             enl1textID: enl1textID, $
             enl2textID: enl2textID, $
             slider_sigID: slider_sigID, $
             toggle: 0, $
             current_view: 0, $
             enl1: 8.0, $
             enl2: 8.0, $
             signif: 0.01, $
             cols1: 0,   $
             rows1: 0,    $ 
             dims1: ptr_new(), $
             map_info: ptr_new(), $
             bnames: ptr_new(), $
             enlim1: ptr_new(), $
             enlim2: ptr_new(), $
             image1: ptr_new(), $
             image2: ptr_new(), $
             Z: ptr_new(), $
             CP: ptr_new() $    
         }
  
  !ORDER = 1
  widget_control, tlb, set_Uvalue=status, /no_copy
  
  XManager, 'wishart_wizard',tlb,/no_block, $
             cleanup='wishart_cleanup'

END
