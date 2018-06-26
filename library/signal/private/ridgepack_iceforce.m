function ridgepack_iceforce(year,name,rec1,rec2,spatial,temporal,grid)

% ridgepack_iceforce - Interpolates and filters forcing fields for ice-ocean models
%
% function ridgepack_iceforce(year,name,rec1,rec2,spatial,temporal,grid)
%
% INPUT:
%
% year     - year for which data is being interpolated
%
% name     - name of variable in the file to be interpolated. For the NCEP
%            dataset, this can take the following values:
%
%            mslp  - which will output geostrophic winds (not MSLP)
%            air   - 2m air temperature
%            shum  - 2m specfic humidity 
%            prate - precipitation rate
%            dswrf - downward shortwave radiation
%            dlwrf - downward longwave radiation
%
%            If multiple fields need to be interpolated, then the field names
%            can be entered as a cell array (e.g. {'mslp','air','shum'})
%
% rec1     - start record (optional). This can be entered in as either an index
%            in the netcdf file being interpolated, or as string specifying
%            time, such as '1-1-1999 06:00:0.0'.
%
% rec2     - end record (optional). Can be either an index or string, as for rec1.
%
% spatial  - if set to true, spatially filters data (optional, default is true)
% temporal - if set to true, time filters data (optional, default is true)
% grid     - chosen interpolation grid
%
%
% OUTPUT:
%
% Output is by way of a new netcdf file written from Matlab
% with the nomenclature icepack_file and the final netcdf structure nc 
% of filtered data.
%
% This routine is currently set up to interpolate and filter NCEP data, but
% can easily be expanded to include the ERA40 dataset. This routine must
% be run in the same directory as the original NCEP data. 
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% Naval Postgraduate School, Department of Oceanography
% 

% check running time
tic

% determine data to be extracted
if nargin<3; rec1=1; end
if nargin<4; rec2=[]; end
if nargin<5; spatial=1 ; end
if nargin<6; temporal=1 ; end
if nargin<7; grid=1 ; end
if iscell(name); 
 namelist=name; 
else
 namelist={name}; 
end

for i=1:length(namelist)

name=char(namelist{i});

 % for mslp, output geostropic wind rather than pressure
 if strcmp(name,'mslp') 

  for gg=1:2

   if gg==1; 
	wind='ugwind'; antiwind='vgwind';
   elseif gg==2; 
	wind='vgwind'; antiwind='ugwind';
   else
	error('wrong wind index')
   end

   % read in the data
   disp('opening the file')
   file=[name,'.',num2str(year),'.nc'];
   nc=ncclone(file,name);

   % regrid the data to the ice-ocean grid, padding the boundaries with 35 extra
   % points for spatial filter boundary conditions twice
   disp('regrid the data')
   if spatial & temporal
    nc=ncregrid(nc,name,'',grid,rec1,rec2,35*2);
   else
    nc=ncregrid(nc,name,'',grid,rec1,rec2,35);
   end

   % filter this field to remove short wavelengths before calculating geostophic
   disp('spatially filter the data pressure')
   nc=ridgepack_filter2(nc,name,{'x','y'},300,'low',70,292,15);

   % calculate the geostrophic velocity based on this data
   nc=ncgeostrophic(nc,name);
   nc=rmfield(nc,antiwind);

   % spatially filter each time slice with a 216=(2*108)th order FIR filter
   % basing the nyquist frequency on the size of the final grid y direction=292 points
   % creating two strutures, one for u and v component filtered geostrophic winds
   % and then combining them into ncguf.
   if spatial
    disp(['spatially filtering ',wind])
    nc=ridgepack_filter2(nc,wind,{'x','y'},300,'low',70,292,15)
   end
 
   % add attribute information to the structure 
   nc.attributes.title=['NCEP2 reanalysis data ',name,' interpolated to ice-tide grid'] ;
   nc.attributes.source='Generated using bi-cubic interpolation and FIR spatial and temporal filters';
   
   % write out the new dataset without time filtering
   disp('writing temporally unfiltered scratch file')
   ncwrite(nc,['icepack_scratch_wind_',num2str(year)]);
 
   % now interpolate and filter in time, pealing off each grid point time series
   disp(['time filtering ',wind])

   % set file name
   fileout=['icepack_',wind,'_',num2str(year)];
   if spatial & temporal
	fileout=[fileout,'_spatial_temporal'];
   elseif spatial
	fileout=[fileout,'_spatial'];
   elseif temporal
	fileout=[fileout,'_temporal'];
   end

   % accelerate filtering taylor made computers from .host file in home directory
   [status,comp]=system('cat ~/.host')
   host=comp(1:length(comp)-1);
   complist={'gemini'}; % add to complist here

   display('Running time filtering in low-memory stripe mode')
 
   jmax=length(nc.y.data);

   for j=1:jmax;
 
     ncnew=ncclone(['icepack_scratch_wind_',num2str(year)],wind,{'y'},{j});
 
     if temporal
      ncnew=ridgepack_timeprocess(ncnew,wind,'','','linear',[0 0 0 3 0 0.0],1.6);
     else
      ncnew=ridgepack_timeprocess(ncnew,wind,'','','linear',[0 0 0 3 0 0.0]);
     end
 
     if j==1 
      nc.time=ncnew.time;
      nc.attributes=ncnew.attributes;
      ncwrite(nc,fileout);
      clear nc
     end

     ncstruct(ncnew)
  
     ncwrite(ncnew,fileout,{'y'},{j});
  
   end
 
  end
 
 else

  % set filename
  if strcmp(name,'air')
   file=[name,'.2m.gauss.',num2str(year),'.nc'];
  elseif strcmp(name,'shum')
   file=[name,'.2m.gauss.',num2str(year),'.nc'];
  elseif strcmp(name,'prate')
   file=[name,'.sfc.gauss.',num2str(year),'.nc'];
  elseif strcmp(name,'dswrf')
   file=[name,'.sfc.gauss.',num2str(year),'.nc'];
  elseif strcmp(name,'dlwrf')
   file=[name,'.sfc.gauss.',num2str(year),'.nc'];
  end

  % read in the data
  disp('Opening the file')
  nc=ncclone(file,name);

  % remove level markers in NCEP2 dataset
  if strcmp(name,'air')
   nc=ncreduce(nc,{'level'});
  elseif strcmp(name,'shum')
   nc=ncreduce(nc,{'level'});
  end

  % regrid the data to the ice-ocean grid, padding the boundaries with 35 extra
  % points for spatial filter boundary conditions twice
  disp('Regrid the data')
  if spatial
   nc=ncregrid(nc,name,'',1,rec1,rec2,35);
  else
   nc=ncregrid(nc,name,'',1,rec1,rec2);
  end

  % spatially filter each time slice with a 216=(2*108)th order FIR filter
  % basing the nyquist frequency on the size of the final grid y direction=292 points
  % creating two strutures, one for u and v component filtered geostrophic winds
  % and then combining them into ncguf.
  if spatial
   disp(['Spatially filtering ',name])
   nc=ridgepack_filter2(nc,name,{'x','y'},300,'low',70,292,15)
  end
 
  % add attribute information to the structure 
  nc.attributes.title=['NCEP2 reanalysis data ',name,' interpolated to ice-tide grid'] ;
  nc.attributes.source='Generated using bi-cubic interpolation.';
  
  % write out the new dataset without time filtering
  disp('Writing temporally unfiltered scratch file')
  ncwrite(nc,'clobber',['icepack_scratch_',name,'_',num2str(year)]);
 
  % now interpolate and filter in time, pealing off each grid point time series
  disp(['Time filtering ',name])

  % set file name
  fileout=['icepack_',name,'_',num2str(year)];
  if spatial & temporal
	fileout=[fileout,'_spatial_temporal'];
  elseif spatial
	fileout=[fileout,'_spatial'];
  elseif temporal
	fileout=[fileout,'_temporal'];
  end

  display('Running time filtering in low-memory stripe mode')

  jmax=length(nc.y.data);

  for j=1:jmax;
 
    ncnew=ncclone(['icepack_scratch_',name,'_',num2str(year)],name,{'y'},{j});

    if temporal
     ncnew=ridgepack_timeprocess(ncnew,name,'','','linear',[0 0 0 3 0 0.0],1.6);
    else
     ncnew=ridgepack_timeprocess(ncnew,name,'','','linear',[0 0 0 3 0 0.0]);
    end
 
    if j==1 
     nc.time=ncnew.time;
     nc.attributes=ncnew.attributes;
     ncwrite(nc,'fill',fileout);
     clear nc
    end
 
    ncwrite(ncnew,'slot',fileout,name,{'y'},{j});
  
  end
 
 end

end

% complete by finding running time
toc
