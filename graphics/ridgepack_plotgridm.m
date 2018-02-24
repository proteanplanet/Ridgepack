function ridgepack_plotgridm(nc,type)

% ridgepack_plotGRIDM - Overlay a model grid on a map
%
% function ridgepack_plotgridm(nc,type)
%
% This function plots a grid on a map either as individual grid
% points or as the boundary lines of the grid.  There are several
% types of plot, including one that provides the mask outline
% and one that simply plots the grid points.  
%
% INPUT:
%
% nc   - nc structure containing the fields 'latitude' and 'longitude'
%        of the model grid.  In addition, the field 'mask' can be
%        used depending on the type of grid plot selected
%
% type - Text string sepecifying one of:
%        'outline'     - line indicating the boundaries of the grid
%                        (this is the default if type is left blank)
%        'maskoutline' - line indicating boundaries with patched mask
%        'points'      - each grid point is plotted as a point
%        'maskpoints'  - unmasked grid point is plotted as a point
%
% Output is graphical only.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
 

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% Check that nc is a structure
if not(isstruct(nc));
 error([inputname(1),' is not a structure']);
elseif nargin<2 | isempty(type)
 type='outline';
end

% check the grid type and check latitude and longitudes
gridtype=ridgepack_gridtype(nc);

% check that a mask exists
if isfield(nc,'mask') & size(nc.mask.data)==size(nc.latitude.data)
 maskexists=true;
else
 maskexists=false;
end

% pass through set of logicals for the type of plot required
if strcmp(type,'outline')
 %color=[0.25 0.25 0.25];
 color=[0 0 1];
 lwidth=1;
 if strcmp(gridtype,'regional')
  plotm(nc.latitude.data(1,:),nc.longitude.data(1,:),'Color',color,'linewidth',lwidth);
  plotm(nc.latitude.data(end,:),nc.longitude.data(end,:),'Color',color,'linewidth',lwidth);
  plotm(nc.latitude.data(:,1),nc.longitude.data(:,1),'Color',color,'linewidth',lwidth);
  plotm(nc.latitude.data(:,end),nc.longitude.data(:,end),'Color',color,'linewidth',lwidth);
 else
  error('Grid is global, so the outline cannot be plotted')
 end
elseif strcmp(type,'maskoutline')
 color=[0.25 0.25 0.25];
 lwidth=0.75;
 if maskexists
  disp('Using ridgepack_maskm');
  mask=nc.mask.data;
  mask(mask==0)=NaN;
  ridgepack_maskm(nc.latitude.data,nc.longitude.data,mask,color,lwidth)
 else
  error('No mask exists in this nc structure');
 end
 gridtype
 if strcmp(gridtype,'regional')
  lwidth=1.5;
  style='.';
  color=[0. 0. 0.];
  lat=nc.latitude.data;
  lon=nc.longitude.data;
  lon(isnan(mask))=NaN;
  lat(isnan(mask))=NaN;
  plotm(lat(1,:),lon(1,:),style,'Color',color,'linewidth',lwidth);
  plotm(lat(end,:),lon(end,:),style,'Color',color,'linewidth',lwidth);
  plotm(lat(:,1),lon(:,1),style,'Color',color,'linewidth',lwidth);
  plotm(lat(:,end),lon(:,end),style,'Color',color,'linewidth',lwidth);
 else
  error('Grid is global, so the outline cannot be plotted')
 end
elseif strcmp(type,'points')
 style='.';
 color=[0.15 0.15 0.15];
 msize=0.75;
 plotm(nc.latitude.data,nc.longitude.data,style,'Color',color,'MarkerSize',msize);
elseif strcmp(type,'maskpoints')
 if maskexists
  style='.';
  color=[0.15 0.15 0.15];
  msize=0.75;
  mask=nc.mask.data;
  mask(mask==0)=NaN;
  lat=nc.latitude.data;
  lon=nc.longitude.data;
  lon(isnan(mask))=NaN;
  lat(isnan(mask))=NaN;
  plotm(lat,lon,style,'Color',color,'MarkerSize',msize);
 else
  error('No mask exists in this nc structure');
 end
else
 error(['The type ',type,'does not exist as an option']);
end

drawnow

if debug; disp(['...Leaving ',mfilename]); end


