function [nc,SCP,STP,STL]=ridgepack_e3smseasaw(ncvert,ncc,varc,threshold,...
                                       centlat,centlon,horizon)

% ridgepack_e3smseasaw - generate a threshold on an unstructured mesh
%
% function [nc]=ridgepack_e3smseasaw(ncc,varc,threshold,ncvert,centlat,centlon,horizon)
%
% This function .... description
%  
% INPUT:
%
% ncc       - netcdf structure including
% varc      - variable in ncc to be used
% threshold - threshold value used as a boundary on unstructured mesh
% ncvert    - vertices on the unstructed mesh
% centlat   - center latitude of plotted satellite view (optional)
% centlon   - center longitude if plotted in satellite view (optional)
% horizon   - horizon, in degrees of satellite view (optional)
%
% Note that centlat, centlon and horizon should not be provided
% if nc is specificied.
% 
% OUTPUT:
% 
% nc - netcdf structure with lats, longs, mesh cell and vertex
%      indices
%
% Ridgepack Version 1.1
% Andrew Roberts, Los Alamos National Laboratory, April 3, 2020 (afroberts@lanl.gov)

if nargin<4
 error('There is no threshold value set, and possibly more')
end

if nargin<5
 centlat=90;
elseif nargin>3 & nargout==1
 error('you are asking to both plot and produce coastal data')
elseif ~isnumeric(centlat)
 error('centlat should be a real number between -90 and 90')
elseif centlat>90 | centlat<-90
 error('centlat should be between -90 and 90')
end

if nargin<6
 centlon=0;
elseif ~isnumeric(centlon)
 error('centlon should be a real number between -90 and 90')
elseif centlat>180 | centlat<-180
 error('centlon should be between -180 and 180')
end

if nargin<7
 horizon=90;
elseif ~isnumeric(horizon)
 error('horizon must be a real number between 0 and 90')
elseif horizon<0 | horizon>90
 error('horizon must be between 0 and 90')
end

% Generate a coastline
cidx=[1:length(ncvert.nCells.data)]'; %coast
[clats,clons,cverts]=ridgepack_e3smperimeter(ncvert,cidx);

% If requested, also generate a threshold
if ~isempty(ncc)

 infill=true;

 cidx=find(ncc.(varc).data>threshold); 
 [tlats,tlons,tverts]=ridgepack_e3smperimeter(ncvert,cidx,infill);

 [xlats,xlons,xverts]=...
                 ridgepack_e3smcontour(tlats,tlons,tverts,cverts);

end

% create netcdf structure
nc.attributes.title='E3SM-MPAS Edge Definition';

nc.attributes.coastal_segments=length(isnan(cverts))+1;

nc.cnpoints.data=[1:length(cverts)];
nc.cnpoints.long_name='number of points on coastal outline';
nc.cnpoints.dimension={'cnpoints'};

nc.clatitude.data=clats;
nc.clatitude.long_name='latitude of coast vertices';
nc.clatitude.dimension={'cnpoints'};

nc.clongitude.data=clons;
nc.clongitude.long_name='longitude of coast vertices';
nc.clongitude.dimension={'cnpoints'};

nc.cvertices.data=cverts;
nc.cvertices.long_name='MPAS vertices on coast';
nc.cvertices.dimension={'cnpoints'};

if ~isempty(ncc)

 nc.attributes.threshold_segments=length(isnan(tverts))+1;

 nc.npoints.data=[1:length(tverts)];
 nc.npoints.long_name='number of points on threshold shapes';
 nc.npoints.dimension={'npoints'};

 nc.latitude.data=tlats;
 nc.latitude.long_name='latitude of threshold shapes';
 nc.latitude.dimension={'npoints'};

 nc.longitude.data=tlons;
 nc.longitude.long_name='longitude of threshold shapes';
 nc.longitude.dimension={'npoints'};

 nc.vertices.data=tverts;
 nc.vertices.long_name='vertex indices on threshold shapes';
 nc.vertices.dimension={'npoints'};

 nc.attributes.contour_segments=length(isnan(xverts))+1;

 nc.xnpoints.data=[1:length(xverts)];
 nc.xnpoints.long_name='number of points on threshold contours';
 nc.xnpoints.dimension={'xnpoints'};

 nc.xlatitude.data=xlats;
 nc.xlatitude.long_name='latitude of threshold contours';
 nc.xlatitude.dimension={'xnpoints'};

 nc.xlongitude.data=xlons;
 nc.xlongitude.long_name='longitude of threshold contours';
 nc.xlongitude.dimension={'xnpoints'};

 nc.xvertices.data=xverts;
 nc.xvertices.long_name='vertex indices of contours';
 nc.xvertices.dimension={'xnpoints'};

end

% create geoshapes
SCP=geoshape(nc.clatitude.data,nc.clongitude.data,...
            'MPAS_Ocean','Coastal Definition',...
            'Geometry','polygon');

if ~isempty(ncc)

 STL=geoshape(nc.xlatitude.data,nc.xlongitude.data,...
             'MPAS_SeaIce',[varc,'>',num2str(threshold),' contour'],...
             'Geometry','line');

 STP=geoshape(nc.latitude.data,nc.longitude.data,...
             'MPAS_SeaIce',[varc,'>',num2str(threshold)],...
             'Geometry','polygon');

else

 STL=[];

 STP=[];

end



