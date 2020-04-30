function [nc,SCP,SCL,STP,STL]=...
             ridgepack_e3smseasaw(ncvert,ncc,varc,threshold)

% ridgepack_e3smseasaw - generate a threshold on an unstructured mesh
%
% function [nc,SCP,SCL,STP,STL]=...
%             ridgepack_e3smseasaw(ncvert,ncc,varc,threshold)
%
% This function generates a coastline as well as a threshold on an 
% unstructured mesh given a scalar field varc in the netcdf 
% structure ncc.
%  
% INPUT:
%
% ncc       - netcdf structure including
% varc      - variable in ncc to be used
% threshold - threshold value used as a boundary on unstructured mesh
% ncvert    - vertices on the unstructed mesh
% 
% OUTPUT:
% 
% nc  - netcdf structure with lats, longs, mesh cell and vertex
%       indices of the coast, contour and shapes of areas
% SCP - Polygons of the E3SM coastline as a geoshape.
% SCL - Line of the E3SM coastline as a geoshape.
% STP - Polygons of the E3SM threshold as a geoshape.
% STL - Line of the E3SM threshold contours as a geoshape.
%
% Ridgepack Version 2.0
% Andrew Roberts, Los Alamos National Laboratory, 2020 

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if nargin==1
 ncc=[];
elseif nargin<4
 error('There is no threshold value set, and possibly more')
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

nc.attributes.coastal_segments=num2str(length(isnan(cverts))+1);

nc.npoints.data=[1:length(cverts)];
nc.npoints.long_name='number of points on coastal outline';
nc.npoints.dimension={'npoints'};

nc.latitude.data=clats;
nc.latitude.long_name='latitude of coast vertices';
nc.latitude.dimension={'npoints'};

nc.longitude.data=clons;
nc.longitude.long_name='longitude of coast vertices';
nc.longitude.dimension={'npoints'};

nc.vertices.data=cverts;
nc.vertices.long_name='MPAS vertices on coast';
nc.vertices.dimension={'npoints'};

if ~isempty(ncc)

 nc.attributes.threshold_segments=num2str(length(isnan(tverts))+1);

 nc.tnpoints.data=[1:length(tverts)];
 nc.tnpoints.long_name='number of points on threshold shapes';
 nc.tnpoints.dimension={'tnpoints'};

 nc.tlatitude.data=tlats;
 nc.tlatitude.long_name='latitude of threshold shapes';
 nc.tlatitude.dimension={'tnpoints'};

 nc.tlongitude.data=tlons;
 nc.tlongitude.long_name='longitude of threshold shapes';
 nc.tlongitude.dimension={'tnpoints'};

 nc.tvertices.data=tverts;
 nc.tvertices.long_name='vertex indices on threshold shapes';
 nc.tvertices.dimension={'tnpoints'};

 nc.attributes.contour_segments=num2str(length(isnan(xverts))+1);

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

% run checks on structure
nc=ridgepack_struct(nc);

% create geoshapes
% coastal line
SCL=geoshape(nc.latitude.data,nc.longitude.data,...
            'MPAS_Ocean','Coastal Definition',...
            'Geometry','line');

% coastal polygon
SCP=geoshape(nc.latitude.data,nc.longitude.data,...
            'MPAS_Ocean','Coastal Definition',...
            'Geometry','polygon');

if ~isempty(ncc)


 % contour line
 STL=geoshape(nc.xlatitude.data,nc.xlongitude.data,...
             'MPAS_SeaIce',[varc,'>',num2str(threshold),' contour'],...
             'Geometry','line');

 % polygon
 STP=geoshape(nc.tlatitude.data,nc.tlongitude.data,...
             'MPAS_SeaIce',[varc,'>',num2str(threshold)],...
             'Geometry','polygon');

else

 STL=[];

 STP=[];

end

% debug stuff
if debug; disp(['...Leaving ',mfilename]); end

