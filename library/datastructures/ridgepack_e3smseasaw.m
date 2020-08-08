function [nc,SCP,SCL,STP,STL]=...
             ridgepack_e3smseasaw(ncvert,ncc,varc,threshold)

% ridgepack_e3smseasaw - Generate a threshold on an unstructured mesh
%
% function [nc,SCP,SCL,STP,STL]=...
%             ridgepack_e3smseasaw(ncvert,ncc,varc,threshold)
%
% This function generates a coastline from an MPAS-Ocean unstructured
% mesh, and, in addition, also generates a region defined by a 
% threshold on a scalar field. Output is in the form of an nc
% structure and as geoshapes for the coast and threshold region.
%  
% INPUT:
%
% ncvert    - netcdf structure of mesh information from MPAS-Ocean.
%             If only ncvert is specified, this generates only a
%             coastline, formed of closed contours.
% ncc       - netcdf structure for variable varc on ncvert 
%             mesh [optional]
% varc      - variable in ncc to be used [character variable, optional]
% threshold - scalar threshold value on an unstructured mesh 
%             [optional]
% 
% OUTPUT:
% 
% nc  - netcdf structure with lats, longs, mesh cell and vertex
%       indices of the coast, contour and shapes of areas.
%       Included in the metadata is the length of the perimeter
%       for each defined region, and the area within.
% SCP - Polygons of the E3SM coastline as a geoshape.
% SCL - Line of the E3SM coastline as a geoshape.
% STP - Polygons of the E3SM threshold as a geoshape.
% STL - Line of the E3SM threshold contours as a geoshape.
%
% Ridgepack Version 2.0
% Andrew Roberts, Los Alamos National Laboratory, 2020 
% afroberts@lanl.gov

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

