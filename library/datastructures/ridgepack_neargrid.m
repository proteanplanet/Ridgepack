function [row,col,lat,lon,rowmat,colmat]=...
           ridgepack_neargrid(nc,lat,lon,switchd,lastrow,lastcol)

% ridgepack_neargrid - Finds the nearest grid point to a latitude and longitude position
%
% This uses the wgs84 ellipsoid for calculations rather than a spherical 
% earth by default..
% 
% function [row,col,lat,lon]=...
%            ridgepack_neargrid(nc,lat,lon,switchd,lastrow,lastcol)
%
% INPUT:
%
% nc      - netcdf structure with latitude and longitude
% lat     - latitude of point
% lon     - longitude of point
% switchd - set to 1 to turn off wgs84 ellipsoid
% lastrow - provide row of previous search to greatly improve efficiency
% lastcol - provide col of previous search to greatly improve efficiency
%
% OUTPUT:
%
% [row,col] - row and column of latmat (or lonmat) matrix
% lat     - latitude of point on grid
% lon     - longitude of point on grid
% ndist - distance quad surrounding point
%
% lastcol and lastrow are used when performing searches along a track
% or path, and the search is related to the previous search. This then
% searches in close proximity to the previous location, thus greatly
% reducing the search area.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

global debug;
%debug=true;
if debug; disp(['Entering ',mfilename,'...']); end

if ~isstruct(nc)
 error('nc is not a structure')
end

if nargin<4
 switchd=0;
end

if nargin<3
 error('incorrect number of inputs')
else
 lon=wrapTo180(lon);
end
 
if isnan(lat) | isnan(lon)
 error(['Input lats and/or longs are incorrect: ',num2str([lat lon])])
end

% assign array
lonmat=nc.longitude.data;
latmat=nc.latitude.data;
rowmat=[];
colmat=[];

% check input
if ndims(latmat)>2
 error('Latitude must be 2D')
elseif ndims(lonmat)>2
 error('Longitude must be 2D')
end

% create grid if needed
if size(latmat,1)==1 | size(latmat,2)==1
 [latgrat,longrat] = meshgrat(latmat,lonmat);
 lonmat=longrat;
 latmat=latgrat;
end

% check arrays again
if size(latmat,1)~=size(lonmat,1) | size(latmat,2)~=size(lonmat,2)  
 error('latmat and lonmat must have the same size')
end


% Keep Longitudes in -180 to 180 convention
lonmat=wrapTo180(lonmat);
lon=wrapTo180(lon);

% set almanac value
if switchd
 alms=almanac('earth','sphere');
else
 alms=almanac('earth','wgs84');
end

if nargin==6 & lastrow>0 & lastcol>0

 lr=lastrow;
 lc=lastcol;

 % set xmin and xmax
 xmin=max(lr-10,1);
 xmax=min(lr+10,size(latmat,1));
 ymin=max(lc-10,1);
 ymax=min(lc+10,size(latmat,2));

 dist=distance(lat,lon,...
               latmat(xmin:xmax,ymin:ymax),...
               lonmat(xmin:xmax,ymin:ymax),...
               alms);

 newdist=sort(dist(:));

 [row,col]=find(dist(:,:)==newdist(1));

 row=max(1,min(xmin+row-1,size(latmat,1)));
 col=max(1,min(ymin+col-1,size(latmat,2)));

 if (row==xmin | row==xmax | col==ymin | col==ymax)
  widesearch=true;
 else
  widesearch=false;
 end

else

 widesearch=true;

end

% grid-wide search
if widesearch

 % search every sl grid point
 dist=distance(lat,lon,...
               latmat(1:1:end,1:1:end),...
               lonmat(1:1:end,1:1:end),alms);
 newdist=sort(dist(:));

 % row and column of nearest point
 [row,col]=find(dist(:,:)==newdist(1));

elseif debug

 disp('No widesearch necessary')

end

lat=latmat(row,col);
lon=lonmat(row,col);

if debug
 disp(['Nearest location on the grid is ',...
                        num2str(latmat(row,col)),', ',...
                        num2str(lonmat(row,col))])
end

if debug; disp(['...Leaving ',mfilename]); end


