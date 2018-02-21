function [row,col,lat,lon,nrow,ncol,ndist]=ncneargrid(nc,lat,lon,switchd)

% NCNEARGRID - Finds the nearest grid point to a latitude and longitude position
%
% This uses the wgs84 ellipsoid for calculations rather than a spherical 
% earth by default..
% 
% function [row,col,lat,lon,nrow,ncol,ndist]=ncneargrid(nc,lat,lon,switchd)
%
% Input:
% nc      - netcdf structure with latitude and longitude
% lat     - latitude of point
% lon     - longitude of point
% switchd - set to 1 to turn off wgs84 ellipsoid
%
% Output:
% [row,col] - row and column of latmat (or longmat) matrix
% lat     - latitude of point on grid
% lon     - longitude of point on grid
% nrow - nine nearest row indices surrounding point
% ncol - nine nearest col indices surrounding point
% ndist - distance quad surrounding point
%
% Written by Andrew Roberts 2012
% Department of Oceanography, Naval Postgraduate School

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
longmat=nc.longitude.data;
latmat=nc.latitude.data;

% check input
if ndims(latmat)>2
 error('Latitude must be 2D')
elseif ndims(longmat)>2
 error('Longitude must be 2D')
end

% create grid if needed
if size(latmat,1)==1 | size(latmat,2)==1
 [latgrat,longrat] = meshgrat(latmat,longmat);
 longmat=longrat;
 latmat=latgrat;
end

% check arrays again
if size(latmat,1)~=size(longmat,1) | size(latmat,2)~=size(longmat,2)  
 error('latmat and longmat must have the same size')
end

% Keep Longitudes in -180 to 180 convention
longmat=wrapTo180(longmat);
lon=wrapTo180(lon);

% find step in latmat and longmat for accelerated search
sl=100;
if size(latmat,1)>sl+1 & size(latmat,2)>sl+1
 el=10;
else
 sl=1;
end

% set almanac value
if switchd
 alms=almanac('earth','sphere');
else
 alms=almanac('earth','wgs84');
end

% Create distance array from the point, using best ellipsoid.
% If stepl is greater than two, only search every second grid point initially.
dist=distance(lat,lon,latmat(1:sl:end,1:sl:end),longmat(1:sl:end,1:sl:end),alms);
newdist=sort(dist(:));

% This loop is for accelerated searches over large grids
if sl>1

 % first pass on sl grid point leap grid
 [row,col]=find(dist(:,:)==newdist(1));

 row=1+sl*(row-1);
 col=1+sl*(col-1);

 xmin=max(1,row-ceil(sl/2));
 xmax=min(row+ceil(sl/2),size(latmat,1));
 ymin=max(1,col-ceil(sl/2));
 ymax=min(col+ceil(sl/2),size(latmat,2));

 % second pass on zoomed-in area
 dist=distance(lat,lon,...
               latmat(xmin:el:xmax,ymin:el:ymax),longmat(xmin:el:xmax,ymin:el:ymax),...
               alms);
 newdist=sort(dist(:));

 % second pass on el grid point leap grid
 [row,col]=find(dist(:,:)==newdist(1));

 row=xmin+el*(row-1);
 col=ymin+el*(col-1);

 xmin=max(1,row-ceil(el/2));
 xmax=min(row+ceil(el/2),size(latmat,1));
 ymin=max(1,col-ceil(el/2));
 ymax=min(col+ceil(el/2),size(latmat,2));

 % third pass on zoomed-in area
 dist=distance(lat,lon,latmat(xmin:xmax,ymin:ymax),longmat(xmin:xmax,ymin:ymax),alms);
 newdist=sort(dist(:));
 
 for i=1:9
  ndist(i)=newdist(i);
  [nrow(i),ncol(i)]=find(dist(:,:)==ndist(i));
 
  nrow(i)=xmin+nrow(i)-1;
  ncol(i)=ymin+ncol(i)-1;
 end

 row=nrow(1);
 col=ncol(1);

 if debug % check that zooming in works

  dist=distance(lat,lon,latmat,longmat,alms);
  [rowcheck,colcheck]=find(dist(:,:)==newdist(1));

  if rowcheck~=row | colcheck~=col
   disp(['rowcheck: ',num2str(rowcheck)]);
   disp(['     row: ',num2str(row)]);
   disp(['colcheck: ',num2str(colcheck)]);
   disp(['     col: ',num2str(col)]);
   error('failed index check')
  end

 end % debug

else

 for i=1:9
  ndist(i)=newdist(i);
  [nrow(i),ncol(i)]=find(dist(:,:)==ndist(i));
 
  nrow(i)=1+sl*(nrow(i)-1);
  ncol(i)=1+sl*(ncol(i)-1);
 end


 row=nrow(1);
 col=ncol(1);

end

lat=latmat(row,col);
lon=longmat(row,col);

if debug
 disp(['Nearest location on the grid is ',num2str(latmat(row,col)),', ',...
                                          num2str(longmat(row,col))])
end

if debug; disp(['...Leaving ',mfilename]); end


