function [row,col,lat,lon]=...
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
% [row,col] - row and column of latmat (or longmat) matrix
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

% set almanac value
if switchd
 alms=almanac('earth','sphere');
else
 alms=almanac('earth','wgs84');
end

% set xmin and xmax
xmin=1;
xmax=size(latmat,1);
ymin=1;
ymax=size(latmat,2);

if nargin==6 & lastrow>0 & lastcol>0

 lr=lastrow;
 lc=lastcol;

 dist=distance(lat,lon,...
               latmat(max(lr-4,xmin):min(lr+4,xmax),...
                      max(lc-4,ymin):min(lc+4,ymax)),...
               longmat(max(lr-4,xmin):min(lr+4,xmax),...
                       max(lc-4,ymin):min(lc+4,ymax)),...
               alms);

 newdist=sort(dist(:));

 [row,col]=find(dist(:,:)==newdist(1));

 if (row==max(lr-4,xmin) | row==min(lr+4,xmax) | ...
     col==max(lc-4,ymin) | col==min(lc+4,ymax))
  widesearch=true;
 else
  widesearch=false;
 end

else

 widesearch=true;

end

% grid-wide search
if widesearch

 % find step in latmat and longmat for accelerated search
 sl=max(1,floor(min(size(longmat,1),size(longmat,2))/10));
 if size(latmat,1)>sl+1 & size(latmat,2)>sl+1
  el=10;
 else
  sl=1;
 end

 % search every sl grid point
 dist=distance(lat,lon,...
               latmat(1:sl:end,1:sl:end),...
               longmat(1:sl:end,1:sl:end),alms);
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

  %el=max(1,floor(min(xmax-xmin,ymax-ymin)/10));
  el=1;

  % second pass on zoomed-in area
  dist=distance(lat,...
                lon,...
                latmat(xmin:el:xmax,ymin:el:ymax),...
                longmat(xmin:el:xmax,ymin:el:ymax),...
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
  dist=distance(lat,lon,...
                latmat(xmin:xmax,ymin:ymax),...
                longmat(xmin:xmax,ymin:ymax),alms);
  newdist=sort(dist(:));
 
  % set final row
  [row,col]=find(dist(:,:)==newdist(1));

  row=xmin+row-1;
  col=ymin+col-1;

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

  % row and column of nearest point
  [row,col]=find(dist(:,:)==newdist(1));

 end

elseif debug

 disp('No widesearch necessary')

end

lat=latmat(row,col);
lon=longmat(row,col);

if debug
 disp(['Nearest location on the grid is ',num2str(latmat(row,col)),', ',...
                                          num2str(longmat(row,col))])
end

if debug; disp(['...Leaving ',mfilename]); end


