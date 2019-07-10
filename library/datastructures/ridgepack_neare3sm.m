function [cell,lat,lon]=ridgepack_neare3sm(nc,searchlat,searchlon,switchd,lastcell)

global debug;
%debug=true;
if debug; disp(['Entering ',mfilename,'...']); end

if ~isstruct(nc)
 error('nc is not a structure')
end

if nargin<3
 error('incorrect number of inputs')
else
 searchlon=wrapTo180(searchlon);
end

if nargin<4
 switchd=0;
end

if nargin<5
 lastcell=0;
end

% set almanac value
if switchd
 alms=almanac('earth','sphere');
else
 alms=almanac('earth','wgs84');
end

if lastcell==0

 dist=distance(searchlat,searchlon,...
               nc.latCell.data(nc.nhsicells.data),...
               nc.lonCell.data(nc.nhsicells.data),...
               alms);

 newdist=sort(dist(:));

 cellidx=find(dist(:)==newdist(1));

 cell=nc.nhsicells.data(cellidx);

 lat=nc.latCell.data(cell);
 lon=nc.lonCell.data(cell);

end


