function [vertices,cells,edge]=ridgepack_e3smnear(ncvert,searchlat,searchlon)

global debug;
%debug=true;
if debug; disp(['Entering ',mfilename,'...']); end

if ~isstruct(nc)
 error('nc is not a structure')
end

if nargin<3
 error('incorrect number of inputs')
elseif 
 searchlon=wrapTo180(searchlon);
end

dist=ridgepack_greatcircle(lat1,lon1,lat2,lon2);

dist=distance(searchlat,searchlon,...
               nc.latCell.data(nc.nhsicells.data),...
               nc.lonCell.data(nc.nhsicells.data),...
               alms);

newdist=sort(dist(:));

cellidx=find(dist(:)==newdist(1));

cell=nc.nhsicells.data(cellidx);

lat=nc.latCell.data(cell);
lon=nc.lonCell.data(cell);

if debug; disp(['Leaving ',mfilename,'...']); end

