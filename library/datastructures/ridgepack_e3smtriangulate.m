function [cell,vert,tvert,incell,cdist,vdist,cidx,vidx,tidx]=...
               ridgepack_e3smtriangulate(ncvert,searchlat,searchlon)

% ridgepack_e3smtriangulate - 
%
% function [cell,vert,tvert,incell,cdist,vdist,cidx,vidx,tidx]=...
%              ridgepack_e3smtriangulate(ncvert,searchlat,searchlon)
%
%

global debug;
%debug=true;
if debug; disp(['Entering ',mfilename,'...']); end

if nargin<3
 error('incorrect number of inputs')
else
 searchlon=wrapTo180(searchlon);
end

if ~isstruct(ncvert)
 error('nccell is not a structure')
elseif ~isfield(ncvert,'nCells')
 error('nc missing nCells')
end

if isfield(ncvert,'latCell')

 [cdist,cangl,cphi]=ridgepack_greatcircle(searchlat,searchlon,...
                                   ncvert.latCell.data*180/pi,...
                                   ncvert.lonCell.data*180/pi);

 [vdist,vangl,vphi]=ridgepack_greatcircle(searchlat,searchlon,...
                                  ncvert.latitude.data*180/pi,...
                                  ncvert.longitude.data*180/pi);

else

 error('There is a problem with ncvert') 

end

% grab the cell and vertex indices, and distances from them
% and calculate phi - the azimuthal angle
[cnewdist,I] = sort(cdist);
cdist=cnewdist(1);
cidx=I(1);
cell=ncvert.nCells.data(cidx);

[vnewdist,I] = sort(vdist);
vdist=vnewdist(1);
vidx=I(1);
vert=ncvert.nVertices.data(vidx);

% Vertex angle: At the vertex from search point to center
slat=searchlat;
slon=searchlon;
clat=ncvert.latCell.data(cidx)*180/pi;
clon=ncvert.lonCell.data(cidx)*180/pi;
vlat=ncvert.latitude.data(vert)*180/pi;
vlon=ncvert.longitude.data(vert)*180/pi;
[x,y,z,phi]=ridgepack_satfwd([slat clat],[slon clon],vlat,vlon);
vphi=wrapTo180((diff(phi))*180/pi);

% Determine the triangulation vertex based on the vertex angle
% In plain language, this determins the triangle relavant to 
% working out if the point falls within cell or not.
maxidx=ncvert.nEdgesOnCell.data(cell);
idx=find(ncvert.verticesOnCell.data(1:maxidx,cell)==vert);
if vphi<=0
 if idx==1
  tvert=ncvert.verticesOnCell.data(maxidx,cell);
 else
  tvert=ncvert.verticesOnCell.data(idx-1,cell);
 end
elseif vphi>0
 if idx==maxidx
  tvert=ncvert.verticesOnCell.data(1,cell);
 else
  tvert=ncvert.verticesOnCell.data(idx+1,cell);
 end
end
tidx=find(ncvert.nVertices.data(tvert));

% Triangulation Vertex angle: At the vertex from tvert to center 
tlat=ncvert.latitude.data(tvert)*180/pi;
tlon=ncvert.longitude.data(tvert)*180/pi;
clat=ncvert.latCell.data(cidx)*180/pi;
clon=ncvert.lonCell.data(cidx)*180/pi;
vlat=ncvert.latitude.data(vert)*180/pi;
vlon=ncvert.longitude.data(vert)*180/pi;
[x,y,z,phi]=ridgepack_satfwd([tlat clat],[tlon clon],vlat,vlon);
tphi=wrapTo180((diff(phi))*180/pi);

% If tphi is greater than vphi, then the point falls within the 
% triangulation between points tvert vert cell.
if abs(tphi)>=abs(vphi)
 incell=true;
else
 incell=false;
end

if debug; disp(['Leaving ',mfilename,'...']); end

