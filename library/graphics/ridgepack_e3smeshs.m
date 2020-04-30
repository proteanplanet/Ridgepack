function [SML]=ridgepack_e3smeshs(ncvert,mask,col)

% ridgepack_e3smeshs - Draws an E3SM scalar mesh 
%
% function ridgepack_e3smeshs(ncvert,mask,col)
%
% This function generates a mesh for a given masked area in a 
% color shaded maps of E3SM fields from MPAS components
% 
% INPUT:
%
% ncvert - Netcdf grid structure from E3SM. It must include the 
%          vectors:
%          nCells, nEdgesOnCell, verticesOnCell, latitude, longitude
% mask   - mask of cell indices to be plotted, given as a vector of 
%          cell indices from nCells that are to be plotted.[optional] 
% col    - This sets the color of the lines [optional]
%
% OUTPUT:
%
% SML - Geoshape of mesh.
%
% Ridgepack Version 2.0
% Andrew Roberts, LANL, 2020 (afroberts@lanl.gov) 

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% check for mask
if (nargin>=2 & isempty(mask)) | nargin<2
 mask=ncvert.nCells.data;
end

% check for contour interval
if nargin>=3 & isempty(col) | nargin<3
 col='b';
end

% get current axes
hmap=get(gcf,'CurrentAxes');
if ~ismap(hmap)
 error('Current axes must be a map')
else
 maphandle=gcm;
 latextrem=maphandle.maplatlimit;
 minlat=max(-90,latextrem(1)-5)*pi/180;
 maxlat=min(90,latextrem(2)+5)*pi/180;
 lidx=find(ncvert.latitude.data>maxlat | ncvert.latitude.data<minlat);
 ncvert.longitude.data(lidx)=NaN;
 ncvert.latitude.data(lidx)=NaN;
end
figout=get(hmap,'OuterPosition');
fontsize=min(11,max(9,10*sqrt((figout(3).^2)+(figout(4).^2))/sqrt(2)));
set(hmap,'fontsize',fontsize);

idx=ncvert.nCells.data;
idxn=intersect(idx,mask);

if length(idxn)>0

      maxsize=ncvert.maxEdges.data(end)+1;

      lat=NaN*zeros(length(idxn),maxsize);
      lon=NaN*zeros(length(idxn),maxsize);

      latitude=[];
      longitude=[];

      for i=1:length(idxn)

       maxidx=ncvert.nEdgesOnCell.data(idxn(i));

       lat(i,1:maxidx)=ncvert.latitude.data(...
               ncvert.verticesOnCell.data(...
                1:maxidx,idxn(i)))*180/pi;
       lon(i,1:maxidx)=ncvert.longitude.data(...
               ncvert.verticesOnCell.data(...
                1:maxidx,idxn(i)))*180/pi;
       
       lon(i,maxidx+1:maxsize)=lon(i,1);
       lat(i,maxidx+1:maxsize)=lat(i,1);
   
       latitude=[latitude NaN lat(i,:)];
       longitude=[longitude NaN lon(i,:)];

      end

      [cc,dd] = mfwdtran(gcm,latitude,longitude);

      plot(cc',dd','Color',col,'Linewidth',0.25)

      clear cc dd lon lat idxn idx

end

drawnow

if nargout>0
 
 % grid line
 SML=geoshape(latitude,longitude,...
             'MPAS_SeaIce','Mesh',...
             'Geometry','line');
 
end

if debug; disp(['...Leaving ',mfilename]); end

