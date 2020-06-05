function [SML]=ridgepack_e3smeshs(ncvert,mask,col,fill)

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
% fill   - Logical set to true if the patch is to be filled 
%          The default is not to fill. [optional]
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

% check for contour interval
if nargin>=4 & isempty(fill) | nargin<4
 fill=false;
end

% get current axes
if nargout==0
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
end

idx=ncvert.nCells.data;
idxn=intersect(idx,mask);

if length(idxn)>0

      maxsize=ncvert.maxEdges.data(end)+1;

      % These lats and longs are for drawing the patch
      lat=NaN*zeros(length(idxn),maxsize);
      lon=NaN*zeros(length(idxn),maxsize);

      % This series is for generating a geoshape
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

      if nargout==0

       % translate las and longs into cartesian coords
       [cc,dd] = mfwdtran(gcm,lat(:,:),lon(:,:));

       if fill
        patch(cc',dd',col,...
             'EdgeColor',[0.5 0.5 1],'LineWidth',0.2)
       else
        patch(cc',dd',col,'FaceColor','none',...
             'EdgeColor',col,'LineWidth',0.2)
       end

       clear cc dd lon lat idxn idx

      end

end


if nargout>0
 
 % grid line
 SML=geoshape(latitude,longitude,...
             'MPAS_SeaIce','Mesh',...
             'Geometry','line');

else

 drawnow
 
end

% debug stuff
if debug; disp(['...Leaving ',mfilename]); end

