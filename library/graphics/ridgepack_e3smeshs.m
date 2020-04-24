function ridgepack_e3smeshs(ncvert,mask,color)

% ridgepack_e3smcolors - Color fills scalar data from E3SM/MPAS on an unstructured mesh 
%
% function ridgepack_e3smcolors(nc,var,ncvert,mask,cont,loglin,ref,horiz,colors,colvals)
%
% This function generates color shaded maps of E3SM fields from MPAS components
% 
% INPUTS:
%
% ncvert - Netcdf grid structure from E3SM. It must include the vectors:
%          nCells, nEdgesOnCell, verticesOnCell, latitude, longitude
%
% mask   - mask of cell indices to be plotted.
%
% color - This sets the color scheme required:
%
% Ridgepack Version 1.2
% Andrew Roberts, LANL, 2019 (afroberts@lanl.gov) 

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% check for mask
if (nargin>=2 & isempty(mask)) | nargin<2
 mask=ncvert.nCells.data;
end

% check for contour interval
if nargin>=3 & isempty(cont) | nargin<3
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

title(['MPAS S-Mesh']);

% debug stuff
if debug; disp(['...Leaving ',mfilename]); end

