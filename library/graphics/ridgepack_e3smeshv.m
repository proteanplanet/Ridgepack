function [STL]=ridgepack_e3smeshv(nccell,mask,col)

% ridgepack_e3smeshv - Draws an E3SM vector mesh 
%
% function ridgepack_e3smeshv(nccell,mask,col)
%
% This function generates a mesh for a given masked area in a 
% color shaded maps of E3SM fields from MPAS components
% 
% INPUT:
%
% nccell - Netcdf grid structure from E3SM. It must include the 
%          vectors:
%          nCells, nEdgesOnCell, verticesOnCell, latitude, longitude
% mask   - mask of cell indices to be plotted, given as a vector of 
%          cell indices from nCells that are to be plotted.[optional] 
% col    - This sets the color of the lines [optional]
%
% OUTPUT:
%
% STL - Geoshape of mesh.
%
% Ridgepack Version 2.0
% Andrew Roberts, LANL, 2020 (afroberts@lanl.gov) 

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% check for mask
if (nargin>=2 & isempty(mask)) | nargin<2
 mask=nccell.nVertices.data;
end

% set color
if nargin<3
 col=0.25*[1 1 1];
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
 lidx=find(nccell.latitude.data>maxlat | nccell.latitude.data<minlat);
 nccell.longitude.data(lidx)=NaN;
 nccell.latitude.data(lidx)=NaN;
end
figout=get(hmap,'OuterPosition');
fontsize=min(11,max(9,10*sqrt((figout(3).^2)+(figout(4).^2))/sqrt(2)));
set(hmap,'fontsize',fontsize);

idx=nccell.nVertices.data;
idxn=intersect(idx,mask);

% pass through vertices within this contour band
if length(idxn)>0

      % grab maximum number of cells per vertex
      maxsize=5;
      vertexDegree=nccell.vertexDegree.data(end);

      lat=NaN*zeros(length(idxn),maxsize+1);
      lon=NaN*zeros(length(idxn),maxsize+1);

      latitude=[];
      longitude=[];

      for i=1:length(idxn)

       xsdf=nccell.cellsOnVertex.data(1:vertexDegree,idxn(i));

       if all(xsdf>0) % triangle falls within mesh

        lat(i,1:length(xsdf))=nccell.latitude.data(xsdf)*180/pi;
        lon(i,1:length(xsdf))=nccell.longitude.data(xsdf)*180/pi;
        lat(i,length(xsdf)+1:maxsize+1)=lat(i,1);
        lon(i,length(xsdf)+1:maxsize+1)=lon(i,1);

       else % coastal-edge vertex 

        xsdf=xsdf(xsdf>0);
        lxsdf=length(xsdf); 
        lat(i,1:lxsdf)=nccell.latitude.data(xsdf)*180/pi;
        lon(i,1:lxsdf)=nccell.longitude.data(xsdf)*180/pi;

        % find coastal edges
        edges=nccell.edgesOnVertex.data(:,idxn(i));
        edges=edges(edges>0);
        for ie=1:length(edges)
         cells=nccell.cellsOnEdge.data(:,edges(ie));
         if ~any(cells==0);  
          edges(ie)=0;
         else
          linkcell(ie)=cells(cells>0);
         end
        end
        exd=find(edges>0);
        linkcell=linkcell(exd);
        edges=edges(exd);

        % make sure cell links from its center to its edge,
        % not to the edge of the adjacent cell.
        if linkcell(2)==xsdf(1)
         edges=edges(end:-1:1);
        end

        % add in mid-edges and vertex to the polygon
        if length(xsdf)<3 & length(edges)==2
         lat(i,lxsdf+1)=nccell.latEdge.data(edges(2))*180/pi;
         lon(i,lxsdf+1)=nccell.lonEdge.data(edges(2))*180/pi;
         lat(i,lxsdf+2)=nccell.latVertex.data(idxn(i))*180/pi;
         lon(i,lxsdf+2)=nccell.lonVertex.data(idxn(i))*180/pi;
         lat(i,lxsdf+3)=nccell.latEdge.data(edges(1))*180/pi;
         lon(i,lxsdf+3)=nccell.lonEdge.data(edges(1))*180/pi;
         lat(i,lxsdf+4:maxsize+1)=lat(i,1);
         lon(i,lxsdf+4:maxsize+1)=lon(i,1);
        else
         error('Too many cell centers/edges for coastal edge')
        end
       
       end

       latitude=[latitude NaN lat(i,:)];
       longitude=[longitude NaN lon(i,:)];

      end

      % translate las and longs in to cartesian coords
      [cc,dd] = mfwdtran(gcm,latitude,longitude);

      % draw one long continuous patch
      plot(cc,dd,'Color',col,'Linewidth',0.1)

 clear cc dd lon lat idxn idx

end

% force the thing to draw
drawnow

if nargout>0

 % grid lines
 STL=geoshape(latitude,longitude,...
             'MPAS_SeaIce','Triangulation',...
             'Geometry','line');
end

% debug stuff
if debug; disp(['...Leaving ',mfilename]); end

