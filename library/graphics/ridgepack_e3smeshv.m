function [STL]=ridgepack_e3smeshv(nccell,mask,col,fill)

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
% mask   - mask of vertex indices to be plotted, given as a vector of 
%          cell indices from nVertices that are to be plotted.
%          [optional] 
% col    - This sets the color of the lines [optional]
% fill   - Logical set to true if the patch is to be filled 
%          The default is not to fill. [optional]
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
if nargin>=3 & isempty(col) | nargin<3
 col=[0.7 0.2 0];
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
  lidx=find(nccell.latitude.data>maxlat | nccell.latitude.data<minlat);
  nccell.longitude.data(lidx)=NaN;
  nccell.latitude.data(lidx)=NaN;
 end
 figout=get(hmap,'OuterPosition');
 fontsize=min(11,max(9,10*sqrt((figout(3).^2)+(figout(4).^2))/sqrt(2)));
 set(hmap,'fontsize',fontsize);
end

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

      end
 
 clear cc dd lon lat idxn idx

end

if nargout>0

 % grid lines
 STL=geoshape(latitude,longitude,...
             'MPAS_SeaIce','Triangulation',...
             'Geometry','line');
else

 drawnow

end

% debug stuff
if debug; disp(['...Leaving ',mfilename]); end

