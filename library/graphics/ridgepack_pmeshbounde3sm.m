function [nc]=ridgepack_pmeshbounde3sm(ncvert,ncc,varc,threshold,...
                                       centlat,centlon,horizon)

% ridgepack_pthresholde3sm - generate a threshold on an unstructured mesh
%
% function [nc]=ridgepack_pthresholde3sm(ncc,varc,threshold,ncvert,centlat,centlon,horizon)
%
% This function .... description
%  
% INPUT:
%
% ncc       - netcdf structure including
% varc      - variable in ncc to be used
% threshold - threshold value used as a boundary on unstructured mesh
% ncvert    - vertices on the unstructed mesh
% centlat   - center latitude of plotted satellite view (optional)
% centlon   - center longitude if plotted in satellite view (optional)
% horizon   - horizon, in degrees of satellite view (optional)
%
% Note that centlat, centlon and horizon should not be provided
% if nc is specificied.
%
% 
% OUTPUT:
% 
% nc - netcdf structure with lats, longs, mesh cell and vertex
%      indices
%
% Ridgepack Version 1.1
% Andrew Roberts, Los Alamos National Laboratory, April 3, 2020 (afroberts@lanl.gov)


if nargin<4
 error('There is no threshold value set, and possibly more')
end

if nargin<5
 centlat=90;
elseif nargin>3 & nargout==1
 error('you are asking to both plot and produce coastal data')
elseif ~isnumeric(centlat)
 error('centlat should be a real number between -90 and 90')
elseif centlat>90 | centlat<-90
 error('centlat should be between -90 and 90')
end

if nargin<6
 centlon=0;
elseif ~isnumeric(centlon)
 error('centlon should be a real number between -90 and 90')
elseif centlat>180 | centlat<-180
 error('centlon should be between -180 and 180')
end

if nargin<7
 horizon=90;
elseif ~isnumeric(horizon)
 error('horizon must be a real number between 0 and 90')
elseif horizon<0 | horizon>90
 error('horizon must be between 0 and 90')
end

% first generate a coastline, then generate a threshold
% line if needed
for ctype=1:2

 % get indices of cells that match the criteria
 if ctype==1
  cidx=[1:length(ncvert.nCells.data)]'; %coast
  [clats,clons,cverts]=ridgepack_findE3SMedges(ncvert,cidx);
 else
  cidx=find(ncc.(varc).data>threshold); %threshold
  [tlats,tlons,tverts]=ridgepack_findE3SMedges(ncvert,cidx);
 end

end % coasts

% Now remove coast from threshold polygons

% Find NaN's seperating the vertices. idx will have a length
% of n-2 where n is the number of closed loops.
idx=find(isnan(tverts));

for i



end





% create netcdf structure

nc.attributes.title='Threshold Edge Definition';

nc.cnpoints.data=[1:length(cverts)];
nc.cnpoints.long_name='number of points on threshold';
nc.cnpoints.dimension={'npoints'};

nc.clatitude.data=clats;
nc.clatitude.long_name='latitude of threshold';
nc.clatitude.dimension={'npoints'};

nc.clongitude.data=clons;
nc.clongitude.long_name='longitude of threshold';
nc.clongitude.dimension={'npoints'};

nc.cvertices.data=cverts;
nc.cvertices.long_name='vertex indices';
nc.cvertices.dimension={'npoints'};

nc.npoints.data=[1:length(tverts)];
nc.npoints.long_name='number of points on threshold';
nc.npoints.dimension={'npoints'};

nc.latitude.data=tlats;
nc.latitude.long_name='latitude of threshold';
nc.latitude.dimension={'npoints'};

nc.longitude.data=tlons;
nc.longitude.long_name='longitude of threshold';
nc.longitude.dimension={'npoints'};

nc.vertices.data=tverts;
nc.vertices.long_name='vertex indices';
nc.vertices.dimension={'npoints'};

% plot data if no output argument is specified
if nargout==0

 [x,y,z,phi,theta]=ridgepack_satfwd(nc.latitude.data,...
                                    nc.longitude.data,...
                                    centlat,centlon,...
                                    2*horizon,1.0001);
 plot3(x,y,z,'Color','m',...
       'LineWidth',0.5-sin(deg2rad(horizon))*0.4)

 hold on

 [x,y,z,phi,theta]=ridgepack_satfwd(nc.clatitude.data,...
                                    nc.clongitude.data,...
                                    centlat,centlon,...
                                    2*horizon,1.0001);
 plot3(x,y,z,'Color','k',...
       'LineWidth',0.5-sin(deg2rad(horizon))*0.4)

 axis off

 axis equal
 view([0 0 0.4])
 axis tight

else

 nc=ridgepack_struct(nc);

end

return

 % mask out areas not in the threshold plotting area
 if nargout==0
  disp('Reducing to just the threshold plotting area')
  maxth=deg2rad(horizon);
  for i=cidx'

   maxidx=ncvert.nEdgesOnCell.data(i);

   la=ncvert.clatitude.data(ncvert.verticesOnCell.data(1:maxidx,i));
   lo=ncvert.clongitude.data(ncvert.verticesOnCell.data(1:maxidx,i));

   [x,y,z,ph,th]=ridgepack_satfwd(rad2deg(la),rad2deg(lo),...
                                  centlat,centlon,horizon,1);

   % filter cells not in frame, and find cropping limit
   if all(isnan(x))
    ncvert.clatitude.data(ncvert.verticesOnCell.data(1:maxidx,i))=NaN;
    ncvert.clongitude.data(ncvert.verticesOnCell.data(1:maxidx,i))=NaN;
   end

  end
 end

