function [nc,gs]=ridgepack_e3smseasaw(ncvert,ncc,varc,threshold,...
                                       centlat,centlon,horizon)

% ridgepack_e3smseasaw - generate a threshold on an unstructured mesh
%
% function [nc]=ridgepack_e3smseasaw(ncc,varc,threshold,ncvert,centlat,centlon,horizon)
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

% First generate a coastline
cidx=[1:length(ncvert.nCells.data)]'; %coast
[clats,clons,cverts]=ridgepack_findE3SMedges(ncvert,cidx);

% If requested, also generate a threshold
if ~isempty(ncc)
 infill=true;
 cidx=find(ncc.(varc).data>threshold); 
 [tlats,tlons,tverts]=ridgepack_findE3SMedges(ncvert,cidx,infill);

 % Now remove coast from threshold polygons

 % Find NaN's seperating the vertices. idx will have a length
 % of n-1 where n is the number of closed loops.
 tvertsnan=[NaN tverts NaN];
 idx=find(isnan(tvertsnan));

 % Also find the number of coastal closed loops
 cvertsnan=[NaN cverts NaN];
 cdx=find(isnan(cvertsnan));

 % create new set of contours 
 xvertsnan=tvertsnan;
 xlatsnan=[NaN tlats NaN];
  xlonsnan=[NaN tlons NaN];

 nbreaks=0; % number of line breaks that need to be inserted
 breakvert=[]; % break of vertex index

 for i=1:length(idx)-1

  % vertices along a close loop
  vertcircleidx=[idx(i)+1:idx(i+1)-1 idx(i)+2:idx(i+1)-1];
  vertcircle=tvertsnan(vertcircleidx);

  for k=1:length(cdx)-1

   coastcirclecdx=[cdx(k)+1:cdx(k+1)-1 cdx(k)+2:cdx(k)+3];
   coastcircle=cvertsnan(coastcirclecdx);

   for j=1:length(vertcircle)-2

    if ~isempty(strfind(coastcircle,vertcircle(j:j+1))) | ...
       ~isempty(strfind(coastcircle,vertcircle(j+1:-1:j)))  
      nbreaks=nbreaks+1;
      breakvert(nbreaks)=vertcircleidx(j);
    end

   end
 
  end
 
 end

 % Now scan for shared sides between the new contours and 
 % the coastline, and place a NaN in between
 sortidx=sort(breakvert);
 for i=length(sortidx):-1:1
  xvertsnan=[xvertsnan(1:sortidx(i)) NaN xvertsnan(sortidx(i)+1:end)];
  xlatsnan=[xlatsnan(1:sortidx(i)) NaN xlatsnan(sortidx(i)+1:end)];
  xlonsnan=[xlonsnan(1:sortidx(i)) NaN xlonsnan(sortidx(i)+1:end)];
 end

 % Now remove sequences of [NaN Number NaN] and replace with NaN
 % Also removing the NaN that padded the start and finish of the 
 % sequence in the previous search, and sequences of NaNs.

 xverts=xvertsnan;
 xlats=xlatsnan;
 xlons=xlonsnan;

 % First, fill in single vertices with NaNs
 for i=2:length(xverts)-1
  if (isnan(xverts(i-1)) & ...
     ~isnan(xverts(i)) & ...
      isnan(xverts(i+1))) 
   xvertsnan(i)=NaN;
   xlatsnan(i)=NaN;
   xlonsnan(i)=NaN;
  end
 end

 xverts=[];
 xlats=[];
 xlons=[];

 % Now, replace repeated NaNs with a single one
 k=0;
 for i=2:length(xvertsnan)
  if ~(isnan(xvertsnan(i-1)) & isnan(xvertsnan(i)))
   k=k+1;
   xverts(k)=xvertsnan(i);
   xlats(k)=xlatsnan(i);
   xlons(k)=xlonsnan(i);
  end
 end

 % quick error check
 for i=2:length(xverts)
  if (isnan(xverts(i-1)) & isnan(xverts(i)))
   error('NaNs repeated')
  end
 end

end

% create netcdf structure
nc.attributes.title='E3SM-MPAS Edge Definition';

nc.attributes.coastal_segments=length(isnan(cverts))+1;

nc.cnpoints.data=[1:length(cverts)];
nc.cnpoints.long_name='number of points on coastal outline';
nc.cnpoints.dimension={'cnpoints'};

nc.clatitude.data=clats;
nc.clatitude.long_name='latitude of coast vertices';
nc.clatitude.dimension={'cnpoints'};

nc.clongitude.data=clons;
nc.clongitude.long_name='longitude of coast vertices';
nc.clongitude.dimension={'cnpoints'};

nc.cvertices.data=cverts;
nc.cvertices.long_name='MPAS vertices on coast';
nc.cvertices.dimension={'cnpoints'};

if isempty(ncc)

 nc.attributes.threshold_segments=length(isnan(tverts))+1;

 nc.npoints.data=[1:length(tverts)];
 nc.npoints.long_name='number of points on threshold shapes';
 nc.npoints.dimension={'npoints'};

 nc.latitude.data=tlats;
 nc.latitude.long_name='latitude of threshold shapes';
 nc.latitude.dimension={'npoints'};

 nc.longitude.data=tlons;
 nc.longitude.long_name='longitude of threshold shapes';
 nc.longitude.dimension={'npoints'};

 nc.vertices.data=tverts;
 nc.vertices.long_name='vertex indices on threshold shapes';
 nc.vertices.dimension={'npoints'};

 nc.attributes.contour_segments=length(isnan(xverts))+1;

 nc.xnpoints.data=[1:length(xverts)];
 nc.xnpoints.long_name='number of points on threshold contours';
 nc.xnpoints.dimension={'xnpoints'};

 nc.xlatitude.data=xlats;
 nc.xlatitude.long_name='latitude of threshold contours';
 nc.xlatitude.dimension={'xnpoints'};

 nc.xlongitude.data=xlons;
 nc.xlongitude.long_name='longitude of threshold contours';
 nc.xlongitude.dimension={'xnpoints'};

 nc.xvertices.data=xverts;
 nc.xvertices.long_name='vertex indices of contours';
 nc.xvertices.dimension={'xnpoints'};

end

% plot data if no output argument is specified
nc=ridgepack_struct(nc);


% create geoshapes
S=geoshape(nc.latitude.data,nc.longitude.data,...
            'MPAS_Threshold','Sea Ice Threshold');

SP=geoshape(nc.latitude.data,nc.longitude.data,...
            'MPAS_Threshold','Sea Ice Threshold',...
            'Geometry','polygon');

STP=geoshape(nc.latitude.data,nc.longitude.data,...
            'MPAS-Seaice',[varc,'>',num2str(threshold)],...
            'Geometry','polygon');

SCL=geoshape(nc.clatitude.data,nc.clongitude.data,...
            'MPAS Coast','Sea Ice Threshold',...
            'Geometry','line');

SCP=geoshape(nc.clatitude.data,nc.clongitude.data,...
            'MPAS-Ocean','Coast',...
            'Geometry','polygon');



