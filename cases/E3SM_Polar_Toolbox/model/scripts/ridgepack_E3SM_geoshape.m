
%regenerate=true;
regenerate=false;

if regenerate

clf
clear

centlat=80; % degrees north
centlon=-150; % degrees east
horizon=60; % degrees of satellite horizon (0-90)

% location of grid file
gridloc='/Users/afroberts/data/MODEL/E3SM/DECK/grid';
gridfile='E3SM_LR_V1_grid.nc';
cgrid=true;

% location of sea ice data
dataloc='/Users/afroberts/data/MODEL/E3SM/DECK/monthly/h1/archive/ice/reduced';
datafile='mpascice.hist.am.timeSeriesStatsMonthly.1980-03-01.nc';
varc='timeMonthly_avg_iceAreaCell';
datatitle='Sea Ice Extent';
threshold=0.15;

% plot location
plotloc='/Users/afroberts/work';

% %%%%%%%%%%%%%%% CHANGE ABOVE THIS LINE %%%%%%%%%%%%%%%%%% %

% obtain grid information (vertices, cell centers)
cd(gridloc)
ncvert=ridgepack_clone(gridfile,{'latVertex','lonVertex','dcEdge',...
                                 'verticesOnCell','indexToCellID',...
                                 'nEdgesOnCell','edgesOnCell',...
                                 'cellsOnEdge'});

nccell=ridgepack_clone(gridfile,{'latCell','lonCell'});

% obtain data
cd(dataloc)
ncc=ridgepack_clone(datafile,varc);

% put up satview
%ridgepack_satview(centlat,centlon,horizon,1,2);

% plot ice edge 
if 1==0
 ridgepack_pthresholde3sm(ncc,varc,threshold,ncvert,...
                                centlat,centlon,horizon);
else
 nc=ridgepack_pthresholde3sm(ncc,varc,threshold,ncvert);
end

end

% find line segments - they are seperated by NaNs
idx=find(isnan(nc.vertices.data));

% account for case where first is a NaN
if idx(1)==1
 idx=idx(2:end);
end

% break into line segments and number them
segment=zeros(size(nc.vertices.data));
segment(1:idx(1))=1;
k=2;
for i=[idx(1)+1:length(nc.vertices.data)];
 segment(i)=k;
 if ~isempty(find(idx(:)==i)) 
  k=k+1;
 end
end

% search check
segmentfound=zeros(size(nc.vertices.data));
segmentfound(1:idx(1))=1;

% now hunt for vertex not in the segment

vert=nc.vertices.data(2:idx(1)-1);

startidx=2;
startvertex=nc.vertices.data(startidx);

%for i=idx

for i=idx(1)

 endidx=i-1;
 endvertex=nc.vertices.data(endidx);

 jdx=find(nc.vertices.data==endvertex);
 jdx=jdx(jdx~=endidx);

 kdx=find(nc.vertices.data==startvertex);
 kdx=kdx(kdx~=startidx);

 while ~isempty(jdx) & ~isempty(kdx)

  % check for intersections
  if length(jdx)>1
   error('Found intersection that should not exist: 1')
  elseif length(kdx)>1
   error('Found intersection that should not exist: 2')
  end

  % case of sequential-order vertices
  if isnan(nc.vertices.data(jdx+1));

   disp('sequential')

   startex=find(isnan(nc.vertices.data(1:jdx)));
   vert=[vert nc.vertices.data(jdx:-1:startex(end)+1)];
   segmentfound(startex(end):jdx+1)=1;
   endidx=startex(end)+1;
   endvertex=vert(end);

  % case of reversed-order vertices
  elseif isnan(nc.vertices.data(jdx-1));

   disp('reversed')

   endx=find(isnan(nc.vertices.data(jdx:end)));
   vert=[vert nc.vertices.data(jdx+1:jdx+endx(1)-2)];
   segmentfound(jdx+1:jdx+endx(1)-2)=1;
   endidx=jdx+endx(1)-2;
   endvertex=vert(end);

  % have found the middle of a segment, which 
  % theoretically should not exist
  else
   error('Found intersection that should not exist: 3')

  end
 
  jdx=find(nc.vertices.data==endvertex);
  jdx=jdx(jdx~=endidx & segmentfound(jdx)~=1);

 end

 if isempty(jdx) & ~isempty(kdx)
  error('Cant find matching start or end index: 4')
 elseif isempty(jdx) & isempty(kdx)
  error('Cant find matching start or end index: 5')
 end

end
 
 
% add title
title([datatitle])

cd(plotloc)
ridgepack_fprint('png','Outfile_Sea_Ice_GeoTif',1,1)




