function [ncnew]=ridgepack_thintrack(nc,resolution,switchd)

if switchd
 alms=almanac('earth','sphere');
else
 alms=almanac('earth','wgs84');
end

i=1; % starting index of old dataset
j=1; % starting index of new dataset

ncnew=nc;
ncnew.attributes.title=[nc.attributes.title,...
                       ' at ~',num2str(resolution),...
                       ' degree resolution'];
ncnew.points.data=j;
ncnew.latitude.data=nc.latitude.data(i);
ncnew.longitude.data=nc.longitude.data(i);
ncnew.distance=nc.latitude;
ncnew.distance.data=0;
ncnew.distance.long_name='distance along track';
ncnew.distance.units='km';

while i<length(nc.latitude.data)-1
 k=i+1;
 while k<length(nc.latitude.data)
  dist=distance(nc.latitude.data(i),nc.longitude.data(i),...
                nc.latitude.data(k),nc.longitude.data(k),alms);
  if dist>resolution; break; end
  k=k+1;
 end
 j=j+1;
 ncnew.points.data(j)=j;
 ncnew.latitude.data(j)=nc.latitude.data(k);
 ncnew.longitude.data(j)=nc.longitude.data(k);
 dist=distance(nc.latitude.data(i),nc.longitude.data(i),...
               nc.latitude.data(k),nc.longitude.data(k),alms);
 ncnew.distance.data(j)=ncnew.distance.data(j-1)+dist;
 i=k;
end

ncnew=ridgepack_struct(ncnew);











