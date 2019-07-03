function [ncnew]=ridgepack_thintrack(nc,resolution,switchd)


if switchd
 alms=almanac('earth','sphere');
else
 alms=almanac('earth','wgs84');
end

ncnew=nc;
ncnew.points.data=1;
ncnew.latitude.data=nc.latitude.data(1);
ncnew.longitude.data=nc.longitude.data(1);

for i=2:length(nc.latitude.data)
 for k=i:length(nc.latitude.data)
  dist=distance(nc.latitude.data(i-1),nc.longitude.data(i-1),...
                nc.latitude.data(k),nc.longitude.data(k),alms);
  if dist>resolution; break; end
 end
 ncnew.points=
 
 
end












