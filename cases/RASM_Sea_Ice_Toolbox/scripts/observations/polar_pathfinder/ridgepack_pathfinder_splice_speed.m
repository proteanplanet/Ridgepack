
% This script splices in a time_bounds function to the SSM/I dataset
%datadir='/Users/aroberts/data/SATELLITE/processed';
datadir='/Volumes/RobertsRaid3/data/SATELLITE/processed';
datafile='Pathfinder_icemotion_monthly_1979_2016_v3_RASM_CICE_time_bounds';
datagrab={'u','v'};

outfile=[datafile,'_speed'];

cd(datadir)

ncu=ridgepack_clone(datafile,'u',1);
ncu.speed=ncu.u;
ncu.speed.long_name='drift speed';

nc=ridgepack_clone(datafile,'time');

for i=1:length(nc.time.data);

 ncr=ridgepack_clone(datafile,datagrab,{'time'},{i});

 ncr.speed=ncu.speed;
 ncr.speed.data=sqrt((ncr.u.data.^2) + (ncr.v.data.^2));

 if i==1
  ridgepack_write(ncr,outfile)
 else
  ridgepack_write(ncr,outfile,{'time'},{0})
 end

end

