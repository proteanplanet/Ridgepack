
clear


% This script splices in a time_bounds function to the SSM/I dataset
%datadir='/Users/aroberts/data/SATELLITE/processed';
datadir='/Volumes/RobertsRaid3/data/SATELLITE/processed';
datafile='Pathfinder_icemotion_monthly_1979_2016_v3_RASM_CICE';
datagrab={'u','v'};

outfile=[datafile,'_time_bounds']

cd(datadir)

nc=ridgepack_clone(datafile,'time');

nc.d2.data=[1 2];
nc.d2.type='NC_INT';
nc.d2.dimension={'d2'}

nc.time_bounds.dimension={'d2','time'};
nc.time_bounds.calendar=nc.time.calendar;
nc.time_bounds.long_name='Time bounds for each sample';
nc.time_bounds.type='NC_DOUBLE';
nc.time_bounds.units=nc.time.units;

for i=1:length(nc.time.data);

 ncr=ridgepack_clone(datafile,datagrab,{'time'},{i});

 ncr.d2=nc.d2;
 ncr.time_bounds=nc.time_bounds;

 datev=datevec(nc.time.data(i));

 datev(2)=datev(2)+1;

 ncr.time_bounds.data(1,1)=nc.time.data(i);
 ncr.time_bounds.data(2,1)=datenum(datev);

 if i==1
  ridgepack_write(ncr,outfile)
 else
  ridgepack_write(ncr,outfile,{'time'},{0})
 end

end

