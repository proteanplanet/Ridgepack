% ridgepack_splice_time_bounds - splice time bounds into RASM-interpolated SSM/I concentration
%
% This script is designed to add a time_bounds field into an existing
% sea ice concentration netcdf file interpolated to the RASM grid
% to enable correct calculation of seasonal means from monthly model output,
% weighted according to the lengths of the months.  It also enables other
% statistics to be calculated.
%
% Andrew Roberts, Naval Postgraduate School, June 2018 (afrobert@nps.edu)

% This script splices in a time_bounds function to the SSM/I dataset
datadir='/Users/aroberts/data/SATELLITE/processed';
datafile='G02202_v3_merged_conc_north_1979_2017_RASM_CICE';
%datagrab={'conc','polehole'};
datagrab={'conc'};

% add time bounds to output file name
outfile=[datafile,'_time_bounds'];

cd(datadir)

% read in data file time array
nc=ridgepack_clone(datafile,'time');

% set the cross-over from SMMR to SSM/I
smmrend=datenum('1987-08-01')

nc.d2.data=[1 2];
nc.d2.type='NC_INT';
nc.d2.dimension={'d2'}

nc.time_bounds.dimension={'d2','time'};
nc.time_bounds.calendar=nc.time.calendar;
nc.time_bounds.long_name='Time bounds for each sample';
nc.time_bounds.type='NC_DOUBLE';
nc.time_bounds.units=nc.time.units;


% add time bounds of 1-day interval for SSM/I and 2-day
% for SMMR, and then write it to output
for i=1:length(nc.time.data);

 ncr=ridgepack_clone(datafile,datagrab,{'time'},{i});

 ncr.d2=nc.d2;
 ncr.time_bounds=nc.time_bounds;

 if nc.time.data(i)<smmrend

  ncr.time_bounds.data(1,1)=nc.time.data(i)-1.0;
  ncr.time_bounds.data(2,1)=nc.time.data(i)+1.0;

 else

  ncr.time_bounds.data(1,1)=nc.time.data(i)-0.5;
  ncr.time_bounds.data(2,1)=nc.time.data(i)+0.5;

 end

 if i==1
  ridgepack_write(ncr,outfile)
 else
  ridgepack_write(ncr,outfile,{'time'},{0})
 end

end

