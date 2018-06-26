% ridgepack_ssmi_to_model - regrid SSM/I concentration to RASM domain
%
% This script regrids SSM/I-derived sea ice concentration on a polar stereographic
% grid to the RASM ice-ocean domain grid, and writes the output in netcdf.
% In order to use this script, you must first peel off the desired field
% from the NOAA CDR dataset using the ridgepack_cdr_ssmi_to_nc script.
%
% Andrew Roberts, Naval Postgraduate School, June 2018 (afrobert@nps.edu)

clear
clf

% set data home and set name here
datahome='~/data/SATELLITE/processed';
dataset='G02202_v3_merged_conc_north_1979_2017';

% specific output file name
outfile=[dataset,'_RASM_CICE'];

cd(datahome)

nctime=ridgepack_clone(dataset,{'time'});

for k=1:length(nctime.time.data)

 ncold=ridgepack_clone(dataset,{'conc'},k)

 ncnew=ridgepack_regrid(ncold,'conc','',7)
 
 ncnew.attributes.title=[ncold.attributes.title,' on RASM POP/CICE mesh'];

 ncnew.conc.data(ncnew.conc.data<10)=0;

 if k==1
  ridgepack_write(ncnew,outfile)
 else
  ridgepack_write(ncnew,outfile,{'time'},{0})
 end

end

