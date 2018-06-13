
clear
clf

datahome='~/data/SATELLITE/processed';

%dataset='nsidc0051_north_1978_2010';
%dataset='nsidc0081_conc_north_2011_2012';
%dataset='NPS_CDR_conc_1978_2010';
%dataset='NPS_CDR_conc_1978_2010_monthly';

dataset='G02202_v3_merged_conc_north_1979_2017';

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

