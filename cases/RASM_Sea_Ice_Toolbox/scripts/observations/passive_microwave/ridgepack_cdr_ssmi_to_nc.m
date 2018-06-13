% ridgepack_cdr_ssmi_to_nc - Generate SSM/I CDR Concentration NetCDF file


clear

%hemisphere='north';
hemisphere='south';

%datahome='~/data/SATELLITE/downloaded/G02202_v2';
datahome='/Volumes/Promise Pegasus/data2/G02202_V3';

dataset='goddard_merged_seaice_conc';

firstyear=1979;
latestyear=2017;

% Set up basic nc structure, using a sample netcdf timeslice from the dataset
if strcmp(hemisphere,'north')
 ncs=ridgepack_clone([datahome,'/',hemisphere,'/daily/2016/seaice_conc_daily_nh_f17_20161231_v03r01.nc']);
elseif strcmp(hemisphere,'south')
 ncs=ridgepack_clone([datahome,'/',hemisphere,'/daily/2016/seaice_conc_daily_sh_f17_20161231_v03r01.nc']);
else
 error('Hemisphere specification is incorrect')
end

nc.attributes.title='CDR Merged SMMR and SSM/I sea ice concentration from G02202_V3';

nc.x=ncs.x;
nc.y=ncs.y;

nc.latitude=ncs.latitude;
nc.longitude=ncs.longitude;

nc.mask.data=squeeze(ncs.goddard_merged_seaice_conc.data*100);
nc.mask.data=~(nc.mask.data==-2 | nc.mask.data==-3 | nc.mask.data==-4);
nc.mask.long_name='land mask';
nc.mask.dimension={'y' 'x'};
nc.mask.type='NC_FLOAT';

nc.time.long_name='time centered at 1200Z on representative day';
nc.time.dimension={'time'};
nc.time.data=[];
nc.time.type='NC_DOUBLE';

nc.polehole.long_name='pole hole mask';
nc.polehole.dimension={'time' 'y' 'x'};
nc.polehole.type='NC_FLOAT';

nc.conc.long_name='CDR Merged SMMR and SSM/I sea ice concentration from G02202_V3';
nc.conc.dimension={'time' 'y' 'x'};
nc.conc.units='%';
nc.conc.type='NC_FLOAT';

clear ncs;

numm=0;

recycle('off')
delete([datahome,'/error.txt'])

for k=firstyear:1:latestyear

 dday=datenum([k 1 1 12 0 0.0]);
 ddayend=datenum([k 12 31 12 0 0.0]);

 m=0;

 for time=dday:1:ddayend

  cd([datahome,'/',hemisphere,'/daily/',num2str(k)]);

  file=dir(['seaice_conc_daily_*_',datestr(time,'yyyymmdd'),'_*.nc']);
  if length(file)>1
    warning(['More than one file in directory for ',datestr(time,'yyyy-mm-dd')])
  end

  if length(file)>0

   disp(['Processing data for ',datestr(time,'yyyymmdd')]);

   try
    ncs=ridgepack_clone(file(1).name,{'goddard_merged_seaice_conc','longitude','latitude'});
   catch
    fid=fopen([datahome,'/error.txt'],'a+');
    fprintf(fid,'%s\n',['Time error for ',datestr(time)]) ;
    fclose(fid);
    continue
   end

   nc.time.data=ncs.time.data+datenum([0 0 0 12 0 0]);

   nc.conc.data=ncs.goddard_merged_seaice_conc.data*100;

   if ~all(isnan(nc.conc.data(:)))

    m=m+1;

    [nc,out]=ridgepack_struct(nc);

    if m==1; nc.mask.data=~(squeeze(nc.conc.data)<0 & ...
                            squeeze(nc.conc.data)>-4.9); end

    nc.polehole.data=(nc.conc.data<-4.9);

    nc.conc.data(nc.conc.data<0)=NaN;

    cd ~/data/SATELLITE/processed

    numm=numm+1;

    outfile=['G02202_v3_merged_conc_',hemisphere,'_',...
             num2str(firstyear),'_',num2str(latestyear)];
    if numm==1;
     ridgepack_write(nc,outfile);
    else
     ridgepack_write(nc,outfile,{'time'},{numm});
    end
   
   end

   clear ncs

  else

   disp(['No data found for ',datestr(time,'yyyymmdd')]);

  end

 end

end

