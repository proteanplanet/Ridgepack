
clear
clf

highres=false;
%highres=true;

stateset=true;
%stateset=false;
%fluxset=true;
fluxset=false;

if highres

 startyear=26;
 endyear=55;
 cd('/Users/afroberts/SIhMSatArray/E3SM/highres/monthly')
 resname='HR_v1'

else

 startyear=26;
 endyear=55;
 cd('/Users/afroberts/SIhMSatArray/E3SM/lrhrequiv/monthly')
 resname='LR_v1'

end

if stateset

 fields={'timeMonthly_avg_iceAreaCell',...
        'timeMonthly_avg_iceVolumeCell',...
        'timeMonthly_avg_snowVolumeCell',...
        'timeMonthly_avg_uVelocityGeo',...
        'timeMonthly_avg_vVelocityGeo',...
        'timeMonthly_avg_uOceanVelocityVertexGeo',...
        'timeMonthly_avg_vOceanVelocityVertexGeo',...
        'timeMonthly_avg_divergence',...
        'timeMonthly_avg_uAirVelocity',...
        'timeMonthly_avg_vAirVelocity'};
 
 setname='state';

else

 fields={'timeMonthly_avg_latentHeatFluxInitialArea',...
        'timeMonthly_avg_sensibleHeatFluxInitialArea',...
        'timeMonthly_avg_longwaveUpInitialArea',...
        'timeMonthly_avg_shortwaveDown'};

 setname='fluxes';

end

%for month=[3,9];
for month=[1,2,3,4,5,6,7,8,9,10,11,12];

 k=0

 for year=startyear:endyear
  k=k+1;
  dataset=['mpascice.hist.am.timeSeriesStatsMonthly.',num2str(year,'%4.4i'),'-',...
           num2str(month,'%2.2i'),'-01.nc'];
  if k==1
   nc=ridgepack_reduce(ridgepack_clone(dataset,fields),{'time'});
  else
   ncnew=ridgepack_reduce(ridgepack_clone(dataset,fields),{'time'});
   for i=1:length(fields)
    nc.(char(fields{i})).data=nc.(char(fields{i})).data+...
                              ncnew.(char(fields{i})).data;
   end
   clear ncnew
  end
 end

 for i=1:length(fields)
  nc.(char(fields{i})).data=nc.(char(fields{i})).data./k;
 end

 nc=rmfield(nc,'attributes');

 nc.attributes.title=['Mean years ',num2str(startyear),' to ',num2str(endyear)];

 outfilename=['mpascice.hist.am.',resname,'.',setname,...
              '.mean.',num2str(month,'%2.2i'),'.',...
               num2str(startyear,'%4.4i'),'_',num2str(endyear,'%4.4i')];

 ridgepack_write(nc,outfilename);

end


