
clear
clf

startyear=1980;
endyear=2014;
%endyear=1980;
%monthset=[1 2 3];
monthset=[4 5 6];
%monthset=[7 8 9];
%monthset=[10 11 12];

generate=true;
%generate=false;

nruns=1;

daynames={'v2.NARRM.historical_0101.mpassi.hist.am.timeSeriesStatsDaily.'};
outnames={'v2.NARRM.historical_0101.mpassi.daily'};
locations={'/Users/afroberts/data/MODEL/E3SM/v2/v2.NARRM.historical_0101/data/ice/hist'};
outlocation='/Users/afroberts/work';

mons={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

fields={'timeDaily_avg_iceAreaCell',...
        'timeDaily_avg_iceVolumeCell',...
        'timeDaily_avg_snowVolumeCell',...
        'timeDaily_avg_uVelocityGeo',...
        'timeDaily_avg_vVelocityGeo'};

for nrun=1:nruns

 % first construct timeseries of only relevant months

 monthtag=['months'];
 monthnum=['months'];
 titlenc=[char(outnames{nrun}),', months:'];
 for months=monthset;
  monthtag=[monthtag,'_',char(mons{months})];
  monthnum=[monthnum,'_',num2str(months,'%2.2i')];
  titlenc=[titlenc,' ',char(mons{months})];
 end
 yeartag=['years_',num2str(startyear,'%4.4i'),'_',num2str(endyear,'%4.4i')];
 titlenc=[titlenc,', years: ',num2str(startyear,'%4.4i'),' to ',num2str(endyear,'%4.4i')];

 dayname=char(daynames{nrun});
 inlocation=char(locations{nrun});
 outfile=[char(outnames{nrun}),'.',monthnum,'.',yeartag];

 if generate

  k=1;

  for year=startyear:endyear

   for month=monthset;

    cd(inlocation)

    infile=[dayname,num2str(year,'%4.4i'),'-',num2str(month,'%2.2i'),'-01.nc'];

    nc=ridgepack_clone(infile,fields);
    nc=rmfield(nc,'attributes'); 
    nc.attributes.title=titlenc;
    nc.time.data=datenum(year,month,nc.time.data);
    nc.time.units=['days since ',num2str(startyear,'%4.4i'),'-01-01 00:00:00.0'];
    nc.time.calendar='noleap';

    nc.timeDaily_avg_Speed.units='m/s';
    nc.timeDaily_avg_Speed.type='NC_FLOAT';
    nc.timeDaily_avg_Speed.long_name='Daily Mean Sea Ice Speed';
    nc.timeDaily_avg_Speed.dimension=nc.timeDaily_avg_uVelocityGeo.dimension;
    nc.timeDaily_avg_Speed.data=sqrt(nc.timeDaily_avg_uVelocityGeo.data.^2 + ...
                                     nc.timeDaily_avg_vVelocityGeo.data.^2);

    nc=ridgepack_struct(nc);

    cd(outlocation)

    if k==1
     ridgepack_write(nc,outfile)
    else
     ridgepack_write(nc,outfile,{'time'},{0})
    end

    k=k+length(nc.time.data);

    clear nc

   end

  end

 end

 % next create means with additional statistics from timeseries

 cd(outlocation)
  
 nc=ridgepack_reduce(ridgepack_clone(outfile),{'time'},{},true);

 nc=rmfield(nc,'attributes');
 titlenc=[char(outnames{nrun}),' mean, months:'];
 for months=monthset;
  titlenc=[titlenc,' ',char(mons{months})];
 end
 titlenc=[titlenc,', years: ',num2str(startyear,'%4.4i'),' to ',num2str(endyear,'%4.4i')];
 nc.attributes.title=titlenc;
 meanoutfile=[char(outnames{nrun}),'.mean.',monthnum,'.',yeartag];

 ridgepack_write(nc,meanoutfile);

end


