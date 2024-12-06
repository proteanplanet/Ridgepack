
clear
clf

startyear=1980;
endyear=2000;

%startyear=2000;
%endyear=2020;

ensemblemember='0051';
%ensemblemember='0101';

%monthsets={[1 2 3]};
%monthsets={[4 5 6]};
%monthsets={[7 8 9]};
%monthsets={[10 11 12]};
%monthsets={[1 2 3],[4 5 6],[7 8 9],[10 11 12]};
monthsets={[4 5 6],[10 11 12]};
%monthsets={[1 2 3],[7 8 9]};

% generate timeseries before creating means
generate=true;
%generate=false;

% non-gridded means for native mesh (keep switched off if running out of memory)
%reducenative=true;
reducenative=false;

nruns=1;

daynames={['v3.LR.historical_',...
           ensemblemember,'.mpassi.hist.am.timeSeriesStatsDaily.']};

outnames={['v3.LR.historical_',...
           ensemblemember,'.mpassi.daily']};

locations={['/Users/afroberts/data/MODEL/E3SM/v3/v3.LR.historical_',...
            ensemblemember,'/hist']};

gridnames={'E3SM_IcoswISC30E3r5.nc'};

gridlocations={'/Users/afroberts/data/MODEL/E3SM/v3/v3.LR.spinup/grid'};

cdrnorthmask='cdr_25km_north_area_mask_sectors.nc';
cdrsouthmask='cdr_25km_south_area_mask_sectors.nc';
cdrlocation='/Users/afroberts/data/data/SATELLITE/processed/G02202_v4';

outlocations={['/Users/afroberts/data/MODEL/E3SM/v3/v3.LR.historical_',...
             ensemblemember,'/processed']};

mons={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

fields={'timeDaily_avg_iceAreaCell',...
        'timeDaily_avg_iceVolumeCell',...
        'timeDaily_avg_snowVolumeCell',...
        'timeDaily_avg_uVelocityGeo',...
        'timeDaily_avg_vVelocityGeo'};

% step through individual simulations

for nrun=1:nruns

 for nm=1:length(monthsets)

  monthset=monthsets{nm};

  % set up tag names and data file names

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
  outfileregridnh=[char(outnames{nrun}),'.nh.regrid.',monthnum,'.',yeartag];
  outfileregridsh=[char(outnames{nrun}),'.sh.regrid.',monthnum,'.',yeartag];
  gridfile=[char(gridlocations{nrun}),'/',char(gridnames{nrun})];

  if generate

   % generate polar stereographic interpolation grid

   [ncnh]=ridgepack_gridgen('',8);
   xnh=ncnh.x.data;
   ynh=ncnh.y.data;
   [xnhm,ynhm]=meshgrid(xnh,ynh);

   [ncsh]=ridgepack_gridgen('',9);
   xsh=ncsh.x.data;
   ysh=ncsh.y.data;
   [xshm,yshm]=meshgrid(xsh,ysh);


   % read in land-sea mask from CDR dataset

   cd(cdrlocation)

   ncnhcdr=ridgepack_clone(cdrnorthmask);
   ncnh.mask=ncnhcdr.mask;

   ncshcdr=ridgepack_clone(cdrsouthmask);
   ncsh.mask=ncshcdr.mask;


   % set up model grid and model polar stereographic projection

   nclatv=ridgepack_clone(gridfile,'latVertex');
   nclonv=ridgepack_clone(gridfile,'lonVertex');
   nclatc=ridgepack_clone(gridfile,'latCell');
   nclonc=ridgepack_clone(gridfile,'lonCell');

   % project onto a polar stereographic grid for the north 
   maskvnh=find(nclatv.latVertex.data>30*pi/180);
   [Xvnh,Yvnh]=ridgepack_geodetictoxy(nclatv.latVertex.data(maskvnh)*180/pi,...
                                    nclonv.lonVertex.data(maskvnh)*180/pi+45,1);
   Xvnh(Xvnh<min(xnh) | Xvnh>max(xnh))=NaN;
   Yvnh(Yvnh<min(ynh) | Yvnh>max(ynh))=NaN;
   maskUVnh=(~isnan(Xvnh) & ~isnan(Yvnh));
   Xvnh=Xvnh(maskUVnh); Yvnh=Yvnh(maskUVnh);

   maskcnh=find(nclatc.latCell.data>30*pi/180);
   [Xcnh,Ycnh]=ridgepack_geodetictoxy(nclatc.latCell.data(maskcnh)*180/pi,...
                                    nclonc.lonCell.data(maskcnh)*180/pi+45,1);
   Xcnh(Xcnh<min(xnh) | Xcnh>max(xnh))=NaN;
   Ycnh(Ycnh<min(ynh) | Ycnh>max(ynh))=NaN;
   maskCCnh=(~isnan(Xcnh) & ~isnan(Ycnh));
   Xcnh=Xcnh(maskCCnh); Ycnh=Ycnh(maskCCnh);

   % project onto a polar stereographic grid for the south
   maskvsh=find(nclatv.latVertex.data<-30*pi/180);
   [Xvsh,Yvsh]=ridgepack_geodetictoxy(nclatv.latVertex.data(maskvsh)*180/pi,...
                                    nclonv.lonVertex.data(maskvsh)*180/pi,-1);
   Xvsh(Xvsh<min(xsh) | Xvsh>max(xsh))=NaN;
   Yvsh(Yvsh<min(ysh) | Yvsh>max(ysh))=NaN;
   maskUVsh=(~isnan(Xvsh) & ~isnan(Yvsh));
   Xvsh=Xvsh(maskUVsh); Yvsh=Yvsh(maskUVsh);

   maskcsh=find(nclatc.latCell.data<-30*pi/180);
   [Xcsh,Ycsh]=ridgepack_geodetictoxy(nclatc.latCell.data(maskcsh)*180/pi,...
                                     nclonc.lonCell.data(maskcsh)*180/pi,-1);
   Xcsh(Xcsh<min(xsh) | Xcsh>max(xsh))=NaN;
   Ycsh(Ycsh<min(ysh) | Ycsh>max(ysh))=NaN;
   maskCCsh=(~isnan(Xcsh) & ~isnan(Ycsh));
   Xcsh=Xcsh(maskCCsh); Ycsh=Ycsh(maskCCsh);
 
   % first construct timeseries of only relevant months


   k=1;

   for year=startyear:endyear

    for month=monthset;

     cd(inlocation)

     % compile all data into a single time series spanning months

     infile=[dayname,num2str(year,'%4.4i'),'-',num2str(month,'%2.2i'),'-01.nc'];

     nc=ridgepack_clone(infile,fields);

     nc=rmfield(nc,'attributes'); 
     nc.attributes.title=titlenc;

     days=nc.time.data;
     nc.time.data=datenum(year,month,days,12,0,0.0);
     nc.time.units=['days since ',num2str(startyear,'%4.4i'),'-01-01 00:00:00.0'];
     nc.time.calendar='noleap';
     nc.time.dimension={'time'};

     nc.d2.data=[1 2];
     nc.d2.type='NC_INT';
     nc.d2.dimension={'d2'};

     nc.time_bounds.dimension={'d2','time'};
     nc.time_bounds.data=[datenum(year,month,days,0,0,0.0); ...
                         datenum(year,month,days,23,59,59.9)];
     nc.time_bounds.units=['days since ',num2str(startyear,'%4.4i'),'-01-01 00:00:00.0'];
     nc.time_bounds.calendar='noleap';

     nc.timeDaily_avg_Speed.units='m/s';
     nc.timeDaily_avg_Speed.type='NC_FLOAT';
     nc.timeDaily_avg_Speed.long_name='Daily Mean Sea Ice Speed';
     nc.timeDaily_avg_Speed.dimension=nc.timeDaily_avg_uVelocityGeo.dimension;
     nc.timeDaily_avg_Speed.data=sqrt(nc.timeDaily_avg_uVelocityGeo.data.^2 + ...
                                     nc.timeDaily_avg_vVelocityGeo.data.^2);

     nc=ridgepack_struct(nc);
     
     % setup output location
     xout=dir(char(outlocations{nrun}));
     if isempty(xout)
      mkdir(char(outlocations{nrun}))
     end
     cd(char(outlocations{nrun}))

     if k==1
      ridgepack_write(nc,outfile)
     else
      ridgepack_write(nc,outfile,{'time'},{0})
     end

     % now interpolate data to SSM/I grid
     %method='nearest'; This does not work properly
     method='cubic';

     if k==1

      ncnh.attributes.title=[titlenc,' on northern polar stereographic mesh'];

      ncnh.time.units=['days since ',num2str(startyear,'%4.4i'),'-01-01 00:00:00.0'];
      ncnh.time.calendar='noleap';
      ncnh.time.dimension={'time'};

      ncnh.d2.data=[1 2];
      ncnh.d2.type='NC_INT';
      ncnh.d2.dimension={'d2'};

      ncnh.time_bounds.dimension={'d2','time'};
      ncnh.time_bounds.units=['days since ',num2str(startyear,'%4.4i'),'-01-01 00:00:00.0'];
      ncnh.time_bounds.calendar='noleap';

      ncnh.u.units='m/s';
      ncnh.u.type='NC_FLOAT';
      ncnh.u.long_name='Sea Ice U-Velocity';
      ncnh.u.dimension={'y','x','time'};

      ncnh.v.units='m/s';
      ncnh.v.type='NC_FLOAT';
      ncnh.v.long_name='Sea Ice U-Velocity';
      ncnh.v.dimension=ncnh.u.dimension;

      ncnh.speed.units='m/s';
      ncnh.speed.type='NC_FLOAT';
      ncnh.speed.long_name='Daily Mean Sea Ice Speed';
      ncnh.speed.dimension=ncnh.u.dimension;

      ncnh.conc.type='NC_FLOAT';
      ncnh.conc.long_name='Sea Ice Concentration';
      ncnh.conc.dimension=ncnh.u.dimension;


      ncsh.attributes.title=[titlenc,' on southern polar stereographic mesh'];

      ncsh.time.units=['days since ',num2str(startyear,'%4.4i'),'-01-01 00:00:00.0'];
      ncsh.time.calendar='noleap';
      ncsh.time.dimension={'time'};

      ncsh.d2.data=[1 2];
      ncsh.d2.type='NC_INT';
      ncsh.d2.dimension={'d2'};

      ncsh.time_bounds.dimension={'d2','time'};
      ncsh.time_bounds.units=['days since ',num2str(startyear,'%4.4i'),'-01-01 00:00:00.0'];
      ncsh.time_bounds.calendar='noleap';

      ncsh.u.units='m/s';
      ncsh.u.type='NC_FLOAT';
      ncsh.u.long_name='Sea Ice U-Velocity';
      ncsh.u.dimension=ncnh.u.dimension;

      ncsh.v.units='m/s';
      ncsh.v.type='NC_FLOAT';
      ncsh.v.long_name='Sea Ice U-Velocity';
      ncsh.v.dimension=ncsh.u.dimension;
 
      ncsh.speed.units='m/s';
      ncsh.speed.type='NC_FLOAT';
      ncsh.speed.long_name='Daily Mean Sea Ice Speed';
      ncsh.speed.dimension=ncsh.u.dimension;

      ncsh.conc.type='NC_FLOAT';
      ncsh.conc.long_name='Sea Ice Concentration';
      ncsh.conc.dimension=ncsh.u.dimension;

     end
    
     for ti=1:length(nc.time.data)

      ncnh.time.data=datenum(year,month,days(ti),12,0,0.0);

      disp(['------- ',datestr(ncnh.time.data),' -------'])

      ncnh.time_bounds.data=[datenum(year,month,days(ti),0,0,0.0); ...
                             datenum(year,month,days(ti),23,59,59.9)];

      Uvnh=nc.timeDaily_avg_uVelocityGeo.data(maskvnh,ti);
      Vvnh=nc.timeDaily_avg_vVelocityGeo.data(maskvnh,ti);
      Ccnh=nc.timeDaily_avg_iceAreaCell.data(maskcnh,ti);
 
      ncnh.u.data(:,:)=griddata(Xvnh(:),Yvnh(:),Uvnh(maskUVnh),xnhm,ynhm,method);
      ncnh.v.data(:,:)=griddata(Xvnh(:),Yvnh(:),Vvnh(maskUVnh),xnhm,ynhm,method);
      ncnh.speed.data(:,:)=sqrt((ncnh.u.data.^2)+(ncnh.v.data.^2));
      ncnh.conc.data(:,:)=griddata(Xcnh(:),Ycnh(:),Ccnh(maskCCnh),xnhm,ynhm,method);
  
      ncnh.u.data(ncnh.mask.data==0)=NaN;
      ncnh.v.data(ncnh.mask.data==0)=NaN;
      ncnh.speed.data(ncnh.mask.data==0)=NaN;
      ncnh.conc.data(ncnh.mask.data==0)=NaN;

      ncsh.time.data=datenum(year,month,days(ti),12,0,0.0);

      ncsh.time_bounds.data=[datenum(year,month,days(ti),0,0,0.0); ...
                             datenum(year,month,days(ti),23,59,59.9)];

      Uvsh=nc.timeDaily_avg_uVelocityGeo.data(maskvsh,ti);
      Vvsh=nc.timeDaily_avg_vVelocityGeo.data(maskvsh,ti);
      Ccsh=nc.timeDaily_avg_iceAreaCell.data(maskcsh,ti);

      ncsh.u.data(:,:)=griddata(Xvsh(:),Yvsh(:),Uvsh(maskUVsh),xshm,yshm,method);
      ncsh.v.data(:,:)=griddata(Xvsh(:),Yvsh(:),Vvsh(maskUVsh),xshm,yshm,method);
      ncsh.speed.data(:,:)=sqrt((ncsh.u.data.^2)+(ncsh.v.data.^2));
      ncsh.conc.data(:,:)=griddata(Xcsh(:),Ycsh(:),Ccsh(maskCCsh),xshm,yshm,method);

      ncsh.u.data(ncsh.mask.data==0)=NaN;
      ncsh.v.data(ncsh.mask.data==0)=NaN;
      ncsh.speed.data(ncsh.mask.data==0)=NaN;
      ncsh.conc.data(ncsh.mask.data==0)=NaN;

      if k==1 & ti==1
       ridgepack_write(ncnh,outfileregridnh)
       ridgepack_write(ncsh,outfileregridsh)
      else
       ridgepack_write(ncnh,outfileregridnh,{'time'},{0})
       ridgepack_write(ncsh,outfileregridsh,{'time'},{0})
      end


     end

     k=k+length(nc.time.data);

     clear nc

    end

   end

   clear Uvnh Vvnh Ccnh ncnh Uvsh Vvsh Ccsh ncsh nclatv nclatc nclonv nclonc

  end


  % next create means with additional statistics from timeseries

  cd(char(outlocations{nrun}))
  

  % non-gridded means for native mesh (keep switched off if running out of memory)

  if reducenative

   ncseries=ridgepack_clone(outfile);
 
   nc=ridgepack_reduce(ncseries,{'time'},{},true);

   nc=rmfield(nc,'attributes');
   titlenc=[char(outnames{nrun}),' mean, months:'];
   for months=monthset;
    titlenc=[titlenc,' ',char(mons{months})];
   end
   titlenc=[titlenc,', years: ',num2str(startyear,'%4.4i'),' to ',num2str(endyear,'%4.4i')];
   nc.attributes.title=titlenc;
   meanoutfile=[char(outnames{nrun}),'.mean.',monthnum,'.',yeartag];
 
   ridgepack_write(nc,meanoutfile);

   clear nc ncseries

  end


  % gridded means

  ncnhm=ridgepack_reduce(ridgepack_clone(outfileregridnh),{'time'},{},true);

  ncnhm=rmfield(ncnhm,'attributes');
  ncnhm.attributes.title=[titlenc,', northern hemisphere'];
  meanoutfilenh=[outfileregridnh,'.mean'];
  ridgepack_write(ncnhm,meanoutfilenh);

  clear ncnhm 


  ncshm=ridgepack_reduce(ridgepack_clone(outfileregridsh),{'time'},{},true);
 
  ncshm=rmfield(ncshm,'attributes');
  ncshm.attributes.title=[titlenc,' southern hemisphere'];
  meanoutfilesh=[outfileregridsh,'.mean'];
  ridgepack_write(ncshm,meanoutfilesh);

  clear ncshm

 end % monthset

end % simulations


