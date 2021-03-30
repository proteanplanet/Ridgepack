%function E3SM_sea_ice_extent_volume(runcell,leg,minthick,mintime,maxtime,pub,pubdir)

clear
close all

% cell array of simulation abbreviation
run='InteRFACE1alphaC';

% cell array of full name of simulation
tag='InteRFACE1alphaC';

% minimum time of flux calculation
mintime=datenum(0030,1,1);

% maxumum time of flux calculation
maxtime=datenum(0031,1,1);

% set hemisphere
arctic=true;
%arctic=false;
if arctic
 hem='Arctic';
else
 hem='Antarctic';
end

% set line width
lw=0.5; 

% set worklocation location
workloc='/Users/afroberts/work';

% extract mesh lat, long, and cell area from the mpassi restart file
cd(['/Users/afroberts/data/MODEL/E3SM/',run,'/grid'])
ncgrid=ridgepack_clone('mpassi.rst.0002-01-01_00000.nc',...
                        {'latCell','lonCell','areaCell'});

% diagnostics
if isempty(run);
  error('Casename is empty');
else
  disp(['CASE: ',run]);
end
 
% set directory of data location
basedir=['/Users/afroberts/data/MODEL/E3SM/',run,'/PI'];

% set model file name nomenclature
filetype1=[tag,'.mpassi.hist.am.timeSeriesStatsMonthly.'];

% set mass fields to be extracted 
fieldm{1}='';

% set energy fields to be extracted 
fielde{1}='timeMonthly_avg_absorbedShortwaveFluxInitialArea';
fielde{2}='timeMonthly_avg_oceanShortwaveFlux';
fielde{3}='timeMonthly_avg_longwaveDown';
fielde{4}='timeMonthly_avg_longwaveUpInitialArea';
fielde{5}='timeMonthly_avg_latentHeatFluxInitialArea';
fielde{6}='timeMonthly_avg_sensibleHeatFluxInitialArea';
fielde{7}='timeMonthly_avg_oceanHeatFlux';
fielde{8}='timeMonthly_avg_shortwaveDown';

% set mass fields to be extracted 
fieldm{1}='timeMonthly_avg_congelation';
fieldm{2}='timeMonthly_avg_frazilFormation';
fieldm{3}='timeMonthly_avg_snowiceFormation';
fieldm{4}='timeMonthly_avg_evaporativeWaterFluxInitialArea';
fieldm{5}='timeMonthly_avg_snowMelt';
fieldm{6}='timeMonthly_avg_surfaceIceMelt';
fieldm{7}='timeMonthly_avg_basalIceMelt';
fieldm{8}='timeMonthly_avg_lateralIceMelt';
fieldm{9}='timeMonthly_avg_oceanFreshWaterFlux';

% generate comma separated field list for energy
for fi=1:length(fielde)
  if fi==1
   fieldes=char(fielde{fi});
  else
   fieldes=[fieldes,',',char(fielde{fi})];
  end
end

% generate comma separated field list for mass
for fi=1:length(fieldm)
  if fi==1
   fieldms=char(fieldm{fi});
  else
   fieldms=[fieldms,',',char(fieldm{fi})];
  end
end

% now get the data
cd([basedir])
startyear=str2num(datestr(mintime,'yyyy'));
endyear=str2num(datestr(maxtime,'yyyy'))-1

% clear existing netcdf structure
clear nce ncm

% find parts of mesh in northern and southern hemisphere and total mesh area
if arctic
  idxe=find(ncgrid.latitude.data>=0);
else
  idxe=find(ncgrid.latitude.data<0);
end
totarea=sum(ncgrid.areaCell.data(idxe));

% grab the data one month at a time to construct a timeseries of integrated quantities
k=1;
for year=startyear:endyear
  for month=1:12
    % grab data
    nctempe=ridgepack_clone([filetype1,...
             num2str(year,'%4.4i'),'-',num2str(month,'%2.2i'),'-01'],fielde);
    nctempm=ridgepack_clone([filetype1,...
             num2str(year,'%4.4i'),'-',num2str(month,'%2.2i'),'-01'],fieldm);

    % start accumulating at first timestep, first building timeseries structure
    if year==startyear & month==1
     nce=nctempe;
     ncm=nctempm;
     for fi=1:length(fielde)
      nce.(char(fielde{fi})).data=10000000000*ones([1 (endyear-startyear+1)*12]);
     end
     for fi=1:length(fieldm)
      ncm.(char(fieldm{fi})).data=10000000000*ones([1 (endyear-startyear+1)*12]);
     end
    end
 
    % accumulate timeseries for subsequent months and years at time step k
    for fi=1:length(fielde)
      nce.(char(fielde{fi})).data(k)=sum(nctempe.(char(fielde{fi})).data(idxe).*ncgrid.areaCell.data(idxe))./totarea;
    end

    for fi=1:length(fieldm)
      ncm.(char(fieldm{fi})).data(k)=sum(nctempm.(char(fieldm{fi})).data(idxe).*ncgrid.areaCell.data(idxe))./totarea;
    end

    nce.time.data(k)=datenum(year,month+1,1);
    ncm.time.data(k)=datenum(year,month+1,1);
    k=k+1;
  end
end
clear nctempe nctempm

if k==1
  error([run,' is not within the time bracket defined'])
end

% clean up netcdf structure for energy
nce=rmfield(nce,{'attributes','nCells'});
nce.attributes.title=[hem,' E3SM MPAS-SI sea ice energy fluxes for ',run];
nce.time.calendar='proleptic_gregorian';
for j=1:length(fielde)
 nce.(char(fielde{j})).dimension={'time'};
end
nce=ridgepack_struct(nce);   

% clean up netcdf structure for mass
ncm=rmfield(ncm,{'attributes','nCells'});
ncm.attributes.title=[hem,' E3SM MPAS-SI sea ice mass fluxes for ',run];
ncm.time.calendar='proleptic_gregorian';
for j=1:length(fieldm)
 ncm.(char(fieldm{j})).dimension={'time'};
end
ncm=ridgepack_struct(ncm);   

% center model values in time (use mean center value for first sample).
nce.time.data(1)=nce.time.data(1)-mean(diff(nce.time.data)/2);
nce.time.data(2:end)=nce.time.data(2:end)-diff(nce.time.data)/2;

% center model values in time (use mean center value for first sample).
ncm.time.data(1)=ncm.time.data(1)-mean(diff(ncm.time.data)/2);
ncm.time.data(2:end)=ncm.time.data(2:end)-diff(ncm.time.data)/2;

% run checks on data structure 
nce=ridgepack_struct(nce);
ncm=ridgepack_struct(ncm);

% set up totals field for energy
nce.total.long_name=[hem,' E3SM MPAS-SI sea ice energy total for ',run];
nce.total.units='km^2';
nce.total.dimension=nce.time.dimension;
nce.total.data=zeros(size(nce.time.data));

% set up totals field for mass
ncm.total.long_name=[hem,' E3SM MPAS-SI sea ice mass total for ',run];
ncm.total.units='km^2';
ncm.total.dimension=nce.time.dimension;
ncm.total.data=zeros(size(nce.time.data));

% plot data for energy in timeseries
ridgepack_multiplot(2,1,1,1,'a'); 
emax=-1000000;
emin=1000000;

% set color scheme
col=colormap(lines(length(fielde))); 

for j=1:length(fielde)
  he(j)=plot(nce.time.data,nce.(char(fielde{j})).data,'Color',col(j,:),'LineWidth',lw);
  hold on
  emax=max([emax nce.(char(fielde{j})).data]);
  emin=min([emin nce.(char(fielde{j})).data]);
end
drawnow
legend(he,fielde,'Location','NorthEast','FontSize',7)
legend('boxoff')
ylim([1.1*emin 1.1*emax])
xlim([nce.time.data(1) nce.time.data(end)])
datetick('x','YY','keeplimits')
ridgepack_clearax('x',1)
ylabel('Energy Flux')
set(gca,'box','on')
hold off

% plot data for mass in timeseries
ridgepack_multiplot(2,1,2,1,'b'); 
mmax=-1000000;
mmin=1000000;

% set color scheme
col=colormap(lines(length(fieldm))); 

for j=1:length(fieldm)
  hm(j)=plot(ncm.time.data,ncm.(char(fieldm{j})).data,'Color',col(j,:),'LineWidth',lw);
  hold on
  mmax=max([mmax ncm.(char(fieldm{j})).data]);
  mmin=min([mmin ncm.(char(fieldm{j})).data]);
end
drawnow
legend(hm,fieldm,'Location','NorthEast','FontSize',7)
legend('boxoff')
ylim([1.1*mmin 1.1*mmax])
xlim([ncm.time.data(1) ncm.time.data(end)])
datetick('x','YY','keeplimits')
ylabel('Mass Flux')
xlabel('Year')
set(gca,'box','on')
hold off

% line up the plot frames and give the whole plot a title
ridgepack_multialign(gcf,[hem,' E3SM MPAS-SeaIce ',...
                  datestr(mintime,'YYYY'),'-',datestr(maxtime,'YYYY')],9)

% printout figure and write energy and mass timeseries
cd(workloc)
ridgepack_fprint('png',[run,'_',hem,'_',datestr(mintime,'YYYY'),'_',datestr(maxtime,'YYYY'),'.png'],1,1)
ridgepack_write(nce,[run,'_',filetype1,'_',hem,'_energy_flux'])
ridgepack_write(ncm,[run,'_',filetype1,'_',hem,'_mass_flux'])


