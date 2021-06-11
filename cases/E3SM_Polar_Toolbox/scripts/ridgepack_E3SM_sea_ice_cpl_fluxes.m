+1%function E3SM_sea_ice_extent_volume(runcell,leg,minthick,mintime,maxtime,pub,pubdir)

clear
close all

% cell array of simulation abbreviation
tag='seaice.20210316.v2beta3GM900.IceReconcile.ne30pg2_EC30to60E2r2.chrysalis';

% cell array of full name of simulation
run='20210316';

% minimum time of flux calculation
mintime=datenum(0061,1,1);

% maxumum time of flux calculation
maxtime=datenum(0065,1,1);

% set line width
lw=0.5; 

% set worklocation location
workloc='/Users/afroberts/work';

% extract mesh lat, long, and cell area from the mpassi restart file
cd(['/Users/afroberts/data/MODEL/E3SM/',run,'/grid'])
ncgrid=ridgepack_clone('mpassi.rst.0006-01-01_00000.nc',...
                        {'latCell','lonCell','areaCell'});
% ncgrid.areaCell.data(ncgrid.latitude.data<0)=0;

% set directory of data location
basedir=['/Users/afroberts/data/MODEL/E3SM/',run,'/flux'];

% set model file name nomenclature
filetype1=['seaice.20210316.v2beta3GM900.IceReconcile.ne30pg2_EC30to60E2r2.chrysalis.cpl.ha.'];

% set state fields to be extracted 
fields{1}='i2xavg_Si_ifrac';

% set energy fields to be extracted 
fielde{1}='i2xavg_Faii_swnet';
fielde{2}='i2xavg_Fioi_swpen';
fielde{3}='x2iavg_Faxa_lwdn';
fielde{4}='i2xavg_Faii_lwup';
fielde{5}='i2xavg_Faii_lat';
fielde{6}='i2xavg_Faii_sen';
fielde{7}='i2xavg_Fioi_melth';
fielde{8}='x2iavg_Fioo_frazil';

% set mass fields to be extracted 
fieldm{1}='x2iavg_Faxa_rain';
fieldm{2}='x2iavg_Faxa_snow';
fieldm{3}='i2xavg_Faii_evap';
fieldm{4}='i2xavg_Fioi_meltw';
fieldm{5}='x2iavg_Fioo_frazil';

% generate comma separated field list for state
for fi=1:length(fields)
  idx=strfind(char(fields{fi}),'_');
  str=char(fields{fi});
  fieldss{fi}=str(idx(2)+1:end);
end

% generate comma separated field list for energy
for fi=1:length(fielde)
  idx=strfind(char(fielde{fi}),'_');
  str=char(fielde{fi});
  fieldes{fi}=str(idx(2)+1:end);
end

% generate comma separated field list for mass
for fi=1:length(fieldm)
  idx=strfind(char(fieldm{fi}),'_');
  str=char(fieldm{fi});
  fieldms{fi}=str(idx(2)+1:end);
end

% now get the data
cd([basedir])
startyear=str2num(datestr(mintime,'yyyy'));
endyear=str2num(datestr(maxtime,'yyyy'))-1

% clear existing netcdf structure
clear nce ncm

% find parts of mesh in northern and southern hemisphere and total mesh area
totarea=sum(ncgrid.areaCell.data(:));
%totarea=1;

% grab the data one month at a time to construct a timeseries of integrated quantities
k=1;
for year=startyear:endyear
  for month=1:12
    % grab data
    nctemps=ridgepack_clone([filetype1,...
             num2str(year,'%4.4i'),'-',num2str(month,'%2.2i')],fields);
    nctempe=ridgepack_clone([filetype1,...
             num2str(year,'%4.4i'),'-',num2str(month,'%2.2i')],fielde);
    nctempm=ridgepack_clone([filetype1,...
             num2str(year,'%4.4i'),'-',num2str(month,'%2.2i')],fieldm);

    % start accumulating at first timestep, first building timeseries structure
    if year==startyear & month==1

     ncs=nctemps;
     nce=nctempe;
     ncm=nctempm;

     for fi=1:length(fields)
      ncs.(char(fields{fi})).data=10000000000*ones([1 (endyear-startyear+1)*12]);
     end
     for fi=1:length(fielde)
      nce.(char(fielde{fi})).data=10000000000*ones([1 (endyear-startyear+1)*12]);
     end
     for fi=1:length(fieldm)
      ncm.(char(fieldm{fi})).data=10000000000*ones([1 (endyear-startyear+1)*12]);
     end

    end
 
    % accumulate timeseries for subsequent months and years at time step k
    for fi=1:length(fields)
      ncs.(char(fields{fi})).data(k)=sum(nctemps.(char(fields{fi})).data.*ncgrid.areaCell.data)./totarea;
    end

    for fi=1:length(fielde)
      %nce.(char(fielde{fi})).data(k)=sum(nctempe.(char(fielde{fi})).data.*ncgrid.areaCell.data.*nctemps.(char(fields{1})).data)./totarea;
      nce.(char(fielde{fi})).data(k)=sum(nctempe.(char(fielde{fi})).data.*ncgrid.areaCell.data)./totarea;
    end

    for fi=1:length(fieldm)
      %ncm.(char(fieldm{fi})).data(k)=sum(nctempm.(char(fieldm{fi})).data.*ncgrid.areaCell.data.*nctemps.(char(fields{1})).data)./totarea;
      ncm.(char(fieldm{fi})).data(k)=sum(nctempm.(char(fieldm{fi})).data.*ncgrid.areaCell.data)./totarea;
    end

    ncs.time.data(k)=datenum(year,month+1,1);
    nce.time.data(k)=datenum(year,month+1,1);
    ncm.time.data(k)=datenum(year,month+1,1);

    k=k+1;

  end
end

clear nctemps nctempe nctempm

if k==1
  error([run,' is not within the time bracket defined'])
end

% clean up netcdf structure for energy
ncs=rmfield(ncs,{'attributes'});
ncs.attributes.title=['E3SM MPAS-SI sea ice state fluxes for ',run];
ncs.time.calendar='proleptic_gregorian';
for j=1:length(fields)
 ncs.(char(fields{j})).dimension={'time'};
end
ncs=ridgepack_struct(ncs);   

nce=rmfield(nce,{'attributes'});
nce.attributes.title=['E3SM MPAS-SI sea ice energy fluxes for ',run];
nce.time.calendar='proleptic_gregorian';
for j=1:length(fielde)
 nce.(char(fielde{j})).dimension={'time'};
end
nce=ridgepack_struct(nce);   

% clean up netcdf structure for mass
ncm=rmfield(ncm,{'attributes'});
ncm.attributes.title=['E3SM MPAS-SI sea ice mass fluxes for ',run];
ncm.time.calendar='proleptic_gregorian';
for j=1:length(fieldm)
 ncm.(char(fieldm{j})).dimension={'time'};
end
ncm=ridgepack_struct(ncm);   

% center model values in time (use mean center value for first sample).
%nce.time.data(1)=nce.time.data(1)-mean(diff(nce.time.data)/2);
%nce.time.data(2:end)=nce.time.data(2:end)-diff(nce.time.data)/2;

% center model values in time (use mean center value for first sample).
%ncm.time.data(1)=ncm.time.data(1)-mean(diff(ncm.time.data)/2);
%ncm.time.data(2:end)=ncm.time.data(2:end)-diff(ncm.time.data)/2;

% run checks on data structure 
ncs=ridgepack_struct(ncs);
nce=ridgepack_struct(nce);
ncm=ridgepack_struct(ncm);

% set up totals field for energy
nce.netEnergyFlux_cpl.long_name=['E3SM MPAS-SI sea ice energy total for ',run];
nce.netEnergyFlux_cpl.units='W';
nce.netEnergyFlux_cpl.dimension=nce.time.dimension;
nce.netEnergyFlux_cpl.data=zeros(size(nce.(char(fielde{1})).data));

%fielde{1}='i2xavg_Faii_swnet';
%fielde{2}='i2xavg_Fioi_swpen';
%fielde{3}='x2iavg_Faxa_lwdn';
%fielde{4}='i2xavg_Faii_lwup';
%fielde{5}='i2xavg_Faii_lat';
%fielde{6}='i2xavg_Faii_sen';
%fielde{7}='i2xavg_Fioi_melth';
%fielde{8}='x2iavg_Fioo_frazil';

nce.netEnergyFlux_cpl.data = nce.(char(fielde{1})).data ...
                            -nce.(char(fielde{2})).data ...
                            +nce.(char(fielde{3})).data ...
                            +nce.(char(fielde{4})).data ...
                            +nce.(char(fielde{5})).data ...
                            +nce.(char(fielde{6})).data ...
                            +nce.(char(fielde{7})).data;

% set up totals field for mass
ncm.netMassFlux_cpl.long_name=['E3SM MPAS-SI sea ice mass total for ',run];
ncm.netMassFlux_cpl.units='kg/s';
ncm.netMassFlux_cpl.dimension=nce.time.dimension;
ncm.netMassFlux_cpl.data=zeros(size(nce.time.data));

%fieldm{1}='x2iavg_Faxa_rain';
%fieldm{2}='x2iavg_Faxa_snow';
%fieldm{3}='i2xavg_Faii_evap';
%fieldm{4}='i2xavg_Fioi_meltw';
%fieldm{5}='x2iavg_Fioo_frazil';

ncm.netMassFlux_cpl.data = ncm.x2iavg_Faxa_rain.data ...
                          +ncm.x2iavg_Faxa_snow.data ...
                          -ncm.i2xavg_Faii_evap.data ...
                          -ncm.i2xavg_Fioi_meltw.data ...
                          +ncm.x2iavg_Fioo_frazil.data;

if 1==1

% plot data for energy in timeseries
ridgepack_multiplot(2,1,1,1,'a'); 
emax=-1000000;
emin=1000000;

% set color scheme
col=colormap(lines(length(fielde)+1)); 

for j=1:length(fielde)
  he(j)=plot(nce.time.data,nce.(char(fielde{j})).data,'Color',col(j,:),'LineWidth',lw);
  hold on
  emax=max([emax nce.(char(fielde{j})).data]);
  emin=min([emin nce.(char(fielde{j})).data]);
end
he(j+1)=plot(nce.time.data,nce.netEnergyFlux_cpl.data,'Color',col(j+1,:),'LineWidth',lw);
drawnow
fieldes{j+1}='Net';
legend(he,fieldes,'Location','NorthEast','FontSize',7)
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
col=colormap(lines(length(fieldm)+1)); 

for j=1:length(fieldm)
  hm(j)=plot(ncm.time.data,ncm.(char(fieldm{j})).data,'Color',col(j,:),'LineWidth',lw);
  hold on
  mmax=max([mmax ncm.(char(fieldm{j})).data]);
  mmin=min([mmin ncm.(char(fieldm{j})).data]);
end
hm(j+1)=plot(ncm.time.data,ncm.netMassFlux_cpl.data,'Color',col(j+1,:),'LineWidth',lw);
fieldms{j+1}='Net';
drawnow
legend(hm,fieldms,'Location','NorthEast','FontSize',7)
legend('boxoff')
ylim([1.1*mmin 1.1*mmax])
xlim([ncm.time.data(1) ncm.time.data(end)])
datetick('x','YY','keeplimits')
ylabel('Mass Flux')
xlabel('Year')
set(gca,'box','on')
hold off

% line up the plot frames and give the whole plot a title
ridgepack_multialign(gcf,['E3SM MPAS-SeaIce ',...
                  datestr(mintime,'YYYY'),'-',datestr(maxtime,'YYYY')],9)

% printout figure and write energy and mass timeseries
cd(workloc)
ridgepack_fprint('png',[run,'_',datestr(mintime,'YYYY'),'_',datestr(maxtime,'YYYY'),'_cpl.png'],1,1)

end

cd('/Users/afroberts/data/MODEL/E3SM/20210316/flux')
ridgepack_write(nce,[run,'_',filetype1,'_energy_flux_cpl'])
ridgepack_write(ncm,[run,'_',filetype1,'_mass_flux_cpl'])


