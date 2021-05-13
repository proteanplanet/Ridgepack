clear
close all

% short name of simulation
%run='20210508';
run='20210427';

% full name of simulation abbreviation
tag='20210427.v2rc1a.CPcheck5.ne30pg2_EC30to60E2r2.chrysalis';

% hemis: 1 = global, 2=NH, 3=SH
hemis=1;

% set mass fields
masstype='NetMass';

% set heat fields
%heattype='NetHeat';
%heattype='Frazil';
%heattype='Sensible';
%heattype='OceanHeat';
%heattype='LongwaveUp';
%heattype='LongwaveDown';
%heattype='NetSW';
heattype='Latent';
%heattype='SnowHeat';

% minimum year of flux calculation
mintime=datenum(0001,1,1);

% maxumum year of flux calculation
maxtime=datenum(0003,1,1);

% set line width in plots
lw=0.5; 

% set plot location
workloc='/Users/afroberts/work';

% diagnostics
if isempty(run);
  error('Casename is empty');
else
  disp(['CASE: ',run]);
end
 
% set directory of data location
basedir=['/Users/afroberts/data/MODEL/E3SM/',run,'/flux'];

% set model file name nomenclature
filetype1=[tag,'.mpassi.hist.am.conservationCheck.'];

%fielde{1}='energyConsSurfaceHeatFlux';

% extract mesh lat, long, and cell area from the mpassi restart file
cd(['/Users/afroberts/data/MODEL/E3SM/',run,'/grid'])
ncgrid=ridgepack_clone('mpassi.rst.0006-01-01_00000.nc',...
                        {'latCell','lonCell','areaCell'});

% hemisphere tags and title
if hemis<=1
 hemisphere='Global';
 hemtag='global';
elseif hemis==2
 hemisphere='Northern Hemisphere';
 hemtag='NH';
elseif hemis==3
 hemisphere='Southern Hemisphere';
 hemtag='SH';
end

% now get the data
cd([basedir])

% clear existing netcdf structure
clear nce

% check for available data, start and end years
xdir=dir([filetype1,'*'])
nc=ridgepack_clone(xdir(1).name);
startyear=max(str2num(nc.xtime.data(1,1:4)),str2num(datestr(mintime,'YYYY')));
nc=ridgepack_clone(xdir(end).name);
endyear=min(str2num(nc.xtime.data(1,1:4)),str2num(datestr(maxtime,'YYYY')));

% grab MPAS-SI data one month at a time to construct a timeseries of integrated quantities
k=0;
for year=startyear:endyear

 % grab data for conservation check
 nctempe=ridgepack_clone([filetype1,num2str(year,'%4.4i')]);

 % check length
 maxlength=length(nctempe.time.data);

 % create time vector for the given year
 for i=1:12
  time([(year-startyear)*12+i])=datenum(year-1,i,1);
 end

 % start accumulating at first timestep, first building timeseries structure
 if year==startyear

  massmpassi=NaN.*ones([1 (endyear-startyear+1)*12]);

  masscpl=NaN.*ones([1 (endyear-startyear+1)*12]);

  heatmpassi=NaN.*ones([1 (endyear-startyear+1)*12]);

  heatcpl=NaN.*ones([1 (endyear-startyear+1)*12]);

 end

 % accumulate timeseries for subsequent months and years at time step k
 datrange=[(year-startyear)*12+1:(year-startyear)*12+maxlength];

 if strcmp(masstype,'NetMass')

  massmpassi(datrange)=nctempe.netMassFlux.data(hemis,:);

 end

 if strcmp(heattype,'NetHeat')

  heatmpassi(datrange)=nctempe.netEnergyFlux.data(hemis,:);

 elseif strcmp(heattype,'Frazil')

  heatmpassi(datrange)=nctempe.energyConsFreezingPotential.data(hemis,:);

 elseif strcmp(heattype,'Sensible')

  heatmpassi(datrange)=nctempe.energyConsSensibleHeatFlux.data(hemis,:);

 elseif strcmp(heattype,'OceanHeat')

  heatmpassi(datrange)=nctempe.energyConsOceanHeatFlux.data(hemis,:);

 elseif strcmp(heattype,'LongwaveUp')

  heatmpassi(datrange)=nctempe.energyConsLongwaveUp.data(hemis,:);

 elseif strcmp(heattype,'LongwaveDown')

  heatmpassi(datrange)=nctempe.energyConsLongwaveDown.data(hemis,:);

 elseif strcmp(heattype,'NetSW')

  heatmpassi(datrange)=nctempe.energyConsAbsorbedShortwaveFlux.data(hemis,:)+...
                       nctempe.energyConsOceanShortwaveFlux.data(hemis,:);

 elseif strcmp(heattype,'Latent')

  heatmpassi(datrange)=nctempe.energyConsLatentHeat2.data(hemis,:)

 elseif strcmp(heattype,'SnowHeat')

  heatmpassi(datrange)=nctempe.energyConsSnowfallHeat.data(hemis,:);

 end

 k=k+1

end

% shift MPAS-SI output by one month to match CPL logs, and set up CPL netcdf field
time=time(2:end);
massmpassi=massmpassi(2:end);
masscpl=NaN*ones(size(massmpassi));
heatmpassi=heatmpassi(2:end);
heatcpl=NaN*ones(size(heatmpassi));

%clear nctempe

% get mass  flux
if strcmp(masstype,'NetMass')
 filename = ['/Users/afroberts/data/MODEL/E3SM/',run,'/flux/mass.txt'];
 [seaicemassNh,seaicemassSh]=ridgepack_E3SM_sea_ice_cpl_heat_check(filename);
end

if hemis<=1
  masscpl=[seaicemassNh(1:length(time))+seaicemassSh(1:length(time))]*1.E-6;
elseif hemis==2
  masscpl=[seaicemassNh(1:length(time))]*1.E-6;
elseif hemis==3
  masscpl=[seaicemassSh(1:length(time))]*1.E-6;
end

% get mass  flux
if strcmp(heattype,'NetHeat')

 filename = ['/Users/afroberts/data/MODEL/E3SM/',run,'/flux/heat.txt'];
 [seaiceheatNh,seaiceheatSh]=ridgepack_E3SM_sea_ice_cpl_heat_check(filename);

elseif strcmp(heattype,'Frazil')

 filename = ['/Users/afroberts/data/MODEL/E3SM/',run,'/flux/hfreeze.txt'];
 [seaiceheatNh,seaiceheatSh]=ridgepack_E3SM_sea_ice_cpl_heat_check(filename);

elseif strcmp(heattype,'Sensible')

 filename = ['/Users/afroberts/data/MODEL/E3SM/',run,'/flux/hsen.txt'];
 [seaiceheatNh,seaiceheatSh]=ridgepack_E3SM_sea_ice_cpl_heat_check(filename);

elseif strcmp(heattype,'OceanHeat')

 filename = ['/Users/afroberts/data/MODEL/E3SM/',run,'/flux/hmelt.txt'];
 [seaiceheatNh,seaiceheatSh]=ridgepack_E3SM_sea_ice_cpl_heat_check(filename);

elseif strcmp(heattype,'LongwaveUp')

 filename = ['/Users/afroberts/data/MODEL/E3SM/',run,'/flux/hlwup.txt'];
 [seaiceheatNh,seaiceheatSh]=ridgepack_E3SM_sea_ice_cpl_heat_check(filename);

elseif strcmp(heattype,'LongwaveDown')

 filename = ['/Users/afroberts/data/MODEL/E3SM/',run,'/flux/hlwdn.txt'];
 [seaiceheatNh,seaiceheatSh]=ridgepack_E3SM_sea_ice_cpl_heat_check(filename);

elseif strcmp(heattype,'NetSW')

 filename = ['/Users/afroberts/data/MODEL/E3SM/',run,'/flux/hnetsw.txt'];
 [seaiceheatNh,seaiceheatSh]=ridgepack_E3SM_sea_ice_cpl_heat_check(filename);

elseif strcmp(heattype,'Latent')

 filename = ['/Users/afroberts/data/MODEL/E3SM/',run,'/flux/hlatvap.txt'];
 [seaiceheatNh,seaiceheatSh]=ridgepack_E3SM_sea_ice_cpl_heat_check(filename);

elseif strcmp(heattype,'SnowHeat')

 filename = ['/Users/afroberts/data/MODEL/E3SM/',run,'/flux/hlatfus.txt'];
 [seaiceheatNh,seaiceheatSh]=ridgepack_E3SM_sea_ice_cpl_heat_check(filename);

end

if hemis<=1
  heatcpl=[seaiceheatNh(1:length(time))+seaiceheatSh(1:length(time))];
elseif hemis==2
  heatcpl=[seaiceheatNh(1:length(time))];
elseif hemis==3
  heatcpl=[seaiceheatSh(1:length(time))];
end

% set color scheme
col=colormap(lines(4)); 

% plot first frame of mass flux
ridgepack_multiplot(3,1,1,1,['a ',masstype]); 

plot(time,massmpassi,'Color',col(1,:))
plot(time,masscpl,'--','Color',col(4,:))
legend({'MPAS-SI','CPL'},'Location','North','Orientation','horizontal')
legend('boxoff')
ylabel('Mass (kg m$^{-2}$ s$^{-1}$)')
xlim([time(1) time(end)])
datetick('x','YY','keeplimits')
grid on

ridgepack_clearax('x',1)

% plot second frame of energy flux
ridgepack_multiplot(3,1,2,1,['b ',heattype]);

plot(time,heatmpassi,'Color',col(1,:))
plot(time,heatcpl,'--','Color',col(4,:))
ylabel('Heat (W m$^{-2}$)')
xlim([time(1) time(end)])
datetick('x','YY','keeplimits')
grid on

ridgepack_clearax('x',1)

% start third frame
ridgepack_multiplot(3,1,3,1,'c'); 

% plot delta mass on left axis
yyaxis left
deltamass=massmpassi-masscpl;
h1=plot(time,deltamass);
ylabel('$\Delta$Mass (kg m$^{-2}$ s$^{-1}$)')
hold on
yl=ylim;
ylim([yl(1)-diff(yl)*0.15 yl(2)])

% plot delta heat on right axis
yyaxis right
deltaheat=heatmpassi-heatcpl;
h2=plot(time,deltaheat);
ylabel('$\Delta$Heat (W m$^{-2}$)')
xlim([time(1) time(end)])
datetick('x','YY','keeplimits')
xlabel('Model Year')
yl=ylim;
ylim([yl(1)-diff(yl)*0.15 yl(2)])
grid on

% add overall diagnostics to plot
summass=sum(deltamass)/length(deltamass);
sumheat=sum(deltaheat)/length(deltaheat);
if summass>0
 legmass=['Net Mass +',num2str(summass),' kg m$^{-2}$ s$^{-1}$'];
else
 legmass=['Net Mass ',num2str(summass),' kg m$^{-2}$ s$^{-1}$'];
end
if sumheat>0
 legheat=['Net Heat +',num2str(sumheat),' W m$^{-2}$'];
else
 legheat=['Net Heat',num2str(sumheat),' W m$^{-2}$'];
end
legend([h1 h2],{legmass,legheat},'Orientation','horizontal','Location','South','Interpreter','Latex')
legend('boxoff')

% align all the plot frames and give it a title
ridgepack_multialign(gcf,['E3SM Sea Ice Conservation Check - ',hemisphere,],9)

% output graphics and netcdf file of fluxes from both model and coupler
cd(workloc)
ridgepack_fprint('png',[run,'_',hemtag,'_',num2str(startyear,'%4.4i'),'_',...
                                num2str(endyear,'%4.4i'),'_',masstype,'_',heattype,'.png'],1,1)

%ridgepack_write(nce,[run,'_',hemtag,'_',num2str(startyear,'%4.4i'),'_',...
%                     num2str(endyear,'%4.4i'),'_',masstype,'_',heattype,'_conservation_check'])


