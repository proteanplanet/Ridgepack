clear
close all

% short name of simulation
run='20210427';

% full name of simulation abbreviation
tag='20210427.v2rc1a.CPcheck3.ne30pg2_EC30to60E2r2.chrysalis';

% hemis: 1 = global, 2=NH, 3=SH
hemis=1;

% minimum year of flux calculation
mintime=datenum(0001,1,1);

% maxumum year of flux calculation
maxtime=datenum(0005,1,1);

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

% set mass fields
massgroup{1}.name='NetMass';
massgroup{1}.mpassi={'netMassFlux'};
massgroup{1}.cpl={'mass'};

% set energy fields to be extracted 
heatgroup{1}.name='NetHeat';
heatgroup{1}.mpassi={'netEnergyFlux'};
heatgroup{1}.cpl={'heat'};

heatgroup{2}.name='FrazilHeat';
heatgroup{2}.mpassi={'energyConsFreezingPotential'};
heatgroup{2}.cpl={'hfreeze'};

heatgroup{3}.name='OceanHeatFlux';
heatgroup{3}.mpassi={'energyConsOceanHeatFlux'};
heatgroup{3}.cpl={'hmelt'};

heatgroup{4}.name='SensibleHeatFlux';
heatgroup{4}.mpassi={'energyConsSensibleHeatFlux'};
heatgroup{4}.cpl={'hsen'};

heatgroup{5}.name='';
heatgroup{5}.mpassi={''};
heatgroup{5}.cpl={''};

heatgroup{6}.name='';
heatgroup{6}.mpassi={''};
heatgroup{6}.cpl={''};

heatgroup{7}.name='';
heatgroup{7}.mpassi={''};
heatgroup{7}.cpl={''};

heatgroup{8}.name='';
heatgroup{8}.mpassi={''};
heatgroup{8}.cpl={''};


%fielde{1}='energyConsAbsorbedShortwaveFlux';
%fielde{1}='energyConsLatentHeat';
%fielde{1}='energyConsOceanShortwaveFlux';
%fielde{1}='energyConsSurfaceHeatFlux';
%fielde{1}='energyConsLongwaveUp';
%fielde{1}='energyConsLongwaveDown';
%fielde{1}='energyConsSnowfallHeat';

%elseif strcmp(char(fielde{1}),'energyConsAbsorbedShortwaveFlux')
% heatname={'hnetsw'};
%elseif strcmp(char(fielde{1}),'energyConsLatentHeat')
% heatname={'hlatfus','hlatvap'};

%heatname='heat';
%heatname='hfreeze';
%heatname='hlatfus';
%heatname='hlwup';
%heatname='hlwdn';
%heatname='hnetsw';
%heatname='hlatvap';
%heatname='hsen';

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

 % start accumulating at first timestep, first building timeseries structure
 if year==startyear

  nce.attributes.title='E3SM MPAS-SI Conservation Check';
  nce.time=nctempe.time;
  nce.time.calendar='Gregorian';
  nce.time.units='days since 0001-01-01';

  nce.(char(massgroup{1}.name))=nctempe.(char(massgroup{1}.mpassi{1}));
  nce.(char(massgroup{1}.name)).data=NaN.*ones([1 (endyear-startyear+1)*12]);
  nce.(char(massgroup{1}.name)).dimension={'time'};
  nce.([char(massgroup{1}.name),'_cpl'])=nce.(char(massgroup{1}.name));
  nce.([char(massgroup{1}.name),'_cpl']).long_name=['coupler ',char(massgroup{1}.name)];

  nce.(char(heatgroup{1}.name))=nctempe.(char(heatgroup{1}.mpassi{1}));
  nce.(char(heatgroup{1}.name)).data=NaN.*ones([1 (endyear-startyear+1)*12]);
  nce.(char(heatgroup{1}.name)).dimension={'time'};
  nce.([char(heatgroup{1}.name),'_cpl'])=nce.(char(heatgroup{1}.name));
  nce.([char(heatgroup{1}.name),'_cpl']).long_name=['coupler ',char(heatgroup{1}.name)];

 end

 % create time vector for the given year
 for i=1:12
  nce.time.data([(year-startyear)*12+i])=datenum(year,i,1);
 end

 % accumulate timeseries for subsequent months and years at time step k
 datrange=[(year-startyear)*12+1:(year-startyear)*12+maxlength];

 nce.(char(massgroup{1}.name)).data(datrange)=nctempe.(char(massgroup{1}.mpassi{1})).data(hemis,:);
 nce.(char(heatgroup{1}.name)).data(datrange)=nctempe.(char(heatgroup{1}.mpassi{1})).data(hemis,:);

 k=k+1

end

% shift MPAS-SI output by one month to match CPL logs, and set up CPL netcdf field
nce.time.data=nce.time.data(2:end);
nce.(char(massgroup{1}.name)).data=nce.(char(massgroup{1}.name)).data(2:end);
nce.([char(massgroup{1}.name),'_cpl']).data=NaN*ones(size(nce.(char(massgroup{1}.name)).data));
nce.(char(heatgroup{1}.name)).data=nce.(char(heatgroup{1}.name)).data(2:end);
nce.([char(heatgroup{1}.name),'_cpl']).data=NaN*ones(size(nce.(char(heatgroup{1}.name)).data));

clear nctempe

%  grab data extracted from CPL logs, or from CPL history files
%cpllogs=false;
cpllogs=true;
if cpllogs

 filename = ['/Users/afroberts/data/MODEL/E3SM/',run,'/flux/monthly_area.txt'];
 [seaiceareaNh,seaiceareaSh]=ridgepack_E3SM_sea_ice_cpl_area_check(filename);

 % get heat flux
 seaiceheatNh=0;
 seaiceheatSh=0;
 for ij=1:length(heatname)
  filename = ['/Users/afroberts/data/MODEL/E3SM/',run,'/flux/monthly_',char(heatname{ij}),'.txt'];
  [Nh,Sh]=ridgepack_E3SM_sea_ice_cpl_heat_check(filename);
  seaiceheatNh=seaiceheatNh+Nh;
  seaiceheatSh=seaiceheatSh+Sh;
 end

 if hemis<=1
  nce.([char(fielde{1}),'_cpl']).data=[seaiceheatNh(1:length(nce.time.data))+...
                                       seaiceheatSh(1:length(nce.time.data))];
 elseif hemis==2
  nce.([char(fielde{1}),'_cpl']).data=[seaiceheatNh(1:length(nce.time.data))];
 elseif hemis==3
  nce.([char(fielde{1}),'_cpl']).data=[seaiceheatSh(1:length(nce.time.data))];
 end

 filename = ['/Users/afroberts/data/MODEL/E3SM/',run,'/flux/monthly_',massname,'.txt'];
 [seaicemassNh,seaicemassSh]=ridgepack_E3SM_sea_ice_cpl_heat_check(filename);
 if hemis<=1
  nce.([char(fielde{2}),'_cpl']).data=[seaicemassNh(1:length(nce.time.data))+...
                                       seaicemassSh(1:length(nce.time.data))]*1.E-6;
 elseif hemis==2
  nce.([char(fielde{2}),'_cpl']).data=[seaicemassNh(1:length(nce.time.data))]*1.E-6;
 elseif hemis==3
  nce.([char(fielde{2}),'_cpl']).data=[seaicemassSh(1:length(nce.time.data))]*1.E-6;
 end

else

 cd(['/Users/afroberts/data/MODEL/E3SM/',run,'/flux'])
 nccple=ridgepack_clone([run,'_seaice.',run,'.v2beta3GM900.IceReconcile.ne30pg2_EC30to60E2r2.chrysalis.cpl.ha._energy_flux_cpl.nc'],'netEnergyFlux_cpl')

 nccplm=ridgepack_clone([run,'_seaice.',run,'.v2beta3GM900.IceReconcile.ne30pg2_EC30to60E2r2.chrysalis.cpl.ha._mass_flux_cpl.nc'],'netMassFlux_cpl')

 for i=1:length(nccple.time.data)

  idx=find(nccple.time.data(i)==nce.time.data);

  if ~isempty(idx)

   nce.([char(fielde{1}),'_cpl']).data(idx)=nccple.netEnergyFlux_cpl.data(i);
   nce.([char(fielde{2}),'_cpl']).data(idx)=nccplm.netMassFlux_cpl.data(i);

  end

 end

end

% run checks on data structure 
nce=ridgepack_struct(nce);

% set color scheme
col=colormap(lines(4)); 

% plot first frame of mass flux
ridgepack_multiplot(3,1,1,1,['a ',char(fielde{2})]); 

plot(nce.time.data,nce.([char(fielde{2})]).data,'Color',col(1,:))
plot(nce.time.data,nce.([char(fielde{2}),'_cpl']).data,'--','Color',col(4,:))
legend({'MPAS-SI','CPL'},'Location','North','Orientation','horizontal')
legend('boxoff')
ylabel('Mass (kg m$^{-2}$ s$^{-1}$)')
xlim([nce.time.data(1) nce.time.data(end)])
datetick('x','YY','keeplimits')
grid on

ridgepack_clearax('x',1)

% plot second frame of energy flux
ridgepack_multiplot(3,1,2,1,['b ',char(fielde{1})]);

plot(nce.time.data,nce.([char(fielde{1})]).data,'Color',col(1,:))
plot(nce.time.data,nce.([char(fielde{1}),'_cpl']).data,'--','Color',col(4,:))
ylabel('Heat (W m$^{-2}$)')
xlim([nce.time.data(1) nce.time.data(end)])
datetick('x','YY','keeplimits')
grid on

ridgepack_clearax('x',1)

% start third frame
ridgepack_multiplot(3,1,3,1,'c'); 

% plot delta mass on left axis
yyaxis left
deltamass=nce.([char(fielde{2})]).data-nce.([char(fielde{2}),'_cpl']).data;
h1=plot(nce.time.data,deltamass)
ylabel('$\Delta$Mass (kg m$^{-2}$ s$^{-1}$)')
hold on
yl=ylim;
ylim([yl(1)-diff(yl)*0.15 yl(2)])

% plot delta heat on right axis
yyaxis right
deltaheat=nce.([char(fielde{1})]).data-nce.([char(fielde{1}),'_cpl']).data;
h2=plot(nce.time.data,deltaheat)
ylabel('$\Delta$Heat (W m$^{-2}$)')
xlim([nce.time.data(1) nce.time.data(end)])
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
                                num2str(endyear,'%4.4i'),'_',char(fielde{1}),'_',char(fielde{2}),'.png'],1,1)

ridgepack_write(nce,[run,'_',hemtag,'_',num2str(startyear,'%4.4i'),'_',...
                     num2str(endyear,'%4.4i'),'_',char(fielde{1}),'_',char(fielde{2}),'_conservation_check'])


