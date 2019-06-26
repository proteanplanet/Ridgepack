%function E3SM_sea_ice_extent_volume(runcell,leg,minthick,mintime,maxtime,pub,pubdir)

clear
close all

runcell={'E3SM-HR-V1','E3SM-HR-V0','E3SM-LR-V1'};
leg={'E3SM-HR-V1','E3SM-HR-V0','E3SM-LR-V1'};

runcell={'E3SM-LR-V1'};
leg={'E3SM-LR-V1'};

runcell={'E3SM-DECK-PI'};
leg={'E3SM-DECK-PI'};

minthick=0.0

% mist start in january and end in december
%mintime=datenum(0054,1,1);
mintime=datenum(0006,1,1);
%mintime=datenum(2000,1,1);
maxtime=datenum(0140,1,1);
%maxtime=datenum(2015,1,1);

%arctic=true;
arctic=false;
cesmle=false;
piomas=false;
%piomas=true;
icesat=false;
%icesat=true;
%runnin=false;
runnin=true;
passive=false;
%passive=true;

clf

if ~arctic
 piomas=false;
 icesat=false;
end

col=colormap(lines(length(runcell))); % colorscheme 
if length(runcell)==2
 col(1,:)=[0 0 1];
 col(2,:)=[1 0 0];
elseif length(runcell)==3
 col(1,:)=[0 0 1];
 col(2,:)=[0 1 0];
 col(3,:)=[1 0 0];
end

%if nargin<6
 pub=false
%end

satcol=0.0*[1 1 1]; % color of Passive Microwave Obs

if cesmle
 cesmcol=0.0*[1 1 1]; % color of CESM results
 cesmdir='/Users/aroberts/data/MODEL/E3SM';
 leg{length(leg)+1}='CESM';
end

if piomas
 piomascol=0.0*[1 1 1]; % color of CESM results
 piomasdir='/Volumes/RobertsRaid3/data_1_NPS_August_2018/MODEL/PIOMAS'
end

lw=0.5; % linewidth

if icesat
 % set ICESat file and mean date settings
 islist={'h_fm04','h_fm05','h_fm06','h_ma07','h_fm08','h_on03','h_on04','h_on05','h_on06','h_on07'};
 ismeanyear=[2004 2005 2006 2007 2008 2003 2004 2005 2006 2007];
 ismeanmonth=[3 3 3 4 3 11 11 11 11 11];
 ismeanday=[1 1 1 1 1 1 1 1 1 1];
 ismeandates=datenum(ismeanyear,ismeanmonth,ismeanday);
 
 isstartmonth=[2 2 2 3 2 10 10 10 10 10];
 isstartday=[1 1 1 1 1 1 1 1 1 1];
 isstartdates=datenum(ismeanyear,isstartmonth,isstartday);

 isendmonth=[3 3 3 4 3 11 11 11 11 11];
 isendday=[31 31 31 30 31 30 30 30 30 30];
 isenddates=datenum(ismeanyear,isendmonth,isendday);

 % load in ICESat data, and add area of grid cells
 ncis=ridgepack_clone('/Volumes/RobertsRaid3/data_1_NPS_August_2018/SATELLITE/processed/kwok_icesat.nc');
end

xemax=0;
xemin=100;
xamax=0;
xamin=100;
xvmax=0;
xvmin=100;
xsvmax=0;
xsvmin=100;

hiso=NaN; % handle for icesat observations

for j=1:length(runcell)

 run=char(runcell{j});

 % pull in grid information
 if strcmp(run,'E3SM-HR-V1')

  cd('/Users/afroberts/SIhMSatArray/E3SM/highres/grid')
  ncgrid=ridgepack_clone('E3SM_hr_grid',{'latVertex','lonVertex',...
           'verticesOnCell','indexToCellID','nEdgesOnCell','areaCell'});
  nclat=ridgepack_clone('E3SM_hr_grid',{'latCell'});
  nclon=ridgepack_clone('E3SM_hr_grid',{'lonCell'});
  basedir=['/Users/afroberts/SIhMSatArray/E3SM/highres/'];

 elseif strcmp(run,'E3SM-LR-V1')

  cd('/Users/afroberts/SIhMSatArray/E3SM/lrhrequiv/grid')
  ncgrid=ridgepack_clone('E3SM_LR_V1_grid',{'latVertex','lonVertex',...
           'verticesOnCell','indexToCellID','nEdgesOnCell','areaCell'});
  nclat=ridgepack_clone('E3SM_LR_V1_grid',{'latCell'});
  nclon=ridgepack_clone('E3SM_LR_V1_grid',{'lonCell'});
  basedir=['/Users/afroberts/SIhMSatArray/E3SM/lrhrequiv/'];

 elseif strcmp(run,'E3SM-DECK-PI')

  cd('/Users/afroberts/SIhMSatArray/E3SM/DECK/grid')
  ncgrid=ridgepack_clone('E3SM_LR_V1_grid',{'latVertex','lonVertex',...
           'verticesOnCell','indexToCellID','nEdgesOnCell','areaCell'});
  nclat=ridgepack_clone('E3SM_LR_V1_grid',{'latCell'});
  nclon=ridgepack_clone('E3SM_LR_V1_grid',{'lonCell'});
  basedir=['/Users/afroberts/SIhMSatArray/E3SM/DECK/monthly/PI/archive/ice/hist'];

 end


 %%%%%%%%%%%%% CONFIGURATION %%%%%%%%%%%%%%

 publish=0; % publish plots to Dropbox

 if isempty(run);
  error('Casename is empty');
 else
  disp(['CASE: ',run]);
 end
 
 rrun=run;

 filetype1='mpascice.hist.am.timeSeriesStatsMonthly.';

 afield='timeMonthly_avg_iceAreaCell';
 sfield='timeMonthly_avg_iceVolumeCell';
 snowfield='timeMonthly_avg_snowVolumeCell';

 satdir=['/Volumes/RobertsRaid3/data_1_NPS_August_2018/SATELLITE/processed'];
 if arctic
  filetype5='G02202_v3_merged_conc_north_1979_2017_extent_area.nc';
 else
  filetype5='G02202_v3_merged_conc_south_1979_2017_extent_area.nc';
 end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if strcmp(run,'E3SM-DECK-PI')
  cd(basedir)
 else 
  cd([basedir,'/monthly'])
 end


 % read in monthly E3SM data
 startyear=str2num(datestr(mintime,'yyyy'));
 endyear=str2num(datestr(maxtime,'yyyy'))-1;

 k=1;
 for year=startyear:endyear
  for month=1:12
   nctemp=ridgepack_clone([filetype1,...
             num2str(year,'%4.4i'),'-',num2str(month,'%2.2i'),'-01'],...
             {afield,sfield,snowfield});
   if year==startyear & month==1
    nc=nctemp;
   else 
    nc.(afield).data(:,k)=nctemp.(afield).data(:);
    nc.(sfield).data(:,k)=nctemp.(sfield).data(:);
    nc.(snowfield).data(:,k)=nctemp.(snowfield).data(:);
   end
   nc.time.data(k)=datenum(year,month+1,1);
   k=k+1;
  end
 end
 clear nctemp
 nc=rmfield(nc,'attributes');
 nc.attributes.title=['E3SM Sea ice area,extent and volume for ',rrun];
 nc.time.calendar='proleptic_gregorian';
 nc.latitude=nclat.latCell;
 nc.longitude=nclon.lonCell;
 nc=ridgepack_struct(nc);   

 % set up timeseries netcdf file
 nce.attributes.title=['E3SM Sea ice area,extent and volume for ',rrun];
 nce.time=nc.time;

 % center model values in time (use mean center value for first sample).
 nce.time.data(1)=nce.time.data(1)-mean(diff(nce.time.data)/2);
 nce.time.data(2:end)=nce.time.data(2:end)-diff(nce.time.data)/2;

 % remove no-leap tag
 %nce.time=rmfield(nce.time,'calendar') ;
 nce=ridgepack_struct(nce);

 % set up extent, area and volume graph
 nce.extent.long_name=['E3SM Sea ice extent for ',rrun];
 nce.extent.units='km^2';
 nce.extent.dimension=nce.time.dimension;
 nce.extent.data=zeros(size(nce.time.data));

 nce.area=nce.extent;
 nce.area.long_name=['E3SM Sea ice area (south of 87^{\circ}N) for ',rrun];

 nce.area84=nce.extent;
 nce.area84.long_name=['E3SM Sea ice area (south of 84^{\circ}N) for ',rrun];

 nce.volume=nce.extent;
 nce.volume.units='km^3';
 nce.volume.long_name=['E3SM Sea ice volume for ',rrun];

 nce.snowvolume=nce.extent;
 nce.snowvolume.units='km^3';
 nce.snowvolume.long_name=['E3SM Sea ice snow volume for ',rrun];

 if icesat
  for k=1:length(ismeanyear);
   tlist{k}=[];

   nce.(char(islist{k})).long_name='ICESat Sea ice volume';
   nce.(char(islist{k})).units='km^3';
   nce.(char(islist{k})).dimension={};
   nce.(char(islist{k})).data=[];
   nce.(char(islist{k})).time=ismeandates(k);

   nce.([char(islist{k}),'_model']).long_name='Sea ice volume on ICESat Grid';
   nce.([char(islist{k}),'_model']).units='km^3';
   nce.([char(islist{k}),'_model']).dimension={};
   nce.([char(islist{k}),'_model']).data=[];
   nce.([char(islist{k}),'_model']).count=0;
   nce.([char(islist{k}),'_model']).time=ismeandates(k);
  end
 end

 for i=1:length(nce.time.data);

  if arctic
   idxe=find(nc.(afield).data(:,i)>0.15 & nclat.latCell.data>30*pi/180);
  else
   idxe=find(nc.(afield).data(:,i)>0.15 & nclat.latCell.data<-30*pi/180);
  end

  cellextent=ncgrid.areaCell.data(idxe); 
  nce.extent.data(i)=sum(cellextent(:))./10^6;
 
  dvec=datevec(nce.time.data(i));
  nce.area.data(i)=NaN;
  nce.area84.data(i)=NaN;

  if arctic
   if dvec(1)<1988 
    idx=find(nc.(afield).data(:,i)>0.15 & ...
             nclat.latCell.data>30*pi/180 & ...
             nclat.latCell.data<84*pi/180);
    cellarea=ncgrid.areaCell.data(idx).*nc.(afield).data(idx,i);
    nce.area84.data(i)=sum(cellarea(:))./10^6;
   else
    idx=find(nc.(afield).data(:,i)>0.15 & ...
             nclat.latCell.data>30*pi/180 & ...
             nclat.latCell.data<87*pi/180);
    cellarea=ncgrid.areaCell.data(idx).*nc.(afield).data(idx,i);
    nce.area.data(i)=sum(cellarea(:))./10^6;
   end
  else
   idx=find(nc.(afield).data(:,i)>0.15 & nclat.latCell.data<-30*pi/180);
   cellarea=ncgrid.areaCell.data(idx).*nc.(afield).data(idx,i);
   nce.area.data(i)=sum(cellarea(:))./10^6;
  end

  cellvolume=ncgrid.areaCell.data(idxe).*squeeze(nc.(sfield).data(idxe,i));
  nce.volume.data(i)=sum(cellvolume(~isnan(cellvolume)))./10^9;

  cellvolume=ncgrid.areaCell.data(idxe).*squeeze(nc.(snowfield).data(idxe,i));
  nce.snowvolume.data(i)=sum(cellvolume(~isnan(cellvolume)))./10^9;

  if icesat
   for k=1:length(ismeanyear);
    if nce.time.data(i)>=isstartdates(k) & (isempty(mintime) | ...
      nce.time.data(i)>=mintime) & ...
      nce.time.data(i)<=isenddates(k) & (isempty(maxtime) | ...
      nce.time.data(i)<=maxtime)
      tlist{k}=[tlist{k} i];
    end
   end
  end
  
 end

 if runnin
  nce.volume.data=nce.volume.data';
  [nce]=ridgepack_runningmean(nce,'volume',12)
  nce.volume.data=nce.volume.data';
 end

 % construct ICESat summaries
 if icesat
  icesat_obs=[];
  icesat_mod=[];
  icesat_tim=[];
  for k=1:length(ismeanyear);
   if ~isempty(tlist{k})
    disp(['ICESat comparison: ',datestr(nce.time.data(min(tlist{k}))),' to ',...
          datestr(nce.time.data(max(tlist{k}))),' for ',char(islist{k})])
 
    ncmodsat=ridgepack_regrid(ridgepack_reduce(nc,{'time'},{[min(tlist{k}) max(tlist{k})]}),sfield,'',ncis);  
    ncmodsat.(sfield).data=ncmodsat.(sfield).data.*ncmodsat.cell_area.data./10^3;
    nce.([char(islist{k}),'_model']).data=sum(ncmodsat.(sfield).data(~isnan(ncis.(char(islist{k})).data)));
 
    icesatdata=ncis.(char(islist{k})).data.*ncis.cell_area.data./10^3;
    nce.(char(islist{k})).data=sum(icesatdata(~isnan(icesatdata)));
 
    icesat_mod=[icesat_mod nce.([char(islist{k}),'_model']).data];
    icesat_obs=[icesat_obs nce.(char(islist{k})).data];
    icesat_tim=[icesat_tim ismeandates(k)];
   else
    nce=rmfield(nce,[char(islist{k}),'_model']);
    nce=rmfield(nce,char(islist{k}));
   end
  end
 end
 
 if j==1; 
  ridgepack_multiplot(4,1,1,1,'a'); 
  h1=gca;
 else
  set(gcf,'CurrentAxes',h1)
 end
 h(j)=plot(nce.time.data,nce.extent.data/(10^6),'Color',col(j,:),'LineWidth',lw)
 xemax=max([xemax nce.extent.data/(10^6)]);
 xemin=min([xemin nce.extent.data/(10^6)]);
 ridgepack_clearax('x',1)
 hold on
 drawnow

 if j==1; 
  ridgepack_multiplot(4,1,2,1,'b'); 
  h2=gca; 
 else
  set(gcf,'CurrentAxes',h2)
 end
 plot(nce.time.data,nce.area84.data/(10^6),'Color',col(j,:),'LineWidth',lw)
 plot(nce.time.data,nce.area.data/(10^6),'Color',col(j,:),'LineWidth',lw)
 xamax=max([xamax nce.area.data/(10^6) nce.area84.data/(10^6)]);
 xamin=min([xamin nce.area.data/(10^6) nce.area84.data/(10^6)]);
 ridgepack_clearax('x',1)
 hold on
 drawnow

 if j==1; 
  ridgepack_multiplot(4,1,3,1,'c'); 
  h3=gca; 
 else
  set(gcf,'CurrentAxes',h3)
 end
 plot(nce.time.data,nce.volume.data/(10^3),'Color',col(j,:),'LineWidth',lw)
 if runnin
  plot(nce.time.data,nce.volume_running.data/(10^3),':','Color',col(j,:),'LineWidth',lw)
 end
 xvmax=max([xvmax nce.volume.data/(10^3)]);
 xvmin=min([xvmin nce.volume.data/(10^3)]);
 ridgepack_clearax('x',1)
 hold on;
 if icesat
  if ~isempty(icesat_mod)
   plot(icesat_tim,icesat_mod/(10^3),'h','Color',col(j,:),'MarkerSize',4)
   xvmax=max([xvmax icesat_mod'/(10^3)]);
   xvmin=min([xvmin icesat_mod'/(10^3)]);
  end
  if ~isempty(icesat_obs) 
   hiso=plot(icesat_tim,icesat_obs/(10^3),'kh','MarkerSize',4);
   xvmax=max([xvmax icesat_obs'/(10^3)]);
   xvmin=min([xvmin icesat_obs'/(10^3)]);
  end
 end
 drawnow

 if j==1;
  ridgepack_multiplot(4,1,4,1,'d');
  h4=gca;
 else
  set(gcf,'CurrentAxes',h4)
 end
 plot(nce.time.data,nce.snowvolume.data/(10^2),'Color',col(j,:),'LineWidth',lw)
 xsvmax=max([xsvmax nce.snowvolume.data/(10^2)]);
 xsvmin=min([xsvmin nce.snowvolume.data/(10^2)]);
 ridgepack_clearax('x',1)
 hold on;
 drawnow

 if j==1
  if isempty(mintime); mintime=nce.time.data(1); end
  if isempty(maxtime); maxtime=nce.time.data(end); end
 else
  if isempty(mintime); mintime=min(mintime,nce.time.data(1)); end
  if isempty(maxtime); maxtime=max(maxtime,nce.time.data(end)); end
 end

 %if strcmp(filetype1(end:end),'.')
 if arctic
  ridgepack_write(nce,[rrun,'_',filetype1,'north_volume_area_extent'])
 else
  ridgepack_write(nce,[rrun,'_',filetype1,'south_volume_area_extent'])
 end

end

if cesmle
 cd([cesmdir])
 ncesm=ridgepack_clone('CESM_perpetual.cice.h.northern_volume_area_extent.nc')
 xvmax=max([xvmax; ncesm.volume.data/(10^3)]);
 xvmin=min([xvmin; ncesm.volume.data/(10^3)]);
end

if piomas
 cd([piomasdir])
 ncpiomas=ridgepack_clone('PIOMAS_Volume.nc',{'volume'},{'time'},...
                  {datestr(mintime)},{datestr(maxtime)})
 xvmax=max([xvmax; ncpiomas.volume.data/(10^3)]);
 xvmin=min([xvmin; ncpiomas.volume.data/(10^3)]);
end

if passive
 cd([satdir])
 if arctic
  ncc=ridgepack_clone(filetype5,{'extent','area','area84'},{'time'},...
             {datestr(mintime)},{datestr(maxtime)},15);
 else
  ncc=ridgepack_clone(filetype5,{'extent','area'},{'time'},...
             {datestr(mintime)},{datestr(maxtime)},15);
 end

 xemax=max([xemax; ncc.extent.data/(10^6)]);
 xemin=min([xemin; ncc.extent.data/(10^6)]);

 if arctic
  xamax=max([xamax; ncc.area.data/(10^6); ncc.area84.data/(10^6)]);
  xamin=min([xamin; ncc.area.data/(10^6); ncc.area84.data/(10^6)]);
 else
  xamax=max([xamax; ncc.area.data/(10^6)]);
  xamin=min([xamin; ncc.area.data/(10^6)]);
 end

 cd([basedir,'/monthly'])
end
 
set(gcf,'CurrentAxes',h1)
if passive
 hs=plot(ncc.time.data,ncc.extent.data/(10^6),'Color',satcol,...
        'LineStyle','--','LineWidth',lw);
end

if cesmle
 hc=plot(ncesm.time.data,ncesm.extent.data/(10^6),'Color',cesmcol,'LineStyle',':');
end
if passive
 legend(hs,{'Passive Microwave Extent'},'Location','NorthEast','FontSize',7)
 legend('boxoff')
end
xlim([mintime maxtime])
%ylim([floor(xemin) ceil(xemax)])
ylim([max(0.,floor(xemin)-1) ceil(xemax)+0.25*(ceil(xemax)-max(0.,floor(xemin)-1))])
datetick('x','YY','keeplimits')
ridgepack_clearax('x',1)
ylabel('$\mathrm{Extent\,(\times\,10^{6}\,km^{2})}$')
set(h1,'box','on')
hold off

set(gcf,'CurrentAxes',h2)
if passive
 hsa=plot(ncc.time.data,ncc.area.data/(10^6),'Color',satcol,...
         'LineStyle','--','LineWidth',lw);
 if arctic
  plot(ncc.time.data,ncc.area84.data/(10^6),'Color',satcol,...
         'LineStyle','--','LineWidth',lw);
 end
end

if cesmle
 hc=plot(ncesm.time.data,ncesm.area.data/(10^6),'Color',cesmcol,'LineStyle',':');
 plot(ncesm.time.data,ncesm.area84.data/(10^6),'Color',cesmcol,'LineStyle',':');
end
if arctic & passive
 if mintime<datenum(1987,12,1)
  legend(hsa,{'NOAA CDR Area South of 87^{\circ} N (84^{\circ} N up to 1988) and model equivalent'},...
              'Location','NorthEast','FontSize',7)
 else
  legend(hsa,{'NOAA CDR Area South of 87^{\circ} N and model equivalent'},...
             'Location','NorthEast','FontSize',7)
 end
elseif passive
 legend(hsa,{'NOAA Climate Data Record'},'Location','NorthEast','FontSize',7)
 legend('boxoff')
end
xlim([mintime maxtime])
%ylim([floor(xamin) ceil(xamax)])
ylim([max(0.,floor(xamin)-1) ceil(xamax)+0.25*(ceil(xamax)-max(0.,floor(xamin)-1))])
datetick('x','YY','keeplimits')
ridgepack_clearax('x',1)
ylabel('$\mathrm{Area\,(\times\,10^{6}\,km^{2})}$')
set(h2,'box','on')
hold off

set(gcf,'CurrentAxes',h3)
if piomas
 hpv=plot(ncpiomas.time.data,ncpiomas.volume.data/(10^3),'Color',piomascol,...
          'LineStyle','--','LineWidth',lw);
end
if ishandle(hiso) & piomas
 legend([hpv hiso],{'PIOMAS','ICESat Arctic Ocean volume and model equivalent'},'Location','NorthEast','FontSize',7,'Orientation','horizontal')
 legend('boxoff')
elseif piomas
 legend([hpv],{'PIOMAS'},'Location','NorthEast','FontSize',7)
 legend('boxoff')
elseif ishandle(hiso) 
 legend(hiso,{'ICESat Arctic Ocean volume and model equivalent'},'Location','NorthEast','FontSize',7)
 legend('boxoff')
end
if cesmle
 hc=plot(ncesm.time.data,ncesm.volume.data/(10^3),'Color',cesmcol,'LineStyle',':');
end
xlim([mintime maxtime])
%ylim([floor(xvmin) ceil(xvmax)])
ylim([max(0.,floor(xvmin)-1) ceil(xvmax)+0.25*(ceil(xvmax)-max(0.,floor(xvmin)-1))])
datetick('x','YY','keeplimits')
ridgepack_clearax('x',1)
ylabel('$\mathrm{Volume\,(\times\,10^{3}\,km^{3})}$')
set(h3,'box','on')
hold off

set(gcf,'CurrentAxes',h4)
if cesmle
 hc=plot(ncesm.time.data,ncesm.snowvolume.data/(10^2),'Color',cesmcol,'LineStyle',':');
end

% Turn off normal legend
legend(h4,leg,'Location','NorthEast','FontSize',7,'Orientation','horizontal')
legend('boxoff')

xlim([mintime maxtime])
%ylim([floor(xsvmin) ceil(xsvmax)])
ylim([max(0.,floor(xsvmin)-1) ceil(xsvmax)+1])
datetick('x','YY','keeplimits')
xlabel(['Year'])
ylabel('$\mathrm{Snow\;Volume\,(\times\,10^{3}\,km^{3})}$')
set(h4,'box','on')
hold off

%if cesmle 
% ridgepack_multilegend([h hc],leg,'North')
%else
% ridgepack_multilegend([h],leg,'North')
%end

if pub
 ridgepack_multialign(gcf)
 cd(pubdir)
else 
 if arctic
  ridgepack_multialign(gcf,['E3SM-HR 1950 Cycled V1 Arctic Sea Ice ',...
                  datestr(mintime,'YYYY'),'-',datestr(maxtime,'YYYY')],9)
 else
  ridgepack_multialign(gcf,['E3SM-HR 1950 Cycled V1 Antarctic Sea Ice ',...
                  datestr(mintime,'YYYY'),'-',datestr(maxtime,'YYYY')],9)
 end
end

if arctic
 nams='_arctic_seaice_volume_graph_';
else
 nams='_antarctic_seaice_volume_graph_';
end

ridgepack_fprint('png',[ridgepack_cellcat(runcell,'_'),nams,datestr(mintime,'YYYY'),'_',datestr(maxtime,'YYYY'),'_',num2str(minthick),'.png'],1,1)
ridgepack_fprint('epsc',[ridgepack_cellcat(runcell,'_'),nams,datestr(mintime,'YYYY'),'_',datestr(maxtime,'YYYY'),'_',num2str(minthick),'.eps'],1,1)



