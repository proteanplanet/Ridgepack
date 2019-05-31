
clear
close all

runcell={'b1850c5_acmev0_highres'};
leg={'E3SM V0'};

config='E3SM V0';

%hemisphere=1;
hemisphere=-1;

startyear=1;
%endyear=60;
endyear=131;

mintime=[];
maxtime=[];

%col=colormap(lines(length(runcell))); % colorscheme 
col(1,:)=[1 0 0]; % colorscheme 
col(2,:)=[0 0 1]; % colorscheme 
col(3,:)=[0 1 0]; % colorscheme 
col(4,:)=[1 1 0]; % colorscheme 
col(5,:)=[0 1 1]; % colorscheme 
col(6,:)=0.5*[1 1 1]; % colorscheme 
col(7,:)=0.0*[1 1 1]; % colorscheme 
colormap(col);
lw=0.5; % linewidth

if hemisphere>0
 hemis='Northern';
else
 hemis='Southern';
end

xemax=0;
xemin=100;
xamax=0;
xamin=100;
xvmax=0;
xvmin=100;
xsvmax=0;
xsvmin=100;

% build timeseries
nct.time.long_name='time';
nct.time.dimension={'time'};
nct.time.units='days since 0000-01-01 00:00:00';
nct.time.calendar='noleap';
nct.d2.long_name='d2';
nct.d2.dimension={'d2'};
nct.d2.units='';
nct.d2.data=[1 2];
nct.time_bounds.long_name='time bounds';
nct.time_bounds.dimension={'d2','time'};
nct.time_bounds.units='days since 0000-01-01 00:00:00';
nct.time_bounds.calendar='noleap';
k=0;
for year=startyear:endyear
 for month=1:12
  k=k+1;
  nct.time.data(k)=...
     (datenum(year,month,1)+datenum(year,month+1,1))/2;
  nct.time_bounds.data(1,k)=datenum(year,month,1);
  nct.time_bounds.data(2,k)=datenum(year,month+1,1);
 end
end

for j=1:length(runcell)

 run=char(runcell{j});

 %%%%%%%%%%%%% CONFIGURATION %%%%%%%%%%%%%%

 if isempty(run);
  error('Casename is empty');
 else
  disp(['CASE: ',run]);
 end
 basedir=['/Users/afroberts/SIhMSatArray/E3SM/v0highres/monthly'];
 %basedir=['/Volumes/MacBookProT3SeaIce/E3SM/highresv0/monthly'];
 
 rrun=run;
 filetype1='.cice.h.';

 afield='aice';
 hfield='hi';
 sfield='hs';

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 cd([basedir])

 % get area
 ncx=ridgepack_clone(...
    'hi.b1850c5_acmev0_highres.cice.h.0130-12',...
    {'tarea','TLAT','TLON'});

 % set up timeseries netcdf file
 nce.attributes.title=[rrun,' ',hemis,...
           ' sea ice area, extent and volume'];
 nce.attributes.source=config;

 nce.time=nct.time;
 nce.d2=nct.d2;
 nce.time_bounds=nct.time_bounds;

 % set up extent, area and volume graph
 nce.extent.long_name=['v0E3SM sea ice extent for ',rrun];
 nce.extent.units='km^2';
 nce.extent.dimension=nce.time.dimension;
 nce.extent.data=zeros(size(nce.time.data));

 nce.area=nce.extent;
 nce.area.long_name=['v0E3SM sea ice area for ',rrun];

 nce.volume=nce.extent;
 nce.volume.units='km^3';
 nce.volume.long_name=['v0E3SM sea ice volume for ',rrun];

 nce.snowvolume=nce.extent;
 nce.snowvolume.units='km^3';
 nce.snowvolume.long_name=['v0E3SM Sea ice snow volume for ',rrun];

 k=0;
 for year=startyear:endyear
  for month=1:12
   k=k+1;
   nca=ridgepack_clone([afield,'.',rrun,filetype1,...
           num2str(year,'%4.4i'),'-',num2str(month,'%2.2i')]);
   nch=ridgepack_clone([hfield,'.',rrun,filetype1,...
           num2str(year,'%4.4i'),'-',num2str(month,'%2.2i')]);
   ncs=ridgepack_clone([sfield,'.',rrun,filetype1,...
           num2str(year,'%4.4i'),'-',num2str(month,'%2.2i')]);

   cellextent=ncx.tarea.data(nca.(afield).data(:,:)>15 ...
                          & hemisphere*ncs.latitude.data(:,:)>0);
   nce.extent.data(k)=sum(cellextent(:))./10^6;
 
   cellarea=ncx.tarea.data.*squeeze(nca.(afield).data(:,:))/100;
   cellarea=cellarea(nca.(afield).data(:,:)>15 ...
                        & hemisphere*nca.latitude.data(:,:)>0); 
   nce.area.data(k)=sum(cellarea(:))./10^6;
 
   cellvolume=ncx.tarea.data.*nch.(hfield).data(:,:);
   cellvolume=cellvolume(nca.(afield).data(:,:)>15 ...
                        & hemisphere*nca.latitude.data(:,:)>0); 
   nce.volume.data(k)=sum(cellvolume(~isnan(cellvolume)))./10^9;

   cellvolume=ncx.tarea.data.*ncs.(sfield).data(:,:);
   cellvolume=cellvolume(nca.(afield).data(:,:)>15 ...
                        & hemisphere*nca.latitude.data(:,:)>0);
   nce.snowvolume.data(k)=sum(cellvolume(~isnan(cellvolume)))./10^9;

  end
 end

 % filter volume timeseries into annual 
 ts=ridgepack_2timeseries(nce,'volume','linear',[0 0 1 0 0 0.0])
 try
  ts=ridgepack_filter1(ts,1/(2*365),'FIR','low')
  voltime=ts.time(365:end-365);
  voldata=ts.data(365:end-365);
 catch
  voltime=NaN;
  voldata=NaN;
 end

 % remove no-leap tag
 nce=ridgepack_struct(nce);
 
 if j==1; 
  ridgepack_multiplot(4,1,1,1); 
  h1=gca;
 else
  set(gcf,'CurrentAxes',h1)
 end
 h(j)=plot(nce.time.data,nce.extent.data/(10^6),...
                'Color',col(j,:),'LineWidth',lw);
 xemax=max([xemax nce.extent.data/(10^6)]);
 xemin=min([xemin nce.extent.data/(10^6)]);
 ridgepack_clearax('x',1)
 hold on
 drawnow

 if j==1; 
  ridgepack_multiplot(4,1,2,1); 
  h2=gca; 
 else
  set(gcf,'CurrentAxes',h2)
 end
 plot(nce.time.data,nce.area.data/(10^6),...
               'Color',col(j,:),'LineWidth',lw);
 xamax=max([xamax nce.area.data/(10^6)]);
 xamin=min([xamin nce.area.data/(10^6)]);
 ridgepack_clearax('x',1)
 hold on
 drawnow

 if j==1; 
  ridgepack_multiplot(4,1,3,1); 
  h3=gca; 
 else
  set(gcf,'CurrentAxes',h3)
 end
 plot(nce.time.data,nce.volume.data/(10^3),...
                'Color',col(j,:),'LineWidth',lw);
 xvmax=max([xvmax nce.volume.data/(10^3)]);
 xvmin=min([xvmin nce.volume.data/(10^3)]);
 ridgepack_clearax('x',1)
 hold on;
 plot(voltime,voldata/(10^3),'Color',col(j,:),...
                'LineWidth',lw,'LineStyle','--');
 drawnow

 if j==1;
  ridgepack_multiplot(4,1,4,1);
  h4=gca;
 else
  set(gcf,'CurrentAxes',h4)
 end
 plot(nce.time.data,nce.snowvolume.data/(10^2),...
                'Color',col(j,:),'LineWidth',lw);
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

 ridgepack_write(nce,[rrun,filetype1,hemis,'_volume_area_extent'])

end

cd([basedir])

set(gcf,'CurrentAxes',h1)
xlim([mintime maxtime])
ylim([max(0.,floor(xemin)-1) ceil(xemax)+1])
datetick('x','YY','keeplimits')
ridgepack_clearax('x',1)
ylabel('Extent ($10^{6}$ km$^2$)')
set(h1,'box','on')
hold off

set(gcf,'CurrentAxes',h2)
xlim([mintime maxtime])
ylim([max(0.,floor(xamin)-1) ceil(xamax)+1])
datetick('x','YY','keeplimits')
ridgepack_clearax('x',1)
ylabel('Area ($10^{6}$ km$^2$)')
set(h2,'box','on')
hold off

set(gcf,'CurrentAxes',h3)
xlim([mintime maxtime])
ylim([max(0.,floor(xvmin)-1) ceil(xvmax)+1])
datetick('x','YY','keeplimits')
xlabel(['Year'])
ylabel('Volume ($10^{3}$ km$^3$)');
set(h3,'box','on')
hold off

set(gcf,'CurrentAxes',h4)
xlim([mintime maxtime])
ylim([max(0.,floor(xsvmin)-1) ceil(xsvmax)+1])
datetick('x','YY','keeplimits')
xlabel(['Year'])
ylabel('Snow Volume ($10^{2}$ km$^3$)');
set(h4,'box','on')
legend(h,leg,'location','NorthEast','Orientation','horizontal')
hold off

ridgepack_multialign(gcf,[hemis,' hemisphere sea ice ',config])

pname=[hemis,'_',run,'_icevol_',datestr(mintime,'YYYY'),'_',datestr(maxtime,'YYYY')];

ridgepack_fprint('png',pname,1,1)
ridgepack_fprint('epsc',pname,1,1)



