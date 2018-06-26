function Ridgepack_RASM_sea_ice_extent_volume(runcell,leg,minthick,mintime,maxtime,pub,pubdir)

% Ridgepack_RASM_sea_ice_extent_volume - produces summary graph of RASM results
%
% This function prepares an analysis plot of RASM sea ice output equivalent
% to the example ridgepack_example_seaice_volume_graph_1980_2010_0.png.
% The script uses a wide array of observational and model data to construct the plot:
%
% 1) ICESat data from JPL
% 2) PIOMAS data from University of Washington
% 3) Ice area and extent from NSIDC NOAA CDR
% 4) RASM ice area (aice), thickness (hi) and snow depth (hs)
%
% The RASM data is assumed to be in the form of a single file for each 
% variable, organized as a timeseries using the 'rasm_cice5_series.bash'
% script on the chosen supercomputer used to integrate RASM. That script
% is located under the model/bash directory of the RASM Sea Ice Toolbox
% in Ridgepack.  The data to produce this plot is located in various places,
% as indicated by the satdir, cesmdir, piomasdir, basedir settings within the 
% script.
%
%
% INPUT:
%
% runcell  - Cell array of RASM or CESM case names to be included on the graph
%            e.g. {'R1009RBRceap01a','R2100aRBRcaaa01a'} for timeseries datafiles
%            with the name, e.g. 'R1009RBRceap01a.cice.h.hi.nc' for ice thickness
%           'hi'.
%
% leg      - Name of simulations to be placed on the plot matching for each of
%            the runs in runcell, e.g. {'RASM 1.1','RASM 2.1'}.
%
% minthick - Minimum thickness of ice to be included in extent and area calculations
% mintime  - Lower time bound on the graph, in MATLAB time format
% maxtime  - Uppder time bound on the graph, in MATLAB time format
% pub      - logical stating whether or not this figure is for publication.
%            If this is set to true, then the header is removed, and the file
%            is written on the pubdir directory instead of the directory of the last
%            case name.
% pubdir   - Directory into which the graphical output is written if pub=true.
%
%
% Andrew Roberts, Naval Postgraduate School, June 2018 (afrobert@nps.edu)

% switch on/off CESM large ensemble overplot
cesmle=false;

% switch on/off piomass data
piomas=true;

clf

col=colormap(lines(length(runcell))); % colorscheme 
if length(runcell)==2
 col(1,:)=[0 0 1];
 col(2,:)=[1 0 0];
elseif length(runcell)==3
 col(1,:)=[0 0 1];
 col(2,:)=[0 1 0];
 col(3,:)=[1 0 0];
end

if nargin<6
 pub=false
end


satcol=0.5*[1 1 1]; % color of Passive Microwave Obs

if cesmle
 cesmcol=0.0*[1 1 1]; % color of CESM results
 cesmdir='/Users/aroberts/data/MODEL/RASM';
 leg{length(leg)+1}='CESM';
end

if piomas
 piomascol=0.0*[1 1 1]; % color of CESM results
 piomasdir='/Users/aroberts/data/MODEL/PIOMAS'
end

lw=0.5; % linewidth

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
ncis=ridgepack_clone([getenv('HOME'),'/data/SATELLITE/processed/kwok_icesat.nc']);

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

 %%%%%%%%%%%%% CONFIGURATION %%%%%%%%%%%%%%

 publish=0; % publish plots to Dropbox

 if isempty(run);
  error('Casename is empty');
 else
  disp(['CASE: ',run]);
 end
 basedir=[getenv('HOME'),'/work/processing/'];
 
 rrun=run;
 %filetype1='.cice.hp.';
 filetype1='.cice.h.';
 filetype1a='.cice.h1.';

 afield='aice';
 sfield='hi';
 snowfield='hs';

 satdir=[getenv('HOME'),'/data/SATELLITE/processed'];
 filetype5='G02202_v3_merged_conc_north_1979_2017_extent_area.nc';

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 try
  cd([basedir,run,'/ice/monthly'])
 catch
  cd([basedir,'r',run(2:end),'/ice/monthly'])
 end

 try
  nca=ridgepack_clone([rrun,filetype1,afield],{afield},{'time'},...
               {datestr(mintime)},{datestr(maxtime)},'year');
  ncs=ridgepack_clone([rrun,filetype1,sfield],{sfield},{'time'},...
               {datestr(mintime)},{datestr(maxtime)},'year');
  ncsnow=ridgepack_clone([rrun,filetype1,snowfield],{snowfield},{'time'},...
               {datestr(mintime)},{datestr(maxtime)},'year');
 catch
  nca=ridgepack_clone([rrun,filetype1a,afield],{afield},{'time'},...
               {datestr(mintime)},{datestr(maxtime)},'year');
  ncs=ridgepack_clone([rrun,filetype1a,sfield],{sfield},{'time'},...
               {datestr(mintime)},{datestr(maxtime)},'year');
  ncsnow=ridgepack_clone([rrun,filetype1a,snowfield],{snowfield},{'time'},...
               {datestr(mintime)},{datestr(maxtime)},'year');
 end

 % set up timeseries netcdf file
 nce.attributes.title=['RACM Sea ice area,extent and volume for ',rrun];
 if length(nca.time.data)<length(ncs.time.data)
  nce.time=nca.time;
 else
  nce.time=ncs.time;
 end

 % center model values in time (use mean center value for first sample).
 nce.time.data(1)=nce.time.data(1)-mean(diff(nce.time.data)/2);
 nce.time.data(2:end)=nce.time.data(2:end)-diff(nce.time.data)/2;

 % remove no-leap tag
 nce.time=rmfield(nce.time,'calendar') ;
 nce=ridgepack_struct(nce);

 % set up extent, area and volume graph
 nce.extent.long_name=['RACM Sea ice extent for ',rrun];
 nce.extent.units='km^2';
 nce.extent.dimension=nce.time.dimension;
 nce.extent.data=zeros(size(nce.time.data));

 nce.area=nce.extent;
 nce.area.long_name=['RACM Sea ice area (south of 87^{\circ}N) for ',rrun];

 nce.area84=nce.extent;
 nce.area84.long_name=['RACM Sea ice area (south of 84^{\circ}N) for ',rrun];

 nce.volume=nce.extent;
 nce.volume.units='km^3';
 nce.volume.long_name=['RACM Sea ice volume for ',rrun];

 nce.snowvolume=nce.extent;
 nce.snowvolume.units='km^3';
 nce.snowvolume.long_name=['RACM Sea ice snow volume for ',rrun];

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

 if max(nca.(afield).data(:))>10; 
  oldconc=true; 
 else
  oldconc=false; 
 end

 for i=1:length(nce.time.data);

  if oldconc
   cellextent=nca.tarea.data(squeeze(nca.(afield).data(i,:,:))>15 & ...
                             squeeze(ncs.(sfield).data(i,:,:))>minthick);
  else
   cellextent=nca.tarea.data(squeeze(nca.(afield).data(i,:,:))>0.15 & ...
                             squeeze(ncs.(sfield).data(i,:,:))>minthick);
  end
  nce.extent.data(i)=sum(cellextent(:))./10^6;
 
  if oldconc
   cellarea=nca.tarea.data.*squeeze(nca.(afield).data(i,:,:))/100;
  else
   cellarea=nca.tarea.data.*squeeze(nca.(afield).data(i,:,:));
  end

  dvec=datevec(nce.time.data(i));
  nce.area.data(i)=NaN;
  nce.area84.data(i)=NaN;

  if dvec(1)<1988 
   if oldconc
    cellarea=cellarea(squeeze(nca.(afield).data(i,:,:))>15 & ...
                      nca.latitude.data(:,:)<84 & ...
                      squeeze(ncs.(sfield).data(i,:,:))>minthick);
   else
    cellarea=cellarea(squeeze(nca.(afield).data(i,:,:))>0.15 & ...
                      nca.latitude.data(:,:)<84 & ...
                      squeeze(ncs.(sfield).data(i,:,:))>minthick);
   end
   nce.area84.data(i)=sum(cellarea(:))./10^6;
  else
   if oldconc
    cellarea=cellarea(squeeze(nca.(afield).data(i,:,:))>15 & ...
                      nca.latitude.data(:,:)<87 & ...
                      squeeze(ncs.(sfield).data(i,:,:))>minthick);
   else
    cellarea=cellarea(squeeze(nca.(afield).data(i,:,:))>0.15 & ...
                      nca.latitude.data(:,:)<87 & ...
                      squeeze(ncs.(sfield).data(i,:,:))>minthick);
   end
   nce.area.data(i)=sum(cellarea(:))./10^6;
  end

  cellvolume=ncs.tarea.data.*squeeze(ncs.(sfield).data(i,:,:));
  if oldconc
   cellvolume=cellvolume(squeeze(nca.(afield).data(i,:,:))>15);
  else
   cellvolume=cellvolume(squeeze(nca.(afield).data(i,:,:))>0.15);
  end
  nce.volume.data(i)=sum(cellvolume(~isnan(cellvolume)))./10^9;


  cellvolume=ncs.tarea.data.*squeeze(ncsnow.(snowfield).data(i,:,:));
  if oldconc
   cellvolume=cellvolume(squeeze(nca.(afield).data(i,:,:))>15);
  else
   cellvolume=cellvolume(squeeze(nca.(afield).data(i,:,:))>0.15);
  end
  nce.snowvolume.data(i)=sum(cellvolume(~isnan(cellvolume)))./10^9;

  for k=1:length(ismeanyear);
   if nce.time.data(i)>=isstartdates(k) & (isempty(mintime) | ...
      nce.time.data(i)>=mintime) & ...
      nce.time.data(i)<=isenddates(k) & (isempty(maxtime) | ...
      nce.time.data(i)<=maxtime)
      tlist{k}=[tlist{k} i];
   end
  end
  
 end

 [nce]=ridgepack_runningmean(nce,'volume',12)

 % construct ICESat summaries
 icesat_obs=[];
 icesat_mod=[];
 icesat_tim=[];
 for k=1:length(ismeanyear);
  if ~isempty(tlist{k})
   disp(['ICESat comparison: ',datestr(nce.time.data(min(tlist{k}))),' to ',...
         datestr(nce.time.data(max(tlist{k}))),' for ',char(islist{k})])

   ncmodsat=ridgepack_regrid(ridgepack_reduce(ncs,{'time'},{[min(tlist{k}) max(tlist{k})]}),sfield,'',ncis);  
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
 
 
 if j==1; 
  ridgepack_multiplot(4,1,1,1,'a'); 
  h1=gca;
 else
  set(gcf,'CurrentAxes',h1)
 end
 h(j)=plot(nce.time.data,nce.extent.data/(10^6),'Color',col(j,:),'LineWidth',lw)
 xemax=max([xemax; nce.extent.data/(10^6)]);
 xemin=min([xemin; nce.extent.data/(10^6)]);
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
 xamax=max([xamax; nce.area.data/(10^6); nce.area84.data/(10^6)]);
 xamin=min([xamin; nce.area.data/(10^6); nce.area84.data/(10^6)]);
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
 plot(nce.time.data,nce.volume_running.data/(10^3),':','Color',col(j,:),'LineWidth',lw)
 xvmax=max([xvmax; nce.volume.data/(10^3)]);
 xvmin=min([xvmin; nce.volume.data/(10^3)]);
 ridgepack_clearax('x',1)
 hold on;
 if ~isempty(icesat_mod)
  plot(icesat_tim,icesat_mod/(10^3),'h','Color',col(j,:),'MarkerSize',4)
  xvmax=max([xvmax; icesat_mod'/(10^3)]);
  xvmin=min([xvmin; icesat_mod'/(10^3)]);
 end
 if ~isempty(icesat_obs) 
  hiso=plot(icesat_tim,icesat_obs/(10^3),'kh','MarkerSize',4);
  xvmax=max([xvmax; icesat_obs'/(10^3)]);
  xvmin=min([xvmin; icesat_obs'/(10^3)]);
 end
 drawnow

 if j==1;
  ridgepack_multiplot(4,1,4,1,'d');
  h4=gca;
 else
  set(gcf,'CurrentAxes',h4)
 end
 plot(nce.time.data,nce.snowvolume.data/(10^2),'Color',col(j,:),'LineWidth',lw)
 xsvmax=max([xsvmax; nce.snowvolume.data/(10^2)]);
 xsvmin=min([xsvmin; nce.snowvolume.data/(10^2)]);
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
 ridgepack_write(nce,[rrun,filetype1,'volume_area_extent'])

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

cd([satdir])
ncc=ridgepack_clone(filetype5,{'extent','area','area84'},{'time'},...
            {datestr(mintime)},{datestr(maxtime)},15);

xemax=max([xemax; ncc.extent.data/(10^6)]);
xemin=min([xemin; ncc.extent.data/(10^6)]);

xamax=max([xamax; ncc.area.data/(10^6); ncc.area84.data/(10^6)]);
xamin=min([xamin; ncc.area.data/(10^6); ncc.area84.data/(10^6)]);

cd([basedir,run,'/ice/monthly/'])

set(gcf,'CurrentAxes',h1)
hs=plot(ncc.time.data,ncc.extent.data/(10^6),'Color',satcol,...
        'LineStyle','--','LineWidth',lw);
if cesmle
 hc=plot(ncesm.time.data,ncesm.extent.data/(10^6),'Color',cesmcol,'LineStyle',':');
end
legend(hs,{'Passive Microwave Extent'},'Location','NorthEast','FontSize',7)
legend('boxoff')
xlim([mintime maxtime])
ylim([max(0.,floor(xemin)-1) ceil(xemax)+0.25*(ceil(xemax)-max(0.,floor(xemin)-1))])
datetick('x','YY','keeplimits')
ridgepack_clearax('x',1)
ylabel('$\mathrm{Extent\,(\times\,10^{6}\,km^{2})}$')
set(h1,'box','on')
hold off

set(gcf,'CurrentAxes',h2)
hsa=plot(ncc.time.data,ncc.area.data/(10^6),'Color',satcol,...
        'LineStyle','--','LineWidth',lw);
plot(ncc.time.data,ncc.area84.data/(10^6),'Color',satcol,...
        'LineStyle','--','LineWidth',lw);
if cesmle
 hc=plot(ncesm.time.data,ncesm.area.data/(10^6),'Color',cesmcol,'LineStyle',':');
 plot(ncesm.time.data,ncesm.area84.data/(10^6),'Color',cesmcol,'LineStyle',':');
end
if mintime<datenum(1987,12,1)
 legend(hsa,{'NOAA CDR Area South of 87^{\circ} N (84^{\circ} N up to 1988) and model equivalent'},...
             'Location','NorthEast','FontSize',7)
else
 legend(hsa,{'NOAA CDR Area South of 87^{\circ} N and model equivalent'},...
             'Location','NorthEast','FontSize',7)
end
legend('boxoff')
xlim([mintime maxtime])
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
legend(h4,leg,'Location','NorthEast','FontSize',7,'Orientation','horizontal')
legend('boxoff')
xlim([mintime maxtime])
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
 ridgepack_multialign(gcf,['Regional Arctic System Model Sea Ice ',...
                  datestr(mintime,'YYYY'),'-',datestr(maxtime,'YYYY')],9)
end

ridgepack_fprint('png',[ridgepack_cellcat(runcell,'_'),'_seaice_volume_graph_',datestr(mintime,'YYYY'),'_',datestr(maxtime,'YYYY'),'_',num2str(minthick),'.png'],1,1)
ridgepack_fprint('epsc',[ridgepack_cellcat(runcell,'_'),'_seaice_volume_graph_',datestr(mintime,'YYYY'),'_',datestr(maxtime,'YYYY'),'_',num2str(minthick),'.eps'],1,1)



