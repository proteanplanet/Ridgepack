
clear
close all

startyear=1980;
endyear=2014;
%endyear=1980;
%monthsets={[1]};
%monthsets={[1 2 3]};
%monthsets={[1 2 3],[4 5 6]};
%monthsets={[1 2 3],[4 5 6],[7 8 9],[10 11 12]};
monthsets={[1 2 3],[7 8 9]};

%hemisphere='nh';
hemisphere='sh';

difference=true;
%difference=false;

%stddev=true;
stddev=false;

%equiv=true;
equiv=false;

density=7;

filetag='NARRM.LR.Talk.';

%ensembles={[1],[2 3 4 5 6],[7 8 9 10 11]};
ensembles={[1],[7 8 9 10 11]};
%ensembles={[1],[2 3 4 5 6]};

filenames={'Pathfinder',...
           'v2.NARRM.historical_0101.mpassi.daily',...
           'v2.NARRM.historical_0151.mpassi.daily',...
           'v2.NARRM.historical_0201.mpassi.daily',...
           'v2.NARRM.historical_0251.mpassi.daily',...
           'v2.NARRM.historical_0301.mpassi.daily',...
           'v2.LR.historical_0101.mpassi.daily',...
           'v2.LR.historical_0151.mpassi.daily',...
           'v2.LR.historical_0201.mpassi.daily',...
           'v2.LR.historical_0251.mpassi.daily',...
           'v2.LR.historical_0301.mpassi.daily'};

locations={'/Users/afroberts/data/data/SATELLITE/processed/pathfinder_v4',...
           '/Users/afroberts/data/MODEL/E3SM/v2/v2.NARRM.historical_0101/data/ice/processed',...
           '/Users/afroberts/data/MODEL/E3SM/v2/v2.NARRM.historical_0151/data/ice/processed',...
           '/Users/afroberts/data/MODEL/E3SM/v2/v2.NARRM.historical_0201/data/ice/processed',...
           '/Users/afroberts/data/MODEL/E3SM/v2/v2.NARRM.historical_0251/data/ice/processed',...
           '/Users/afroberts/data/MODEL/E3SM/v2/v2.NARRM.historical_0301/data/ice/processed',...
           '/Users/afroberts/data/MODEL/E3SM/v2/v2.LR.historical_0101/data/ice/processed',...
           '/Users/afroberts/data/MODEL/E3SM/v2/v2.LR.historical_0151/data/ice/processed',...
           '/Users/afroberts/data/MODEL/E3SM/v2/v2.LR.historical_0201/data/ice/processed',...
           '/Users/afroberts/data/MODEL/E3SM/v2/v2.LR.historical_0251/data/ice/processed',...
           '/Users/afroberts/data/MODEL/E3SM/v2/v2.LR.historical_0301/data/ice/processed'};

cdrlocation='/Users/afroberts/data/data/SATELLITE/processed/G02202_v4'

rnametags={'Pathfinder','NARRM Ensemble','LR Ensemble'};

mons={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

if strcmp(hemisphere,'nh')
 hem=1;
elseif strcmp(hemisphere,'sh')
 hem=-1;
end

alphatag=['abcdefghijklmnopqrstuvwxyz'];

% proceed in this order so that differences are taken from observations in each column
for cols=1:length(monthsets)
 for rows=1:length(ensembles)

  monthnum=['months'];
  for months=monthsets{cols};
   monthnum=[monthnum,'_',num2str(months,'%2.2i')];
  end
  yeartag=['years_',num2str(startyear,'%4.4i'),'_',num2str(endyear,'%4.4i')];

  % get mean data
  for ei=1:length(ensembles{rows})

    % get sensemle file
    if strcmp(char(filenames{ensembles{rows}(ei)}),'Pathfinder')
     infileregrid=['Pathfinder.',hemisphere,'.regrid.daily.icemotion.',monthnum,'.',yeartag,'.mean'];
     modeldata=false;
    else
     infileregrid=[char(filenames{ensembles{rows}(ei)}),'.',hemisphere,...
                   '.regrid.',monthnum,'.',yeartag,'.mean'];
     modeldata=true;
    end

    cd(char(locations{ensembles{rows}(ei)}))

    % get mean data from the entire ensemble
    if ei==1

     nc=ridgepack_clone(infileregrid);

     nc.umean=nc.u;
     nc.vmean=nc.v;
     nc.smean=nc.speed;
     if modeldata
      nc.cmean=nc.conc;
      nc.concensemble=nc.conc;
      nc.concensemble.dimension{end+1}='ensembles';
      nc.ensembles.data=[1:length(ensembles{rows})];
      nc.ensembles.long_name='ensemble index';
     end

    else

     ncnext=ridgepack_clone(infileregrid);

     nc.umean.data=nc.umean.data+ncnext.u.data;
     nc.vmean.data=nc.vmean.data+ncnext.v.data;
     nc.smean.data=nc.smean.data+ncnext.speed.data;
     nc.cmean.data=nc.cmean.data+ncnext.conc.data;
     nc.concensemble.data(:,:,ei)=ncnext.conc.data;

    end

    clear ncnext

  end

  nc.umean.data=nc.umean.data./length(ensembles{rows});
  nc.vmean.data=nc.vmean.data./length(ensembles{rows});
  nc.smean.data=nc.smean.data./length(ensembles{rows});
  if modeldata
    nc.cmean.data=nc.cmean.data./length(ensembles{rows});
  end

  % set multiplot
  ridgepack_multiplot(length(ensembles),length(monthsets),rows,cols,...
                       alphatag((rows-1)*length(monthsets)+cols));

  % generate underlying map
  if hem==1
   if rows==1 & cols==1 
    ridgepack_polarm('seaice','grid','label','noland')
   else
    ridgepack_polarm('seaice','grid','noland')
   end
  elseif hem==-1
   if rows==1 & cols==1 
    ridgepack_polarm('antarctic','grid','label','noland')
   else
    ridgepack_polarm('antarctic','grid','noland')
   end
  else
   error('hem needs to be either 1 (north) or -1 (south)')
  end

  % remove speed data for less than 15% concentration
  if strcmp(char(filenames{ensembles{rows}(1)}),'Pathfinder')
   nco=nc;
   if stddev
    nco.speed_std.data(isnan(nco.speed.data))=NaN;
   else
    nco.speed.data(nco.mask.data==0)=NaN;
   end
  else
   nc.umean.data(nc.cmean.data<0.15)=NaN;
   nc.vmean.data(nc.cmean.data<0.15)=NaN;
   if difference
    nc.smean.data(nc.cmean.data<0.15)=NaN;
    if ~isstruct(nco)
     error('No observed speed loaded')
    elseif length(ensembles{rows})==1
     ncdiff=ridgepack_ttest2(nc,'smean',nco,'speed',1);

     cs=1; % cs sets the size of the checks (8=eight grid cells a side)
     for ai=1:cs:size(ncdiff.diff.data,1)-cs;
     for bi=(mod(ai,2*cs)+1):2*cs:size(ncdiff.diff.data,2)-cs;
      if mean(mean(ncdiff.h.data(ai:ai+cs-1,bi:bi+cs-1)))<0.5 & ...
         mean(mean(nc.mask.data(ai:ai+cs-1,bi:bi+cs-1)))>0.5
         ncdiff.diff.data(ai:ai+cs-1,bi:bi+cs-1)=NaN;
      end
     end
     end
    else
     ncdiff=nc;
     ncdiff.diff=nc.smean;
     ncdiff.diff.long_name='Speed Difference';
     ncdiff.diff.data=nc.smean.data-nco.smean.data;
    end
   elseif stddev
    nc.speed_std.data(nc.conc.data<0.15 | nc.mask.data==0)=NaN;
   else
    nc.smean.data(nc.cmean.data<0.15 | nc.mask.data==0)=NaN;
   end
  end

  % plot speed data
  if difference
   obsconts=[0:0.02:0.16];
   obsfield='smean';
   spconts=[-0.06:0.02:-0.02 0.02:0.02:0.16];
   colfield='diff';
  elseif stddev
   spconts=[0:0.02:0.12];
   colfield='speed_std';
  elseif equiv
   spconts=[0:100:2000];
   colfield='speed_equiv';
  else
   spconts=[0:0.02:0.22];
   colfield='smean';
  end

  if rows==1 & cols==length(monthsets)
   if difference
    ridgepack_pcolorm(nco,obsfield,{},{},obsconts,'linear',0,'vertical','parula')
   else
    ridgepack_pcolorm(nc,colfield,{},{},spconts,'linear',0,'vertical')
   end
  else
   if difference & rows==1 
    ridgepack_pcolorm(nc,obsfield,{},{},obsconts,'linear',0,'none','parula')
   elseif difference & rows==2 & cols==length(monthsets)
    ridgepack_pcolorm(ncdiff,colfield,{},{},spconts,'linear',0,'vertical')
    if length(ensembles)>2
     ridgepack_cbshare(gca)
    end
   elseif difference
    ridgepack_pcolorm(ncdiff,colfield,{},{},spconts,'linear',0,'none')
   else
    ridgepack_pcolorm(nc,colfield,{},{},spconts,'linear',0,'none')
   end
  end

  % plot mask
  ridgepack_maskm(nc.latitude.data,nc.longitude.data,nc.mask.data)  

  % turn vectors from lat-lon to x-y grid if they are in geo coordinates
  if strcmp(char(filenames{ensembles{rows}(1)}),'Pathfinder')
   ui=nc.umean.data;
   vi=nc.vmean.data;
  else
   [th,z]=cart2pol(nc.umean.data,nc.vmean.data);
   if strcmp(filenames{ensembles{rows}(1)},'v2.NARRM.historical_0101.mpassi.daily')
    [ui,vi]=pol2cart(th+deg2rad(nc.turn.data),z);
   else
    [ui,vi]=pol2cart(th-deg2rad(nc.turn.data),z);
   end
  end


  % generate streamlines on x-y grid
  [vertices arrowvertices]=streamslice(nc.x.data,nc.y.data,ui,vi,density);

  % plot streamlines on the map
  latitude=[];
  longitude=[];
  for i=1:length(vertices)
   if ~isempty(vertices{i})
    xy=vertices{i};
    [lat,lon]=ridgepack_xytogeodetic(xy(:,1),xy(:,2),hem);
    if hem==1; lon=lon-45; end
    [c,d]=mfwdtran(gcm,lat,lon,gca,'surface');
    plot(c,d,'Color',0*[1 1 1],'LineWidth',0.5)
    latitude=[latitude NaN lat'];
    longitude=[longitude NaN lon'];
   end
  end

  % plot streamline arrows on the map
  latitudev=[];
  longitudev=[];
  for i=1:length(arrowvertices)
   if ~isempty(arrowvertices{i})
    xy=arrowvertices{i};
    [lat,lon]=ridgepack_xytogeodetic(xy(:,1),xy(:,2),hem);
    if hem==1; lon=lon-45; end
    [c,d]=mfwdtran(gcm,lat,lon,gca,'surface');
    plot(c,d,'Color',0*[1 1 1],'LineWidth',0.5)
    latitudev=[latitudev NaN lat'];
    longitudev=[longitudev NaN lon'];
   end
  end

  % add extent to the plot
  if ~strcmp(char(filenames{ensembles{rows}(1)}),'Pathfinder')

   for ei=1:length(ensembles{rows})
    concmask=zeros(size(nc.conc.data));
    concmask(squeeze(nc.concensemble.data(:,:,ei))>0.15)=1;
    concmask(nc.mask.data==0)=NaN;
    %ridgepack_maskm(nc.latitude.data,nc.longitude.data,concmask,[0 0 0.75],0.25)
    ridgepack_maskm(nc.latitude.data,nc.longitude.data,concmask,[0 0 0.75],0.25,':')
   end

   monthnum=['months'];
   for months=monthsets{cols};
    monthnum=[monthnum,'_',num2str(months,'%2.2i')];
   end
   yeartag=[num2str(startyear,'%4.4i'),'_',num2str(endyear,'%4.4i')];

   cd(cdrlocation)

   if strcmp(hemisphere,'nh')
    infileconc=['G02202_v4_merged_north_r00_1979_2022.',monthnum,'_',yeartag,'.conc.mean'];
   else
    infileconc=['G02202_v4_merged_south_r00_1979_2022.',monthnum,'_',yeartag,'.conc.mean'];
   end

   nccdr=ridgepack_clone(infileconc);
   concmask=zeros(size(nccdr.conc.data));
   concmask(squeeze(nccdr.conc.data(:,:))>15)=1;
   concmask(nccdr.mask.data==0)=NaN;

   ridgepack_maskm(nccdr.latitude.data,nccdr.longitude.data,concmask,'m',0.5)

  end

  % add mean speed to plot
  if strcmp(char(filenames{ensembles{rows}(1)}),'Pathfinder') 
   meanspeed=mean(nc.smean.data(~isnan(nc.smean.data(:))));
   infomean=['$|\bar{u}|$=',num2str(meanspeed,'%3.3f')];
  else
   meanspeed=mean(ncdiff.diff.data(~isnan(ncdiff.diff.data(:))));
   if meanspeed>0
    infomean=['$\Delta|\bar{u}|$=+',num2str(meanspeed,'%3.3f')];
   elseif meanspeed<0
    infomean=['$\Delta|\bar{u}|$=-',num2str(meanspeed,'%3.3f')];
   else
    infomean=['$\Delta|\bar{u}|$=0'];
   end
  end
  if strcmp(hemisphere,'nh')
   if strcmp(char(filenames{ensembles{rows}(1)}),'Pathfinder')
    latmean=40;
    lonmean=-10;
   else
    latmean=52;
    lonmean=85;
   end
  elseif strcmp(hemisphere,'sh')
   latmean=-48;
   lonmean=130;
  end
  if rows==1 & cols==1
  else
  end
  textm(latmean,lonmean,infomean,'Color',[0 0 0.75],...
        'Interpreter','latex','Fontsize',7,'HorizontalAlignment','right')

  % fill in titles
  if rows==1
   if length(monthsets{cols})==1
    title(mons(monthsets{cols}))
   else
    idx=monthsets{cols};
    title([char(mons(idx(1))),' - ',char(mons(idx(end)))])
   end
  else
   title('')
  end
 
  if cols==1
   ylabel(char(rnametags{rows}))
  end

 end
end

ridgepack_multialign(gcf)

cd('/Users/afroberts/work')

if difference
 outfileregrid=[filetag,hemisphere,'.stream.',yeartag,'.diff.png'];
elseif stddev
 outfileregrid=[filetag,hemisphere,'.stream.',yeartag,'.std.png'];
elseif equiv
 outfileregrid=[filetag,hemisphere,'.stream.',yeartag,'.equiv.png'];
else
 outfileregrid=[filetag,hemisphere,'.stream.',yeartag,'.mean.png'];
end

ridgepack_fprint('png',outfileregrid,1,3);




