
clear
close all

%generate=true;
generate=false;

%plotfullellipse=true;
plotfullellipse=false;

als='abcdefghijklmnopqrstuvwxyz';

%grid=false;
grid=true;

%plottimeseries=true;
plottimeseries=false;

%itqrange=false;
itqrange=true;

%label=true;
label=false;

%plotcross=true;
plotcross=false;

titletab='Sea Ice E3SM V2 PI Control';
filetabe='test';

cases=[1 2];
ccols=lines(length(cases));
casename={'v2.LR.piControl','v2.NARRM.piControl'};

yearst=1;
yearen=500;

vars={'totalIceExtent','totalIceVolume','totalSnowVolume'};
fact=[10^6 10^3 10^2];

lquantile=0.25;
uquantile=0.75;

maxcols=3;

% data acquisition
if generate 

 timevars=vars;
 timevars{end+1}='xtime';

 % create processed files for each case
 for i=cases

  k=1;

  cd(['/Users/afroberts/data/MODEL/E3SM/v2/',char(casename{i}),'/data/ice/hist']);

  for yearx=yearst:1:yearen

   files=dir([char(casename{i}),'.mpassi.hist.am.regionalStatistics.',...
              num2str(yearx,'%4.4i'),'.*.nc']);
 
   if length(files)<12
    error(['Missing files for year ',num2str(yearx,'%4.4i')])
   end

   for j=1:length(files)
    if k==1
     nc=ridgepack_clone(files(j).name,timevars);
     nc=rmfield(nc,'time');
     nc.time.data=datenum(nc.xtime.data,'yyyy-mm-dd_HH:MM:SS')';
     nc.time.dimension={'time'};
     for k=1:length(vars)
      nc.(char(vars{k})).data=nc.(char(vars{k})).data./fact(k);
      nc.(char(vars{k})).dimension={'nRegion','time'};
     end
     nc=rmfield(nc,{'xtime','StrLen','attributes','nRegions'})
     nc.attributes.title=[char(casename{i}),' lemnisc statistics years ',...
                         num2str(yearst,'%4.4i'),'-',num2str(yearen,'%4.4i')];
     clear nctime
    else
     ncadd=ridgepack_clone(files(j).name,timevars);
     nc.time.data=[nc.time.data datenum(ncadd.xtime.data,'yyyy-mm-dd_HH:MM:SS')'];
     for k=1:length(vars)
      nc.(char(vars{k})).data=[nc.(char(vars{k})).data ncadd.(char(vars{k})).data./fact(k)];
     end
     clear ncadd 
    end
    k=k+1;
   end

  end

  for k=1:365
   space=[k:365:365*floor(length(nc.time.data)./365)];
   annualtimeaverage(k)=nc.time.data(k);
   for j=1:3

    extentannualmean(j,k)=mean(nc.totalIceExtent.data(j,space));
    extentannualmedian(j,k)=median(nc.totalIceExtent.data(j,space));
    extentannualstd(j,k)=std(nc.totalIceExtent.data(j,space));
    extentannualequiv(j,k)=ridgepack_equiv(nc.totalIceExtent.data(j,space));
    extentlowerquantile(j,k)=quantile(nc.totalIceExtent.data(j,space),lquantile);
    extentupperquantile(j,k)=quantile(nc.totalIceExtent.data(j,space),uquantile);

    volumeannualmean(j,k)=mean(nc.totalIceVolume.data(j,space));
    volumeannualmedian(j,k)=median(nc.totalIceVolume.data(j,space));
    volumeannualstd(j,k)=std(nc.totalIceVolume.data(j,space));
    volumeannualequiv(j,k)=ridgepack_equiv(nc.totalIceVolume.data(j,space));
    volumelowerquantile(j,k)=quantile(nc.totalIceVolume.data(j,space),lquantile);
    volumeupperquantile(j,k)=quantile(nc.totalIceVolume.data(j,space),uquantile);

    snowannualmean(j,k)=mean(nc.totalSnowVolume.data(j,space));
    snowannualmedian(j,k)=median(nc.totalSnowVolume.data(j,space));
    snowannualstd(j,k)=std(nc.totalSnowVolume.data(j,space));
    snowannualequiv(j,k)=ridgepack_equiv(nc.totalSnowVolume.data(j,space));
    snowlowerquantile(j,k)=quantile(nc.totalSnowVolume.data(j,space),lquantile);
    snowupperquantile(j,k)=quantile(nc.totalSnowVolume.data(j,space),uquantile);

   end
  end

  nc.attributes.casename=char(casename{i});
  nc.attributes.regions='nRegions: 1=global, 2=Northern Hem, 3=Southern Hem';

  nc.dayofyear.data=[1:365];
  nc.dayofyear.long_name='day of year';
  nc.dayofyear.dimension={'dayofyear'};
  nc.dayofyear.units='day of year';
  nc.dayofyear.type='NC_INT';

  nc.nRegion.data=[1 2 3];
  nc.nRegion.long_name='Regions: 1=global, 2=Northern Hem, 3=Southern Hem';
  nc.nRegion.dimension={'nRegion'};
  nc.nRegion.type='NC_INT';

  nc.nsamp.data=length(space).*ones(365,1);
  nc.nsamp.long_name='Sample Size';
  nc.nsamp.dimension={'dayofyear'};
  nc.nsamp.type='NC_INT';

  nc.extentdaymean.data=extentannualmean;
  nc.extentdaymean.long_name=[char(casename{i}),'sea ice extent mean'];
  nc.extentdaymean.dimension={'nRegion','dayofyear'};
  nc.extentdaymean.units='km^2';

  nc.extentdaymedian.data=extentannualmedian;
  nc.extentdaymedian.long_name=[char(casename{i}),'sea ice extent median'];
  nc.extentdaymedian.dimension={'nRegion','dayofyear'};
  nc.extentdaymedian.units='km^2';

  nc.extentdaystd.data=extentannualstd;
  nc.extentdaystd.long_name=[char(casename{i}),'sea ice extent standard deviation'];
  nc.extentdaystd.dimension={'nRegion','dayofyear'};
  nc.extentdaystd.units='km^2';

  nc.extentdayequiv.data=extentannualequiv;
  nc.extentdayequiv.long_name=[char(casename{i}),'sea ice extent equivalent sample size'];
  nc.extentdayequiv.dimension={'nRegion','dayofyear'};
  nc.extentdayequiv.units='km^2';

  nc.extentlowerquantile.data=extentlowerquantile;
  nc.extentlowerquantile.long_name=[num2str(lquantile),' quantile of ice extent'];
  nc.extentlowerquantile.dimension={'nRegion','dayofyear'};
  nc.extentlowerquantile.units='km^2';

  nc.extentupperquantile.data=extentupperquantile;
  nc.extentupperquantile.long_name=[num2str(uquantile),' quantile of ice extent'];
  nc.extentupperquantile.dimension={'nRegion','dayofyear'};
  nc.extentupperquantile.units='km^2';

  % ---

  nc.volumedaymean.data=volumeannualmean;
  nc.volumedaymean.long_name=[char(casename{i}),'sea ice volume mean'];
  nc.volumedaymean.dimension={'nRegion','dayofyear'};
  nc.volumedaymean.units='km^3';

  nc.volumedaymedian.data=volumeannualmedian;
  nc.volumedaymedian.long_name=[char(casename{i}),'sea ice volume median'];
  nc.volumedaymedian.dimension={'nRegion','dayofyear'};
  nc.volumedaymedian.units='km^3';

  nc.volumedaystd.data=volumeannualstd;
  nc.volumedaystd.long_name=[char(casename{i}),'sea ice volume standard deviation'];
  nc.volumedaystd.dimension={'nRegion','dayofyear'};
  nc.volumedaystd.units='km^3';

  nc.volumedayequiv.data=volumeannualequiv;
  nc.volumedayequiv.long_name=[char(casename{i}),'sea ice volume equivalent sample size'];
  nc.volumedayequiv.dimension={'nRegion','dayofyear'};
  nc.volumedayequiv.units='km^3';

  nc.volumelowerquantile.data=volumelowerquantile;
  nc.volumelowerquantile.long_name=[num2str(lquantile),' quantile of ice volume'];
  nc.volumelowerquantile.dimension={'nRegion','dayofyear'};
  nc.volumelowerquantile.units='km^3';

  nc.volumeupperquantile.data=volumeupperquantile;
  nc.volumeupperquantile.long_name=[num2str(uquantile),' quantile of ice volume'];
  nc.volumeupperquantile.dimension={'nRegion','dayofyear'};
  nc.volumeupperquantile.units='km^3';

  % ---

  nc.snowdaymean.data=snowannualmean;
  nc.snowdaymean.long_name=[char(casename{i}),'sea ice snow volume mean'];
  nc.snowdaymean.dimension={'nRegion','dayofyear'};
  nc.snowdaymean.units='km^3';

  nc.snowdaymedian.data=snowannualmedian;
  nc.snowdaymedian.long_name=[char(casename{i}),'sea ice snow volume median'];
  nc.snowdaymedian.dimension={'nRegion','dayofyear'};
  nc.snowdaymedian.units='km^3';

  nc.snowdaystd.data=snowannualstd;
  nc.snowdaystd.long_name=[char(casename{i}),'sea ice snow volume standard deviation'];
  nc.snowdaystd.dimension={'nRegion','dayofyear'};
  nc.snowdaystd.units='km^3';

  nc.snowdayequiv.data=snowannualequiv;
  nc.snowdayequiv.long_name=[char(casename{i}),'sea ice snow volume equivalent sample size'];
  nc.snowdayequiv.dimension={'nRegion','dayofyear'};
  nc.snowdayequiv.units='km^3';

  nc.snowlowerquantile.data=snowlowerquantile;
  nc.snowlowerquantile.long_name=[num2str(lquantile),' quantile of ice snow'];
  nc.snowlowerquantile.dimension={'nRegion','dayofyear'};
  nc.snowlowerquantile.units='km^3';

  nc.snowupperquantile.data=snowupperquantile;
  nc.snowupperquantile.long_name=[num2str(uquantile),' quantile of ice snow'];
  nc.snowupperquantile.dimension={'nRegion','dayofyear'};
  nc.snowupperquantile.units='km^3';

  nc=ridgepack_struct(nc);

  cd(['/Users/afroberts/data/MODEL/E3SM/v2/',char(casename{i}),'/data/ice/processed']);

  outfile=[char(casename{i}),'.lemnisc.',...
           num2str(yearst,'%4.4i'),'-',num2str(yearen,'%4.4i'),'.',filetabe];

  nc=ridgepack_struct(nc);

  ridgepack_write(nc,outfile);

%  clear nc

 end

end

%return


for kcols=1:maxcols

 ridgepack_multiplot(1,maxcols,1,kcols,als(kcols))

 % set plot type
 if kcols==1
  titl2=['Sea Ice Extent \times{10^6} km^2'];
  xlab2=['Northern Hemisphere'];
  ylab2=['Southern Hemisphere'];
  %xlims=[0 20];
  %ylims=[0 20];
  xlims=[4 20];
  ylims=[0 20];
  xyticks=5;
  globalconts=30;
  globlab=['Global Extent'];
 elseif kcols==2
  titl2=['Sea Ice Volume \times{10^3} km^3'];
  xlab2=['Northern Hemisphere'];
  %xlims=[0 40];
  %ylims=[0 25];
  xlims=[5 40];
  ylims=[0 23];
  xyticks=5;
  globalconts=45;
  globlab=['Global Volume'];
 elseif kcols==3
  titl2=['Snow Volume \times{10^2} km^3'];
  xlab2=['Northern Hemisphere'];
  xlims=[0 35];
  ylims=[0 50];
  xyticks=10;
  globalconts=50;
  globlab=['Global Snow Volume'];
 end
  
 % create background axis
 axis square
 xlim(xlims)
 ylim(ylims)
 set(gca,'XTick',[0:xyticks:max([xlims ylims])],'YTick',[0:xyticks:max([xlims ylims])])
 ar=diff(xlims)./diff(ylims);

 % underlay global grid
 if grid

  graytype=0.7;
  x=0:max(xlims);
  y=0:max(ylims);
  [xx,yy]=meshgrid(x,y);
  z=xx+yy;
  conts=[0:xyticks:xlims(end)+ylims(end)];
  [C,h]=contour(xx,yy,z,conts,'Color',graytype*[1 1 1]);
  hold on

  for i=1:length(conts)

   theta=atan(1./ar);
   xpos=xlims(1)+(conts(i)-xlims(1))./(tan(theta)+1);
   ypos=xlims(1)+(conts(i)-xlims(1))-xpos;

   if xpos>xlims(1) & ypos>ylims(1) & xpos<xlims(end) & ypos<ylims(end)
    if conts(i)==globalconts
      text(xpos,ypos,globlab,...
         'Rotation',-180*atan(tan(pi*45./180)*ar)/pi,...
         'FontSize',8,'HorizontalAlignment','center',...
         'VerticalAlignment','bottom',...
         'Interpreter','Tex',...
         'Margin',1,'Color',graytype*[1 1 1]);
    else
     text(xpos,ypos,num2str(conts(i)),...
        'Rotation',-180*atan(tan(pi*45./180)*ar)/pi,...
        'FontSize',8,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom',...
        'Interpreter','Tex',...
        'Margin',1,'Color',graytype*[1 1 1]);
    end
   end

  end

 end

 % plot data from each case
 for i=cases

  cd(['/Users/afroberts/data/MODEL/E3SM/v2/',char(casename{i}),'/data/ice/processed']);
  infile=[char(casename{i}),'.lemnisc.',...
          num2str(yearst,'%4.4i'),'-',num2str(yearen,'%4.4i'),'.',filetabe];
  nc=ridgepack_clone(infile);

  if i==cases(1)
   if kcols==1
    mu1=nc.extentdaymean.data(1,:);
    std1=nc.extentdaystd.data(1,:);
    equiv1=nc.extentdayequiv.data(1,:);
   elseif kcols==2
    mu1=nc.volumedaymean.data(1,:); 
    std1=nc.volumedaystd.data(1,:);
    equiv1=nc.volumedayequiv.data(1,:);
   elseif kcols==3
    mu1=nc.snowdaymean.data(1,:);
    std1=nc.snowdaystd.data(1,:);
    equiv1=nc.snowdayequiv.data(1,:);
   end
  else
   if kcols==1
    mu2=nc.extentdaymean.data(1,:);
    std2=nc.extentdaystd.data(1,:);
    equiv2=nc.extentdayequiv.data(1,:);
   elseif kcols==2
    mu2=nc.volumedaymean.data(1,:); 
    std2=nc.volumedaystd.data(1,:);
    equiv2=nc.volumedayequiv.data(1,:);
   elseif kcols==3
    mu2=nc.snowdaymean.data(1,:);
    std2=nc.snowdaystd.data(1,:);
    equiv2=nc.snowdayequiv.data(1,:);
   end
   t = (mu1-mu2)./sqrt(((std1.^2)./equiv1) + ((std2.^2)./equiv2));
   df = (equiv1+equiv2-2);
   tcrit=tinv(0.995,df); % change 0.995 99%, 0.975 95%
   hcrit=ones(size(t));
   hcrit(isnan(t))=0;
   hcrit(-tcrit<t & t<tcrit)=0;
  end

  if kcols==1
   upperquantile=nc.extentupperquantile.data;
   lowerquantile=nc.extentlowerquantile.data;
   daymedian=nc.extentdaymedian.data;
   daymean=nc.extentdaymean.data;
   totalseries=nc.totalIceExtent.data;
  elseif kcols==2
   upperquantile=nc.volumeupperquantile.data;
   lowerquantile=nc.volumelowerquantile.data;
   daymedian=nc.volumedaymedian.data;
   daymean=nc.volumedaymean.data;
   totalseries=nc.totalIceVolume.data;
  elseif kcols==3
   upperquantile=nc.snowupperquantile.data;
   lowerquantile=nc.snowlowerquantile.data;
   daymedian=nc.snowdaymedian.data;
   daymean=nc.snowdaymean.data;
   totalseries=nc.totalSnowVolume.data;
  end

  if itqrange

   ocol=0.90*[1 1 1];
   alpha=0.1;

   for k=1:1:365

    theta=[0:10:90]*pi/180;

    ae=abs(upperquantile(2,k)-daymedian(2,k));
    be=abs(upperquantile(3,k)-daymedian(3,k));

    xe=ae*cos(theta);
    ye=be*sin(theta);

    xs=[daymedian(2,k) ...
        xe+daymedian(2,k) ...
        daymedian(2,k)];
    ys=[daymedian(3,k) ...
        ye+daymedian(3,k) ...
        daymedian(3,k)];
    patch(xs,ys,ocol,'EdgeColor','none','FaceAlpha',alpha)

    theta=[90:10:180]*pi/180;

    ae=abs(lowerquantile(2,k)-daymedian(2,k));
    be=abs(upperquantile(3,k)-daymedian(3,k));

    xe=ae*cos(theta);
    ye=be*sin(theta);

    ys=[daymedian(3,k) ...
        ye+daymedian(3,k) ...
        daymedian(3,k)];
    xs=[daymedian(2,k) ...
        xe+daymedian(2,k) ...
        daymedian(2,k)];
    patch(xs,ys,ocol,'EdgeColor','none','FaceAlpha',alpha)

    theta=[180:10:270]*pi/180;

    ae=abs(lowerquantile(2,k)-daymedian(2,k));
    be=abs(lowerquantile(3,k)-daymedian(3,k));

    xe=ae*cos(theta);
    ye=be*sin(theta);

    xs=[daymedian(2,k) ...
        xe+daymedian(2,k) ...
        daymedian(2,k)];
    ys=[daymedian(3,k) ...
        ye+daymedian(3,k) ...
        daymedian(3,k)];
    patch(xs,ys,ocol,'EdgeColor','none','FaceAlpha',alpha)

    theta=[270:10:360]*pi/180;

    ae=abs(upperquantile(2,k)-daymedian(2,k));
    be=abs(lowerquantile(3,k)-daymedian(3,k));

    xe=ae*cos(theta);
    ye=be*sin(theta);

    xs=[daymedian(2,k) ...
        xe+daymedian(2,k) ...
        daymedian(2,k)];
    ys=[daymedian(3,k) ...
        ye+daymedian(3,k) ...
        daymedian(3,k)];
    patch(xs,ys,ocol,'EdgeColor','none','FaceAlpha',alpha)

   end

  end

  if plottimeseries

   plot(totalseries(2,:),totalseries(3,:),'Color',0.9*[1 1 1]);

   if plotcross

    ocol=[0.8500 0.3250 0.0980]

    k=datenum(0000,5,21);

    space=[k:365:365*floor(length(nc.time.data)./365)];

    plot(totalseries(2,space),totalseries(3,space),'.',...
        'Color',0.5*[1 1 1]);

    plot([lowerquantile(2,k) upperquantile(2,k)],...
         [daymedian(3,k) daymedian(3,k)],...
        'Color',ocol,'LineWidth',1);

    plot([daymedian(2,k) daymedian(2,k)],...
         [lowerquantile(3,k) upperquantile(3,k)],...
         'Color',ocol,'LineWidth',1);

   end

  end

  if plotfullellipse

    theta=[0:10:90]*pi/180;

    ae=abs(upperquantile(2,k)-daymedian(2,k));
    be=abs(upperquantile(3,k)-daymedian(3,k));

    xe=ae*cos(theta);
    ye=be*sin(theta);
 
    plot(xe+daymedian(2,k),ye+daymedian(3,k),...
          'Color',ocol,'LineWidth',1)

    theta=[90:10:180]*pi/180;
 
    ae=abs(lowerquantile(2,k)-daymedian(2,k));
    be=abs(upperquantile(3,k)-daymedian(3,k));

    xe=ae*cos(theta);
    ye=be*sin(theta);

    plot(xe+daymedian(2,k),ye+daymedian(3,k),...
         'Color',ocol,'LineWidth',1)

    theta=[180:10:270]*pi/180;

    ae=abs(lowerquantile(2,k)-daymedian(2,k));
    be=abs(lowerquantile(3,k)-daymedian(3,k));

    xe=ae*cos(theta);
    ye=be*sin(theta);

    plot(xe+daymedian(2,k),ye+daymedian(3,k),...
          'Color',ocol,'LineWidth',1)

    theta=[270:10:360]*pi/180;

    ae=abs(upperquantile(2,k)-daymedian(2,k));
    be=abs(lowerquantile(3,k)-daymedian(3,k));

    xe=ae*cos(theta);
    ye=be*sin(theta);
 
    plot(xe+daymedian(2,k),ye+daymedian(3,k),...
         'Color',ocol,'LineWidth',1)

  end

 end

 for i=cases

  cd(['/Users/afroberts/data/MODEL/E3SM/v2/',char(casename{i}),'/data/ice/processed']);
  infile=[char(casename{i}),'.lemnisc.',...
          num2str(yearst,'%4.4i'),'-',num2str(yearen,'%4.4i'),'.',filetabe];
  nc=ridgepack_clone(infile);

  % grab daily mean  
  if kcols==1
   daymean=nc.extentdaymean.data;
  elseif kcols==2
   daymean=nc.volumedaymean.data;
  elseif kcols==3
   daymean=nc.snowdaymean.data;
  end

  % grab t-test information
  if i==cases(1)
   if kcols==1
    mu1=nc.extentdaymean.data(1,:);
    std1=nc.extentdaystd.data(1,:);
    equiv1=nc.extentdayequiv.data(1,:);
   elseif kcols==2
    mu1=nc.volumedaymean.data(1,:);
    std1=nc.volumedaystd.data(1,:);
    equiv1=nc.volumedayequiv.data(1,:);
   elseif kcols==3
    mu1=nc.snowdaymean.data(1,:);
    std1=nc.snowdaystd.data(1,:);
    equiv1=nc.snowdayequiv.data(1,:);
   end
  else
   if kcols==1
    mu2=nc.extentdaymean.data(1,:);
    std2=nc.extentdaystd.data(1,:);
    equiv2=nc.extentdayequiv.data(1,:);
   elseif kcols==2
    mu2=nc.volumedaymean.data(1,:);
    std2=nc.volumedaystd.data(1,:);
    equiv2=nc.volumedayequiv.data(1,:);
   elseif kcols==3
    mu2=nc.snowdaymean.data(1,:);
    std2=nc.snowdaystd.data(1,:);
    equiv2=nc.snowdayequiv.data(1,:);
   end
   t = (mu1-mu2)./sqrt(((std1.^2)./equiv1) + ((std2.^2)./equiv2));
   df = (equiv1+equiv2-2);
   tcrit=tinv(0.995,df); % change 0.995 99%, 0.975 95%
   hcrit=ones(size(t));
   hcrit(isnan(t))=0;
   hcrit(-tcrit<t & t<tcrit)=0;
  end

  % plot mean that is statistically significant
  if i==cases(1)
   plot(daymean(2,:),daymean(3,:),'Color',ccols(i,:));
  else
   plotmean=daymean;
   plotmean(:,hcrit==0)=NaN;
   plot(plotmean(2,:),plotmean(3,:),'-','Color',ccols(i,:));
   plotmean=daymean;
   plotmean(:,hcrit==1)=NaN;
   plot(plotmean(2,:),plotmean(3,:),'--','Color',ccols(i,:));
  end

  if label

    theta=atan(ar.*diff(daymean(3,[2 3]))./diff(daymean(2,[2 3])));

    text(daymean(2,1),daymean(3,1),'January',...
        'Rotation',theta*180/pi,...
        'FontSize',6,'HorizontalAlignment','left',...
        'VerticalAlignment','bottom',...
        'Margin',1,'Color','b','Interpreter','Tex');

  end

  idx=datenum(0000,3,21);
  ha=plot(daymean(2,idx),daymean(3,idx),'.','Color',ccols(i,:));

  if label

    theta=atan(ar*diff(daymean(3,[idx-2 idx+2]))./...
                  diff(daymean(2,[idx-2 idx+2])));

    text(daymean(2,idx),daymean(3,idx),{'Equinox'},...
        'Rotation',theta*180/pi,...
        'FontSize',6,'HorizontalAlignment','center',...
        'VerticalAlignment','top',...
        'Margin',1,'Color','b','Interpreter','Tex');

  end

  idx=datenum(0000,6,21);
  ha=plot(daymean(2,idx),daymean(3,idx),'o',...
           'MarkerFaceColor','w','MarkerSize',2,'Color',ccols(i,:));
 
  if label

    theta=atan(ar*diff(daymean(3,[idx-2 idx+2]))./...
                  diff(daymean(2,[idx-2 idx+2])));

    text(daymean(2,idx),daymean(3,idx),{'Solstice'},...
        'Rotation',theta*180/pi,...
        'FontSize',5,'HorizontalAlignment','center',...
        'VerticalAlignment','top',...
        'Margin',3,'Color','b','Interpreter','Tex');

  end
 
  idx=datenum(0000,9,21);
  plot(daymean(2,idx),daymean(3,idx),'.','Color',ccols(i,:));

  if label

    theta=atan(ar*diff(daymean(3,[idx-2 idx+2]))./...
                  diff(daymean(2,[idx-2 idx+2])));

    text(daymean(2,idx),daymean(3,idx),{'Equinox'},...
        'Rotation',theta*180/pi,...
        'FontSize',6,'HorizontalAlignment','center',...
        'VerticalAlignment','top',...
        'Margin',1,'Color','b','Interpreter','Tex');

  end
 
  idx=datenum(0000,12,21);
  ha=plot(daymean(2,idx),daymean(3,idx),'o',...
           'MarkerFaceColor','w','MarkerSize',2,'Color',ccols(i,:));

  if label

    theta=atan(ar*diff(daymean(3,[idx-2 idx+2]))./...
                  diff(daymean(2,[idx-2 idx+2])));

    text(daymean(2,idx),daymean(3,idx),{'Solstice'},...
        'Rotation',theta*180/pi,...
        'FontSize',5,'HorizontalAlignment','center',...
        'VerticalAlignment','top',...
        'Margin',3,'Color','b','Interpreter','Tex');

  end
 
  % plot arrow
  theta=atan(ar*diff(daymean(3,[end end-2]))./...
                diff(daymean(2,[end end-2])));
  z=sqrt(diff(daymean(3,[end end-2])).^2+...
         diff(daymean(2,[end end-2])).^2);

  xarrow=[daymean(2,end)+z.*cos(theta-150*pi/180) ...
          daymean(2,end) ...
          daymean(2,end)+z.*cos(theta+150*pi/180)];
  yarrow=[daymean(3,end)+z.*sin(theta-150*pi/180) ...
          daymean(3,end) ...
          daymean(3,end)+z.*sin(theta+150*pi/180)];

  plot(xarrow,yarrow,'Color',ccols(i,:))

 end

 xlabel(xlab2,'Interpreter','Tex','Fontsize',10)
 if kcols==1
  ylabel(ylab2,'Interpreter','Tex','Fontsize',10)
 end
 title(titl2,'Interpreter','Tex','Fontsize',10,'FontWeight','normal')

 %  legend([hp],{'V2 PI Control Reference'});
 %  legend('boxoff');

end

ridgepack_multialign(gcf,'',12)

cd('/Users/afroberts/work')

ridgepack_fprint('png',['E3SM_sea_ice_lemnisc.',...
                 num2str(yearst,'%4.4i'),'-',num2str(yearen,'%4.4i'),'.',filetabe,'.png'],1,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%

function equivn=ridgepack_equiv(seriesdata)

% Calculate equivalent sample size based on Wilks (2006)
% "Statistical Methods in the Atmospheric Sciences" with further 
% information from the National Institute of Standards. Note
% that when all weights are equal but not 1, the lag-1 autocorrelation
% coefficient is identical to the off-diagonal terms of
% r=corrcoef(nc.(name).data(j,1:end-1),nc.(name).data(j,2:end)); 
% which is a verification of the code is correctly approximating
% unbiased weighted variance and covariance for lag-1.

% equally weight all samples, but set up for potential change later
weight=ones(size(seriesdata));
coeff=sum(abs(weight(:))>0);

% calculate lag-1 weight
w1=weight(1:end-1);
w2=weight(2:end);

% calculate lag-1 timeseries (offset by 1)
x1=seriesdata(1:end-1);
x2=seriesdata(2:end);

% get lag-1 length
N=length(x1);

% calculate reduced means
mean1=nansum(w1.*x1)./nansum(w1);
mean2=nansum(w2.*x2)./nansum(w2);

% calculated weighted covariance
covar1=w1.*(x1-mean1);
covar2=w2.*(x2-mean2);
covariance=nansum(covar1.*covar2)./((N-1)*nansum(w1.*w2)./N);

% calculate variance 
var1=nansum(w1.*(x1-mean1).^2)./((N-1)*nansum(w1)./N);
var2=nansum(w2.*(x2-mean2).^2)./((N-1)*nansum(w2)./N);

% lag-1 autocorrelation coefficient
r1=covariance/sqrt(var1.*var2);

equivn=max(0,coeff.*((1-r1)./(1+r1)));

end

