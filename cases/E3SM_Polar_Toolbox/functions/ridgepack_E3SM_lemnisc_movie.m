
clear
close all

%generate=true;
generate=false;

%generateobs=true;
generateobs=false;

% plotmean=true;
plotmean=false;

%plotfullellipse=true;
plotfullellipse=false;

als='abcdefghijklmnopqrstuvwxyz';

grid=false;
%grid=true;

plottimeseries=true;
%plottimeseries=false;

%itqrange=true;
itqrange=false;

%label=true;
label=false;

%yearlabel=true;
yearlabel=false;

%plotcross=true;
plotcross=false;

%plotequinoxtrend=true;
plotequinoxtrend=false;

%observations=true;
observations=false;

titletab='Sea Ice E3SM Preindustrial';

filetabe='control';
%filetabe='industrial3';
%filetabe='PI';

%ensemblenames={'PSLV','Icedge'};
%ensemblenames={'LR','NARRM'};
%ensemblenames={'NARRM'};
ensemblenames={'LR'};

%legnames={'LR 5-member','NARRM 5-member'};
%legnames={'LR PI Control','NARRM PI Control'};
%legnames={'NARRM PI Control'};
legnames={'LR PI Control'}
%legnames={'PSLV','Icedge'};

%ensemblecases={[1 2 3 4 5],[6 7 8 9 10]};
%ensemblecases={[1 2 3 4 5]};
%ensemblecases={[11],[12]};
%ensemblecases={[1],[2]};
ensemblecases={[1]};

%yearrange={[1980 1999],[2000 2014]};
%yearrange={[1980 2014]};
yearrange={[0001 0500]};
%yearrange={[51 100]};

yearsto=1980;
yeareno=2014;

maxcols=2;

% case names
%casenames={'20231014.v3alpha04_trigrid_pslv.piControl.chrysalis',...
%           '20231025.v3alpha04_trigrid_vslim.piControl.chrysalis'};

% case directories where the regional sea ice data is held
%dirnames={'/Users/afroberts/data/MODEL/E3SM/pslv',...
%          '/Users/afroberts/data/MODEL/E3SM/Icedge'};

% directories where the processed lemnisc data is or will be written
%eprnames={'/Users/afroberts/data/MODEL/E3SM/pslv',...
%          '/Users/afroberts/data/MODEL/E3SM/Icedge'};

casenames={'v2.NARRM.piControl'};

%casenames={'v2.LR.piControl',...
%           'v2.NARRM.piControl'};

%casenames={'v2.LR.historical_0101',...
%           'v2.LR.historical_0151',...
%           'v2.LR.historical_0201',...
%           'v2.LR.historical_0251',...
%           'v2.LR.historical_0301',...
%           'v2.NARRM.historical_0101',...
%           'v2.NARRM.historical_0151',...
%           'v2.NARRM.historical_0201',...
%           'v2.NARRM.historical_0251',...
%           'v2.NARRM.historical_0301',...
%           'v2.LR.piControl',...
%           'v2.NARRM.piControl'};

% find ensemble index of case
caseensembleindex=zeros([1 length(casenames)]);
for j=1:length(ensemblecases)
 for i=ensemblecases{j}
  caseensembleindex(i)=j;
 end
end

% get case and ensemble locations
for i=1:length(casenames)
 dirnames{i}=['/Users/afroberts/data/MODEL/E3SM/v2/',char(casenames{i}),'/data/ice/hist'];
 eprnames{i}=['/Users/afroberts/data/MODEL/E3SM/v2/v2.',...
              char(ensemblenames{caseensembleindex(i)}),'/processed'];
end

if observations
 nlemniscs=length(ensemblenames)+1;
 legnames{nlemniscs}='NOAA CDR'; % observational
else
 nlemniscs=length(ensemblenames);
end

ccols=lines(nlemniscs);

vars={'totalIceExtent','totalIceVolume','totalSnowVolume','averageAlbedo','totalKineticEnergy'};
fact=[10^6 10^3 10^2 1 10^12];

lquantile=0.25;
uquantile=0.75;

if length(yearrange)>2
 error('Year range must only include two windows')
end

% data acquisition
if generate 

 timevars=vars;
 timevars{end+1}='xtime';

 % years cases
 for yi=1:length(yearrange)

 yearst=yearrange{yi}(1);
 yearen=yearrange{yi}(2);

 % assemble ensemble statistics
 for l=1:length(ensemblenames)

  if length(ensemblecases)>length(casenames)
   error('Requested cases greater than number of casenames')
  end

  k=1;

  % create processed files for each case
  for i=ensemblecases{l}

   cd(char(dirnames{i}));
 
   for yearx=yearst:1:yearen

    files=dir([char(casenames{i}),'.mpassi.hist.am.regionalStatistics.',...
               num2str(yearx,'%4.4i'),'.*.nc']);
 
    if length(files)<12
     error(['Missing files for year ',num2str(yearx,'%4.4i'),' case ',char(casenames{i})])
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
      nc.attributes.title=[char(ensemblenames{l}),' lemnisc ensemble statistics years ',...
                           num2str(yearst,'%4.4i'),'-',num2str(yearen,'%4.4i')];
      nc.attributes.cases=[''];
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

   nc.attributes.cases=[nc.attributes.cases,char(casenames{i}),' '];

  end % cases in ensemble loop

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

    kineticannualmean(j,k)=mean(nc.totalKineticEnergy.data(j,space));
    kineticannualmedian(j,k)=median(nc.totalKineticEnergy.data(j,space));
    kineticannualstd(j,k)=std(nc.totalKineticEnergy.data(j,space));
    kineticannualequiv(j,k)=ridgepack_equiv(nc.totalKineticEnergy.data(j,space));
    kineticlowerquantile(j,k)=quantile(nc.totalKineticEnergy.data(j,space),lquantile);
    kineticupperquantile(j,k)=quantile(nc.totalKineticEnergy.data(j,space),uquantile);

    albedoannualmean(j,k)=mean(nc.averageAlbedo.data(j,space));
    albedoannualmedian(j,k)=median(nc.averageAlbedo.data(j,space));
    albedoannualstd(j,k)=std(nc.averageAlbedo.data(j,space));
    albedoannualequiv(j,k)=ridgepack_equiv(nc.averageAlbedo.data(j,space));
    albedolowerquantile(j,k)=quantile(nc.averageAlbedo.data(j,space),lquantile);
    albedoupperquantile(j,k)=quantile(nc.averageAlbedo.data(j,space),uquantile);

   end
  end

  nc.attributes.ensemblename=char(ensemblenames{l});
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
  nc.extentdaymean.long_name=[char(ensemblenames{l}),'sea ice extent mean'];
  nc.extentdaymean.dimension={'nRegion','dayofyear'};
  nc.extentdaymean.units='km^2';

  nc.extentdaymedian.data=extentannualmedian;
  nc.extentdaymedian.long_name=[char(ensemblenames{l}),'sea ice extent median'];
  nc.extentdaymedian.dimension={'nRegion','dayofyear'};
  nc.extentdaymedian.units='km^2';

  nc.extentdaystd.data=extentannualstd;
  nc.extentdaystd.long_name=[char(ensemblenames{l}),'sea ice extent standard deviation'];
  nc.extentdaystd.dimension={'nRegion','dayofyear'};
  nc.extentdaystd.units='km^2';

  nc.extentdayequiv.data=extentannualequiv;
  nc.extentdayequiv.long_name=[char(ensemblenames{l}),'sea ice extent equivalent sample size'];
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
  nc.volumedaymean.long_name=[char(ensemblenames{l}),'sea ice volume mean'];
  nc.volumedaymean.dimension={'nRegion','dayofyear'};
  nc.volumedaymean.units='km^3';

  nc.volumedaymedian.data=volumeannualmedian;
  nc.volumedaymedian.long_name=[char(ensemblenames{l}),'sea ice volume median'];
  nc.volumedaymedian.dimension={'nRegion','dayofyear'};
  nc.volumedaymedian.units='km^3';

  nc.volumedaystd.data=volumeannualstd;
  nc.volumedaystd.long_name=[char(ensemblenames{l}),'sea ice volume standard deviation'];
  nc.volumedaystd.dimension={'nRegion','dayofyear'};
  nc.volumedaystd.units='km^3'

  nc.volumedayequiv.data=volumeannualequiv;
  nc.volumedayequiv.long_name=[char(ensemblenames{l}),'sea ice volume equivalent sample size'];
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
  nc.snowdaymean.long_name=[char(ensemblenames{l}),'sea ice snow volume mean'];
  nc.snowdaymean.dimension={'nRegion','dayofyear'};
  nc.snowdaymean.units='km^3';

  nc.snowdaymedian.data=snowannualmedian;
  nc.snowdaymedian.long_name=[char(ensemblenames{l}),'sea ice snow volume median'];
  nc.snowdaymedian.dimension={'nRegion','dayofyear'};
  nc.snowdaymedian.units='km^3';

  nc.snowdaystd.data=snowannualstd;
  nc.snowdaystd.long_name=[char(ensemblenames{l}),'sea ice snow volume standard deviation'];
  nc.snowdaystd.dimension={'nRegion','dayofyear'};
  nc.snowdaystd.units='km^3';

  nc.snowdayequiv.data=snowannualequiv;
  nc.snowdayequiv.long_name=[char(ensemblenames{l}),'sea ice snow volume equivalent sample size'];
  nc.snowdayequiv.dimension={'nRegion','dayofyear'};
  nc.snowdayequiv.units='km^3';

  nc.snowlowerquantile.data=snowlowerquantile;
  nc.snowlowerquantile.long_name=[num2str(lquantile),' quantile of ice snow volume'];
  nc.snowlowerquantile.dimension={'nRegion','dayofyear'};
  nc.snowlowerquantile.units='km^3';

  nc.snowupperquantile.data=snowupperquantile;
  nc.snowupperquantile.long_name=[num2str(uquantile),' quantile of ice snow volume'];
  nc.snowupperquantile.dimension={'nRegion','dayofyear'};
  nc.snowupperquantile.units='km^3';

  % ---

  nc.kineticdaymean.data=kineticannualmean;
  nc.kineticdaymean.long_name=[char(ensemblenames{l}),'sea ice kinetic energy mean'];
  nc.kineticdaymean.dimension={'nRegion','dayofyear'};
  nc.kineticdaymean.units='TJ';

  nc.kineticdaymedian.data=kineticannualmedian;
  nc.kineticdaymedian.long_name=[char(ensemblenames{l}),'sea ice kinetic energy median'];
  nc.kineticdaymedian.dimension={'nRegion','dayofyear'};
  nc.kineticdaymedian.units='TJ';

  nc.kineticdaystd.data=kineticannualstd;
  nc.kineticdaystd.long_name=[char(ensemblenames{l}),'sea ice kinetic energy standard deviation'];
  nc.kineticdaystd.dimension={'nRegion','dayofyear'};
  nc.kineticdaystd.units='TJ';

  nc.kineticdayequiv.data=kineticannualequiv;
  nc.kineticdayequiv.long_name=[char(ensemblenames{l}),'sea ice kinetic energy equiv sample size'];
  nc.kineticdayequiv.dimension={'nRegion','dayofyear'};
  nc.kineticdayequiv.units='TJ';

  nc.kineticlowerquantile.data=kineticlowerquantile;
  nc.kineticlowerquantile.long_name=[num2str(lquantile),' quantile of ice kinetic energy'];
  nc.kineticlowerquantile.dimension={'nRegion','dayofyear'};
  nc.kineticlowerquantile.units='TJ';

  nc.kineticupperquantile.data=kineticupperquantile;
  nc.kineticupperquantile.long_name=[num2str(uquantile),' quantile of ice kinetic energy'];
  nc.kineticupperquantile.dimension={'nRegion','dayofyear'};
  nc.kineticupperquantile.units='TJ';

  % ---

  nc.albedodaymean.data=albedoannualmean;
  nc.albedodaymean.long_name=[char(ensemblenames{l}),'sea ice albedo mean'];
  nc.albedodaymean.dimension={'nRegion','dayofyear'};

  nc.albedodaymedian.data=albedoannualmedian;
  nc.albedodaymedian.long_name=[char(ensemblenames{l}),'sea ice albedo median'];
  nc.albedodaymedian.dimension={'nRegion','dayofyear'};

  nc.albedodaystd.data=albedoannualstd;
  nc.albedodaystd.long_name=[char(ensemblenames{l}),'sea ice albedo standard deviation'];
  nc.albedodaystd.dimension={'nRegion','dayofyear'};

  nc.albedodayequiv.data=albedoannualequiv;
  nc.albedodayequiv.long_name=[char(ensemblenames{l}),'sea ice albedo equiv sample size'];
  nc.albedodayequiv.dimension={'nRegion','dayofyear'};

  nc.albedolowerquantile.data=albedolowerquantile;
  nc.albedolowerquantile.long_name=[num2str(lquantile),' quantile of ice albedo'];
  nc.albedolowerquantile.dimension={'nRegion','dayofyear'};

  nc.albedoupperquantile.data=albedoupperquantile;
  nc.albedoupperquantile.long_name=[num2str(uquantile),' quantile of ice albedo'];
  nc.albedoupperquantile.dimension={'nRegion','dayofyear'};

  nc=ridgepack_struct(nc);

  cd(char(eprnames{l}));

  outfile=[char(ensemblenames{l}),'.ensemble.lemnisc.',...
           num2str(yearst,'%4.4i'),'-',num2str(yearen,'%4.4i'),'.',filetabe];

  nc=ridgepack_struct(nc);

  ridgepack_write(nc,outfile);

 end

 end % years

end

if generateobs & observations

 for yi=1:length(yearrange)

  yearst=yearrange{yi}(1);
  yearen=yearrange{yi}(2);

  clear nc 
  clear extentannualmean extentannualmedian 
  clear extentannualstd extentannualequiv 
  clear extentlowerquantile extentupperquantile

  k=1;

  cd(['~/data/data/SATELLITE/processed/G02202_v4']);

  ncobs=ridgepack_clone('G02202_v4_merged_global_r00_1979_2022');

  years=str2num(datestr(ncobs.time.data,'YYYY'));

  idx=find(years>=yearsto & years<=yeareno);

  k=1
  ncobs.time.data=ncobs.time.data(idx);
  ncobs.extent.data=ncobs.extent.data(:,idx)/fact(k);
  ncobs=rmfield(ncobs,{'attributes'})
  ncobs.attributes.title=['G02202_v4 extent lemnisc statistics years ',...
                       num2str(yearsto,'%4.4i'),'-',num2str(yeareno,'%4.4i')];
  for k=1:365
   space=[k:365:365*floor(length(ncobs.time.data)./365)];
   annualtimeaverage(k)=ncobs.time.data(k);
   for j=1:3

    extentannualmean(j,k)=mean(ncobs.extent.data(j,space));
    extentannualmedian(j,k)=median(ncobs.extent.data(j,space));
    extentannualstd(j,k)=std(ncobs.extent.data(j,space));
    extentannualequiv(j,k)=ridgepack_equiv(ncobs.extent.data(j,space));
    extentlowerquantile(j,k)=quantile(ncobs.extent.data(j,space),lquantile);
    extentupperquantile(j,k)=quantile(ncobs.extent.data(j,space),uquantile);

   end
  end

  ncobs.attributes.dataset='NOAA/NSIDC G02202 V4 Sea Ice Climate Data Record from SMMR, SSM/I, SSMI/S'
  ncobs.attributes.regions='nRegions: 1=global, 2=Northern Hem, 3=Southern Hem';

  ncobs.dayofyear.data=[1:365];
  ncobs.dayofyear.long_name='day of year';
  ncobs.dayofyear.dimension={'dayofyear'};
  ncobs.dayofyear.units='day of year';
  ncobs.dayofyear.type='NC_INT';

  ncobs.nsamp.data=length(space).*ones(365,1);
  ncobs.nsamp.long_name='Sample Size';
  ncobs.nsamp.dimension={'dayofyear'};
  ncobs.nsamp.type='NC_INT';

  ncobs.extentdaymean.data=extentannualmean;
  ncobs.extentdaymean.long_name=['NOAA/NSIDC G02202 V4 CDR sea ice extent mean'];
  ncobs.extentdaymean.dimension={'nRegion','dayofyear'};
  ncobs.extentdaymean.units='km^2';

  ncobs.extentdaymedian.data=extentannualmedian;
  ncobs.extentdaymedian.long_name=['NOAA/NSIDC G02202 V4 CDR sea ice extent median'];
  ncobs.extentdaymedian.dimension={'nRegion','dayofyear'};
  ncobs.extentdaymedian.units='km^2';

  ncobs.extentdaystd.data=extentannualstd;
  ncobs.extentdaystd.long_name=['NOAA/NSIDC G02202 V4 CDR sea ice extent standard deviation'];
  ncobs.extentdaystd.dimension={'nRegion','dayofyear'};
  ncobs.extentdaystd.units='km^2';

  ncobs.extentdayequiv.data=extentannualequiv;
  ncobs.extentdayequiv.long_name=['NOAA/NSIDC G02202 V4 CDR sea ice extent equivalent sample size'];
  ncobs.extentdayequiv.dimension={'nRegion','dayofyear'};
  ncobs.extentdayequiv.units='km^2';

  ncobs.extentlowerquantile.data=extentlowerquantile;
  ncobs.extentlowerquantile.long_name=[num2str(lquantile),' quantile of ice extent'];
  ncobs.extentlowerquantile.dimension={'nRegion','dayofyear'};
  ncobs.extentlowerquantile.units='km^2';

  ncobs.extentupperquantile.data=extentupperquantile;
  ncobs.extentupperquantile.long_name=[num2str(uquantile),' quantile of ice extent'];
  ncobs.extentupperquantile.dimension={'nRegion','dayofyear'};
  ncobs.extentupperquantile.units='km^2';

  ncobs=ridgepack_struct(ncobs);

  outfile=['G02202_v4_merged.lemnisc.',...
           num2str(yearsto,'%4.4i'),'-',num2str(yeareno,'%4.4i'),'.',filetabe];

  ncobs=ridgepack_struct(ncobs);

  ridgepack_write(ncobs,outfile);

 end % years

end

%return

framerateindays=365;
lengthtimeseries=365*500;
%lengthtimeseries=365*(diff(yearrange)+1);

frameidx=0;

for iframe=1:framerateindays:lengthtimeseries-framerateindays

frameidx=frameidx+1;

filenamemodifier=[filetabe,'.',num2str(framerateindays),'.',...
                  num2str(iframe,'%6.6i'),'.',num2str(frameidx,'%6.6i'),'.'];

for kcols=1:maxcols

 ridgepack_multiplot(1,maxcols,1,kcols,als(kcols))

 % set plot type
 if kcols==1
  titl2=['Sea Ice Extent \times{10^6} km^2'];
  xlab2=['Northern Hemisphere'];
  ylab2=['Southern Hemisphere'];
  %xlims=[0 20];
  %ylims=[0 20];
  xlims=[4 22];
  ylims=[0 21];
  xyticks=5;
  globalconts=30;
  globlab=['Global Extent'];
  filenamemodifier=[filenamemodifier,'extent.'];
 elseif kcols==2
  titl2=['Sea Ice Volume \times{10^3} km^3'];
  xlab2=['Northern Hemisphere'];
  %xlims=[0 40];
  %ylims=[0 25];
  xlims=[5 40];
  ylims=[0 25];
  xyticks=5;
  globalconts=45;
  globlab=['Global Volume'];
  filenamemodifier=[filenamemodifier,'volume.'];
 elseif kcols==3
  titl2=['Snow Volume \times{10^2} km^3'];
  xlab2=['Northern Hemisphere'];
  xlims=[0 33];
  ylims=[0 65];
  xyticks=10;
  globalconts=50;
  globlab=['Global Snow Volume'];
  filenamemodifier=[filenamemodifier,'snow.'];
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

  for i=2:length(conts)-1

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

 % plot data for given year ranges
 for yi=1:length(yearrange)

 yearst=yearrange{yi}(1);
 yearen=yearrange{yi}(2);

 % plot data from each lemnisc
 if kcols==1 
  nlemn=nlemniscs;
 elseif observations
  nlemn=nlemniscs-1;
 end

 for l=1:nlemn

  if l==nlemniscs & observations % observed extent
   cd(['~/data/data/SATELLITE/processed/G02202_v4']);
   obsfile=['G02202_v4_merged.lemnisc.',...
             num2str(yearsto,'%4.4i'),'-',num2str(yeareno,'%4.4i'),'.',filetabe];
   nc=ridgepack_clone(obsfile);
  else
   cd(char(eprnames{l}));
   infile=[char(ensemblenames{l}),'.ensemble.lemnisc.',...
           num2str(yearst,'%4.4i'),'-',num2str(yearen,'%4.4i'),'.',filetabe];
   nc=ridgepack_clone(infile);
  end

  if l==1
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
  elseif l<nlemniscs % not compared with observations
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
   tcrit=tinv(0.975,df); % change 0.9995 99.9%, 0.995 99%, 0.975 95%
   hcrit=ones(size(t));
   hcrit(isnan(t))=0;
   hcrit(-tcrit<t & t<tcrit)=0;
  end

  if kcols==1
   upperquantile=nc.extentupperquantile.data;
   lowerquantile=nc.extentlowerquantile.data;
   daymedian=nc.extentdaymedian.data;
   daymean=nc.extentdaymean.data;
   if l==nlemniscs & observations
    totalseries=nc.extent.data;
   else
    totalseries=nc.totalIceExtent.data;
   end
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

   %plot(totalseries(2,:),totalseries(3,:),'Color',0.9*[1 1 1]);

   plot(totalseries(2,1:iframe),totalseries(3,1:iframe),'Color',0.9*[1 1 1]);
   plot(totalseries(2,iframe:iframe+365),totalseries(3,iframe:iframe+365),'Color','b');

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

 end % years


 % plot data for given year ranges
 for yi=1:length(yearrange)

 yearst=yearrange{yi}(1);
 yearen=yearrange{yi}(2);

 for l=1:nlemn

  if l==nlemniscs & observations % observed extent
   cd(['~/data/data/SATELLITE/processed/G02202_v4']);
   obsfile=['G02202_v4_merged.lemnisc.',...
             num2str(yearsto,'%4.4i'),'-',num2str(yeareno,'%4.4i'),'.',filetabe];
   nc=ridgepack_clone(obsfile);
  else
   cd(char(eprnames{l}));
   infile=[char(ensemblenames{l}),'.ensemble.lemnisc.',...
           num2str(yearst,'%4.4i'),'-',num2str(yearen,'%4.4i'),'.',filetabe];
   nc=ridgepack_clone(infile);
  end

  % grab daily mean  
  if kcols==1
   daymean=nc.extentdaymean.data;
  elseif kcols==2
   daymean=nc.volumedaymean.data;
  elseif kcols==3
   daymean=nc.snowdaymean.data;
  end

  if kcols==1
   upperquantile=nc.extentupperquantile.data;
   lowerquantile=nc.extentlowerquantile.data;
   daymedian=nc.extentdaymedian.data;
   daymean=nc.extentdaymean.data;
   if l==nlemniscs & observations
    totalseries=nc.extent.data;
   else
    totalseries=nc.totalIceExtent.data;
   end
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


  % grab t-test information
  if l==1
   if kcols==1
    mu1=nc.extentdaymean.data(1,:);
    std1=nc.extentdaystd.data(1,:);
    equiv1=nc.extentdayequiv.data(1,:);
    if l==nlemniscs & observations
     totalseries=nc.extent.data;
    else
     totalseries=nc.totalIceExtent.data;
    end
   elseif kcols==2
    mu1=nc.volumedaymean.data(1,:);
    std1=nc.volumedaystd.data(1,:);
    equiv1=nc.volumedayequiv.data(1,:);
    totalseries=nc.totalIceVolume.data;
   elseif kcols==3
    mu1=nc.snowdaymean.data(1,:);
    std1=nc.snowdaystd.data(1,:);
    equiv1=nc.snowdayequiv.data(1,:);
    totalseries=nc.totalSnowVolume.data;
   end
  elseif l<nlemniscs | ~observations
   if kcols==1
    mu2=nc.extentdaymean.data(1,:);
    std2=nc.extentdaystd.data(1,:);
    equiv2=nc.extentdayequiv.data(1,:);
    if l==nlemniscs & observations
     totalseries=nc.extent.data;
    else
     totalseries=nc.totalIceExtent.data;
    end
   elseif kcols==2
    mu2=nc.volumedaymean.data(1,:);
    std2=nc.volumedaystd.data(1,:);
    equiv2=nc.volumedayequiv.data(1,:);
    totalseries=nc.totalIceVolume.data;
   elseif kcols==3
    mu2=nc.snowdaymean.data(1,:);
    std2=nc.snowdaystd.data(1,:);
    equiv2=nc.snowdayequiv.data(1,:);
    totalseries=nc.totalSnowVolume.data;
   end
   t = (mu1-mu2)./sqrt(((std1.^2)./equiv1) + ((std2.^2)./equiv2));
   df = (equiv1+equiv2-2);
   tcrit=tinv(0.995,df); % change 0.995 99%, 0.975 95%
   hcrit=ones(size(t));
   hcrit(isnan(t))=0;
   hcrit(-tcrit<t & t<tcrit)=0;
  end

  % plot mean that is statistically significant
  if plotmean
   if l==nlemniscs & observations
    h(l)=plot(daymean(2,:),daymean(3,:),'-','Color',ccols(l,:));
   elseif l==1 
    h(l)=plot(daymean(2,:),daymean(3,:),'-','Color',ccols(l,:));
   elseif l<nlemniscs | ~observations
    plotmean=daymean;
    %plotmean(:,hcrit==1)=NaN;
    plot(plotmean(2,:),plotmean(3,:),':','Color',ccols(l,:));
    %plotmean=daymean;
    plotmean(:,hcrit==0)=NaN;
    h(l)=plot(plotmean(2,:),plotmean(3,:),'-','Color',ccols(l,:));
   end

   if plotequinoxtrend

    dayste=str2num(datestr(nc.time.data,'dd'));
    monthste=str2num(datestr(nc.time.data,'mm'));
    yearste=str2num(datestr(nc.time.data,'yyyy'));

    yearseries=min(yearste):1:max(yearste);

    % equinox (March and September, equ=3 and 9)
    for equ=[3 9]
     for yse=1:length(yearseries)
      idxs=find(dayste==21 & monthste==equ & yearste==yearseries(yse));
      xseries(yse)=mean(totalseries(2,idxs));
      yseries(yse)=mean(totalseries(3,idxs));
      ztime(yse)=nc.time.data(idxs(1));
     end

     % information
     disp('---------------------------')
     if equ==3
      disp(['March Equinox Ensemble Mean Lemnisc: ',char(legnames(l))])
     elseif equ==9
      disp(['September Equinox Ensemble Mean Lemnisc: ',char(legnames(l))])
     end
     disp(titl2)
     disp(['    NH ',num2str(yearseries(1)),': ',num2str(xseries(1))])     
     disp(['    SH ',num2str(yearseries(1)),': ',num2str(yseries(1))])     
     disp(['Global ',num2str(yearseries(1)),': ',num2str(xseries(1)+yseries(1))])     
     disp(['    NH ',num2str(yearseries(end)),': ',num2str(xseries(end))])     
     disp(['    SH ',num2str(yearseries(end)),': ',num2str(yseries(end))])     
     disp(['Global ',num2str(yearseries(end)),': ',num2str(xseries(end)+yseries(end))])     
     disp(['    NH Delta: ',num2str(diff(xseries([1 end])))])     
     disp(['    SH Delta: ',num2str(diff(yseries([1 end])))])     
     disp(['Global Delta: ',num2str(diff(xseries([1 end]))+diff(yseries([1 end])))])     
     disp('---------------------------')

     p2x=polyfit(ztime,xseries,1);
     f2x=polyval(p2x,ztime);
     decadetrend2x=10*100*(f2x(end)-f2x(1))/(length(ztime).*f2x(1));

     p2y=polyfit(ztime,yseries,1);
     f2y=polyval(p2y,ztime);
     decadetrend2y=10*100*(f2y(end)-f2y(1))/(length(ztime).*f2x(1));
 
     plot(f2x,f2y,':','Color',ccols(l,:),'linewidth',0.4)
     plot(f2x([1 end]),f2y([1 end]),'.','Color',ccols(l,:),'linewidth',0.4)
 
     theta=atan(ar*diff(f2y([1 end]))./diff(f2x([1 end])));

     text(f2x(end),f2y(end),[num2str(decadetrend2x,'%3.1f')],...
         'Rotation',theta*180/pi,...
         'FontSize',5,'HorizontalAlignment','Left',...
         'VerticalAlignment','bottom',...
         'Margin',1,'Color',ccols(l,:),'Interpreter','Tex');

     text(f2x(end),f2y(end),[num2str(decadetrend2y,'%3.1f')],...
         'Rotation',theta*180/pi,...
         'FontSize',5,'HorizontalAlignment','Left',...
         'VerticalAlignment','top',...
         'Margin',1,'Color',ccols(l,:),'Interpreter','Tex');

    end

   end

   if label & l==1 & kcols==2

    theta=atan(ar.*diff(daymean(3,[2 3]))./diff(daymean(2,[2 3])));

    text(daymean(2,1),daymean(3,1),'January',...
        'Rotation',theta*180/pi,...
        'FontSize',6,'HorizontalAlignment','left',...
        'VerticalAlignment','bottom',...
        'Margin',1,'Color',ccols(l,:),'Interpreter','Tex');

   elseif yearlabel & l==1 

    theta=atan(ar.*diff(daymean(3,[2 3]))./diff(daymean(2,[2 3])));

    text(daymean(2,1),daymean(3,1),[num2str(yearst,'%4.4i'),'-',num2str(yearen,'%4.4i')],...
        'Rotation',theta*180/pi,...
        'FontSize',6,'HorizontalAlignment','left',...
        'VerticalAlignment','bottom',...
        'Margin',1,'Color',ccols(l,:),'Interpreter','Tex');


   end

   idx=datenum(0000,3,21);
   ha=plot(daymean(2,idx),daymean(3,idx),'.','Color',ccols(l,:));

   if label & l==1 & kcols==2

    theta=atan(ar*diff(daymean(3,[idx-2 idx+2]))./...
                  diff(daymean(2,[idx-2 idx+2])));

    text(daymean(2,idx),daymean(3,idx),{'Equinox'},...
        'Rotation',theta*180/pi,...
        'FontSize',6,'HorizontalAlignment','center',...
        'VerticalAlignment','top',...
        'Margin',1,'Color',ccols(l,:),'Interpreter','Tex');

   end

   idx=datenum(0000,6,21);
   ha=plot(daymean(2,idx),daymean(3,idx),'o',...
           'MarkerFaceColor','w','MarkerSize',2,'Color',ccols(l,:));
 
   if label & l==1 & kcols==2

    theta=atan(ar*diff(daymean(3,[idx-2 idx+2]))./...
                  diff(daymean(2,[idx-2 idx+2])));

    text(daymean(2,idx),daymean(3,idx),{'Solstice'},...
        'Rotation',theta*180/pi,...
        'FontSize',5,'HorizontalAlignment','center',...
        'VerticalAlignment','top',...
        'Margin',3,'Color',ccols(l,:),'Interpreter','Tex');

   end
 
   idx=datenum(0000,9,21);
   plot(daymean(2,idx),daymean(3,idx),'.','Color',ccols(l,:));

   if label & l==1 & kcols==2

    theta=atan(ar*diff(daymean(3,[idx-2 idx+2]))./...
                  diff(daymean(2,[idx-2 idx+2])));

    text(daymean(2,idx),daymean(3,idx),{'Equinox'},...
        'Rotation',theta*180/pi,...
        'FontSize',6,'HorizontalAlignment','center',...
        'VerticalAlignment','top',...
        'Margin',1,'Color',ccols(l,:),'Interpreter','Tex');

   end
 
   idx=datenum(0000,12,21);
   ha=plot(daymean(2,idx),daymean(3,idx),'o',...
           'MarkerFaceColor','w','MarkerSize',2,'Color',ccols(l,:));

   if label & l==1 & kcols==2

    theta=atan(ar*diff(daymean(3,[idx-2 idx+2]))./...
                  diff(daymean(2,[idx-2 idx+2])));

    text(daymean(2,idx),daymean(3,idx),{'Solstice'},...
        'Rotation',theta*180/pi,...
        'FontSize',5,'HorizontalAlignment','center',...
        'VerticalAlignment','top',...
        'Margin',3,'Color',ccols(l,:),'Interpreter','Tex');

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

   plot(xarrow,yarrow,'Color',ccols(l,:))

  end

 end

 end % years

 xlabel(xlab2,'Interpreter','Tex','Fontsize',10)
 if kcols==1
  ylabel(ylab2,'Interpreter','Tex','Fontsize',10)
  legendnames=legnames;
  if plotmean
   if observations & ~yearlabel
    legendnames{nlemniscs}=[num2str(yearsto,'%4.4i'),'-',num2str(yeareno,'%4.4i'),' ',...
                            char(legnames{nlemniscs})];
   end
   legend(h([1:nlemniscs]),legendnames{1:nlemniscs},'location','southwest','FontSize',8);
   legend('boxoff');
  end
 end
 title(titl2,'Interpreter','Tex','Fontsize',10,'FontWeight','normal')

end

ridgepack_multialign(gcf,['Year = ',num2str(floor((iframe+365)./365),'%3.3i')],12)

cd('/Users/afroberts/work')

yeartag=[];
for yi=1:length(yearrange)
 yearst=yearrange{yi}(1);
 yearen=yearrange{yi}(2);
 yeartag=[yeartag,num2str(yearst,'%4.4i'),'-',num2str(yearen,'%4.4i'),'.'];
end

ensembletag=[];
for ei=1:length(ensemblenames)
 ensembletag=[ensembletag,char(ensemblenames{ei}),'.'];
end

filenamemodifier=[filenamemodifier,num2str(maxcols),'.'];

ridgepack_fprint('png',['E3SM_sea_ice_lemnisc.',yeartag,ensembletag,filenamemodifier,'png'],1,2);

clf

end % iframe


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

