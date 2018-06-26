function RASM_sea_ice_summary_climatology(rasmcases,quicknames,startyear,endyear,minthick,format,pub,pubdir,monthset)

delete(gcf);

alpha='abcdefghijklmnopqrstuvwxyz';

% set pub to false
if nargin<7
 pub=false
end

% number of categories in ice model and min plotted thickness
ncats=5;

% month averages for each column
if nargin<9
 monthset={[01 02 03],[04 05 06],[07 08 09],[10 11 12]};
elseif ~iscell(monthset)
 error('monthset must be a cell array')
end

% multiplot specifications
nrows=length(rasmcases);
ncols=length(monthset);

% check number of columns
if ncols>12
 error('Number of columns greater than 12')
end

% set basemap for format==1 based on number of columns
if format==1 & ncols<5
 basemap='seaicerasm';
elseif format==1
 basemap='seaice';
elseif format==2
 basemap='centralarctic2';
else
 error('format is incorrect')
end

if length(rasmcases)~=length(quicknames)
 error('quicknames not same length as rasmcases')
end

if ~pub; disp('This is not for a publication'); end

% read in mask
ncm=ridgepack_clone('/Users/aroberts/data/MODEL/RASM/RASM_POPCICE_GRID_MASKS_AND_METRICS',...
               {'mask_centralarctic','latitude','longitude'});

count=0;

for j=1:nrows

 % set rasm cases
 rasmcase=char(rasmcases{j});
 quickname=char(quicknames{j});

 % set cases
 home=getenv('HOME');
 dirdata=[home,'/data'];
 dircase=[home,'/work/processing/',rasmcase,'/ice/monthly'];
 
 % do all work in case directory
 cd(dircase)
 

 %%%%%%%% MODEL DATA %%%%%%%%%%%%%%%%%%%%%%%%%

 % use sea ice area to get formar right
 try

  % delimiter and mean ending of variables
  delim='.cice.h.';
  meanend='';
 
  % Mean sea ice concentration
  fielda='aice'; filea=[dircase,'/',rasmcase,delim,fielda];

  % run a test to check for file format
  nctest=ridgepack_clone(filea,fielda,1);

 catch

  % delimiter and mean ending of variables
  delim='.cice.h1.';
  meanend='_m';
 
  % Mean sea ice concentration
  fielda='aice'; filea=[dircase,'/',rasmcase,delim,fielda];

  % run a test to check for file format
  nctest=ridgepack_clone(filea,fielda,1);

 end

 clear nctest

 % Category Volume 
 fieldvn=['vicen',meanend];
 filevn=[dircase,'/',rasmcase,delim,char(fieldvn)];
 fieldan=['aicen',meanend];
 filean=[dircase,'/',rasmcase,delim,char(fieldan)];

 % Velocity descriptors
 if format==2
  fieldu=['uvel',meanend]; fileu=[dircase,'/',rasmcase,delim,fieldu];
  fieldv=['vvel',meanend]; filev=[dircase,'/',rasmcase,delim,fieldv];
 end

 % Mean thickness
 fieldh='hi'; fileh=[dircase,'/',rasmcase,delim,fieldh];
 

 %%%%%%%% OBSERVATIONS %%%%%%%%%%%%%%%%%%%%%%%
 
 % Observed sea ice concentration
 fieldoa='conc'; 
 fileoa=[dirdata,...
  '/SATELLITE/processed/G02202_v3_merged_conc_north_1979_2017_RASM_CICE_time_bounds'];

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 m=1;
 
 for k=1:ncols

  % set months
  months=monthset{k};

  % counter for adding letter to frame
  count=count+1

  % get velocity data
  if format==2
   ncum=ridgepack_timesubset(fileu,fieldu,months,startyear,endyear);
   ncvm=ridgepack_timesubset(filev,fieldv,months,startyear,endyear);
  end

  % get total area and thickness 
  ncam=ridgepack_timesubset(filea,fielda,months,startyear,endyear);
  nchm=ridgepack_timesubset(fileh,fieldh,months,startyear,endyear);
 
  % calculate adjusted thickness and concentration
  ncvn=ridgepack_timesubset(filevn,fieldvn,months,startyear,endyear);
  ncan=ridgepack_timesubset(filean,fieldan,months,startyear,endyear);

  ncvn.(fieldvn).data=ncvn.(fieldvn).data./ncan.(fieldan).data;
  ncvn.(fieldvn).data(ncvn.(fieldvn).data<minthick)=0;
  ncan.(fieldan).data(ncvn.(fieldvn).data==0)=0;
  ncvn.(fieldvn).data=ncvn.(fieldvn).data.*ncan.(fieldan).data;

  nchm.(fieldh).data=squeeze(sum(ncvn.(fieldvn).data,1));
  ncam.(fielda).data=squeeze(sum(ncan.(fieldan).data,1));

  ncoam=ridgepack_timesubset(fileoa,fieldoa,months,startyear,endyear);

  % fix line of extent in observations and model 
  ncam.(fielda).data(ncam.(fielda).data<0.15)=0;
  ncam.(fielda).data(ncam.(fielda).data>0.14)=1;

  nchm.(fieldh).data(ncam.(fielda).data<0.5)=NaN;

  if format==2
   ncum.(fieldu).data(ncam.(fielda).data<0.5)=NaN;
   ncvm.(fieldv).data(ncam.(fielda).data<0.5)=NaN;
  end

  ncoam.(fieldoa).data(ncoam.(fieldoa).data<15)=0;
  ncoam.(fieldoa).data(ncoam.(fieldoa).data>14)=1;
  if startyear<1988
   ncoam.(fieldoa).data(ncoam.latitude.data>=84)=1;
  else
   ncoam.(fieldoa).data(ncoam.latitude.data>87)=1;
  end
 
  if pub
   ridgepack_multiplot(nrows,ncols,j,k,alpha(count))
  else
   ridgepack_multiplot(nrows,ncols,j,k)
  end

  if j==1 & k==1
   ridgepack_polarm(basemap,'grid','label')
  else
   ridgepack_polarm(basemap)
  end

  if j==1 & k==1
   ridgepack_pcolorm(nchm,fieldh,{},{},[0 0.1 0.5 1:0.5:6],'linear',0,'vertical');
   ridgepack_multicb;
  else
   ridgepack_pcolorm(nchm,fieldh,{},{},[0 0.1 0.5 1:0.5:6],'linear',0,'none');
  end

  if format==2
   vecthin=3;
   if j==1 & k==ncols
    ridgepack_quiverm(ncum,fieldu,ncvm,fieldv,{},{},vecthin,'k',0.1)
   else
    ridgepack_quiverm(ncum,fieldu,ncvm,fieldv,{},{},vecthin,'k',0.1)
    ridgepack_vecdelete
   end
  elseif format==3
   disp('NOT IMPLEMENTED')
  end
 
  if j==1;
   if ncols<5
    title([datestr(datenum(0,months(1),1),'mmmm'),'-',...
           datestr(datenum(0,months(end),1),'mmmm')],...
           'FontWeight','normal','FontSize',10)
   else
    title([datestr(datenum(0,months(1),1),'mmm'),'-',...
           datestr(datenum(0,months(end),1),'mmm')],...
           'FontWeight','normal','FontSize',10)
   end
  else
   title('');
  end
 
  if k==1;
   ylabel(quickname,'FontSize',10)
  end
 
  h1=ridgepack_maskm(ncoam.latitude.data,ncoam.longitude.data,ncoam.(fieldoa).data,'m',0.4);
 
  drawnow

  m=m+3;
 
 end

end

if ~pub
 ridgepack_multilegend([h1],{['Satellite']},'South');
end

if pub
 ridgepack_multialign(gcf,'',11,[0 0 0],3)
 if nargin>=8
  cd(pubdir)
 end
else
 ridgepack_multialign(gcf,['RASM mean state ',num2str(startyear),'-',num2str(endyear)],10,[0 0 0],3);
end

if format==1
 ridgepack_fprint('png',[ridgepack_cellcat(rasmcases,'_'),'_state_summary_',num2str(startyear),...
          '_',num2str(endyear),'_',num2str(minthick),'.png'],1,1)
 ridgepack_fprint('epsc',[ridgepack_cellcat(rasmcases,'_'),'_state_summary_',num2str(startyear),...
          '_',num2str(endyear),'_',num2str(minthick),'.eps'],1,1)
elseif format==2
 ridgepack_fprint('png',[ridgepack_cellcat(rasmcases,'_'),'_state_summary_central_vel_',num2str(startyear),...
          '_',num2str(endyear),'_',num2str(minthick),'.png'],1,1)
 ridgepack_fprint('epsc',[ridgepack_cellcat(rasmcases,'_'),'_state_summary_central_vel_',num2str(startyear),...
          '_',num2str(endyear),'_',num2str(minthick),'.eps'],1,1)
elseif format==3
 ridgepack_fprint('png',[ridgepack_cellcat(rasmcases,'_'),'_state_summary_central_sig_',num2str(startyear),...
          '_',num2str(endyear),'_',num2str(minthick),'.png'],1,1)
 ridgepack_fprint('epsc',[ridgepack_cellcat(rasmcases,'_'),'_state_summary_central_sig_',num2str(startyear),...
          '_',num2str(endyear),'_',num2str(minthick),'.eps'],1,1)
end


