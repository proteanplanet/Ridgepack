function Ridgepack_RASM_sea_ice_pathfinder_stream_diff(rasmcases,quicknames,...
                               startyear,endyear,format,pub,pubdir,monthset)

delete(gcf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear
%clf %rasmcases={'r26RBRCICE5gx','r26RBRCICE5g8','r26RBRCICE5g0'};
%quicknames={'EVP','Revised-EVP','Anisotropic'};
%startyear=1990;
%endyear=1992;

alpha='abcdefghijklmnopqrstuvwxyz';

% month averages for each column
if nargin<8
 monthset={[01 02 03],[04 05 06],[07 08 09],[10 11 12]};
 %monthset={[01 02 03]};
elseif ~iscell(monthset)
 error('monthset must be a cell array')
end

% multiplot specifications
nrows=length(rasmcases)+1;
ncols=length(monthset);

if length(rasmcases)~=length(quicknames)
 error('quicknames not same length as rasmcases')
end

% read in mask
ridgepack_mask=ridgepack_clone('/Volumes/RobertsRaid3/data/MODEL/RASM/RASM_POPCICE_GRID_MASKS_AND_METRICS',...
               {'mask_centralarctic','latitude','longitude'});

% set field names
meanend='';
fieldu=['uvel',meanend]; 
fieldv=['vvel',meanend]; 

% generate speed files if they don't yet exist
Ridgepack_RASM_sea_ice_make_speed(rasmcases,fieldu,fieldv);

count=0;

for j=1:nrows

 % set rasm cases
 if j>1
  rasmcase=char(rasmcases{j-1});
  quickname=char(quicknames{j-1});
 else
  rasmcase='Pathfinder';
  quickname='Pathfinder';
 end

 % set cases
 dirdata=['/Volumes/RobertsRaid3/data'];
 if j>1
  dircase=['/Volumes/RobertsRaid3/work/processing/',rasmcase,'/ice/monthly'];
  cd(dircase)
 end

 if j>1

  %%%%%%%% MODEL DATA %%%%%%%%%%%%%%%%%%%%%%%%%

  % use sea ice area to get formar right
   % delimiter and mean ending of variables
   delim='.cice.h.';
 
   % field names
   fieldu=['uvel',meanend]; 
   fieldv=['vvel',meanend]; 
   fielda=['aice',meanend]; 
   fields=['speed',meanend]; 

   % set file locations and names
   fileu=[dircase,'/',rasmcase,delim,fieldu];
   filev=[dircase,'/',rasmcase,delim,fieldv];
   filea=[dircase,'/',rasmcase,delim,fielda];
   files=[dircase,'/',rasmcase,delim,fields];
  
   % run a test to check for file format
   nctest=ridgepack_clone(fileu,fieldu,1);

 else 

  clear nctest

  % Observed sea ice parthfinder ice motion
  fieldu='u'; 
  fileou=[dirdata,...
   '/SATELLITE/processed/Pathfinder_icemotion_monthly_1979_2016_v3_RASM_CICE_time_bounds'];

  fieldv='v'; 
  fileov=fileou;

  fields='speed'; 
  fileos=[dirdata,...
   '/SATELLITE/processed/Pathfinder_icemotion_monthly_1979_2016_v3_RASM_CICE_time_bounds_speed'];

 end

 for k=1:ncols

  count=count+1;

  % set months
  months=monthset{k};

  if j>1

    ncu=ridgepack_timesubset(fileu,fieldu,months,startyear,endyear);
    ncv=ridgepack_timesubset(filev,fieldv,months,startyear,endyear);
    nca=ridgepack_timesubset(filea,fielda,months,startyear,endyear);
    ncs=ridgepack_timesubset(files,fields,months,startyear,endyear);
 
    ncu.(fieldu).data(nca.(fielda).data<0.15)=0;
    ncv.(fieldv).data(nca.(fielda).data<0.15)=0;
    ncs.(fields).data(nca.(fielda).data<0.15)=0;

    ncu.(fieldu).data=ncu.(fieldu).data*100;
    ncv.(fieldv).data=ncv.(fieldv).data*100;
    ncs.(fields).data=ncs.(fields).data*100;
    ncs.([fields,'_std']).data=ncs.([fields,'_std']).data*100;

    ncu.(fieldu).units='\times\,10^{-2}';
    ncv.(fieldv).units='\times\,10^{-2}';
    ncs.(fields).units={'\times\,10^{-2}','m\;s^{-1}'};
 
    ncso=ridgepack_timesubset(fileos,fields,months,startyear,endyear);

  else
 
    % extract pathfinder data
    ncu=ridgepack_timesubset(fileou,fieldu,months,startyear,endyear);
    ncv=ridgepack_timesubset(fileov,fieldv,months,startyear,endyear);
    ncso=ridgepack_timesubset(fileos,fields,months,startyear,endyear);

    ncu.(fieldu).units='\times\,10^{-2}';
    ncv.(fieldv).units='\times\,10^{-2}';
    ncso.(fields).units='\times\,10^{-2}';
 
  end

  % mask for central arctic
  if format==2
   ncu.(fieldu).data=ncu.(fieldu).data.*ridgepack_mask.mask_centralarctic.data;
   ncv.(fieldv).data=ncv.(fieldv).data.*ridgepack_mask.mask_centralarctic.data;
   if j>1
    ncs.(fields).data(ridgepack_mask.mask_centralarctic.data==0)=NaN;
   end
   ncso.(fields).data(ridgepack_mask.mask_centralarctic.data==0)=NaN;
  end

  if pub
   ridgepack_multiplot(nrows,ncols,j,k,alpha(count))
  else
   ridgepack_multiplot(nrows,ncols,j,k)
  end

  if j==1
   contcol=[0:1:6];
   if k==ncols
    ridgepack_image(ncso,'y','x','speed',{},{},contcol,'linear',0,'vertical','parula')
   else
    ridgepack_image(ncso,'y','x','speed',{},{},contcol,'linear',0,'none','parula')
   end
  else
   contcol=[-5:1:-1 1:1:5];
   ncdiff=ridgepack_ttest2(ncs,'speed',ncso,'speed',2);

   % create checker board for values that are not statistically significant
   cs=8; % cs sets the size of the checks (8=eight grid cells a side)
   for ai=1:cs:size(ncdiff.diff.data,1)-cs;
   for bi=(mod(ai,2*cs)+1):2*cs:size(ncdiff.diff.data,2)-cs;
    if mean(mean(ncdiff.h.data(ai:ai+cs-1,bi:bi+cs-1)))<0.5 & ...
       mean(mean(ridgepack_mask.mask_centralarctic.data(ai:ai+cs-1,bi:bi+cs-1)))>0.5
       ncdiff.diff.data(ai:ai+cs-1,bi:bi+cs-1)=NaN;
    end
   end
   end
   if k==ncols & j==max(2,1+floor(length(rasmcases)/2))
    ridgepack_image(ncdiff,'y','x','diff',{},{},contcol,'linear',0,'vertical')
    if nrows>2
     ridgepack_cbshare(gca);
    end
   else
    ridgepack_image(ncdiff,'y','x','diff',{},{},contcol,'linear',0,'none')
   end
  end
 
  % plot speed (swap u and v and take transpose since x and y axes are switched)
  h=streamslice(ncv.(fieldv).data',ncu.(fieldu).data',5);
  if j==1
   set(h(:),'Color','w','Linewidth',0.4);
  else
   set(h(:),'Color','k','Linewidth',0.4);
  end

  if format==1
   %ylim([5 1275]);
   %xlim([20 710]);
   ylim([100 1250]);
   xlim([20 710]);
  elseif format==2
   ylim([450 850]);
   xlim([160 540]);
  else
   error('format not represented')
  end
  ridgepack_clearax;
  axis ij
 
  if j==1;
   title([datestr(datenum(0,months(1),1),'mmm'),'-',...
          datestr(datenum(0,months(end),1),'mmm')],...
         'FontWeight','normal','FontSize',10)
  else
   title('');
  end
 
  if k==1;
   ylabel(quickname,'FontWeight','normal','FontSize',10)
  end
 
  drawnow

 end

end

if pub
 ridgepack_multialign(gcf,'',9,[0 0 0],5)
 if nargin==7
  cd(pubdir)
 end
else
 ridgepack_multialign(gcf,['RASM comparison with Pathfinder motion ',num2str(startyear),'-',num2str(endyear)]);
end

if format==1
 ridgepack_fprint('png',[ridgepack_cellcat(rasmcases,'_'),'_pathfinder_climatology_diff_',...
           num2str(startyear),'_',num2str(endyear),'_panarctic','.png'],1,1);
 if pub
  ridgepack_fprint('epsc',[ridgepack_cellcat(rasmcases,'_'),'_pathfinder_climatology_diff_',...
           num2str(startyear),'_',num2str(endyear),'_panarctic','.eps'],1,1);
 end
elseif format==2
 ridgepack_fprint('png',[ridgepack_cellcat(rasmcases,'_'),'_pathfinder_climatology_diff_',...
           num2str(startyear),'_',num2str(endyear),'_centralarctic','.png'],1,1);
 if pub
  ridgepack_fprint('epsc',[ridgepack_cellcat(rasmcases,'_'),'_pathfinder_climatology_diff_',...
           num2str(startyear),'_',num2str(endyear),'_centralarctic','.eps'],1,1);
 end
end
 
%ridgepack_fprint('epsc',[char(rasmcases{nrows-1}),'_pathfinder_climatology_',num2str(startyear),'_',num2str(endyear),'.eps'],1,1);


