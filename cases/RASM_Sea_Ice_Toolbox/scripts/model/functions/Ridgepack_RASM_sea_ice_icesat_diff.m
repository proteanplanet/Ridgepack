function Ridgepack_RASM_sea_ice_icesat(rasmcases,springfall,quicknames,columns,pub,pubdir)

delete(gcf);

alpha='abcdefghijklmnopqrstuvwxyz';

if nargin<3; quicknames=rasmcases; end

if nargin<5; pub=false; end

if columns>3
 titletext='thickness difference from Kwok and Cunningham (2008) on SSM/I mesh';
else
 titletext='thickness difference against ICESat';
end


home=getenv('HOME');
dirdata=['/Volumes/RobertsRaid3/data'];

% multiplot specifications
nrows=length(rasmcases)+1;
ncols=min(columns,5);

%%%%%%%% MODEL DATA %%%%%%%%%%%%%%%%%%%%%%%%%

% Mean thickness information for first case
dircase=['/Volumes/RobertsRaid3/work/processing/',char(rasmcases{1}),'/ice/monthly'];
fieldh='hi'; fileh=[dircase,'/',char(rasmcases{1}),'.cice.h.',fieldh];


%%%%%%%% OBSERVATIONS %%%%%%%%%%%%%%%%%%%%%%%

% Observed sea ice pathfinder ice motion 
fileo=[dirdata,'/SATELLITE/processed/kwok_icesat'];

% ICESat field list and timing arrays
fieldlist={'h_fm04','h_fm05','h_fm06','h_ma07','h_fm08','h_on03','h_on04','h_on05','h_on06','h_on07'};
taglist={'February-March 2004','February-March 2005','February-March 2006','March-April 2007','Februry-March 2008','October-November 2003','October-November 2004','October-November 2005','October-November 2006','October-November 2007'};
startmonth=[2,2,2,3,2,10,10,10,10,10];
startyear=[2004 2003];


%%%%%%%% TIME SPECIFICATION %%%%%%%%%%%%%%%%%

% time required by row and columns
for k=1:ncols
 if strcmp(springfall,'spring')
  colyear{k}=startyear(1)+k-1;
  colyear{k}=startyear(1)+k-1;
  collabel{k}=char(taglist{k});
 elseif strcmp(springfall,'fall')
  colyear{k}=startyear(2)+k-1;
  colyear{k}=startyear(2)+k-1;
  collabel{k}=char(taglist{k+5});
 else
  error('Neither spring nor fall are set')
 end
end

% Model time offset months
modeloffset=+1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do all work in case directory
cd(dircase)

% get times required in model data and set the appropriate field name from 
% ICESat netcdf file for each time period.
nctime=ridgepack_clone(fileh,{'time'});
lengthtime=length(nctime.time.data);
chartime=datestr(nctime.time.data,'yyyy-mm-dd HH:MM:SS');
slistm=zeros(1,ncols);
elistm=zeros(1,ncols);
for k=1:ncols
 if strcmp(springfall,'spring')
  month=[startmonth(k)+modeloffset startmonth(k)+1+modeloffset];
  year=startyear(1)+k-1;
 elseif strcmp(springfall,'fall')
  month=[startmonth(k+5)+modeloffset startmonth(k+5)+1+modeloffset];
  year=startyear(2)+k-1;
 end
 placement=[];
 for m=1:length(month)
 for n=1:length(year)
 for l=1:lengthtime
   index=findstr(chartime(l,:),[num2str(year(n),'%4.4i'),'-',num2str(month(m),'%2.2i'),'-']);
   if ~isempty(index); placement=[placement l]; end
  end
 end
 end
 if isempty(placement)
  slistm(k)=NaN;
  elistm(k)=NaN;
 elseif max(placement)-min(placement)~=length(placement)-1
  error('The time indices in this file are not incremental')
 else
  slistm(k)=min(placement);
  elistm(k)=max(placement);
 end
end

slistm
elistm

if all(isnan(slistm(:)));
 disp('No model data found')
 return;
end

% Work through time loop to construct diagram
for k=1:ncols

 for j=1:nrows

   ridgepack_multiplot(nrows,ncols,j,k,alpha(k+(j-1)*ncols))

   if strcmp(springfall,'spring')
    fieldo=char(fieldlist{k});
   elseif strcmp(springfall,'fall')
    fieldo=char(fieldlist{k+5});
   end

   if (j==1) 

    ncicesat=ridgepack_clone(fileo,{fieldo,'latitude','longitude','mask'});

    colscale=[0 0.25 0.5:0.5:2.0 3.0:1:6];

    ridgepack_image(ncicesat,'x','y',fieldo,{},{},colscale,'linear',3.0,'vertical','parula');

    if k<ncols 
     ridgepack_cbdelete(gca);
    end

    ridgepack_clearax

    if k==1 
     ylabel('ICESat')
    else
     ylabel('')
    end

    title(char(collabel{k}),'FontWeight','normal')
   
   else 

    dircase=['/Volumes/RobertsRaid3/work/processing/',char(rasmcases{j-1}),'/ice/monthly'];
    cd(dircase)

    fileh=[dircase,'/',char(rasmcases{j-1}),'.cice.h.',fieldh];

    nch=ridgepack_reduce(ridgepack_clone(fileh,fieldh,{'time'},{slistm(k)},{elistm(k)}),{'time'});
    ncmodel=ridgepack_regrid(nch,fieldh,'',ncicesat);
    
    ncmodel.diff=ncmodel.(fieldh);
    ncmodel.diff.data=ncmodel.(fieldh).data-ncicesat.(fieldo).data;

    %colscale=[-3:0.5:-0.5 -0.25 0.25 0.5:0.5:3.0];
    colscale=[-3 -2 -1 -0.5 0.5 1 2 3];

    ridgepack_image(ncmodel,'x','y','diff',{},{},colscale);

    if k==ncols & j==max(2,1+floor(length(rasmcases)/2))
     if nrows>2
      ridgepack_cbshare(gca);
     end
    else
     ridgepack_cbdelete(gca);
    end
 
    ridgepack_clearax

    if k==1 
     ylabel(char(quicknames{j-1}))
    else
     ylabel('')
    end

    title('')

   end

   axis tight
   axis equal

   drawnow

 end

end

if strcmp(springfall,'spring')
 if pub
  ridgepack_multialign(gcf,'',11,[0 0 0],2)
  if nargin>=6
   cd(pubdir)
  end
 else
  ridgepack_multialign(gcf,['Spring ',titletext],11,[0 0 0],2);
 end
 drawnow
 ridgepack_fprint('png',[ridgepack_cellcat(rasmcases,'_'),'_icesat_spring_diff_',num2str(startyear(1))],1,1)
 ridgepack_fprint('epsc',[ridgepack_cellcat(rasmcases,'_'),'_icesat_spring_diff_',num2str(startyear(1))],1,1)
elseif strcmp(springfall,'fall')
 if pub
  ridgepack_multialign(gcf,'',11,[0 0 0],2)
  if nargin>=6
   cd(pubdir)
  end
 else
  ridgepack_multialign(gcf,['Fall ',titletext]);
 end
 drawnow
 ridgepack_fprint('png',[ridgepack_cellcat(rasmcases,'_'),'_icesat_fall_diff_',num2str(startyear(1))],1,1)
 ridgepack_fprint('epsc',[ridgepack_cellcat(rasmcases,'_'),'_icesat_fall_diff_',num2str(startyear(1))],1,1)
end




