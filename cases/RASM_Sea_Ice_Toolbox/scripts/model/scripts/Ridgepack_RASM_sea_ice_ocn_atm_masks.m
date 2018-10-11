% RASM script to be used with Ridgepack
%
% This MATLAB script generates masks on the 9km RASM Ocean/Ice
% grid and then calculates the appropriate masks, with area 
% fractions on the associated WRF grid.  In this case, the WRF
% grid is assumed to be a 25km grid.  The actual resolution is 
% determined automatically from the history files provided for 
% the atmosphere (9km is assumed for ice and ocean). The files 
% used and their location are set by the user.  
% 
% Written by Andrew Roberts
% Naval Postgraduate School, April 2015
% Los Alamos National Laboratory, Updated October 2018


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clf

titletext='RASM Mask';

% set location of datafiles
home=getenv('HOME');
dirmode=[home,'/data/MODEL/RASM/'];

% set names of datafiles to be used to construct masks
cicefile='R2100lRBRca.cice.h.2005-01.nc'; % CICE
popfile='R2100lRBRca.pop.h.1997-09.nc'; % POP
cplfile='R2100lRBRca.cpl.ha.2007-01.nc'; % CPL
wrffile='R2100lRBRca.wrf.ha.1997-09.nc'; % WRF

popprocess=1;
wrfprocess=1;

% do all work in model directory
cd(dirmode)

% add regional masks
region={'mask_centralarctic',...
        'mask_bering','mask_okhotsk','mask_japan','mask_pacific',...
        'mask_kara','mask_barents','mask_norwegian','mask_greenland',...
        'mask_denmarkstrait','mask_labrador','mask_baffin',...
        'mask_hudsonbay','mask_archipelago','mask_baltic',...
        'mask_northsea','mask_atlantic','mask_canada','mask_nansen',...
        'mask_shelf','mask_chukchi','mask_eastsiberian','mask_laptev'};

longname={'All of Central Arctic',...
          'Bering Sea','Sea of Okhotsk','Sea of Japan','Pacific Boundary',...
          'Kara Sea','Barents Sea','Norwegian Sea','Greenland Sea',...
          'Denmark Strait and South','Labrador Sea','Baffin Bay',...
          'Hudson Bay','Canadian Archipelago','Baltic Sea',...
          'North Sea','Atlantic Boundary','Canada-Makarov Basins',...
          'Amundsen-Nansen Basins','Shelf Area','Chukchi Sea',...
          'East Siberian Sea','Laptev Sea'}

includedregions=[1:1:20];

x{1}=[833 847 835 811 801 775 769 767 753 739 721 712 707 697 703 691 689 671 671 ...
      559 507 489 485 479 481 453 449 623 635 647 675 679 735 767 811 833 829];
y{1}=[422 392 356 356 346 318 304 298 278 268 251 254 244 220 204 196 196 178 140 ...
      156 208 236 252 280 328 344 700 692 528 518 514 502 516 516 500 491 484];
 
x{2}=[529 472 439 389 373 355 344 335 327 321 317 308 305 296 293 294 293 297 ...
      311 314 325 338 375 476 487 527 551 529];
y{2}=[181 136 130 148 155 165 172 179 186 192 200 212 220 226 241 255 282 289 ...
      330 341 363 389 366 296 256 212 174 181];

x{3}=[415 399 377 317 265 253 235 215 201 191 183 179 191 202 203 211 ...
      223 224 241 287 311 465 425];
y{3}=[456 346 364 398 402 412 424 438 462 482 502 516 529 529 526 526 ...
      522 522 516 510 610 534 342];

x{4}=[289 283 251 215 179 177 169 149 113  91  97  93  71 111 301 291];
y{4}=[512 506 510 522 530 544 546 548 558 600 628 648 658 670 678 508];

x{5}=[ 43  79 611 625 40  40];
y{5}=[662 662 658  40 40 666];

x{6}=[731 704 675 665 657 642 620 643 710 752 760 767 777 750 731 722 731];
y{6}=[519 478 503 514 517 518 543 709 701 652 644 638 629 585 577 570 520];

x{7}=[365 619 911 891 849 845 839 828 824 825 763 730 673 643];
y{7}=[530 696 698 582 513 510 507 502 499 488 516 468 514 518];

x{8}=[893 1031 1042 996 923 887 843];
y{8}=[592  592  433 476 498 500 510];

x{9}=[819 809 957 1037 1024 937];
y{9}=[534 390 382  434  460 528];

x{10}=[991 1069 1100 1089 1069 1057 1049 1037 957];
y{10}=[294  290  290  328  356  382  410  434 382];

x{11}=[1100 1129 1033 1015 1009 983 980 976 969 945 929 955 975];
y{11}=[ 295  182  164  178  184 184 186 186 186 184 218 234 334];

x{12}=[975 955 983 945 929 861 851 837 829 821 811 795 827];
y{12}=[306 234 184 184 218 232 218 234 258 276 294 352 360];

x{13}=[1015 979 933 925 869 857 869 799 1101];
y{13}=[ 180 204 200 218 228 212 174  16   24];

x{14}=[803 595 785 871 857 831];
y{14}=[346 192  18 178 216 348];

x{15}=[1059 967 913 911 1129 1091];
y{15}=[ 646 598 602 706  702  656];

x{16}=[1011 1021 1061 1091 1105 1123 1129 1165 1165 1173 1091 1013];
y{16}=[ 580  542  542  534  530  540  554  582  590  614  692  618];

x{17}=[1240 1240 943 943];
y{17}=[  40  670 670  40];

x{18}=[573 582 584 584 579 572 566 565 568 584 598 605 607 ...
       605 599 596 594 598 607 624 639 641 861 896 693 529];
y{18}=[244 274 294 305 313 323 337 341 346 371 397 411 424 ...
       439 452 470 483 492 500 502 504 516 574 373 159 160];

x{19}=x{18};
y{19}=y{18};

x{20}=x{1};
y{20}=y{1};

x{21}=x{1};
y{21}=y{1};

x{22}=x{1};
y{22}=y{1};

x{23}=x{1};
y{23}=y{1};

legcount=0;
legtext='';


if popprocess

 % model history files
 filem=[dirmode,cicefile];
 filep=[dirmode,popfile];

 % get basic model file
 ncnew=ridgepack_clone(filem,{'TLAT','TLON','tmask','tarea','ANGLET'});
 ncnew=rmfield(ncnew,{'time','time_bounds','d2'});
 ncnew.mask.data(isnan(ncnew.mask.data))=0;

 % load in depths from POP file
 nca=ridgepack_clone(filep,{'HT','TLONG','TLAT'});
 ncnew.HT=nca.HT;
 ncnew.HT.data=nca.HT.data';
 %ncnew.latitude.data=nca.latitude.data';
 %ncnew.longitude.data=nca.longitude.data';
 %nca=ridgepack_clone(filep,{'ULONG','ULAT'});
 %ncnew.ULAT.data=nca.latitude.data';
 %ncnew.ULON.data=nca.longitude.data';
 clear nca;

 % clean up areas, lat/long and attributes
 ncnew.area=ncnew.tarea;
 ncnew=rmfield(ncnew,{'ULAT','ULON','uarea','tarea','blkmask','ANGLE','turn'});
 ncnew.attributes.title='ice-ocean mask for RASM ocn/atm 9/25km';
 ncnew.attributes.source='Regional Arctic System Model';
 ncnew.attributes=rmfield(ncnew.attributes,{'comment'});

 % add mask key array
 ncnew.mask_key=ncnew.HT;
 ncnew.mask_key.data=NaN*zeros(size(ncnew.mask_key.data));
 ncnew.mask_key.data(ncnew.mask.data==0)=NaN;
 ncnew.mask_key.long_name='Reference Masks';
 ncnew.mask_key.units='';

 for i=includedregions

  regions=char(region{i});
  ncnew.(regions)=ncnew.mask;
  ncnew.(regions).long_name=[num2str(i,'%2.2i'),' ',char(longname{i})];

  [X,Y]=meshgrid(ncnew.x.data,ncnew.y.data);
  [IN ON]=inpolygon(X,Y,x{i},y{i});

  % smooth out the points
  zi=filter2(fspecial('disk',10),IN,'same');
  IN=zi; IN(zi<0.5)=0; IN(zi>=0.5)=1;

  ncnew.(regions).data=IN;
  ncnew.(regions).data(ncnew.mask.data<1)=0;

  if i==2 | i==3
   ncnew.(regions).data(ncnew.mask_centralarctic.data==1 & IN==1)=0;
  elseif i==4
   ncnew.(regions).data(ncnew.(char(region{3})).data==1 & IN==1)=0;
  elseif i==5
   ncnew.(regions).data((ncnew.(char(region{2})).data==1 | ...
                       ncnew.(char(region{3})).data==1 | ...
                       ncnew.(char(region{4})).data==1 | ...
                       ncnew.mask_centralarctic.data==1) & ...
                       IN==1)=0;
  elseif i==6
   ncnew.(regions).data(ncnew.mask_centralarctic.data==1 & ...
                        IN==1)=0;
  elseif i==7
   ncnew.(regions).data((ncnew.(char(region{6})).data==1 | ...
                         ncnew.mask_centralarctic.data==1 ) & ...
                        IN==1)=0;
  elseif i==8
   ncnew.(regions).data((ncnew.(char(region{7})).data==1 | ...
                         ncnew.mask_centralarctic.data==1) & ...
                         IN==1)=0;
  elseif i==9
   ncnew.(regions).data((ncnew.(char(region{8})).data==1 | ...
   			ncnew.(char(region{7})).data==1 | ...
                         ncnew.mask_centralarctic.data==1) & ...
                        IN==1)=0;
  elseif i==10
   ncnew.(regions).data((ncnew.(char(region{9})).data==1 | ...
                        ncnew.mask_centralarctic.data==1) & ...
                        IN==1)=0;
  elseif i==11
   ncnew.(regions).data((ncnew.(char(region{10})).data==1 | ...
			ncnew.(char(region{7})).data==1 | ...
                        ncnew.mask_centralarctic.data==1) & ...
                        IN==1)=0;
  elseif i==12
   ncnew.(regions).data((ncnew.(char(region{11})).data==1 | ...
                        ncnew.mask_centralarctic.data==1) & ...
                        IN==1)=0;
  elseif i==13
   ncnew.(regions).data((ncnew.(char(region{11})).data==1 | ...
			ncnew.(char(region{12})).data==1 | ...
                        ncnew.mask_centralarctic.data==1) & ...
                        IN==1)=0;
  elseif i==14
   ncnew.(regions).data((ncnew.(char(region{12})).data==1 | ...
			ncnew.(char(region{13})).data==1 | ...
                        ncnew.mask_centralarctic.data==1) & ...
                        IN==1)=0;
  elseif i==16
   ncnew.(regions).data((ncnew.(char(region{8})).data==1 | ...
			ncnew.(char(region{15})).data==1 | ...
                        ncnew.mask_centralarctic.data==1) & ...
                        IN==1)=0;
  elseif i==17
   ncnew.(regions).data((ncnew.(char(region{16})).data==1 | ...
			ncnew.(char(region{11})).data==1 | ...
			ncnew.(char(region{10})).data==1 | ...
			ncnew.(char(region{8})).data==1 | ...
			ncnew.(char(region{13})).data==1 | ...
			ncnew.(char(region{12})).data==1 | ...
			ncnew.(char(region{15})).data==1 | ...
			ncnew.(char(region{9})).data==1 | ...
                        ncnew.mask_centralarctic.data==1) & ...
                        IN==1)=0;
  elseif i==18
   ncnew.(regions).data((ncnew.mask_centralarctic.data==0 | ...
		        wrapTo360(ncnew.longitude.data)<145 | ...
		        wrapTo360(ncnew.longitude.data)>305) & ...
                       IN==1)=0;
  elseif i==19
   ncnew.(regions).data((ncnew.mask_centralarctic.data==0 | ...
                        ncnew.(char(region{18})).data==1) & ...
                        IN==1)=0;
  elseif i==20
   ncnew.(regions).data((ncnew.(char(region{18})).data==1 | ...
		        ncnew.(char(region{19})).data==1 ) & ...
                        IN==1)=0;
  elseif i==21
   ncnew.(regions).data((ncnew.(char(region{18})).data==1 | ...
		        ncnew.(char(region{19})).data==1) & ...
		        wrapTo360(ncnew.longitude.data)<182 & ...
                        IN==1)=0;
  elseif i==22
   ncnew.(regions).data((ncnew.(char(region{18})).data==1 | ...
		        ncnew.(char(region{19})).data==1 | ...
		        ncnew.(char(region{21})).data==1 | ...
		        wrapTo360(ncnew.longitude.data)>190 | ...
		        wrapTo360(ncnew.longitude.data)<142) & ...
                        IN==1)=0;
  elseif i==23
   ncnew.(regions).data((ncnew.(char(region{18})).data==1 | ...
		        ncnew.(char(region{19})).data==1 | ...
		        ncnew.(char(region{22})).data==1 | ...
		        ncnew.(char(region{6})).data==1 | ...
		        wrapTo360(ncnew.longitude.data)>150 | ...
		        wrapTo360(ncnew.longitude.data)<90) & ...
                        IN==1)=0;
  end


  % labeling for the figure
  legcount=legcount+1;
  ncnew.mask_key.data(ncnew.(regions).data==1)=legcount;
  [Y,X]=find(ncnew.(regions).data==1);
  xpos{legcount}=median(X);
  ypos{legcount}=median(Y);
  labe{legcount}=num2str(legcount,'%2.2i');
  legtext=[legtext,num2str(i,'%2.2i'),'-',char(longname{i}),'  '];
  if mod(legcount,4)==0
   leg{legcount/4}=legtext;
   legtext='';
  elseif i==includedregions(end)
   leg{ceil(legcount/4)}=legtext;
   legtext='';
  end

 end

 ncnew=ridgepack_struct(ncnew);

 clf 

 ridgepack_image(ncnew,'x','y','mask_key',{},{},[1:legcount+1],'linear',0,'vertical')
 ridgepack_clearax
 ridgepack_cbdelete

 if isfield(ncnew,'mask_canada') & ...
    isfield(ncnew,'mask_shelf') & ...
    isfield(ncnew,'mask_nansen')
  maskcentral=ncnew.mask_shelf.data;
  maskcentral(ncnew.mask_shelf.data==0)=NaN;
  ridgepack_mask(ncnew.x.data,ncnew.y.data,maskcentral','k')
  maskcentral=ncnew.mask_nansen.data;
  maskcentral(ncnew.mask_nansen.data==0)=NaN;
  ridgepack_mask(ncnew.x.data,ncnew.y.data,maskcentral','k')
  maskcentral=ncnew.mask_canada.data;
  maskcentral(ncnew.mask_canada.data==0)=NaN;
  ridgepack_mask(ncnew.x.data,ncnew.y.data,maskcentral','k')
 end

 if isfield(ncnew,'mask_centralarctic')
  text(xpos{1},ypos{1},labe{1},'Color','w','HorizontalAlignment','center')
  startcount=2;
 else
  startcount=1;
 end

 for i=startcount:legcount
  text(xpos{i},ypos{i},labe{i},'Color','k','HorizontalAlignment','center')
 end
 xlabel(leg)
 title('RASM regional masks on the ice-ocean domain')
 
 ridgepack_fprint('png',['RASM_ocean_mask_key'],1,1);

 ncnew=ridgepack_shuffle(ncnew,{'x' 'y'});

 ridgepack_write(ncnew,'RASM_POPCICE_GRID_MASKS_AND_METRICS_25km')

else 

 ncnew=ridgepack_clone('RASM_POPCICE_GRID_MASKS_AND_METRICS_25km')

end

if wrfprocess

 ncwrfland=ridgepack_reduce(ridgepack_clone(cplfile,...
             {'x2aavg_Sf_lfrac','doma_lat','doma_lon'}),{'time'});

 ncwrfgrid=ncwrfland;

 ncwrfheight=ridgepack_reduce(ridgepack_clone(wrffile,'HGT'),{'time'});

 ncnew=ridgepack_shuffle(ncnew,{'y' 'x'});

 ncwrfnew=ridgepack_regrid(ncnew,'mask_centralarctic','',ncwrfgrid);
 
 ncwrfnew.attributes.title='RASM WRF/VIC Model Metrics';
 ncwrfnew.attributes.source='';

 ncwrfnew.mask=ncwrfgrid.latitude;
 ncwrfnew.mask.data=1.-ncwrfland.x2aavg_Sf_lfrac.data;
 ncwrfnew.mask.long_name='Ocean Mask in RASM WRF/VIC';
 ncwrfnew.mask.units='';
 ncwrfnew.mask=rmfield(ncwrfnew.mask,'standard_name');

 ncwrfnew.HGT=ncwrfnew.mask;
 ncwrfnew.HGT.data=ncwrfheight.HGT.data';
 ncwrfnew.HGT.long_name='Terrain Height In WRF';
 ncwrfnew.HGT.units='m';

 ncwrfland=ridgepack_shuffle(ncwrfland,{'y' 'x'});
 ncwrfgrid=ridgepack_shuffle(ncwrfgrid,{'y' 'x'});
 ncwrfnew=ridgepack_shuffle(ncwrfnew,{'y' 'x'});

 ncnew.mask_key.data(~isnan(ncnew.mask_key.data))=1;
 ncnew.mask_key.data(1,:)=0;
 ncnew.mask_key.data(end,:)=0;
 ncnew.mask_key.data(:,end)=0;
 ncnew.mask_key.data(:,1)=0;
 nctemp=ridgepack_regrid(ncnew,'mask_key','',ncwrfgrid);
 broadareamask=nctemp.mask_key.data;
 broadareamask(broadareamask<0.95)=0;
 broadareamask(broadareamask>0)=1;
 broadareamask(isnan(broadareamask))=0;

 totalmask1=zeros(size(ncwrfnew.mask.data));
 totalmask2=zeros(size(ncwrfnew.mask.data));

 % interpolate to regions
 for i=includedregions

  regions=char(region{i})
  ncwrfnew.(regions)=ncwrfnew.mask;
  ncwrfnew.(regions).long_name=[num2str(i,'%2.2i'),' ',char(longname{i})];
 
  ncnew.(regions).data(ncnew.mask.data==0)=NaN;
  nctemp=ridgepack_regrid(ncnew,regions,'',ncwrfgrid);
  ncwrfnew.(regions).data=broadareamask.*nctemp.(regions).data;

  if i<18
   totalmask1=totalmask1+ncwrfnew.(regions).data;
  end

  if i>1
   totalmask2=totalmask2+ncwrfnew.(regions).data;
  end

 end

 for k=1,2

  % normalize the mask first with arctic subdivisions
  diff1mask=ncwrfnew.mask.data./totalmask1;
  diff2mask=ncwrfnew.mask.data./totalmask2;

  % interpolate to regions
  totalmask1=zeros(size(ncwrfnew.mask.data));
  totalmask2=zeros(size(ncwrfnew.mask.data));

  for i=includedregions

   regions=char(region{i});

   if i<18
    ncwrfnew.(regions).data=ncwrfnew.(regions).data.*diff1mask;
   else
    ncwrfnew.(regions).data=ncwrfnew.(regions).data.*diff2mask;
   end
   ncwrfnew.(regions).data(isnan(ncwrfnew.(regions).data))=0;
   if k==1; ncwrfnew.(regions).data(ncwrfnew.(regions).data<0.01)=0; end

   clf
   ridgepack_image(ncwrfnew,'x','y',regions,{},{},[0:0.01:1])
%   ridgepack_fprint('png',[num2str(i,'%2.2i'),'_wrf_',regions],1,1); 

   if i<18
    totalmask1=totalmask1+ncwrfnew.(regions).data;
   end

   if i>1
    totalmask2=totalmask2+ncwrfnew.(regions).data;
   end

  end

 end

 ncwrfnew.ocean_fraction=ncwrfnew.mask;
 ncwrfnew.ocean_fraction.long_name='POP/CICE fraction of grid cell';

 ncwrfnew.mask.data(ncwrfnew.mask.data>=0.5)=1;
 ncwrfnew.mask.data(ncwrfnew.mask.data<1)=0;

 ncwrfnew=ridgepack_struct(ncwrfnew);

 ncwrfnew=ridgepack_shuffle(ncwrfnew,{'x' 'y'});

 ridgepack_write(ncwrfnew,'RASM_WRFVIC_GRID_MASKS_AND_METRICS_25km')

end


