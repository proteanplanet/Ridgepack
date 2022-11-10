
clf
clear

for col=1:2
for row=1:2


gridchoice=1;

if col==1 & row==1
 plotchoice=1;
 largescale=true;
 bathymetry=false;
 zoomedareas=true;
 coast=true;
 compositemesh=false;
 alphannotation='a';
elseif col==2 & row==1
 plotchoice=2;
 largescale=true;
 bathymetry=false;
 zoomedareas=true;
 coast=true;
 compositemesh=false;
 alphannotation='b';
elseif col==1 & row==2
 plotchoice=1;
 largescale=false;
 bathymetry=false;
 zoomedareas=false;
 coast=true;
 compositemesh=false;
 alphannotation='c';
elseif col==2 & row==2
 plotchoice=2;
 largescale=false;
 bathymetry=false;
 zoomedareas=false;
 coast=true;
 compositemesh=false;
 alphannotation='d';
else
 error('Plotchoice does not exist')
end

ridgepack_multiplot(2,2,row,col,alphannotation)

fileg{1}.name='WC14r03';
fileg{1}.outname='WC14r03';
fileg{1}.title='WC 14-60~km mesh';

fileg{2}.name='EC30to60E2r2';
fileg{2}.outname='EC30to60E2r2';
fileg{2}.title='EC 30-60km E2 r2';

sector{1}.centlat=60; % degrees north
sector{1}.centlon=-90; % degrees east
sector{1}.horizon=60; % degrees of satellite horizon (0-90)
sector{1}.altitude=1; % Mean Earth radius multiple
sector{1}.annotation=1; % no annotation

sector{2}.centlat=-40; % degrees north
sector{2}.centlon=-90; % degrees east
sector{2}.horizon=60; % degrees of satellite horizon (0-90)
sector{2}.altitude=1; % Mean Earth radius multiple
sector{2}.annotation=0; % no annotation

zoom{1}.centlat=75;
zoom{1}.centlon=-94;
zoom{1}.horizon=7;
zoom{1}.altitude=1;
zoom{1}.name='Canadian Archipelago';
zoom{1}.annotation=1; % Arctic Ship Tracks
zoom{1}.polar=1; % plot 20m bathygraph

zoom{2}.centlat=28;
zoom{2}.centlon=-72;
zoom{2}.horizon=9;
zoom{2}.altitude=1;
zoom{2}.name='Atlantic Coast';
zoom{2}.annotation=0; % no annotation
zoom{2}.polar=0;

% plot location
plotloc='/Users/afroberts/Science/publications/2022_E3SM_V2_NARRM';

% grid location
if strcmp(char(fileg{gridchoice}.name),'WC14r03')
 gridloc=['/Users/afroberts/data/MODEL/E3SM/v2/v2.NARRM.piControl/grid'];
 gridfile='mpaso.rst.0002-01-01_00000.nc';
 %shiplocs=[1 2 3 4 6 20];
 shiplocs=[];
elseif strcmp(char(fileg{gridchoice}.name),'EC30to60E2r2')
 gridloc=['/Users/afroberts/data/MODEL/E3SM/v2/v2.LR.piControl/grid'];
 gridfile='mpaso.rst.0002-01-01_00000.nc';
 %shiplocs=[1 2 3 4 6 20];
 shiplocs=[];
else
 error('No grid found')
end

% obtain grid information
cd(gridloc)
ncvert=ridgepack_clone(gridfile,{'latVertex','lonVertex',...
                                'verticesOnCell','indexToCellID',...
                                'nEdgesOnCell','edgesOnCell',...
                                'cellsOnEdge','cellsOnVertex',...
                                'edgesOnVertex','bottomDepth'});

nccell=ridgepack_clone(gridfile,{'latCell','lonCell',...
                                'verticesOnCell','indexToCellID',...
                                'nEdgesOnCell','edgesOnCell',...
                                'cellsOnEdge','cellsOnVertex',...
                                'edgesOnVertex','bottomDepth'});

ncedge=ridgepack_clone(gridfile,{'latEdge','lonEdge'});

ncvert.latCell=nccell.latitude;
ncvert.lonCell=nccell.longitude;
ncvert.latEdge=ncedge.latitude;
ncvert.lonEdge=ncedge.longitude;
ncvert.latVertex=ncvert.latitude;
ncvert.lonVertex=ncvert.longitude;

nccell.latVertex=ncvert.latitude;
nccell.lonVertex=ncvert.longitude;
nccell.latEdge=ncedge.latitude;
nccell.lonEdge=ncedge.longitude;
nccell.latCell=nccell.latitude;
nccell.lonCell=nccell.longitude;

% read in coast or else generate coast and write it out
if coast
 coastname=[char(fileg{gridchoice}.outname),'_Coast.nc'];
 x=dir(coastname);
 if isempty(x)
  nccoast=ridgepack_e3smcoastm(ncvert);
  ridgepack_write(nccoast,coastname)
 else
  nccoast=ridgepack_clone(coastname);
 end
end

% create 20m isobath
if bathymetry
 bathname=[char(fileg{gridchoice}.outname),'_20mIsobath.nc'];
 x=dir(bathname);
 if isempty(x)
  ncisobath20=ridgepack_e3smseasaw(ncvert,ncvert,'bottomDepth',100);
  ridgepack_write(ncisobath20,bathname)
 else
  ncisobath20=ridgepack_clone(bathname);
 end
end

% invert bathymetry
ncvert.bottomDepth.data=-ncvert.bottomDepth.data;

% load low-resolution coast
gridloclr='/Users/afroberts/data/MODEL/E3SM/v2/v2.LR.piControl/grid';
cd(gridloclr)
coastnamelr='EC30to60E2r2_Coast';
nccoastlr=ridgepack_clone(coastnamelr);

% load in shipping data
shiploc='/Users/afroberts/data/SHIPPING';
cd(shiploc)
ncship1=ridgepack_clone('InteRFACE_Shiptracks',...
               {'track01_longitude','track01_latitude'});
ncship2=ridgepack_clone('InteRFACE_Shiptracks',...
               {'track02_longitude','track02_latitude'});
ncship3=ridgepack_clone('InteRFACE_Shiptracks',...
               {'track03_longitude','track03_latitude'});
ncship4=ridgepack_clone('InteRFACE_Shiptracks',...
               {'track04_longitude','track04_latitude'});
ncship6=ridgepack_clone('InteRFACE_Shiptracks',...
               {'track06_longitude','track06_latitude'});
ncship20=ridgepack_clone('InteRFACE_Shiptracks',...
               {'track20_longitude','track20_latitude'});

% move to plot location
%mkdir(plotloc)
cd(plotloc)

for setting=plotchoice

 if largescale

  centlat=sector{setting}.centlat;   % degrees north
  centlon=sector{setting}.centlon;   % degrees east
  horizon=sector{setting}.horizon;   % degrees of satellite horizon (0-90)
  altitude=sector{setting}.altitude; % Mean Earth radius multiple

  ridgepack_satview(centlat,centlon,horizon)

  if bathymetry

   % reverse colorbar
   cont=[-5000:500:-1500 -1000:250:-250 -100 -50:10:10]; 
   colbarcont{1}='\downarrow';
   for i=2:length(cont)-1;
    colbarcont{i}=num2str(-cont(i));
   end
   colbarcont{length(cont)}='\uparrow';

   % render colors
   ridgepack_e3smsatcol(ncvert,'bottomDepth',ncvert,cont,0,...
                        centlat,centlon,horizon,altitude,...
                        true,false,'linear','bluered');

   % make top color grey
   cmap=colormap;
   cmap(end,:)=0.95*[1 1 1];
   colormap(cmap)

   % add colormap
   ridgepack_colorbar(cont,'m','linear','vertical',0,colbarcont)
   clear colbarcont

   % add 20m isobath
   if bathymetry
    hb=ridgepack_e3smsatthreshold(ncisobath20,centlat,centlon,horizon,...
                                 [0.9290 0.6940 0.1250]);
   end

   % plot coast
   if coast
    ridgepack_e3smsatcoast(nccoast,centlat,centlon,horizon)
   end

   if sector{setting}.annotation==1 
    for shipi=shiplocs
     [x,y,z,phi,theta]=...
      ridgepack_satfwd(eval(['ncship',num2str(shipi),'.latitude.data']),...
                       eval(['ncship',num2str(shipi),'.longitude.data']),...
                       centlat,centlon,horizon,1.001*altitude);
     hship=plot3(x,y,z,'r-');
    end
   end

  else

   ridgepack_e3smsatmeshs(ncvert,centlat,centlon,horizon,altitude);

   if sector{setting}.annotation==1
    for shipi=shiplocs
     [x,y,z,phi,theta]=...
      ridgepack_satfwd(eval(['ncship',num2str(shipi),'.latitude.data']),...
                       eval(['ncship',num2str(shipi),'.longitude.data']),...
                       centlat,centlon,horizon,1.001*altitude);
     hship=plot3(x,y,z,'r-');
    end
   end
  
   if coast
    ridgepack_e3smsatcoast(nccoast,centlat,centlon,horizon)
   end


   alpha='cd';
   for j=1:length(zoom) 
    satlat=zoom{j}.centlat;   % degrees north
    satlon=zoom{j}.centlon;   % degrees east
    sathor=zoom{j}.horizon;   % degrees of satellite horizon (0-90)
    ridgepack_sathorizon(centlat,centlon,horizon,...
                  satlat,satlon,sathor,[1 0 0],alpha(j));
   end
 
  end

 else

  centlat=zoom{setting}.centlat;   % degrees north
  centlon=zoom{setting}.centlon;   % degrees east
  horizon=zoom{setting}.horizon;   % degrees of satellite horizon (0-90)
  altitude=zoom{setting}.altitude; % Mean Earth radius multiple

  ridgepack_satview(centlat,centlon,horizon)
  ridgepack_e3smsatmeshs(ncvert,centlat,centlon,horizon,altitude);

  % add 20m isobath
  if bathymetry
   hk=ridgepack_e3smsatthreshold(ncisobath20,centlat,centlon,horizon,[0 0 1]);
   set(hk,'LineWidth',0.5)
  end

  % add in coastline
  if coast
   hc=ridgepack_e3smsatcoast(nccoast,centlat,centlon,horizon)
  
   % add in low resolution coastline
   hl=ridgepack_e3smsatcoast(nccoastlr,centlat,centlon,horizon,[0.83 0.5 0])
  end

  if bathymetry

   ridgepack_satview(centlat,centlon,horizon,1,1)

   % reverse colorbar
   cont=[-5000:1000:-2000 -1500 -1000:250:-250 -100 -50:10:10];
   colbarcont{1}='\downarrow';
   for i=2:length(cont)-1;
    colbarcont{i}=num2str(-cont(i));
   end
   colbarcont{length(cont)}='\uparrow';
   cmap=colormap;

   % render colors
   ridgepack_e3smsatcol(ncvert,'bottomDepth',ncvert,cont,0,...
                        centlat,centlon,horizon,altitude,...
                        true,false,'linear','bluered');

   % make top color grey
   cmap=colormap;
   cmap(end,:)=0.95*[1 1 1];
   colormap(cmap)

   % add colormap
   ridgepack_colorbar(cont,'m','linear','vertical',0,colbarcont)
   clear cmap colbarcont

   % add 20m isobath
   hb=ridgepack_e3smsatthreshold(ncisobath20,centlat,centlon,horizon,...
                                 [0.9290 0.6940 0.1250]);

  elseif compositemesh

   ridgepack_satview(centlat,centlon,horizon)
   ridgepack_e3smsatmeshv(nccell,centlat,centlon,horizon,altitude);

   ridgepack_e3smsatcoast(nccoast,centlat,centlon,horizon)

   if zoom{setting}.annotation==1 
    for shipi=[1 2 3 4 6 20]
     [x,y,z,phi,theta]=...
      ridgepack_satfwd(eval(['ncship',num2str(shipi),'.latitude.data']),...
                       eval(['ncship',num2str(shipi),'.longitude.data']),...
                       centlat,centlon,horizon,1.001*altitude);
     hship=plot3(x,y,z,'-','Color','k');
    end
   end

  end

  ridgepack_multilegend([hc hl],{'NARRM Coastline','Standard E3SM Coastline'},'South')

 end

end

end
end

%ridgepack_multilegend([hship],{'Arctic Coast Shipping Channels'},'South') 

ridgepack_multialign(gcf)

ridgepack_fprint('png','E3SM_v2_mesh_NARRM',1,2)

