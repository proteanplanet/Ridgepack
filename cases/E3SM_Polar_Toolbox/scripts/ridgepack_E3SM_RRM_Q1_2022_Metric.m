
clf
clear

%largescale=true;
largescale=false;

%bathymetry=true;
bathymetry=false;

zoomedareas=true;
%zoomedareas=false;

%coast=false;
coast=true;

gridchoice=4;

fileg{1}.name='WC12r01';
fileg{1}.outname='WC12';
fileg{1}.title=' WC 12-60~km mesh';

fileg{2}.name='WC14r03';
fileg{2}.outname='WC14r03';
fileg{2}.title=' WC 14-60~km mesh';

fileg{3}.name='DECK';
fileg{3}.outname='DECK';
fileg{3}.title=' DECK 30-60~km standard mesh';

fileg{4}.name='ECwISC30to60E1r02';
fileg{4}.outname='ECwISC30to60E1r02';
fileg{4}.title='EC wISC 30-60km E1 r02';

fileg{5}.name='EC30to60E2r2';
fileg{5}.outname='EC30to60E2r2';
fileg{5}.title='EC 30-60km E2 r2';

fileg{6}.name='EC15to60E2r4';
fileg{6}.outname='EC15to60E2r4';
fileg{6}.title='EC 15-60km E2 r4';

fileg{7}.name='SOwISC12to60E2r4';
fileg{7}.outname='SOwISC12to60E2r4';
fileg{7}.title='Southern Ocean w/Ice Shelves 12-60km E2 r4';

fileg{8}.name='icosahedron7';
fileg{8}.outname='icosahedron7';
fileg{8}.title='icosahedron7 tides mesh';

fileg{9}.name='WC14';
fileg{9}.outname='WC14';
fileg{9}.title='WC14 Mesh';

sector{1}.centlat=90; % degrees north
sector{1}.centlon=0; % degrees east
sector{1}.horizon=60; % degrees of satellite horizon (0-90)
sector{1}.altitude=1; % Mean Earth radius multiple
sector{1}.annotation=1; % add Arctic Shipping

sector{2}.centlat=-60; % degrees north
sector{2}.centlon=-20; % degrees east
sector{2}.horizon=60; % degrees of satellite horizon (0-90)
sector{2}.altitude=1; % Mean Earth radius multiple
sector{2}.annotation=0; % no annotation

sector{3}.centlat=-90; % degrees north
sector{3}.centlon=0; % degrees east
sector{3}.horizon=60; % degrees of satellite horizon (0-90)
sector{3}.altitude=1; % Mean Earth radius multiple
sector{3}.annotation=0; % no annotation

sector{4}.centlat=60; % degrees north
sector{4}.centlon=-90; % degrees east
sector{4}.horizon=60; % degrees of satellite horizon (0-90)
sector{4}.altitude=1; % Mean Earth radius multiple
sector{4}.annotation=1; % no annotation

zoom{1}.centlat=75;
zoom{1}.centlon=-94;
zoom{1}.horizon=10;
zoom{1}.altitude=1;
zoom{1}.name='Canadian Archipelago';
zoom{1}.annotation=1; % Arctic Ship Tracks
zoom{1}.polar=1; % plot 20m bathygraph

zoom{2}.centlat=-70;
zoom{2}.centlon=-43;
zoom{2}.horizon=15;
zoom{2}.altitude=1;
zoom{2}.name='Weddell Sea';
zoom{2}.annotation=0; % no annotation
zoom{2}.polar=1; 

if largescale
 plotchoice=[1:length(sector)];
 %plotchoice=7;
else
 plotchoice=[1:length(zoom)];
 %plotchoice=[1 3 4 5 7 37 38];
 %plotchoice=[21];
end

% plot location
plotloc='/Users/afroberts/Science/publications/2022_Q1_RRM';

% grid location
if strcmp(char(fileg{gridchoice}.name),'DECK')
 gridloc=['/Users/afroberts/data/MODEL/E3SM/EC_60_30_Old/grid'];
 gridfile='init.nc';
 shiplocs=[2 6];
elseif strcmp(char(fileg{gridchoice}.name),'WC12')
 gridloc=['/Users/afroberts/data/MODEL/E3SM/WC12/',...
          char(fileg{gridchoice}.name)];
 gridfile='initial_state.nc';
 shiplocs=[1 2 3 4 6 20];
elseif strcmp(char(fileg{gridchoice}.name),'WC14')
 gridloc=['/Users/afroberts/data/MODEL/E3SM/20210527I/grid'];
 gridfile='InteRFACE1alphaA.mpaso.rst.0002-01-01_00000.nc';
 shiplocs=[1 2 3 4 6 20];
elseif strcmp(char(fileg{gridchoice}.name),'WC14r03')
 gridloc=['/Users/afroberts/data/MODEL/E3SM/WC14/r03'];
 gridfile='initial_state.nc';
 shiplocs=[1 2 3 4 6 20];
elseif strcmp(char(fileg{gridchoice}.name),'ECwISC30to60E1r02')
 gridloc=['/Users/afroberts/data/MODEL/E3SM/ECwISC30to60E1r02'];
 gridfile='ocean.ECwISC30to60E1r02.200408.nc';
 shiplocs=[1 2 3 4 6 20];
elseif strcmp(char(fileg{gridchoice}.name),'EC30to60E2r2')
 gridloc=['/Users/afroberts/data/MODEL/E3SM/PIControlSI/grid'];
 gridfile='mpaso.rst.0002-01-01_00000.nc';
 shiplocs=[1 2 3 4 6 20];
elseif strcmp(char(fileg{gridchoice}.name),'EC15to60E2r4')
 gridloc=['/Users/afroberts/data/MODEL/E3SM/EC15to60E2r4/grid'];
 gridfile='initial_state.nc';
 shiplocs=[1 2 3 4 6 20];
elseif strcmp(char(fileg{gridchoice}.name),'SOwISC12to60E2r4')
 gridloc=['/Users/afroberts/data/MODEL/E3SM/SOwISC12to60E2r4/grid'];
 gridfile='initial_state.nc';
 shiplocs=[1 2 3 4 6 20];
elseif strcmp(char(fileg{gridchoice}.name),'icosahedron7')
 gridloc=['/Users/afroberts/data/MODEL/E3SM/icosahedron_7/grid'];
 gridfile='initial_state.nc';
 shiplocs=[1 2 3 4 6 20];
end

%gridlochr=['/Users/afroberts/data/MODEL/E3SM/highres/grid'];
%gridfilehr='E3SM_hr_grid.nc';

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

% load high-resolution coast
%cd(gridlochr)
%coastnamehr='E3SM_HR_V1_Coast';
%nccoasthr=ridgepack_clone(coastnamehr);

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
cd(plotloc)

for setting=plotchoice

 clf

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

   if zoomedareas;

    for j=1:length(zoom)
     satlat=zoom{j}.centlat;   % degrees north
     satlon=zoom{j}.centlon;   % degrees east
     sathor=zoom{j}.horizon;   % degrees of satellite horizon (0-90)
     ridgepack_sathorizon(centlat,centlon,horizon,...
                          satlat,satlon,sathor,[0.83 0.5 0],num2str(j));
    end

   end
 
  end

  title(['Sector ',num2str(setting),' ',char(fileg{gridchoice}.title)])

  if bathymetry
   ridgepack_fprint('png',[fileg{gridchoice}.outname,...
                          '_sector_bathymetry_',num2str(setting)],1,2)
  else
   ridgepack_fprint('png',[fileg{gridchoice}.outname,...
                          '_sector_mesh_',num2str(setting)],1,2)
  end

 else

  centlat=zoom{setting}.centlat;   % degrees north
  centlon=zoom{setting}.centlon;   % degrees east
  horizon=zoom{setting}.horizon;   % degrees of satellite horizon (0-90)
  altitude=zoom{setting}.altitude; % Mean Earth radius multiple

  ridgepack_multiplot(1,2,1,1)
  ridgepack_satview(centlat,centlon,horizon)
  ridgepack_e3smsatmeshs(ncvert,centlat,centlon,horizon,altitude);

  % add 20m isobath
  if bathymetry
   hk=ridgepack_e3smsatthreshold(ncisobath20,centlat,centlon,horizon,...
                                 [0 0 1]);
   set(hk,'LineWidth',0.5)
  end

  % add in coastline
  %if coast
  % ridgepack_e3smsatcoast(nccoast,centlat,centlon,horizon)
  %
  % % add in high resolution coastline
  % h=ridgepack_e3smsatcoast(nccoasthr,centlat,centlon,horizon,[0 0.8 0])
  %end

  ridgepack_multiplot(1,2,1,2)

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

   % add coast
   if coast
    ridgepack_e3smsatcoast(nccoast,centlat,centlon,horizon)
    ridgepack_multilegend([h hb],...
         {'E3SM-HR V1 coastline','20 m isobath'},'South')
   end

  else

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
    ridgepack_multilegend([hship],{'Arctic Coast Shipping Channels'},'South') 
   end

  end

  ridgepack_multialign(gcf,...
             [num2str(setting),': ',...
              zoom{setting}.name,' ',char(fileg{gridchoice}.title)]);

  if bathymetry
   ridgepack_fprint('png',[fileg{gridchoice}.outname,...
                          '_zoom_bathymetry_',num2str(setting)],1,2)
  else
   ridgepack_fprint('png',[fileg{gridchoice}.outname,...
                          '_zoom_mesh_',num2str(setting)],1,2)
  end

 end

end

