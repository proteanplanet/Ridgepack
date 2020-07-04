
clf
clear

largescale=true;
%largescale=false;

gridchoice=2;

fileg{1}.name='CUSP12';
fileg{1}.outname='WC12';
fileg{1}.title=' WC 12-60~km mesh';
fileg{2}.name='CUSP14';
fileg{2}.outname='WC14';
fileg{2}.title=' WC 14-60~km mesh';
fileg{3}.name='DECK';
fileg{3}.outname='DECK';
fileg{3}.title=' DECK 30-60~km standard mesh';

%sector{1}.centlat=90; % degrees north
%sector{1}.centlon=0; % degrees east
%sector{1}.horizon=60; % degrees of satellite horizon (0-90)
%sector{1}.altitude=1; % Mean Earth radius multiple

sector{1}.centlat=75; % degrees north
sector{1}.centlon=-100; % degrees east
sector{1}.horizon=50; % degrees of satellite horizon (0-90)
sector{1}.altitude=1; % Mean Earth radius multiple

sector{2}.centlat=0; % degrees north
sector{2}.centlon=0; % degrees east
sector{2}.horizon=60; % degrees of satellite horizon (0-90)
sector{2}.altitude=1; % Mean Earth radius multiple

sector{3}.centlat=0; % degrees north
sector{3}.centlon=90; % degrees east
sector{3}.horizon=60; % degrees of satellite horizon (0-90)
sector{3}.altitude=1; % Mean Earth radius multiple

sector{4}.centlat=0; % degrees north
sector{4}.centlon=180; % degrees east
sector{4}.horizon=60; % degrees of satellite horizon (0-90)
sector{4}.altitude=1; % Mean Earth radius multiple

sector{5}.centlat=0; % degrees north
sector{5}.centlon=-90; % degrees east
sector{5}.horizon=60; % degrees of satellite horizon (0-90)
sector{5}.altitude=1; % Mean Earth radius multiple

sector{6}.centlat=-90; % degrees north
sector{6}.centlon=0; % degrees east
sector{6}.horizon=60; % degrees of satellite horizon (0-90)
sector{6}.altitude=1; % Mean Earth radius multiple

zoom{1}.centlat=75;
zoom{1}.centlon=-94;
zoom{1}.horizon=10;
zoom{1}.altitude=1;
zoom{1}.name='Canadian Archipelago';

zoom{2}.centlat=62;
zoom{2}.centlon=23;
zoom{2}.horizon=10.5;
zoom{2}.altitude=1;
zoom{2}.name='North Sea, Norwegian Sea, Baltic Sea, and White Sea';

zoom{3}.centlat=74;
zoom{3}.centlon=71;
zoom{3}.horizon=9;
zoom{3}.altitude=1;
zoom{3}.name='Kara Sea and Gulf of Ob''';

zoom{4}.centlat=72;
zoom{4}.centlon=-145;
zoom{4}.horizon=7;
zoom{4}.altitude=1;
zoom{4}.name='Beaufort Sea';

zoom{5}.centlat=64;
zoom{5}.centlon=-172;
zoom{5}.horizon=6;
zoom{5}.altitude=1;
zoom{5}.name='Bering Sea and Chukchi Sea';

zoom{6}.centlat=55;
zoom{6}.centlon=148;
zoom{6}.horizon=11.5;
zoom{6}.altitude=1;
zoom{6}.name='Sea of Okhotsk';

zoom{7}.centlat=75;
zoom{7}.centlon=135;
zoom{7}.horizon=8;
zoom{7}.altitude=1;
zoom{7}.name='Laptev Sea';

zoom{8}.centlat=43.5;
zoom{8}.centlon=135;
zoom{8}.horizon=12.5;
zoom{8}.altitude=1;
zoom{8}.name='Sea of Japan';

zoom{9}.centlat=60.75;
zoom{9}.centlon=-76.5;
zoom{9}.horizon=10;
zoom{9}.altitude=1;
zoom{9}.name='Hudson Bay';

zoom{10}.centlat=56;
zoom{10}.centlon=-145;
zoom{10}.horizon=7;
zoom{10}.altitude=1;
zoom{10}.name='Gulf of Alaska';

zoom{11}.centlat=49;
zoom{11}.centlon=-62;
zoom{11}.horizon=7;
zoom{11}.altitude=1;
zoom{11}.name='Gulf of St.Lawrence';

zoom{12}.centlat=53;
zoom{12}.centlon=-175;
zoom{12}.horizon=12;
zoom{12}.altitude=1;
zoom{12}.name='North Pacific';

zoom{13}.centlat=32;
zoom{13}.centlon=-120;
zoom{13}.horizon=8;
zoom{13}.altitude=1;
zoom{13}.name='California Coast';

zoom{14}.centlat=47;
zoom{14}.centlon=-128;
zoom{14}.horizon=8;
zoom{14}.altitude=1;
zoom{14}.name='North Pacific Coast';

zoom{15}.centlat=65;
zoom{15}.centlon=-40;
zoom{15}.horizon=6.5;
zoom{15}.altitude=1;
zoom{15}.name='South Greenland';

zoom{16}.centlat=76;
zoom{16}.centlon=-40;
zoom{16}.horizon=8;
zoom{16}.altitude=1;
zoom{16}.name='North Greenland';

zoom{17}.centlat=54;
zoom{17}.centlon=-4;
zoom{17}.horizon=5.5;
zoom{17}.altitude=1;
zoom{17}.name='British Isles';

zoom{18}.centlat=40;
zoom{18}.centlon=0;
zoom{18}.horizon=9;
zoom{18}.altitude=1;
zoom{18}.name='West Mediterranean and Bay of Biscay'; 

zoom{19}.centlat=36;
zoom{19}.centlon=22;
zoom{19}.horizon=12.5;
zoom{19}.altitude=1;
zoom{19}.name='East Mediterranean';

zoom{20}.centlat=37.5;
zoom{20}.centlon=-75.5;
zoom{20}.horizon=2.5;
zoom{20}.altitude=1;
zoom{20}.name='Chesapeake Bay and Delaware Bay';

zoom{21}.centlat=34.5;
zoom{21}.centlon=-73;
zoom{21}.horizon=9.8;
zoom{21}.altitude=1;
zoom{21}.name='Atlantic Coast';

zoom{22}.centlat=22;
zoom{22}.centlon=-107;
zoom{22}.horizon=8;
zoom{22}.altitude=1;
zoom{22}.name='Mexican Pacific';

zoom{23}.centlat=24;
zoom{23}.centlon=-90;
zoom{23}.horizon=8.5;
zoom{23}.altitude=1;
zoom{23}.name='Gulf of Mexico';

zoom{24}.centlat=12;
zoom{24}.centlon=47;
zoom{24}.horizon=8;
zoom{24}.altitude=1;
zoom{24}.name='Gulf of Aden';

zoom{25}.centlat=27;
zoom{25}.centlon=55;
zoom{25}.horizon=8;
zoom{25}.altitude=1;
zoom{25}.name='Persian Gulf';

zoom{26}.centlat=9;
zoom{26}.centlon=100;
zoom{26}.horizon=8;
zoom{26}.altitude=1;
zoom{26}.name='Andaman Sea, Strait of Malacca, and Gulf of Thailand';

zoom{27}.centlat=-4;
zoom{27}.centlon=107;
zoom{27}.horizon=7;
zoom{27}.altitude=1;
zoom{27}.name='Java Sea';

zoom{28}.centlat=-5;
zoom{28}.centlon=124;
zoom{28}.horizon=12;
zoom{28}.altitude=1;
zoom{28}.name='Indonesian Throughflow: Celebes Sea, Banda Sea and Timor Sea';

zoom{29}.centlat=-33;
zoom{29}.centlon=28;
zoom{29}.horizon=10;
zoom{29}.altitude=1;
zoom{29}.name='Agulhas Bank and Plateau';

zoom{30}.centlat=-58;
zoom{30}.centlon=-63;
zoom{30}.horizon=8;
zoom{30}.altitude=1;
zoom{30}.name='Drake Passage';

zoom{31}.centlat=-70;
zoom{31}.centlon=-43;
zoom{31}.horizon=10;
zoom{31}.altitude=1;
zoom{31}.name='Weddell Sea';

zoom{32}.centlat=-67;
zoom{32}.centlon=-85;
zoom{32}.horizon=8;
zoom{32}.altitude=1;
zoom{32}.name='Bellingshausen Sea';

zoom{33}.centlat=-70;
zoom{33}.centlon=-120;
zoom{33}.horizon=8;
zoom{33}.altitude=1;
zoom{33}.name='Amundsen Sea';

zoom{34}.centlat=-75;
zoom{34}.centlon=-170;
zoom{34}.horizon=8;
zoom{34}.altitude=1;
zoom{34}.name='Ross Sea';

zoom{35}.centlat=-66;
zoom{35}.centlon=2;
zoom{35}.horizon=8;
zoom{35}.altitude=1;
zoom{35}.name='Maud Rise';

zoom{36}.centlat=-79;
zoom{36}.centlon=94;
zoom{36}.horizon=20;
zoom{36}.altitude=1;
zoom{36}.name='East Antarctic Coast';

zoom{37}.centlat=71;
zoom{37}.centlon=-114;
zoom{37}.horizon=5;
zoom{37}.altitude=1;
zoom{37}.name='Northwest Passsage: Victoria Island';

zoom{38}.centlat=71;
zoom{38}.centlon=-156.5;
zoom{38}.horizon=2;
zoom{38}.altitude=1;
zoom{38}.name='Utqiagvik';

zoom{39}.centlat=20;
zoom{39}.centlon=-75;
zoom{39}.horizon=8;
zoom{39}.altitude=1;
zoom{39}.name='Caribbean Sea';

zoom{40}.centlat=16;
zoom{40}.centlon=115;
zoom{40}.horizon=11;
zoom{40}.altitude=1;
zoom{40}.name='South China Sea';

zoom{41}.centlat=21;
zoom{41}.centlon=-158;
zoom{41}.horizon=4;
zoom{41}.altitude=1;
zoom{41}.name='Hawaii';

zoom{42}.centlat=12;
zoom{42}.centlon=165;
zoom{42}.horizon=5;
zoom{42}.altitude=1;
zoom{42}.name='Marshall Islands';

% grab shipping data
load('/Users/afroberts/data/SHIPPING/total_tracks')

if largescale
 % plotchoice=[1:length(sector)];
 plotchoice=1;
else
 % plotchoice=[1:length(zoom)];
 plotchoice=[1 37];
end


% plot location
plotloc='/Users/afroberts/work';

% grid location
if strcmp(char(fileg{gridchoice}.name),'DECK')
 gridloc=['/Users/afroberts/data/MODEL/E3SM/',...
          char(fileg{gridchoice}.name),'/grid'];
 gridfile='E3SM_LR_V1_grid.nc';
elseif strcmp(char(fileg{gridchoice}.name),'CUSP12')
 gridloc=['/Users/afroberts/data/MODEL/E3SM/CUSP/',...
          char(fileg{gridchoice}.name)];
 gridfile='initial_state.nc';
elseif strcmp(char(fileg{gridchoice}.name),'CUSP14')
 gridloc=['/Users/afroberts/data/MODEL/E3SM/CUSP/',...
          char(fileg{gridchoice}.name)];
 gridfile='initial_state.nc';
end

gridlochr=['/Users/afroberts/data/MODEL/E3SM/highres/grid'];
gridfilehr='E3SM_hr_grid.nc';

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
coastname=[char(fileg{1}.outname),'_Coast'];
x=dir(coastname);
if isempty(x)
 nccoast=ridgepack_e3smcoastm(ncvert);
 ridgepack_write(nccoast,coastname)
else
 nccoast=ridgepack_clone(coastname);
end

% load high-resolution coast
cd(gridlochr)
coastnamehr='E3SM_HR_V1_Coast';
nccoasthr=ridgepack_clone(coastnamehr);

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
  ridgepack_e3smsatmeshs(ncvert,centlat,centlon,horizon,altitude);
  ridgepack_e3smsatcoast(nccoast,centlat,centlon,horizon)

  %zoomidx=[1:length(zoom)];
  zoomidx=[1 37];
  %zoomidx=[1];

  for j=zoomidx
   satlat=zoom{j}.centlat;   % degrees north
   satlon=zoom{j}.centlon;   % degrees east
   sathor=zoom{j}.horizon;   % degrees of satellite horizon (0-90)
   %ridgepack_sathorizon(centlat,centlon,horizon,...
   %                     satlat,satlon,sathor,[0.83 0.5 0],num2str(j));
   ridgepack_sathorizon(centlat,centlon,horizon,...
                        satlat,satlon,sathor,[0.83 0.5 0]);
  end

  %title(['Sector ',num2str(setting),' ',char(fileg{gridchoice}.title)])

  for ij=1:length(shiptrack)
   [x,y,z]=ridgepack_satfwd(shiptrack{ij}.lat,shiptrack{ij}.lon,...
                            centlat,centlon,horizon,1.001*altitude);
   if ij==4
    plot3(x,y,z,':','Color','c')
   else
    plot3(x,y,z,'-','Color','c')
   end
  end

  ridgepack_fprint('png',['InteRFACE_',fileg{gridchoice}.outname,...
                         '_sector_',num2str(setting)],1,2)
 else

  centlat=zoom{setting}.centlat;   % degrees north
  centlon=zoom{setting}.centlon;   % degrees east
  horizon=zoom{setting}.horizon;   % degrees of satellite horizon (0-90)
  altitude=zoom{setting}.altitude; % Mean Earth radius multiple

  ridgepack_satview(centlat,centlon,horizon)
  ridgepack_e3smsatmeshs(ncvert,centlat,centlon,horizon,altitude);
  ridgepack_e3smsatcoast(nccoast,centlat,centlon,horizon)
  h=ridgepack_e3smsatcoast(nccoasthr,centlat,centlon,horizon,[0 0.8 0])

  for ij=1:length(shiptrack)
   [x,y,z]=ridgepack_satfwd(shiptrack{ij}.lat,shiptrack{ij}.lon,...
                            centlat,centlon,horizon,1.001*altitude);
   if ij==4
    plot3(x,y,z,':','Color',[0.83 0.5 0],'Linewidth',1.0)
   else
    plot3(x,y,z,'-','Color',[0.83 0.5 0],'Linewidth',1.0)
   end
  end


% ridgepack_multilegend(h,{'E3SM-HR V1 coastline'},'South')

%  ridgepack_multialign(gcf,...
%             [num2str(setting),': ',...
%              zoom{setting}.name,' ',char(fileg{gridchoice}.title)]);

  ridgepack_fprint('png',['InteRFACE_',fileg{gridchoice}.outname,...
                          '_zoom_',num2str(setting)],1,2)
 end

end

