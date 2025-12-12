
close all
clear

% Configs:
% 1 - E3SM V3 Coupled Paper Marine Mesh (original)
% 2 - E3SM v3 Coupled Paper Northern Hemisphere (revision)
% 3 - E3SM v3 Coupled Paper Southern Hemisphere (revision)
% 4 - Trial Atmospheric/Land Figure

config=3;

% title
%maintitle='MPAS E3SM V3 Mesh';
%if confg==1; maintitle='';

% configurations, order left to right, top to bottom
if config==1

 % number of cols and rows
 rows=2;
 cols=2;

 % grids to plot
 grids=[4 4 4 4];

 % projections 
 projections=[65 65 67 67];

 % large scale projection, or not
 largescales={true,true,true,true};
 %largescales={true,true};

 % if large scale include zoomed areas
 zoomedareas={false,false,true,false};
 %zoomedareas={false,false}

 % add bathymetry 
 bathymetrys={false,true,false,true};
 %bathymetrys={false,true};

 % use rossby metric
 rossbymetrics={false,false,false,false};
 %rossbymetrics={false,false};

 % apply resolution shading
 resolutionmetrics={true,false,true,false};
 %resolutionmetrics={true,false};

 % apply mesh distortion shading
 distortions={false,false,false,false};
 %distortions={false,false};

 % overlay mesh 
 overlaymeshs={false,false,false,false};
 %overlaymeshs={false,false};

elseif config==2

 % number of cols and rows
 rows=2;
 cols=1;

 % grids to plot
 grids=[4 4];

 % projections 
 projections=[65 65];

 % large scale projection, or not
 largescales={true,true};

 % if large scale include zoomed areas
 zoomedareas={false,false};

 % add bathymetry 
 bathymetrys={false,true};

 % use rossby metric
 rossbymetrics={false,false};

 % apply resolution shading
 resolutionmetrics={true,false};

 % apply mesh distortion shading
 distortions={false,false};

 % overlay mesh 
 overlaymeshs={false,false};

elseif config==3

 % number of cols and rows
 rows=2;
 cols=1;

 % grids to plot
 grids=[4 4];

 % projections 
 projections=[67 67];

 % large scale projection, or not
 largescales={true,true};

 % if large scale include zoomed areas
 zoomedareas={false,false};

 % add bathymetry 
 bathymetrys={false,true};

 % use rossby metric
 rossbymetrics={false,false};

 % apply resolution shading
 resolutionmetrics={true,false};

 % apply mesh distortion shading
 distortions={false,false};

 % overlay mesh    
 overlaymeshs={false,false};

end

% plot location
plotloc='/Users/afroberts/work';

% add in meshtitle
meshtitle=false;

% outfile names
outfilename=['_config',num2str(config)];
grifilename='test';

% set legend flags
%shiplegflag=true;
shiplegflag=false;
iceshelflegflag=true;
bathlegflag=true;
legendh=[];
legendcell={};
drawcb=true;

% labeling
alphass='abcd';

% set up ocean mesh grid file choices
fileo{1}.name='oQU480';
fileo{1}.outname='oQU480';
fileo{1}.title='MPAS-SI oQU480 Column Test Grid';
  
fileo{2}.name='SOwISC12to30E3r3';
fileo{2}.outname='SOwISC12to30E3r3';
fileo{2}.title='MPAS E3SM V3 SORRM R3';
  
fileo{3}.name='RRSwISC6to18E3r5';
fileo{3}.outname='RRSwISC6to18E3r5';
fileo{3}.title='MPAS E3SM V3 High Resolution R5';
  
fileo{4}.name='IcoswISC30E3r5';
fileo{4}.outname='IcoswISC30E3r5';
fileo{4}.title='MPAS E3SM V3 30km Equal Area Mesh r5';
  
% define projection
zoom{1}.centlat=75;
zoom{1}.centlat=75;
zoom{1}.centlat=75;
zoom{1}.centlon=-94;
zoom{1}.horizon=10;
zoom{1}.altitude=1;
zoom{1}.name='Canadian Archipelago';
zoom{1}.annotation=1; % Arctic Ship Tracks
zoom{1}.polar=1; % plot min bathygraph
  
zoom{2}.centlat=62;
zoom{2}.centlon=23;
zoom{2}.horizon=10.5;
zoom{2}.altitude=1;
zoom{2}.name='North Sea, Norwegian Sea, Baltic Sea, and White Sea';
zoom{2}.annotation=0; % no annotation
zoom{2}.polar=0; 
  
zoom{3}.centlat=74;
zoom{3}.centlon=71;
zoom{3}.horizon=9;
zoom{3}.altitude=1;
zoom{3}.name='Kara Sea and Gulf of Ob''';
zoom{3}.annotation=1; % Arctic Ship Tracks
zoom{3}.polar=1; 
  
zoom{4}.centlat=72;
zoom{4}.centlon=-145;
zoom{4}.horizon=7;
zoom{4}.altitude=1;
zoom{4}.name='Beaufort Sea';
zoom{4}.annotation=1; % Arctic Ship Tracks
zoom{4}.polar=1; 
 
zoom{5}.centlat=64;
zoom{5}.centlon=-172;
zoom{5}.horizon=6;
zoom{5}.altitude=1;
zoom{5}.name='Bering Sea and Chukchi Sea';
zoom{5}.annotation=1; % Arctic Ship Tracks
zoom{5}.polar=1; 
  
zoom{6}.centlat=55;
zoom{6}.centlon=148;
zoom{6}.horizon=11.5;
zoom{6}.altitude=1;
zoom{6}.name='Sea of Okhotsk';
zoom{6}.annotation=0; % no annotation
zoom{6}.polar=1; 
  
zoom{7}.centlat=75;
zoom{7}.centlon=135;
zoom{7}.horizon=8;
zoom{7}.altitude=1;
zoom{7}.name='Laptev Sea';
zoom{7}.annotation=1; % Arctic Ship Tracks
zoom{7}.polar=1; 
  
zoom{8}.centlat=43.5;
zoom{8}.centlon=135;
zoom{8}.horizon=12.5;
zoom{8}.altitude=1;
zoom{8}.name='Sea of Japan';
zoom{8}.annotation=0; % no annotation
zoom{8}.polar=0; 
  
zoom{9}.centlat=60.75;
zoom{9}.centlon=-76.5;
zoom{9}.horizon=10;
zoom{9}.altitude=1;
zoom{9}.name='Hudson Bay';
zoom{9}.annotation=0; % no annotation
zoom{9}.polar=1; 
  
zoom{10}.centlat=56;
zoom{10}.centlon=-150;
zoom{10}.horizon=7;
zoom{10}.altitude=1;
zoom{10}.name='Gulf of Alaska';
zoom{10}.annotation=0; % no annotation
zoom{10}.polar=0; 
  
zoom{11}.centlat=49;
zoom{11}.centlon=-62;
zoom{11}.horizon=7;
zoom{11}.altitude=1;
zoom{11}.name='Gulf of St.Lawrence';
zoom{11}.annotation=0; % no annotation
zoom{11}.polar=0; 
  
zoom{12}.centlat=53;
zoom{12}.centlon=-175;
zoom{12}.horizon=12;
zoom{12}.altitude=1;
zoom{12}.name='North Pacific';
zoom{12}.annotation=0; % no annotation
zoom{12}.polar=0; 
  
zoom{13}.centlat=32;
zoom{13}.centlon=-120;
zoom{13}.horizon=8;
zoom{13}.altitude=1;
zoom{13}.name='California Coast';
zoom{13}.annotation=0; % no annotation
zoom{13}.polar=0; 
  
zoom{14}.centlat=47;
zoom{14}.centlon=-128;
zoom{14}.horizon=8;
zoom{14}.altitude=1;
zoom{14}.name='North Pacific Coast';
zoom{14}.annotation=0; % no annotation
zoom{14}.polar=0; 
 
zoom{15}.centlat=65;
zoom{15}.centlon=-40;
zoom{15}.horizon=6.5;
zoom{15}.altitude=1;
zoom{15}.name='South Greenland';
zoom{15}.annotation=0; % no annotation
zoom{15}.polar=0; 
 
zoom{16}.centlat=76;
zoom{16}.centlon=-40;
zoom{16}.horizon=8;
zoom{16}.altitude=1;
zoom{16}.name='North Greenland';
zoom{16}.annotation=0; % no annotation
zoom{16}.polar=1; 
 
zoom{17}.centlat=54;
zoom{17}.centlon=-4;
zoom{17}.horizon=5.5;
zoom{17}.altitude=1;
zoom{17}.name='British Isles';
zoom{17}.annotation=0; % no annotation
zoom{17}.polar=0; 
 
zoom{18}.centlat=40;
zoom{18}.centlon=0;
zoom{18}.horizon=9;
zoom{18}.altitude=1;
zoom{18}.name='West Mediterranean and Bay of Biscay'; 
zoom{18}.annotation=0; % no annotation
zoom{18}.polar=0; 
 
zoom{19}.centlat=36;
zoom{19}.centlon=22;
zoom{19}.horizon=12.5;
zoom{19}.altitude=1;
zoom{19}.name='East Mediterranean';
zoom{19}.annotation=0; % no annotation
zoom{19}.polar=0; 
 
zoom{20}.centlat=37.5;
zoom{20}.centlon=-75.5;
zoom{20}.horizon=2.5;
zoom{20}.altitude=1;
zoom{20}.name='Chesapeake Bay and Delaware Bay';
zoom{20}.annotation=0; % no annotation
zoom{20}.polar=0; 
 
zoom{21}.centlat=34.5;
zoom{21}.centlon=-73;
zoom{21}.horizon=9.8;
zoom{21}.altitude=1;
zoom{21}.name='Atlantic Coast';
zoom{21}.annotation=0; % no annotation
zoom{21}.polar=0; 
 
zoom{22}.centlat=22;
zoom{22}.centlon=-107;
zoom{22}.horizon=8;
zoom{22}.altitude=1;
zoom{22}.name='Mexican Pacific';
zoom{22}.annotation=0; % no annotation
zoom{22}.polar=0; 
 
zoom{23}.centlat=24;
zoom{23}.centlon=-90;
zoom{23}.horizon=8.5;
zoom{23}.altitude=1;
zoom{23}.name='Gulf of Mexico';
zoom{23}.annotation=0; % no annotation
zoom{23}.polar=0; 
 
zoom{24}.centlat=12;
zoom{24}.centlon=47;
zoom{24}.horizon=8;
zoom{24}.altitude=1;
zoom{24}.name='Gulf of Aden';
zoom{24}.annotation=0; % no annotation
zoom{24}.polar=0; 
 
zoom{25}.centlat=27;
zoom{25}.centlon=55;
zoom{25}.horizon=8;
zoom{25}.altitude=1;
zoom{25}.name='Persian Gulf';
zoom{25}.annotation=0; % no annotation
zoom{25}.polar=0; 
 
zoom{26}.centlat=9;
zoom{26}.centlon=100;
zoom{26}.horizon=8;
zoom{26}.altitude=1;
zoom{26}.name='Andaman Sea, Strait of Malacca, and Gulf of Thailand';
zoom{26}.annotation=0; % no annotation
zoom{26}.polar=0; 
 
zoom{27}.centlat=-4;
zoom{27}.centlon=107;
zoom{27}.horizon=7;
zoom{27}.altitude=1;
zoom{27}.name='Java Sea';
zoom{27}.annotation=0; % no annotation
zoom{27}.polar=0; 
 
zoom{28}.centlat=-5;
zoom{28}.centlon=124;
zoom{28}.horizon=12;
zoom{28}.altitude=1;
zoom{28}.name='Indonesian Throughflow: Celebes Sea, Banda Sea and Timor Sea';
zoom{28}.annotation=0; % no annotation
zoom{28}.polar=0; 
 
zoom{29}.centlat=-33;
zoom{29}.centlon=28;
zoom{29}.horizon=10;
zoom{29}.altitude=1;
zoom{29}.name='Agulhas Bank and Plateau';
zoom{29}.annotation=0; % no annotation
zoom{29}.polar=0; 
 
zoom{30}.centlat=-58;
zoom{30}.centlon=-63;
zoom{30}.horizon=8;
zoom{30}.altitude=1;
zoom{30}.name='Drake Passage';
zoom{30}.annotation=0; % no annotation
zoom{30}.polar=0; 
 
zoom{31}.centlat=-72;
zoom{31}.centlon=-43;
zoom{31}.horizon=12.4;
zoom{31}.altitude=1;
zoom{31}.name='Weddell Sea';
zoom{31}.annotation=2; % ice shelf annotation
zoom{31}.polar=1; 
 
zoom{32}.centlat=-69;
zoom{32}.centlon=-85;
zoom{32}.horizon=8;
zoom{32}.altitude=1;
zoom{32}.name='Bellingshausen Sea';
zoom{32}.annotation=2; % ice shelf annotation
zoom{32}.polar=1; 
 
zoom{33}.centlat=-73;
zoom{33}.centlon=-120;
zoom{33}.horizon=8;
zoom{33}.altitude=1;
zoom{33}.name='Amundsen Sea';
zoom{33}.annotation=2; % ice shelf annotation
zoom{33}.polar=1; 
 
zoom{34}.centlat=-78;
zoom{34}.centlon=-170;
zoom{34}.horizon=8;
zoom{34}.altitude=1;
zoom{34}.name='Ross Sea';
zoom{34}.annotation=2; % ice shelf annotation
zoom{34}.polar=1; 
 
zoom{35}.centlat=-66;
zoom{35}.centlon=2;
zoom{35}.horizon=8;
zoom{35}.altitude=1;
zoom{35}.name='Maud Rise';
zoom{35}.annotation=2; % ice shelf annotation
zoom{35}.polar=1; 
 
zoom{36}.centlat=-88;
zoom{36}.centlon=0;
zoom{36}.horizon=27;
zoom{36}.altitude=1;
zoom{36}.name='Antarctica';
zoom{36}.annotation=2; % ice shelf annotation
zoom{36}.polar=1; 
 
zoom{37}.centlat=71;
zoom{37}.centlon=-114;
zoom{37}.horizon=5;
zoom{37}.altitude=1;
zoom{37}.name='Northwest Passsage: Victoria Island';
zoom{37}.annotation=1; % Arctic Ship Tracks
zoom{37}.polar=1; 
 
zoom{38}.centlat=71;
zoom{38}.centlon=-156.5;
zoom{38}.horizon=2;
zoom{38}.altitude=1;
zoom{38}.name='Utqiagvik';
zoom{38}.annotation=1; % Arctic Ship Tracks
zoom{38}.polar=1; 
 
zoom{39}.centlat=20;
zoom{39}.centlon=-75;
zoom{39}.horizon=8;
zoom{39}.altitude=1;
zoom{39}.name='Caribbean Sea';
zoom{39}.annotation=0; % no annotation
zoom{39}.polar=0; 
 
zoom{40}.centlat=16;
zoom{40}.centlon=115;
zoom{40}.horizon=11;
zoom{40}.altitude=1;
zoom{40}.name='South China Sea';
zoom{40}.annotation=0; % no annotation
zoom{40}.polar=0; 
 
zoom{41}.centlat=21;
zoom{41}.centlon=-158;
zoom{41}.horizon=4;
zoom{41}.altitude=1;
zoom{41}.name='Hawaii';
zoom{41}.annotation=0; % no annotation
zoom{41}.polar=0; 
 
zoom{42}.centlat=12;
zoom{42}.centlon=165;
zoom{42}.horizon=5;
zoom{42}.altitude=1;
zoom{42}.name='Marshall Islands';
zoom{42}.annotation=0; % no annotation
zoom{42}.polar=0; 
 
zoom{43}.centlat=-41;
zoom{43}.centlon=173;
zoom{43}.horizon=15;
zoom{43}.altitude=1;
zoom{43}.name='New Zealand';
zoom{43}.annotation=0; % no annotation
zoom{43}.polar=0; 
 
zoom{44}.centlat=-40;
zoom{44}.centlon=146;
zoom{44}.horizon=10;
zoom{44}.altitude=1;
zoom{44}.name='Tasmania';
zoom{44}.annotation=0; % no annotation
zoom{44}.polar=0; 
 
zoom{45}.centlat=-72;
zoom{45}.centlon=-30;
zoom{45}.horizon=8;
zoom{45}.altitude=1;
zoom{45}.name='Filchner-Ronne Coupling';
zoom{45}.annotation=2; % ice shelf annotation
zoom{45}.polar=0; 
 
zoom{46}.centlat=-66;
zoom{46}.centlon=144;
zoom{46}.horizon=5;
zoom{46}.altitude=1;
zoom{46}.name='Mertz Polynya';
zoom{46}.annotation=2; % ice shelf annotation
zoom{46}.polar=0; 
 
zoom{47}.centlat=-72;
zoom{47}.centlon=174;
zoom{47}.horizon=5;
zoom{47}.altitude=1;
zoom{47}.name='Terra Nova Polynya';
zoom{47}.annotation=2; % ice shelf annotation
zoom{47}.polar=0; 
 
zoom{48}.centlat=3;
zoom{48}.centlon=-47;
zoom{48}.horizon=8;
zoom{48}.altitude=1;
zoom{48}.name='Amazon Outflow';
zoom{48}.annotation=0; % no annotation
zoom{48}.polar=0; 
 
zoom{49}.centlat=-12;
zoom{49}.centlon=138;
zoom{49}.horizon=8;
zoom{49}.altitude=1;
zoom{49}.name='Arafura Sea';
zoom{49}.annotation=0; % no annotation
zoom{49}.polar=0; 
 
zoom{50}.centlat=-49;
zoom{50}.centlon=-70;
zoom{50}.horizon=8;
zoom{50}.altitude=1;
zoom{50}.name='Patagonia';
zoom{50}.annotation=0; % no annotation
zoom{50}.polar=0; 

zoom{51}.centlat=14.5;
zoom{51}.centlon=166.4;
zoom{51}.horizon=0.75;
zoom{51}.altitude=1;
zoom{51}.name='Pacific Detail on Sphere';
zoom{51}.annotation=0; % no annotation
zoom{51}.polar=0; 
 
zoom{52}.centlat=26.5;
zoom{52}.centlon=144.0;
zoom{52}.horizon=20;
zoom{52}.altitude=1;
zoom{52}.name='Pacific 4xCO2 Crash Point';
zoom{52}.annotation=0; % no annotation
zoom{52}.polar=0; 
 
zoom{53}.centlat=26.5;
zoom{53}.centlon=-144.0;
zoom{53}.horizon=20;
zoom{53}.altitude=1;
zoom{53}.name='East Pacific Pentagon';
zoom{53}.annotation=0; % no annotation
zoom{53}.polar=0; 

zoom{54}.centlat=-26.5;
zoom{54}.centlon=180.0;
zoom{54}.horizon=20;
zoom{54}.altitude=1;
zoom{54}.name='South Pacific Pentagon';
zoom{54}.annotation=0; % no annotation
zoom{54}.polar=0; 

zoom{55}.centlat=26.5;
zoom{55}.centlon=-72;
zoom{55}.horizon=20;
zoom{55}.altitude=1;
zoom{55}.name='Atlantic Pentagon';
zoom{55}.annotation=0; % no annotation
zoom{55}.polar=0; 

zoom{56}.centlat=90; % degrees north
zoom{56}.centlon=-100; % degrees east
zoom{56}.horizon=30; % degrees of satellite horizon (0-90)
zoom{56}.altitude=1; % Mean Earth radius multiple
zoom{56}.name='Arctic';
zoom{56}.annotation=1; % add Arctic Shipping
zoom{56}.polar=0; 

zoom{57}.centlat=0; % degrees north
zoom{57}.centlon=0; % degrees east
zoom{57}.horizon=60; % degrees of satellite horizon (0-90)
zoom{57}.altitude=1; % Mean Earth radius multiple
zoom{57}.name='Greewhich Equator';
zoom{57}.annotation=0; % no annotation
zoom{57}.polar=0; 

zoom{58}.centlat=0; % degrees north
zoom{58}.centlon=90; % degrees east
zoom{58}.horizon=60; % degrees of satellite horizon (0-90)
zoom{58}.altitude=1; % Mean Earth radius multiple
zoom{58}.name='Mid East Equator';
zoom{58}.annotation=0; % no annotation
zoom{58}.polar=0; 

zoom{59}.centlat=0; % degrees north
zoom{59}.centlon=180; % degrees east
zoom{59}.horizon=60; % degrees of satellite horizon (0-90)
zoom{59}.altitude=1; % Mean Earth radius multiple
zoom{59}.name='Date Line Equator';
zoom{59}.annotation=0; % no annotation
zoom{59}.polar=0; 

zoom{60}.centlat=0; % degrees north
zoom{60}.centlon=-90; % degrees east
zoom{60}.horizon=60; % degrees of satellite horizon (0-90)
zoom{60}.altitude=1; % Mean Earth radius multiple
zoom{60}.name='Mid West Equator';
zoom{60}.annotation=0; % no annotation
zoom{60}.polar=0; 

zoom{61}.centlat=-90; % degrees north
zoom{61}.centlon=0; % degrees east
zoom{61}.horizon=60; % degrees of satellite horizon (0-90)
zoom{61}.altitude=1; % Mean Earth radius multiple
zoom{61}.name='Antarctic';
zoom{61}.annotation=2; % ice shelves
zoom{61}.polar=0; 

zoom{62}.centlat=60; % degrees north
zoom{62}.centlon=-90; % degrees east
zoom{62}.horizon=60; % degrees of satellite horizon (0-90)
zoom{62}.altitude=1; % Mean Earth radius multiple
zoom{62}.name='Arctic Combination 1';
zoom{62}.annotation=1; % add Arctic Shipping
zoom{62}.polar=0; 

zoom{63}.centlat=70; % degrees north
zoom{63}.centlon=-100; % degrees east
zoom{63}.horizon=50; % degrees of satellite horizon (0-90)
zoom{63}.altitude=1; % Mean Earth radius multiple
zoom{63}.name='Arctic Combination 2';
zoom{63}.annotation=1; % add Arctic Shipping
zoom{63}.polar=0; 

zoom{64}.centlat=90; % degrees north
zoom{64}.centlon=-100; % degrees east
zoom{64}.horizon=30; % degrees of satellite horizon (0-90)
zoom{64}.altitude=1; % Mean Earth radius multiple
zoom{64}.name='Arctic Combination 3';
zoom{64}.annotation=1; % add Arctic Shipping
zoom{64}.polar=0; 

zoom{65}.centlat=60; % degrees north
zoom{65}.centlon=-100; % degrees east
zoom{65}.horizon=60; % degrees of satellite horizon (0-90)
zoom{65}.altitude=1; % Mean Earth radius multiple
zoom{65}.name='United States Domain';
zoom{65}.annotation=1; % add Arctic Shipping
zoom{65}.polar=0; 

zoom{66}.centlat=-60; % degrees north
zoom{66}.centlon=135; % degrees east
zoom{66}.horizon=60; % degrees of satellite horizon (0-90)
zoom{66}.altitude=1; % Mean Earth radius multiple
zoom{66}.name='Antarctic Pan View Including Australia';
zoom{66}.annotation=2; % ice shelf annotation
zoom{66}.polar=1; 

zoom{67}.centlat=-88; % degrees north
zoom{67}.centlon=0; % degrees east
zoom{67}.horizon=28; % degrees of satellite horizon (0-90)
zoom{67}.altitude=1; % Mean Earth radius multiple
zoom{67}.name='Antarctic Pan View Including Atlantic';
zoom{67}.annotation=2; % ice shelf annotation
zoom{67}.polar=1; 


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
ncship5=ridgepack_clone('InteRFACE_Shiptracks',...
               {'track05_longitude','track05_latitude'});
ncship6=ridgepack_clone('InteRFACE_Shiptracks',...
               {'track06_longitude','track06_latitude'});
ncship7=ridgepack_clone('InteRFACE_Shiptracks',...
               {'track07_longitude','track07_latitude'});
ncship8=ridgepack_clone('InteRFACE_Shiptracks',...
               {'track08_longitude','track08_latitude'});
ncship9=ridgepack_clone('InteRFACE_Shiptracks',...
               {'track09_longitude','track09_latitude'});
ncship10=ridgepack_clone('InteRFACE_Shiptracks',...
               {'track10_longitude','track10_latitude'});
ncship11=ridgepack_clone('InteRFACE_Shiptracks',...
               {'track11_longitude','track11_latitude'});

lk=0;

for nrow=1:rows
for ncol=1:cols

 lk=lk+1;

 % multiplot only if multiple grids
 ridgepack_multiplot(rows,cols,nrow,ncol,alphass(lk))

 % grid information
 if strcmp(char(fileo{grids(lk)}.name),'oQU480')
  gridloc='/Users/afroberts/work';
  gridfile='mpassi.rst.0001-04-01_00000.nc';
  shiplocs=[];
 elseif strcmp(char(fileo{grids(lk)}.name),'SOwISC12to30E3r3')
  gridloc='/Users/afroberts/data/MODEL/E3SM/v3/SOwISC12to30E3r3';
  gridfile='mpaso.SOwISC12to30E3r3.20240829.nc';
  shiplocs=[1 2 3 4 6];
 elseif strcmp(char(fileo{grids(lk)}.name),'RRSwISC6to18E3r5')
  gridloc='/Users/afroberts/data/MODEL/E3SM/v3/RRSwISC6to18E3r5';
  gridfile='initial_state.nc';
  shiplocs=[1 2 3 4 6];
 elseif strcmp(char(fileo{grids(lk)}.name),'IcoswISC30E3r5')
  gridloc='/Users/afroberts/data/MODEL/E3SM/v3/IcoswISC30E3r5';
  gridfile='mpaso.IcoswISC30E3r5.20231120.nc';
  shiplocs=[1 2 3 4 6];
 else
  error('Unable to find mesh')
 end
  
 % extract grid data
 cd(gridloc)
 ncvert=ridgepack_clone(gridfile,{'latVertex','lonVertex',...
                                  'verticesOnCell','indexToCellID',...
                                  'nEdgesOnCell','edgesOnCell',...
                                  'cellsOnEdge','cellsOnVertex',...
                                  'edgesOnVertex','bottomDepth',...
                                  'landIceMask','landIceDraft',...
                                  'landIceFraction','restingThickness',...
                                  'maxLevelCell','fCell',...
                                  'areaCell','dcEdge'});
 nccell=ridgepack_clone(gridfile,{'latCell','lonCell',...
                                  'verticesOnCell','indexToCellID',...
                                  'nEdgesOnCell','edgesOnCell',...
                                  'cellsOnEdge','cellsOnVertex',...
                                  'edgesOnVertex','bottomDepth',...
                                  'landIceMask','landIceDraft',...
                                  'landIceFraction','fCell',...
                                  'areaCell','dcEdge'});
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
  
 % read in coast or generate coast and write it out
 coastname=[char(fileo{grids(lk)}.outname),'_Coast.nc'];
 x=dir(coastname);
 if isempty(x)
  nccoast=ridgepack_e3smcoastm(ncvert);
  ridgepack_write(nccoast,coastname)
 else
  nccoast=ridgepack_clone(coastname);
 end

 % check bounds on landIceFraction if it exists;
 iceshelfplot=false;
 if isfield(ncvert,'landIceFraction')
    if any(ncvert.landIceFraction.data>1+10^-12) 
     error('landIceFraction is out of bounds: greater than 1')
    elseif any(ncvert.landIceFraction.data<0)
     error('landIceFraction is out of bounds: less than zero')
    elseif any(ncvert.landIceFraction.data>0)
     iceshelf=[char(fileo{grids(lk)}.outname),'_iceshelf.nc'];
     x=dir(iceshelf);
     if isempty(x) 
      nciceshelf=ridgepack_e3smseasaw(ncvert,ncvert,'landIceFraction',0.5);
      ridgepack_write(nciceshelf,iceshelf);
     else
      nciceshelf=ridgepack_clone(iceshelf);
     end
     iceshelfplot=true;
    end
    if any(ncvert.landIceDraft.data>0)
     error('landIceDraft is out greater than zero')
    elseif any(ncvert.landIceDraft.data<-ncvert.bottomDepth.data)
     error('landIceDraft is less than the ocean draft')
    else
     ncvert.cavityThickness=ncvert.bottomDepth;
     ncvert.cavityThickness.long_name='Cavity thickness under ice shelves';
     ncvert.cavityThickness.data=-ncvert.bottomDepth.data-ncvert.landIceDraft.data;
     ncvert.cavityThickness.data(ncvert.landIceMask.data==0)=NaN;
     ncvert.iceShelfMask=ncvert.landIceMask;
     ncvert.iceShelfMask.data(ncvert.landIceMask.data==1)=NaN;
     ncvert.iceShelfMask.data(ncvert.landIceMask.data==0)=1;
    end
    if any(ncvert.landIceDraft.data==0 & ncvert.landIceFraction.data>0)
     xidx=find(ncvert.landIceDraft.data(:)==0 & ncvert.landIceFraction.data(:)>0);
     disp('land ice draft    land ice fraction')
     disp([num2str(ncvert.landIceDraft.data(xidx))])
     disp([num2str(ncvert.landIceFraction.data(xidx))])
     %error('cavity thickness equal to draft')
    end
 else
    disp('No land ice fraction in this mesh')
 end

 % create isobath contours
 % minimum isobath
 mbatname=[char(fileo{grids(lk)}.outname),'_minIsobath.nc'];
 x=dir(mbatname);
 minIso=ceil(10*min(ncvert.bottomDepth.data(:)))./10;
 mdepthc=[num2str(minIso),'m'];
 if isempty(x) 
  disp([mdepthc,' minimum isobath...'])
  ncisobathmin=ridgepack_e3smseasaw(ncvert,ncvert,'bottomDepth',minIso);
  ridgepack_write(ncisobathmin,mbatname)
 else
  ncisobathmin=ridgepack_clone(mbatname);
 end

 % 50m isobath
 greename=[char(fileo{grids(lk)}.outname),'_50mIsobath.nc'];
 x=dir(greename);
 if isempty(x)                    
  disp('50m isobath...')
  ncisobath50=ridgepack_e3smseasaw(ncvert,ncvert,'bottomDepth',50);
  ridgepack_write(ncisobath50,greename)
 else                             
  ncisobath50=ridgepack_clone(greename);
 end                              

 % 500m isobath
 deepname=[char(fileo{grids(lk)}.outname),'_500mIsobath.nc'];
 x=dir(deepname);
 if isempty(x)                    
  disp('500m isobath...')         
  ncisobath500=ridgepack_e3smseasaw(ncvert,ncvert,'bottomDepth',500);
  ridgepack_write(ncisobath500,deepname)
 else                             
  ncisobath500=ridgepack_clone(deepname);
 end                              

 % calculate selective metrics if desired
 if rossbymetrics{lk}
  ncvert.rossby=ncvert.bottomDepth;
  ncvert.rossby.long_name='Model Resolution to Rossby Radius Ratio';
  ncvert.rossby.units='';
  ncvert.rossby.data=abs(sqrt(ncvert.areaCell.data).*ncvert.fCell.data)./...
                     sqrt(ncvert.bottomDepth.data.*9.80665);

 elseif resolutionmetrics{lk}
  ncvert.resolution=ncvert.bottomDepth;
  ncvert.resolution.long_name='Model Resolution';
  ncvert.resolution.units='km';
  ncvert.resolution.data=sqrt(ncvert.areaCell.data)./1000;

 elseif distortions{lk}
  ncvert.distortion=ncvert.bottomDepth;
  ncvert.distortion.long_name='Mesh Distortion';
  ncvert.distortion.units='';

  for iCell=1:length(ncvert.bottomDepth.data)
   edges=ncvert.edgesOnCell.data(:,iCell);
   maxedge=max(ncvert.dcEdge.data(edges(edges>0)));
   minedge=min(ncvert.dcEdge.data(edges(edges>0)));
   ncvert.distortion.data(iCell)=minedge./maxedge;
  end

 end

 % invert bathymetry
 ncvert.bottomDepth.data=-ncvert.bottomDepth.data;
  
 if largescales{lk}

    centlat=zoom{projections(lk)}.centlat;   % degrees north
    centlon=zoom{projections(lk)}.centlon;   % degrees east
    horizon=zoom{projections(lk)}.horizon;   % degrees of satellite horizon (0-90)
    altitude=zoom{projections(lk)}.altitude; % Mean Earth radius multiple
  
    ridgepack_satview(centlat,centlon,horizon)

    if bathymetrys{lk}

     if rossbymetrics{lk}

      % outfilename modifier
      outfilename=[outfilename,'rossbymetric_'];

      % reverse colorbar
      cont=[0 0.001 0.002:0.002:0.008 0.01:0.01:0.09 0.1 0.2:0.2:0.8 1.0 2.0 3.0 4.0];

      % render colors
      ridgepack_e3smsatcol(ncvert,'rossby',ncvert,cont,0,...
                           centlat,centlon,horizon,altitude,...
                           true,false,'linear','bluered');

      % add colormap
      if ncol==cols & nrow==rows
       ridgepack_colorbar(cont,'\Delta{x}/L_R','linear','vertical',0)
      end

     elseif resolutionmetrics{lk}

      % outfilename modifier
      outfilename=[outfilename,'resolution_'];

      % reverse colorbar
      minres=floor(min(ncvert.resolution.data(:)));
      maxres=ceil(max(ncvert.resolution.data(:)));
      increment=floor(10*(maxres-minres)./18)./10;
      maxres=minres+18.*increment;
      cont=[minres-increment:increment:maxres+increment];

      % render colors
      ridgepack_e3smsatcol(ncvert,'resolution',ncvert,cont,0,...
                           centlat,centlon,horizon,altitude,...
                           true,false,'linear','parula');

      % add colormap
      if ncol==cols & nrow==rows
       ridgepack_colorbar(cont,'\Delta{x} (km)','linear','vertical',0)
      end

     elseif distortions{lk}

      % outfilename modifier
      outfilename=[outfilename,'distortion_'];

      % reverse colorbar
      %mindistortion=floor(10*min(ncvert.distortion.data))./10;
      %maxdistortion=ceil(10*max(ncvert.distortion.data))./10;
      cont=[0.85:0.01:1];

      % render colors
      ridgepack_e3smsatcol(ncvert,'distortion',ncvert,cont,0,...
                           centlat,centlon,horizon,altitude,...
                           true,false,'linear','bluered');

      % add colormap
      if ncol==cols & nrow==rows
       ridgepack_colorbar(cont,'Distortion','linear','vertical',0)
      end

     else
 
      % outfilename modifier
      outfilename=[outfilename,'bathymetry_'];
  
      % reverse colorbar
      cont=[-5000:1000:-2000 -1500 -1000:250:-250 -100 -50:10:-20 0 10]; 
      colbarcont{1}='\downarrow';
      for i=2:length(cont)-1
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
      if ncol==cols & drawcb
       ridgepack_colorbar(cont,'m','linear','vertical',0,colbarcont)
       clear cmap colbarcont
       drawcb=false;
       if rows>1 & nrow==1
        ridgepack_cbshare(gca)
       end
      end

     end

     % add min, 50m, and 500m isobath
     hd=ridgepack_e3smsatthreshold(ncisobath500,centlat,centlon,horizon,...
                                  [0.9290 0.6940 0.1250]);
     if projections(lk)~=67
      hg=ridgepack_e3smsatthreshold(ncisobath50,centlat,centlon,horizon,...
                                  [0.4660 0.6740 0.1880]);
      ht=ridgepack_e3smsatthreshold(ncisobathmin,centlat,centlon,horizon,...
                                  [0.8500 0.3250 0.0980]);
     end
  
     % plot coast
     ridgepack_e3smsatcoast(nccoast,centlat,centlon,horizon)

     % plot ice shelf edge
     if iceshelfplot
      hs=ridgepack_e3smsatthreshold(nciceshelf,centlat,centlon,horizon,0.45*[1 0 0]);
      if projections(lk)==6 & iceshelflegflag
       legendh(end+1)=hs;
       legendcell{end+1}='ice shelf edge';
       iceshelflegflag=false;
      end
     end

     % add minimum depth
     if any(diff(grids))~=0
      xlims=xlim;
      ylims=ylim;
      text((xlims(2)-diff(xlims)/2)+cosd(45)*diff(xlims)/2,...
           (ylims(2)-diff(ylims)/2)-sind(45)*diff(ylims)/2,...
           [mdepthc,' min'],...
           'Color',[0.8500 0.3250 0.0980],...
           'FontSize',5,...
           'Rotation',45,...
           'VerticalAlignment','top',...
           'HorizontalAlignment','center')
     end

    else % not bathymetry
  
     % outfilename modifier
     outfilename=[outfilename,'_'];

     if zoom{projections(lk)}.annotation==2 & iceshelfplot % render ice shelf 
      ridgepack_e3smsatmeshs(ncvert,centlat,centlon,horizon,altitude,ncvert,'iceShelfMask',true);
     else
      ridgepack_e3smsatmeshs(ncvert,centlat,centlon,horizon,altitude);
     end

     ridgepack_e3smsatcoast(nccoast,centlat,centlon,horizon)
  
     if zoomedareas{lk} & lk==1;
 
      % outfilename modifier
      outfilename=[outfilename,'zoomedareas_'];
  
      for j=1:length(zoom)
       satlat=zoom{j}.centlat;   % degrees north
       satlon=zoom{j}.centlon;   % degrees east
       sathor=zoom{j}.horizon;   % degrees of satellite horizon (0-90)
       ridgepack_sathorizon(centlat,centlon,horizon,...
                            satlat,satlon,sathor,[0 0 0.5],num2str(j));
       %ridgepack_sathorizon(centlat,centlon,horizon,...
       %                     satlat,satlon,sathor,[0 0.5 0]);
      end
  
     end

     if resolutionmetrics{lk}

      ncvert.resolution=ncvert.bottomDepth;
      ncvert.resolution.long_name='Model Resolution';
      ncvert.resolution.units='km';
      ncvert.resolution.data=sqrt(ncvert.areaCell.data)./1000;

      ncisores25=ridgepack_e3smseasaw(ncvert,ncvert,'resolution',25);
      ncisores27=ridgepack_e3smseasaw(ncvert,ncvert,'resolution',27);
      ncisores28=ridgepack_e3smseasaw(ncvert,ncvert,'resolution',28);

      if config==3
       vcols(:,1)=[0 0.0 0];
       vcols(:,2)=[0 0.5 0];
       vcols(:,3)=[0 0.5 0];
      else
       vcols(:,1)=[0.0 0 0];
       vcols(:,2)=[0.5 0 0];
       vcols(:,3)=[0.9 0 0];
      end

      hay=ridgepack_e3smsatthreshold(ncisores25,centlat,centlon,horizon,...
                                   vcols(:,1),'-');
      hby=ridgepack_e3smsatthreshold(ncisores27,centlat,centlon,horizon,...
                                   vcols(:,2),'-');
      hby=ridgepack_e3smsatthreshold(ncisores28,centlat,centlon,horizon,...
                                   vcols(:,3),'-');


      % make isolines of resolution
      for ixl=1:3

       % arctic labeling
       if ixl==1
        latmark=[max(ncisores25.xlatitude.data) max(ncisores25.xlatitude.data)];
        reslabel='25';
        colx=vcols(:,1);
       elseif ixl==2
        latmark=[max(ncisores27.xlatitude.data) max(ncisores27.xlatitude.data)];
        reslabel='27';
        colx=vcols(:,2);
       elseif ixl==3
        latmark=[max(ncisores28.xlatitude.data) max(ncisores28.xlatitude.data)];
        reslabel='28';
        colx=vcols(:,3);
       end

       lonmark=[-2 2];
      
       % find local angle of labels
       [xl,yl,zl,phl,thl]=ridgepack_satfwd(latmark,lonmark,centlat,centlon,horizon,1);

       % find grid label angle & rotate to readable angle
       rotation=atan2d(diff(yl),diff(xl));
       if isnan(rotation)
        rotation=90;
       end

       if ~isnan(rotation)

        % plot labels inside frame
        rotation=wrapTo180(rotation);
        if rotation>91 | rotation<-91
         rotation=mod(rotation+180,0);
        end

        gridheight=1.01; % height of grid superimposed on plot
        text(xl(1),yl(1),2*gridheight*zl(1),reslabel,...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom',...
               'FontSize',5,...
               'Color',colx,...
               'Rotation',rotation);

        if ixl==1
         text(xl(1),yl(1),2*gridheight*zl(1),'km',...
               'HorizontalAlignment','center',...
               'VerticalAlignment','top',...
               'FontSize',4,...
               'Color',colx,...
               'Rotation',rotation);
        end

       end

       % antarctic labeling
       if ixl>1
        if ixl==2
         latmark=[min(ncisores27.xlatitude.data) min(ncisores27.xlatitude.data)];
         reslabel='27';
         colx=vcols(:,2);
        elseif ixl==3
         latmark=[min(ncisores28.xlatitude.data) min(ncisores28.xlatitude.data)];
         reslabel='28';
         colx=vcols(:,3);
        end
 
        lonmark=[-176 -178];
 
        % find local angle of labels
        [xl,yl,zl,phl,thl]=ridgepack_satfwd(latmark,lonmark,centlat,centlon,horizon,1);
 
        % find grid label angle & rotate to readable angle
        rotation=atan2d(diff(yl),diff(xl));
        if isnan(rotation)
         rotation=90;
        end

        if ~isnan(rotation)

         % plot labels inside frame
         rotation=wrapTo180(rotation);
         if rotation>91 | rotation<-91
          rotation=mod(rotation+180,0);
         end
 
         gridheight=1.01; % height of grid superimposed on plot
         text(xl(1),yl(1),2*gridheight*zl(1),reslabel,...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom',...
               'FontSize',5,...
               'Color',colx,...
               'Rotation',rotation);

         if ixl==3
          text(xl(1),yl(1),2*gridheight*zl(1),'km',...
               'HorizontalAlignment','center',...
               'VerticalAlignment','top',...
               'FontSize',5,...
               'Color',colx,...
               'Rotation',rotation);
         end

        end
 
       end

      end

     end
  
     if zoom{projections(lk)}.annotation==1
  
      for shipi=shiplocs
       latship=eval(['ncship',num2str(shipi),'.latitude.data']);
       lonship=eval(['ncship',num2str(shipi),'.longitude.data']);
       [x,y,z,phi,theta]=...
        ridgepack_satfwd(latship,lonship,centlat,centlon,horizon,1.001*altitude);
        %hship=plot3(x,y,z,'Color','#FF8800');
        hship=plot3(x,y,z,'Color',[0 0.5 0.1]);

       if shipi==4 | shipi==6

        % label text
        if shipi==4 
         shiplabel={'Northwest','Passage'};
        elseif shipi==6
         shiplabel={'Northern','Sea Route'};
        end

        % find local angle of labels
        [xl,yl,zl,phl,thl]=ridgepack_satfwd(latship([end-5 end]),lonship([end-5 end]),...
                                            centlat,centlon,horizon,1);
                          
        % find grid label angle & rotate to readable angle
        rotation=atan2d(diff(yl),diff(xl));
        if isnan(rotation)
         rotation=90;
        end

        if ~isnan(rotation)

         % plot labels inside frame
         rotation=wrapTo180(rotation);
         if rotation>91 | rotation<-91
          rotation=mod(rotation+180,0);
         end

         gridheight=1.01; % height of grid superimposed on plot
         text(xl(1),yl(1),2*gridheight*zl(1),shiplabel,...
               'HorizontalAlignment','left',...
               'VerticalAlignment','middle',...
               'FontSize',4,...
               'Color',0.99*[1 1 1],...
               'Rotation',rotation);

        end
 
       end

      end
  
      if ~isempty(shiplocs) & shiplegflag
       legendh(end+1)=hship
       legendcell{end+1}='Arctic coastal routes';
       shiplegflag=false;
      end

     elseif zoom{projections(lk)}.annotation==2 & iceshelfplot % render ice shelf 

      % plot ice shelf edge
      hs=ridgepack_e3smsatthreshold(nciceshelf,centlat,centlon,horizon,0.45*[0.9 0 0]);
      hsleg='ice shelf edge';

      % render ice shelf cavity thickness
      cont=[-5000:1000:-2000 -1500 -1000:250:-250 -100 -50:10:-20 0 10]; 
      colbarcont{1}='\downarrow';
      for i=2:length(cont)-1;
       colbarcont{i}=num2str(-cont(i));
      end
      colbarcont{length(cont)}='\uparrow';

      ridgepack_e3smsatcol(ncvert,'cavityThickness',ncvert,cont,0,...
                           centlat,centlon,horizon,altitude,...
                           true,false,'linear','bluered',false);

      % make top color grey
      cmap=colormap;
      cmap(end,:)=0.95*[1 1 1];
      colormap(cmap)

      % add colormap
      if ncol==cols & drawcb
       ridgepack_colorbar(cont,'m','linear','vertical',0,colbarcont)
       clear cmap colbarcont
       drawcb=false;
       if rows>1 & nrow==1
        ridgepack_cbshare(gca)
       end
      end

      if meshtitle
       title('Cavity Thickness','FontWeight','normal')
      end

      % plot ice shelf edge
      hs=ridgepack_e3smsatthreshold(nciceshelf,centlat,centlon,horizon,0.45*[1 0 0]);
      if zoom{projections(lk)}.annotation==2 & iceshelflegflag
       legendh(end+1)=hs;
       legendcell{end+1}='ice shelf edge';
       iceshelflegflag=false;
      end


      % add min, 50m, and 500m isobath
      hd=ridgepack_e3smsatthreshold(ncisobath500,centlat,centlon,horizon,...
                                    [0.9290 0.6940 0.1250]);
      if projections(lk)~=67
       hg=ridgepack_e3smsatthreshold(ncisobath50,centlat,centlon,horizon,...
                                    [0.4660 0.6740 0.1880]);
       ht=ridgepack_e3smsatthreshold(ncisobathmin,centlat,centlon,horizon,...
                                    [0.8500 0.3250 0.0980]);
      end

     end

    end  % bathymetry conditional

    % bathymetry legend
    if bathymetrys{lk} | (zoom{projections(lk)}.annotation==2 & iceshelfplot)
     if any(diff(grids))~=0 & bathlegflag
      if projections(lk)==67
       legendh(end+1)=[hd];
       legendcell{end+1}='500m isobath';
      else
       legendh(end+1:end+3)=[ht hg hd];
       legendcell{end+1}='minimum depth';
       legendcell{end+1}='50m';
       legendcell{end+1}='500m isobaths';
      end
      bathlegflag=false;
     elseif bathlegflag
      if projections(lk)==67
       legendh(end+1)=[hd];
       legendcell{end+1}='500m isobath';
      else
       legendh(end+1:end+3)=[ht hg hd];
       legendcell{end+1}=[mdepthc,' minimum'];
       legendcell{end+1}='50m';
       legendcell{end+1}='500m isobaths';
      end
      bathlegflag=false;
     end
    end

    % label mesh if desired
    if meshtitle
     title([char(fileo{grids(lk)}.outname)],'FontWeight','normal')
    end

 else % not largescale
 
    % outfilename modifier
    outfilename='_zoom_';
  
    centlat=zoom{projections(lk)}.centlat;   % degrees north
    centlon=zoom{projections(lk)}.centlon;   % degrees east
    horizon=zoom{projections(lk)}.horizon;   % degrees of satellite horizon (0-90)
    altitude=zoom{projections(lk)}.altitude; % Mean Earth radius multiple
  
    ridgepack_satview(centlat,centlon,horizon)

    if bathymetrys{lk} & zoom{projections(lk)}.annotation==2 % render ice shelf 
     ridgepack_e3smsatmeshs(ncvert,centlat,centlon,horizon,altitude,ncvert,'iceShelfMask');
    else
     ridgepack_e3smsatmeshs(ncvert,centlat,centlon,horizon,altitude);
    end
  
    % add isobaths
    if bathymetrys{lk}
     hk=ridgepack_e3smsatthreshold(ncisobath500,centlat,centlon,horizon,...
        [0.9290 0.6940 0.1250]);
     if projections(lk)~=67
      hk=ridgepack_e3smsatthreshold(ncisobath50,centlat,centlon,horizon,...
        [0.4660 0.6740 0.1880]);
      hk=ridgepack_e3smsatthreshold(ncisobathmin,centlat,centlon,horizon,...
        [0.8500 0.3250 0.0980]);
     end
     set(hk,'LineWidth',0.5)
    end
  
    % add in coastline
    ridgepack_e3smsatcoast(nccoast,centlat,centlon,horizon)

    if zoom{projections(lk)}.annotation==2 & iceshelfplot % render ice shelf 

     % plot ice shelf edge
     hs=ridgepack_e3smsatthreshold(nciceshelf,centlat,centlon,horizon,0.45*[0.9 0 0]);
     hsleg='ice shelf edge';

     % render ice shelf cavity thickness
     cont=[-5000:1000:-2000 -1500 -1000:250:-250 -100 -50:10:-20 0 10]; 
     ridgepack_e3smsatcol(ncvert,'cavityThickness',ncvert,cont,0,...
                          centlat,centlon,horizon,altitude,...
                          true,false,'linear','bluered',false);

     if meshtitle
      title('Cavity Thickness','FontWeight','normal')
     end

    elseif zoom{projections(lk)}.annotation==-1
  
     % add in high resolution coastline
     hs=ridgepack_e3smsatcoast(nccoasthr,centlat,centlon,horizon,'m')
     hsleg='E3SM-HR V1 coastline';

    end
 
    xlims=xlim;
 
    if bathymetrys{lk}
 
     % outfilename modifier
     outfilename=[outfilename,'bathymetry_'];
  
     ridgepack_satview(centlat,centlon,horizon,1,1)
  
     % reverse colorbar
     cont=[-5000:1000:-2000 -1500 -1000:250:-250 -100 -50:10:-20 0 10]; 
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

     % add colorbar
     if ncol==cols & drawcb
      ridgepack_colorbar(cont,'m','linear','vertical',0,colbarcont)
      clear colbarcont
      drawcb=false;
      if rows>1 & nrow==1
       ridgepack_cbshare(gca)
      end
     end
  
     % add min, 50m and 500m isobath
     hd=ridgepack_e3smsatthreshold(ncisobath500,centlat,centlon,horizon,...
                                   [0.9290 0.6940 0.1250]);
     if projections(lk)~=67
      hg=ridgepack_e3smsatthreshold(ncisobath50,centlat,centlon,horizon,...
                                   [0.4660 0.6740 0.1880]);
      ht=ridgepack_e3smsatthreshold(ncisobathmin,centlat,centlon,horizon,...
                                   [0.8500 0.3250 0.0980]);
     end

     % plot ice shelf edge
     if zoom{projections(lk)}.annotation==2 & iceshelfplot
      hs=ridgepack_e3smsatthreshold(nciceshelf,centlat,centlon,horizon,0.45*[0.9 0 0]);
      if iceshelflegflag
       legendh(end+1)=hs;
       legendcell{end+1}='ice shelf edge';
       iceshelflegflag=false;
      end
     end
  
     % add coast
     ridgepack_e3smsatcoast(nccoast,centlat,centlon,horizon)

     % add legend
     if any(diff(grids))~=0 & bathlegflag
      if projections(lk)==67
       legendh(end+1)=[hd];
       legendcell{end+1}='500m isobath';
      else
       legendh(end+1:end+3)=[ht hg hd];
       legendcell{end+1}='minimum depth';
       legendcell{end+1}='50m';
       legendcell{end+1}='500m isobaths';
      end
      bathlegflag=false;
     elseif bathlegflag
      if projections(lk)==67
       legendh(end+1)=[hd];
       legendcell{end+1}='500m isobath';
      else
       legendh(end+1:end+3)=[ht hg hd];
       legendcell{end+1}='minimum depth';
       legendcell{end+1}='50m';
       legendcell{end+1}='500m isobaths';
      end
      bathlegflag=false;
     end

     % add minimum depth
     if any(diff(grids))~=0
      xlims=xlim;
      ylims=ylim;
      text((xlims(2)-diff(xlims)/2)+cosd(45)*diff(xlims)/2,...
           (ylims(2)-diff(ylims)/2)+sind(45)*diff(ylims)/2,...
           [mdepthc,' min'],...
           'Color',[0.8500 0.3250 0.0980],...
           'FontSize',5,...
           'Rotation',-45,...
           'VerticalAlignment','bottom',...
           'HorizontalAlignment','center')
     end

    else
  
     ridgepack_satview(centlat,centlon,horizon)
     ridgepack_e3smsatmeshv(nccell,centlat,centlon,horizon,altitude);
     ridgepack_e3smsatcoast(nccoast,centlat,centlon,horizon)

     if overlaymeshs{lk}

      % add edge points
      %[x,y,z,phi,theta]=...
      % ridgepack_satfwd(ncedge.latitude.data*180/pi,ncedge.longitude.data*180/pi,...
      %                  centlat,centlon,horizon,1.001*altitude);
      %plot3(x,y,z,'.','Color','b');

      % plot atan(lat)=0.5
      tlont=[0:0.001:360];
      tlatt=atand(0.5)*ones(size(tlont));
      [x,y,z,phi,theta]=...
       ridgepack_satfwd(tlatt,tlont,centlat,centlon,horizon,1.001*altitude);
      plot3(x,y,z,'Color','c');

      % outfilename modifier
      outfilename=[outfilename,'mesh+composite_'];

      % add atmospheric model points
      ncat=ridgepack_clone('/Users/afroberts/data/MODEL/E3SM/atmosgrid/ne30np4_latlon.091226.nc')

      [x,y,z,phi,theta]=...
       ridgepack_satfwd(ncat.latitude.data,ncat.longitude.data,...
                        centlat,centlon,horizon,1.001*altitude);
      plot3(x,y,z,'.','Color',[0 0.5 0]);

     else

      % outfilename modifier
      outfilename=[outfilename,'mesh_'];

     end

     % plot ice shelf edge
     if zoom{projections(lk)}.annotation==1 

      for shipi=[1 2 3 4 6]
       [x,y,z,phi,theta]=...
        ridgepack_satfwd(eval(['ncship',num2str(shipi),'.latitude.data']),...
                         eval(['ncship',num2str(shipi),'.longitude.data']),...
                         centlat,centlon,horizon,1.001*altitude);
       hship=plot3(x,y,z,'-','Color','k');
      end

      if ~isempty(shiplocs) & shiplegflag
       legendh(end+1)=hship;
       legendcell{end+1}='Arctic coastal routes';
       shiplegflag=false;
      end

     elseif zoom{projections(lk)}.annotation==2 & iceshelfplot

      hx=ridgepack_e3smsatthreshold(nciceshelf,centlat,centlon,horizon,0.45*[0 0 1]);
      if projections(lk)==6 & iceshelflegflag
       legendh(end+1)=hx;
       legendcell{end+1}='ice shelf edge';
       iceshelflegflag=false;
      end

     end
  
    end % bathymetrys{lk}

    if meshtitle
     text(1.15*xlims(1),0,fileo{grids(lk)}.outname,...
         'Rotation',90,'HorizontalAlignment','center')
    end

 end % largescale or not

 if lk==1
  grifilename=[fileo{grids(lk)}.outname];
 else
  grifilename=[grifilename,'_',fileo{grids(lk)}.outname];
 end

end
end

% multlegend
if lk>1
 lh=ridgepack_multilegend(legendh,legendcell,'South')
 set(lh,'FontSize',7)
 set(lh,'Position',[(1-lh.Position(3))/2 lh.Position(2) ...
                     lh.Position(3) lh.Position(4)])
else
 legend(legendh,legendcell,'Location','SouthOutside','Orientation','Horizontal')
 legend('boxoff')
end

% multialign
ridgepack_multialign(gcf);

% move to plot location
cd(plotloc)
  
% print png
ridgepack_fprint('png',[grifilename,outfilename,num2str(projections(lk))],1,2)

