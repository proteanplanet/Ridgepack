
clf
clear

gridchoice=1;

fileg{1}.name='CUSP12';
fileg{1}.title=' CUSP 12-60~km';
fileg{2}.name='CUSP14';
fileg{2}.title=' CUSP 14-60~km';

meshtype='tesselation';
%meshtype='triangulation';

largescale=true;
%largescale=false;

if strcmp(meshtype,'tesselation')
 plotchoice=[1:7]
elseif strcmp(meshtype,'triangulation')
 plotchoice=[3 5]
end

plotchoice=[1 2 3 4 5];

plotchoice=[1];

sector{1}.centlat=90; % degrees north
sector{1}.centlon=0; % degrees east
sector{1}.horizon=60; % degrees of satellite horizon (0-90)
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
zoom{2}.horizon=11;
zoom{2}.altitude=1;
zoom{2}.name='Baltic Sea and White Sea';

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

zoom{5}.centlat=63;
zoom{5}.centlon=-176;
zoom{5}.horizon=8;
zoom{5}.altitude=1;
zoom{5}.name='Bering Sea and Chukchi Sea';

zoom{6}.centlat=55;
zoom{6}.centlon=150;
zoom{6}.horizon=11;
zoom{6}.altitude=1;
zoom{6}.name='Sea of Okhotsk';

zoom{7}.centlat=75;
zoom{7}.centlon=135;
zoom{7}.horizon=8;
zoom{7}.altitude=1;
zoom{7}.name='Laptev Sea';

zoom{8}.centlat=40;
zoom{8}.centlon=135;
zoom{8}.horizon=8;
zoom{8}.altitude=1;
zoom{8}.name='Sea of Japan';

zoom{9}.centlat=61;
zoom{9}.centlon=-76.5;
zoom{9}.horizon=10;
zoom{9}.altitude=1;
zoom{9}.name='Hudson Bay';

zoom{10}.centlat=55;
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

zoom{13}.centlat=33;
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
zoom{15}.horizon=8;
zoom{15}.altitude=1;
zoom{15}.name='South Greenland';

zoom{16}.centlat=76;
zoom{16}.centlon=-40;
zoom{16}.horizon=9;
zoom{16}.altitude=1;
zoom{16}.name='North Greenland';

zoom{17}.centlat=55;
zoom{17}.centlon=-5;
zoom{17}.horizon=7;
zoom{17}.altitude=1;
zoom{17}.name='British Isles';

zoom{18}.centlat=40;
zoom{18}.centlon=0;
zoom{18}.horizon=9;
zoom{18}.altitude=1;
zoom{18}.name='West Mediterranean';

zoom{19}.centlat=35;
zoom{19}.centlon=22;
zoom{19}.horizon=14;
zoom{19}.altitude=1;
zoom{19}.name='East Mediterranean';

zoom{20}.centlat=35;
zoom{20}.centlon=20;
zoom{20}.horizon=10;
zoom{20}.altitude=1;
zoom{20}.name='East Mediterranean';

zoom{20}.centlat=34;
zoom{20}.centlon=-73;
zoom{20}.horizon=10;
zoom{20}.altitude=1;
zoom{20}.name='Atlantic Coast';









% define contours
% define contours
cont=[5 10:1:31];
ref=20;

% plot location
plotloc='/Users/afroberts/work';

% grid location
gridloc=['/Users/afroberts/data/MODEL/E3SM/CUSP/',...
         char(fileg{gridchoice}.name)];
gridfile='initial_state.nc';
coastname='CUSP Version 2 14km'; % grid name

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

cd(plotloc)

for setting=plotchoice

 clf

 if largescale
  centlat=sector{setting}.centlat;   % degrees north
  centlon=sector{setting}.centlon;   % degrees east
  horizon=sector{setting}.horizon;   % degrees of satellite horizon (0-90)
  altitude=sector{setting}.altitude; % Mean Earth radius multiple
 else
  centlat=zoom{setting}.centlat;   % degrees north
  centlon=zoom{setting}.centlon;   % degrees east
  horizon=zoom{setting}.horizon;   % degrees of satellite horizon (0-90)
  altitude=zoom{setting}.altitude; % Mean Earth radius multiple
 end

 ridgepack_satview(centlat,centlon,horizon)

 if strcmp(meshtype,'tesselation')
  ridgepack_e3smsatmeshs(ncvert,centlat,centlon,horizon,altitude);
 elseif strcmp(meshtype,'triangulation')
  ridgepack_e3smsatmeshv(nccell,centlat,centlon,horizon,altitude);
 else
  error('nothing to see here')
 end

 ridgepack_e3smsatcoast(ncvert,centlat,centlon,horizon)

 if largescale
  for j=1:length(zoom)
   satlat=zoom{j}.centlat;   % degrees north
   satlon=zoom{j}.centlon;   % degrees east
   sathor=zoom{j}.horizon;   % degrees of satellite horizon (0-90)
   ridgepack_sathorizon(centlat,centlon,horizon,...
                        satlat,satlon,sathor,[0.83 0.5 0])
  end
 end

 if largescale
  title(['Sector ',num2str(setting),' ',char(fileg{gridchoice}.title)])
 else
  title([zoom{setting}.name,' ',char(fileg{gridchoice}.title)])
 end

 ridgepack_fprint('png',[fileg{gridchoice}.name,'_',...
                          meshtype,'_',num2str(setting)],1,2)

end












