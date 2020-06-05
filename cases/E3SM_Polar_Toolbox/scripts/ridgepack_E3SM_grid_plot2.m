
clf
clear

fileg='CUSP14'

%tesselation=true;
tesselation=false;

if tesselation
 plotchoice=[1:7]
else
 plotchoice=[3 5]
end

plotchoice=1;

for setting=plotchoice

clf

if setting==1

 % Arctic Center
 centlat=80; % degrees north
 centlon=-100; % degrees east
 horizon=60; % degrees of satellite horizon (0-90)
 altitude=1; % Mean Earth radius multiple
 satlat=72; % degrees north
 satlon=-102; % degrees east
 sathorizon=20; % degrees of satellite horizon (0-90)
 titled='Central Arctic';

elseif setting==2

 % Arctic Center
 centlat=75; % degrees north
 centlon=-102; % degrees east
 horizon=20; % degrees of satellite horizon (0-90)
 altitude=1; % Mean Earth radius multiple
 satlat=72; % degrees north
 satlon=-102; % degrees east
 sathorizon=8; % degrees of satellite horizon (0-90)
 titled='Canadian Archipelago';

elseif setting==3

 % Arctic Center
 centlat=72; % degrees north
 centlon=-102; % degrees east
 horizon=8; % degrees of satellite horizon (0-90)
 altitude=1; % Mean Earth radius multiple
 titled='Northwest Passage';

elseif setting==4

 centlat=0; % degrees north
 centlon=-30; % degrees east
 horizon=60; % degrees of satellite horizon (0-90)
 altitude=1; % Mean Earth radius multiple
 titled='Atlantic';

elseif setting==5

 % Arctic Center
 centlat=60; % degrees north
 centlon=-180; % degrees east
 horizon=10; % degrees of satellite horizon (0-90)
 altitude=1; % Mean Earth radius multiple
 titled='North Pacific';

elseif setting==6

 % Arctic Center
 centlat=-90; % degrees north
 centlon=0; % degrees east
 horizon=60; % degrees of satellite horizon (0-90)
 altitude=1; % Mean Earth radius multiple
 titled='Antarctic';

end

% define contours
cont=[5 10:1:31];
ref=20;

% plot location
plotloc='/Users/afroberts/work';

% grid location
gridloc=['/Users/afroberts/data/MODEL/E3SM/CUSP/',fileg];
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

if setting<7

% ridgepack_satview(centlat,centlon,horizon)

% if tesselation

%  ridgepack_e3smsatmeshs(ncvert,centlat,centlon,horizon,altitude);

% else

%  ridgepack_e3smsatmeshv(nccell,centlat,centlon,horizon,altitude);

% end

% if setting==2 | setting==3

%   ridgepack_sathorizon(centlat,centlon,horizon,...
%                        satlat,satlon,sathorizon,[0.83 0.5 0])

% end

 %cont=[-10 0 10 25 50 100 250 500:500:5000];
 cont=[-8000:500:-1500 -1000:100:-200 -100:10:-10 0 10];
 ref=0;
 lighting=true;
 ncvert.bottomDepth.data=-ncvert.bottomDepth.data;

% ridgepack_e3smsatcol(ncvert,'bottomDepth',ncvert,cont,0,...
%                              centlat,centlon,horizon,altitude,...
%                              lighting);

 ridgepack_e3smsatrough(ncvert,'bottomDepth',ncvert,cont,0,...
                              centlat,centlon,horizon,altitude,...
                              lighting);

 ridgepack_e3smsatcoast(ncvert,centlat,centlon,horizon)

 

elseif setting==7

 ridgepack_polarm('centralarctic','grid','label','noland')

 ridgepack_e3smeshs(ncvert)

 ridgepack_e3smcoastm(ncvert)

 titled='Arctic Ocean';

end

if strcmp(fileg,'CUSP14')
 title([titled,' CUSP 14-60km',])
else
 title([titled,' CUSP 12-60km',])
end


if tesselation

 ridgepack_fprint('png',[fileg,'_tesselation_',num2str(setting)],1,2)

else

 ridgepack_fprint('png',[fileg,'_triangulation_',num2str(setting)],1,2)

end

end












