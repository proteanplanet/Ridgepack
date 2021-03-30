
clf
clear

gridchoice=2;

fileg{1}.name='WC12r01';
fileg{1}.outname='WC12';
fileg{1}.title=' WC 12-60~km mesh';
fileg{2}.name='WC14r03';
fileg{2}.outname='WC14r03';
fileg{2}.title=' WC 14-60~km mesh';
fileg{3}.name='DECK';
fileg{3}.outname='DECK';
fileg{3}.title=' DECK 30-60~km standard mesh';


sector{1}.centlat=90; % degrees north
sector{1}.centlon=-90; % degrees east
sector{1}.horizon=40; % degrees of satellite horizon (0-90)
sector{1}.altitude=1; % Mean Earth radius multiple
sector{1}.annotation=1; % add Arctic Shipping

% plot location
plotloc='/Users/afroberts/work';

% grid location
if strcmp(char(fileg{gridchoice}.name),'DECK')
 gridloc=['/Users/afroberts/data/MODEL/E3SM/',...
          char(fileg{gridchoice}.name),'/grid'];
 gridfile='E3SM_LR_V1_grid.nc';
elseif strcmp(char(fileg{gridchoice}.name),'WC12')
 gridloc=['/Users/afroberts/data/MODEL/E3SM/WC12/',...
          char(fileg{gridchoice}.name)];
 gridfile='initial_state.nc';
elseif strcmp(char(fileg{gridchoice}.name),'WC14r03')
 gridloc=['/Users/afroberts/data/MODEL/E3SM/WC14/r03'];
 gridfile='initial_state.nc';
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

% invert bathymetry
ncvert.bottomDepth.data=-ncvert.bottomDepth.data;

% read in coast or else generate coast and write it out
coastname=[char(fileg{1}.outname),'_Coast'];
x=dir(coastname);
if isempty(x)
 nccoast=ridgepack_e3smcoastm(ncvert);
 ridgepack_write(nccoast,coastname)
else
 nccoast=ridgepack_clone(coastname);
end

shiploc='/Users/afroberts/data/SHIPPING';
cd(shiploc)

% get the tracks
%shipl=zeros([1 19]);
shipl=zeros([1 1]);
for shipi=1:length(shipl)
     trackname=['track',num2str(shipi,'%2.2i')];
     ncship=ridgepack_clone('InteRFACE_Shiptracks',...
               {[trackname,'_longitude'],[trackname,'_latitude']});
     [cell,oncell,celllist,oncelllist,celldist]=...
              ridgepack_e3smontrack(ncvert,...
                       ncship.latitude.data,ncship.longitude.data);
end

%legend(hship,shipleg,'Location','north')
%legend('boxoff')

%title(titletext)

%cd(plotloc)

%ridgepack_fprint('png',['PMs_Shipping_Figure_',num2str(sets)],1,2)

