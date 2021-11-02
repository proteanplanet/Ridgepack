
close all
clear

% legend
%leg=false;
leg=true;

tag='draft2';

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
sector{1}.centlon=-100; % degrees east
sector{1}.horizon=39; % degrees of satellite horizon (0-90)
sector{1}.altitude=1; % Mean Earth radius multiple
sector{1}.annotation=1; % add Arctic Shipping

nsr061=[1603 2166 3012];
nsr061string={'Sannikov$\rightarrow$','Vilkitsky$\rightarrow$','Kara$\rightarrow$'}
nsr061angle=[125 95 95]

nsr062=[1603 2166 3012 870];
nsr062string={'Sannikov$\rightarrow$','Vilkitsky$\rightarrow$','Kara$\rightarrow$','Long$\rightarrow$'}
nsr062angle=[125 95 95 145]

nsr063=[2166 3012];
nsr063string={'Vilkitsky$\rightarrow$','Kara$\rightarrow$'}
nsr063angle=[95 95]

% plot location
%plotloc='/Users/afroberts/work';
plotloc='/Users/afroberts/data/SHIPPING';

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
 ridgepack_write(nccoast,coastname);
else
 nccoast=ridgepack_clone(coastname);
end

setting=1;

notation='abc';

maxsets=3;
minsets=1;

for sets=minsets:maxsets

if sets==1
 centlat=90;   % degrees north
 centlon=-80;   % degrees east
 horizon=35;
% centlat=65;   % degrees north
% centlon=-45   % degrees east
% horizon=30;
elseif sets==2
 centlat=90;   % degrees north
 centlon=-80;   % degrees east
 horizon=35;
% centlat=65;   % degrees north
% centlon=-45;   % degrees east
% horizon=30;
else
 centlat=90;   % degrees north
 centlon=-80;   % degrees east
 horizon=35;
end
altitude=sector{setting}.altitude; % Mean Earth radius multiple

if maxsets>minsets
 ridgepack_multiplot(1,maxsets-minsets+1,1,sets-minsets+1,notation(sets-minsets+1));
else
 ridgepack_multiplot(1,maxsets-minsets+1,1,sets-minsets+1);
end

ridgepack_satview(centlat,centlon,horizon);

% reverse colorbar
%cont=[-5000:500:-1500 -1000:250:-250 -100 -50:10:10];
cont=[-5000 0 10];
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
%   ridgepack_colorbar(cont,'m','linear','vertical',0,colbarcont)
%   clear colbarcont

ridgepack_e3smsatcoast(nccoast,centlat,centlon,horizon,0.85*[1 1 1],0.1);

if sets==1
    ships=[5 12 13 14 17 18];
    complement=[6 11 15 16];
    mainroute=[];
    cols=colormap(lines(length(ships))); 
    titletext='Wrangel Island Routes';
elseif sets==2
    ships=[6 11 15 16];
    complement=[5 12 13 14];
    mainroute=[];
    cols=colormap(lines(length(ships))); 
    titletext='Proliv Longa Routes';
elseif sets==3
    ships=[20];
    complement=[];
    mainroute=[];
    cols=colormap(lines(length(ships))); 
    titletext='Optimal Kara Sea Route';
end

shiploc='/Users/afroberts/data/SHIPPING';
cd(shiploc)

for shipi=1:length(ships)
     ncship=ridgepack_clone('InteRFACE_Shiptracks',...
               {['track',num2str(ships(shipi),'%2.2i'),'_longitude'],...
                ['track',num2str(ships(shipi),'%2.2i'),'_latitude']});
     shipl(shipi)=length(ncship.latitude.data);
end

[shipz,Iorder]=sort(shipl,'descend');

% plot all ship tracks
k=0;
for shipi=Iorder
     k=k+1;
     ncship=ridgepack_clone('InteRFACE_Shiptracks',...
               {['track',num2str(ships(shipi),'%2.2i'),'_longitude'],...
                ['track',num2str(ships(shipi),'%2.2i'),'_latitude'],...
                ['track',num2str(ships(shipi),'%2.2i'),'_distance']});
     [x,y,z,phi,theta]=...
      ridgepack_satfwd(ncship.latitude.data,ncship.longitude.data,...
                       centlat,centlon,horizon,1.001*altitude);
     hship(k)=plot3(x,y,z,':','Color',cols(k,:));
     if shipi<=length(complement)
      shipleg{k}=['Track ',num2str(ships(shipi),'%2.2i'),': ',...
                 num2str(length(ncship.latitude.data)),' NM [',...
                 num2str(complement(shipi),'%2.2i'),']'];
     else
      shipleg{k}=['Track ',num2str(ships(shipi),'%2.2i'),': ',...
                 num2str(length(ncship.latitude.data)),' NM'];
     end
end


% add annotations
shiptrack=6;
ncship=ridgepack_clone('InteRFACE_Shiptracks',...
               {['track',num2str(shiptrack,'%2.2i'),'_longitude'],...
                ['track',num2str(shiptrack,'%2.2i'),'_latitude'],...
                ['track',num2str(shiptrack,'%2.2i'),'_distance']});

[x,y,z,phi,theta]=ridgepack_satfwd(ncship.latitude.data,ncship.longitude.data,...
                       centlat,centlon,horizon,1.001*altitude);
if sets==1
 nsr06=nsr061;
 nsr06string=nsr061string;
 nsr06angle=nsr061angle;
elseif sets==2
 nsr06=nsr062;
 nsr06string=nsr062string;
 nsr06angle=nsr062angle;
elseif sets==3
 nsr06=nsr063;
 nsr06string=nsr063string;
 nsr06angle=nsr063angle;
end

for nsk=1:length(nsr06)
       text(x(nsr06(nsk)),y(nsr06(nsk)),1.05*z(nsr06(nsk)),...
          char(nsr06string{nsk}),...
          'FontSize',5,...
          'Margin',1,...
          'Color',0.25*[1 1 1],...
          'Rotation',nsr06angle(nsk)+180,...
          'HorizontalAlignment','right',...
          'VerticalAlignment','middle');
       disti=ncship.(['track',num2str(shiptrack,'%2.2i'),'_distance']).data(nsr06(nsk))
       disti*1.852
end

% mark in sabetta
[x,y,z,phi,theta]=ridgepack_satfwd(71.2733,72.0725,...
                       centlat,centlon,horizon,1.001*altitude);
plot3(x,y,z,'b.')

% plot most likely ship track
if ~isempty(mainroute)
     ncship=ridgepack_clone('InteRFACE_Shiptracks',...
               {['track',num2str(mainroute,'%2.2i'),'_longitude'],...
                ['track',num2str(mainroute,'%2.2i'),'_latitude']});
     [x,y,z,phi,theta]=...
      ridgepack_satfwd(ncship.latitude.data,ncship.longitude.data,...
                       centlat,centlon,horizon,1.001*altitude);
     plot3(x,y,z,'r-')
     disp('-------')
     disp(ncship.latitude.long_name);
     disp([num2str(length(ncship.latitude.data)),' NM']);
     disp('-------')
end


if leg
 legend(hship,shipleg,'Location','SouthOutside')
 legend('boxoff')
else
 for shipi=Iorder
  field=['track',num2str(ships(shipi),'%2.2i'),'_latitude'];
  ncship=ridgepack_clone('InteRFACE_Shiptracks',{field});
  idx=strfind(ncship.(field).long_name,',');
  disp([ncship.(field).long_name(idx(1)+2:end),': ',char(shipleg{k})]);
 end 
end

title(titletext)

clear shipz Iorder shipl ships mainroute cols titletext shipleg hship

end

ridgepack_multialign(gcf)

cd(plotloc)

ridgepack_fprint('png',['PMs_Shipping_Figure_',tag],1,2)

