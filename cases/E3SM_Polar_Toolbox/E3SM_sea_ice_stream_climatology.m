
clear

%cases={'v0-HR','v1-HR','v1-LR'}
%labels={'v0 HR','v1 HR','v1 LR'}

cases={'v1-HR','v1-LR'}
labels={'v1 HR','v1 LR'}

%cases={'v0-HR'}
%labels={'v0 HR'}

%cases={'v1-HR'}
%labels={'v1 HR'}

%cases={'v1-LR'}
%labels={'v1 LR'}

% set this to change location of data
laptop=true;
%laptop=false;

% density of streamlines
density=7;

% switch on ice speed colors
speedcolor=true;
%speedcolor=false;

% set up size of plot
nrows=length(cases);

% set columns for each quarter
ncols=4;

% alphabetic counter
alpha='abcdefghijklmnopqrstuvwxyz';

% progress through arctic and antarctic
for lk=2:-1:1

% close figure window
close all

% set contour colormap, with grey between 10 and 15
cont=[0:2:10 15 20:10:60];
[cmap]=ridgepack_colormap(cont,12.5);

% alaphabetic counter
ma=0;

if lk==1
 hemisphere=1;
else
 hemisphere=-1;
end

for rows=1:nrows

casenow=char(cases{rows});

% get latitude and longitude
if strcmp(casenow,'v1-LR')

 if laptop
  cd('/Volumes/MacBookProT3SeaIce/E3SM/lrhrequiv/grid')
 else
  cd('/Users/afroberts/SIhMSatArray/E3SM/DECK/grid')
 end
 gridfile='oEC60to30v3_60layer.restartFrom_anvil0926.171101.nc';
 landm = shaperead('E3SM_V1_C_grid_Coast.shp',...
                   'UseGeoCoords',true);
elseif strcmp(casenow,'v1-HR')

 if laptop
  cd('/Volumes/MacBookProT3SeaIce/E3SM/highres/grid')
 else 
  cd('/Users/afroberts/SIhMSatArray/E3SM/highres/grid')
 end
 gridfile='gridFieldsRRS18to6v3.nc';
 landm = shaperead('E3SM_HR_V1_C_grid_Coast.shp',...
                   'UseGeoCoords',true);

elseif strcmp(casenow,'v0-HR')

 if laptop
  cd('/Volumes/MacBookProT3SeaIce/E3SM/highresv0/grid')
 else 
  cd('/Users/afroberts/SIhMSatArray/E3SM/highresv0/grid')
 end
 gridfile='hi.b1850c5_acmev0_highres.cice.h.0130-12.nc';

end

if strcmp(casenow,'v0-HR')

 nclatv=ridgepack_clone(gridfile,'latt_bounds');
 nclatv.latVertex=nclatv.latt_bounds;
 nclatv=rmfield(nclatv,'latt_bounds');
 nclonv=ridgepack_clone(gridfile,'lont_bounds');
 nclonv.lonVertex=nclonv.lont_bounds;
 nclonv=rmfield(nclonv,'lont_bounds');
 nclatc=ridgepack_clone(gridfile,'TLAT');
 nclatc.latCell=nclatc.TLAT;
 nclatc=rmfield(nclatc,'TLAT');
 nclonc=ridgepack_clone(gridfile,'TLON');
 nclonc.lonCell=nclonc.TLON;
 nclonc=rmfield(nclonc,'TLON');

else

 nclatv=ridgepack_clone(gridfile,'latVertex');
 nclonv=ridgepack_clone(gridfile,'lonVertex');
 nclatc=ridgepack_clone(gridfile,'latCell');
 nclonc=ridgepack_clone(gridfile,'lonCell');
 ncvonc=ridgepack_clone(gridfile,{'verticesOnCell','nEdgesOnCell'});

end

for cols=1:ncols

% counter for alphabetic labels
ma=ma+1;

% set seasonal boundaries
headers={'Jan-Mar','Apr-Jun','Jul-Sep','Oct-Dec'};
if cols==1
 minmonth=1;maxmonth=3;
elseif cols==2
 minmonth=4;maxmonth=6;
elseif cols==3
 minmonth=7;maxmonth=9;
elseif cols==4
 minmonth=10;maxmonth=12;
end

% get velocity on native grid
for i=minmonth:maxmonth
 if strcmp(casenow,'v1-LR')
  if laptop
   cd('/Volumes/MacBookProT3SeaIce/E3SM/lrhrequiv/monthly')
  else
   cd('/Users/afroberts/SIhMSatArray/E3SM/DECK/monthly/h1/archive/ice/reduced')
  end
  infile=['mpascice.hist.am.LR_V1.state.mean.',...
           num2str(i,'%2.2i'),'.0026_0055.nc'];
 elseif strcmp(casenow,'v1-HR')
  if laptop
   cd('/Volumes/MacBookProT3SeaIce/E3SM/highres/monthly')
  else
   cd('/Users/afroberts/SIhMSatArray/E3SM/highres/monthly')
  end
  infile=['mpascice.hist.am.HR_V1.state.mean.',...
           num2str(i,'%2.2i'),'.0026_0055.nc'];
 elseif strcmp(casenow,'v0-HR')
  if laptop
   cd('/Volumes/MacBookProT3SeaIce/E3SM/highresv0/monthly')
  else
   cd('/Users/afroberts/SIhMSatArray/E3SM/highresv0/monthly')
  end
  infile=['b1850c5_acmev0_highres.HR_v0.state.mean.',...
           num2str(i,'%2.2i'),'.0031_0060.nc'];
 end

 if i==minmonth
  if strcmp(casenow,'v0-HR')
   nc=ridgepack_clone(infile,...
          {'uvel','vvel','aice','latitude','longitude'});
   nc.timeMonthly_avg_uVelocityGeo=nc.uvel;
   nc.timeMonthly_avg_vVelocityGeo=nc.vvel;
   nc.aice.data=nc.aice.data/100;
   nc.timeMonthly_avg_iceAreaCell=nc.aice;
   nc=rmfield(nc,{'uvel','vvel','aice'});
  else
   nc=ridgepack_clone(infile,...
                   {'timeMonthly_avg_uVelocityGeo',...
                    'timeMonthly_avg_vVelocityGeo',...
                    'timeMonthly_avg_iceAreaCell'});
  end
 else
  if strcmp(casenow,'v0-HR')
   ncnew=ridgepack_clone(infile,...
          {'uvel','vvel','aice','latitude','longitude'});
   nc.timeMonthly_avg_uVelocityGeo.data=...
         nc.timeMonthly_avg_uVelocityGeo.data+ncnew.uvel.data;
   nc.timeMonthly_avg_vVelocityGeo.data=...
         nc.timeMonthly_avg_vVelocityGeo.data+ncnew.vvel.data;
   ncnew.aice.data=ncnew.aice.data/100;
   nc.timeMonthly_avg_iceAreaCell.data=...
          nc.timeMonthly_avg_iceAreaCell.data+ncnew.aice.data;
   clear ncnew
  else
   ncnew=ridgepack_clone(infile,...
                   {'timeMonthly_avg_uVelocityGeo',...
                    'timeMonthly_avg_vVelocityGeo',...
                    'timeMonthly_avg_iceAreaCell'});
   nc.timeMonthly_avg_uVelocityGeo.data=...
         nc.timeMonthly_avg_uVelocityGeo.data+...
         ncnew.timeMonthly_avg_uVelocityGeo.data;
   nc.timeMonthly_avg_vVelocityGeo.data=...
         nc.timeMonthly_avg_vVelocityGeo.data+...
         ncnew.timeMonthly_avg_vVelocityGeo.data;
   nc.timeMonthly_avg_iceAreaCell.data=...
         nc.timeMonthly_avg_iceAreaCell.data+...
         ncnew.timeMonthly_avg_iceAreaCell.data;
   clear ncnew
  end
 end

end

nc.timeMonthly_avg_uVelocityGeo.data=...
         nc.timeMonthly_avg_uVelocityGeo.data/...
         (maxmonth-minmonth+1);
nc.timeMonthly_avg_vVelocityGeo.data=...
         nc.timeMonthly_avg_vVelocityGeo.data/...
         (maxmonth-minmonth+1);
nc.timeMonthly_avg_iceAreaCell.data=...
         nc.timeMonthly_avg_iceAreaCell.data/...
         (maxmonth-minmonth+1);

% find (approximately correct) turning angle of on tripolar grid
% for graphics purposes only
if strcmp(casenow,'v0-HR')
 u=nc.timeMonthly_avg_uVelocityGeo.data;
 v=nc.timeMonthly_avg_vVelocityGeo.data;
 lat=nclatc.latCell.data;
 lon=nclonc.lonCell.data;

 clear az
 az = azimuth(lat(1:end-1,1:end-1),lon(1:end-1,1:end-1),...
              lat(2:end,1:end-1),lon(2:end,1:end-1));
 az(end+1,end+1)=az(end,end);
 az(end-1,end)=az(end,end);
 az(end,end-1)=az(end,end);

 [th,z]=cart2pol(u,v);
 [uu,vv]=pol2cart(th-deg2rad(az),z);

 uu(isnan(lat))=0;
 vv(isnan(lat))=0;

 nc.timeMonthly_avg_uVelocityGeo.data=uu;
 nc.timeMonthly_avg_vVelocityGeo.data=vv;
end

% create masks to use only all area north of 40N (or south of 40S)
% mask at 10% concentration for calculating streamlines
if strcmp(casenow,'v0-HR')
 maskc=find(hemisphere*nclatc.latCell.data>40);
 maskv=maskc;
else
 maskc=find(hemisphere*nclatc.latCell.data>40*pi/180 & ...
            nc.timeMonthly_avg_iceAreaCell.data>0.01);
 maskv=NaN*zeros([7 length(maskc)]);
 for i=1:length(maskc)
  maskv(1:ncvonc.nEdgesOnCell.data(maskc(i)),i)=...
           ncvonc.verticesOnCell.data(...
           1:ncvonc.nEdgesOnCell.data(maskc(i)),maskc(i));
 end
 maskv=unique(sort(maskv(~isnan(maskv(:)))));
end

% project onto a polar stereographic projection
if strcmp(casenow,'v0-HR')
 [X,Y]=ridgepack_geodetictoxy(nclatc.latCell.data(maskv),...
                             nclonc.lonCell.data(maskv),...
                             hemisphere);
 Xc=X; Yc=Y;
 Xi=X; Yi=Y;
else
 [X,Y]=ridgepack_geodetictoxy(nclatv.latVertex.data(maskv)*180/pi,...
                             nclonv.lonVertex.data(maskv)*180/pi,...
                             hemisphere);
 [Xc,Yc]=ridgepack_geodetictoxy(nclatc.latCell.data(maskc)*180/pi,...
                               nclonc.lonCell.data(maskc)*180/pi,...
                               hemisphere);
 Xi=Xc; Yi=Yc;
end

U=nc.timeMonthly_avg_uVelocityGeo.data(maskv);
V=nc.timeMonthly_avg_vVelocityGeo.data(maskv);
CONC=nc.timeMonthly_avg_iceAreaCell.data(maskc);
Speed=sqrt(nc.timeMonthly_avg_uVelocityGeo.data.^2+...
           nc.timeMonthly_avg_vVelocityGeo.data.^2);

% change mask to 15% concentration for speed 
if strcmp(casenow,'v0-HR')
 maskc=find(hemisphere*nclatc.latCell.data>40 & ...
           nc.timeMonthly_avg_iceAreaCell.data>0.15);
else
 maskc=find(hemisphere*nclatc.latCell.data>40*pi/180 & ...
           nc.timeMonthly_avg_iceAreaCell.data>0.15);
end

% generate speed on native grid, interpolating from vertices
% (N.B. this needs to be improved with intra-cell weightings)
if strcmp(casenow,'v0-HR')
 Speedc=Speed(maskc);
 Latc=NaN*zeros([length(maskc) 5]);
 Lonc=NaN*zeros([length(maskc) 5]);
 for i=1:length(maskc)
  Latc(i,1:4)=nclatv.latVertex.data(1:4,maskc(i))';
  Lonc(i,1:4)=nclonv.lonVertex.data(1:4,maskc(i))';
  Latc(i,5)=Latc(i,1);
  Lonc(i,5)=Lonc(i,1);
 end
else
 Speedc=NaN*zeros([length(maskc) 1]);
 Latc=NaN*zeros([length(maskc) 8]);
 Lonc=NaN*zeros([length(maskc) 8]);
 for i=1:length(maskc)
  maxidx=ncvonc.nEdgesOnCell.data(maskc(i));
  maxv=ncvonc.verticesOnCell.data(1:maxidx,maskc(i));
  Speedc(i)=100*sum(Speed(maxv))./maxidx;
  Latc(i,1:maxidx)=nclatv.latVertex.data(maxv)*180/pi;
  Lonc(i,1:maxidx)=nclonv.lonVertex.data(maxv)*180/pi;
  Latc(i,maxidx+1:8)=Latc(i,1);
  Lonc(i,maxidx+1:8)=Lonc(i,1);
 end
end

% generate polar stereographic x-y grid
if hemisphere==1
 [ncgrid]=ridgepack_gridgen('',10);
else
 [ncgrid]=ridgepack_gridgen('',11);
end
x=ncgrid.x.data;
y=ncgrid.y.data;
X(X<min(x) | X>max(x))=NaN;
Y(Y<min(y) | Y>max(y))=NaN;
maskUVs=(~isnan(X) & ~isnan(Y));
X=X(maskUVs); Y=Y(maskUVs);
U=U(maskUVs); V=V(maskUVs);

Xc(Xc<min(x) | Xc>max(x))=NaN;
Yc(Yc<min(y) | Yc>max(y))=NaN;
maskCs=(~isnan(Xc) & ~isnan(Yc));
Xc=Xc(maskCs); Yc=Yc(maskCs);
CONC=CONC(maskCs);

[x,y]=meshgrid(x,y);

% interpolate data to polar stereographic
method='nearest';
u=griddata(X,Y,U,x,y,method);
v=griddata(X,Y,V,x,y,method);
conc=griddata(Xc,Yc,CONC,x,y,method);
u(conc<0.15)=0;
v(conc<0.15)=0;

% turn vectors from lat-lon to x-y grid
[th,z]=cart2pol(u,v);
[ui,vi]=pol2cart(th+deg2rad(ncgrid.turn.data),z);

% generate streamlines on x-y grid
[vertices arrowvertices]=streamslice(x,y,ui,vi,density);

% set of multiplot
ridgepack_multiplot(nrows,ncols,rows,cols,alpha(ma));

% create map, and get the structure
if hemisphere==1
 if cols==1 & rows==1
  ridgepack_polarm('seaice','grid','label','noland')
 else
  ridgepack_polarm('seaice','noland')
 end
 textm(42,-10,[num2str(0.01*median(Speedc(~isnan(Speedc))),'%4.2f'),...
               '~m~s$^{-1}$'],'FontSize',7,'Color','b',...
               'HorizontalAlignment','right')
else
 if cols==1 & rows==1
  ridgepack_polarm('antarctic','grid','label','noland')
 else
  ridgepack_polarm('antarctic','noland')
 end
 textm(-87,0,[num2str(0.01*median(Speedc(~isnan(Speedc))),'%4.2f'),...
               '~m~s$^{-1}$'],'FontSize',8,'Color','b',...
               'HorizontalAlignment','center')
end

% color the speed of the ice drift
if speedcolor
 for j=1:length(cont)
  if j<length(cont)
   idxn=find(Speedc(:)>cont(j) & Speedc(:)<=cont(j+1));
  else
   idxn=find(Speedc(:)>cont(j));
  end
  if length(idxn)>0
   [zindex,truecolor]=ridgepack_colorindex(Speedc(idxn),cont,5);
   [c,d] = mfwdtran(gcm,Latc(idxn,:),Lonc(idxn,:));
   patch(c',d',truecolor(1,:),'EdgeColor','none')
   drawnow
  else
   disp(['No contours for ',num2str(cont(j))])
  end
 end
end

% draw the native model coast
if strcmp(casenow,'v0-HR')
  nc.mask=nc.timeMonthly_avg_iceAreaCell;
  nc.mask.data(~isnan(nc.timeMonthly_avg_iceAreaCell.data))=1;
  nc.mask.data(isnan(nc.timeMonthly_avg_iceAreaCell.data))=0;

  ridgepack_maskm(nc.latitude.data,...
                  nc.longitude.data,...
                  nc.mask.data)
  drawnow
else
 [c,d] = mfwdtran(gcm,[landm.Lat],[landm.Lon]);
 h1=line(c,d,'Color',0.5*[1 1 1]);
 clear c d 
end

% plot streamlines on the map
for i=1:length(vertices)
 if ~isempty(vertices{i})
  xy=vertices{i};
  [lat,lon]=ridgepack_xytogeodetic(xy(:,1),xy(:,2),hemisphere);
  [c,d]=mfwdtran(gcm,lat,lon,gca,'surface');
  plot(c,d,'Color',0*[1 1 1])
 end
end

% plot streamline arrow on the map
for i=1:length(arrowvertices)
 if ~isempty(arrowvertices{i})
  xy=arrowvertices{i};
  [lat,lon]=ridgepack_xytogeodetic(xy(:,1),xy(:,2),hemisphere);
  [c,d]=mfwdtran(gcm,lat,lon,gca,'surface');
  plot(c,d,'Color',0*[1 1 1])
 end
end

% colorbar
if speedcolor & rows==1 & cols==1
 ridgepack_colorbar(cont,{'\times 10^{-2}','m~s^{-1}'});
 ridgepack_multicb(gca);
end

% titles
if rows==1 
 title(char(headers{cols}),'FontSize',12);
end
if cols==1
 ylabel(char(labels{rows}),'FontSize',12);
end

end
end

% line up the plots
ridgepack_multialign(gcf)

if laptop
 cd('/Volumes/MacBookProT3SeaIce/E3SM/highres/plots')
else 
 cd('/Users/afroberts/SIhMSatArray/E3SM/highres/plots')
end

plotfile='streamline';

if hemisphere==1
 hem='north.';
else
 hem='south.';
end

ridgepack_fprint('png',[hem,plotfile,'.png'],1,2)
ridgepack_fprint('epsc',[hem,plotfile,'.eps'],1,2)

end % lk





