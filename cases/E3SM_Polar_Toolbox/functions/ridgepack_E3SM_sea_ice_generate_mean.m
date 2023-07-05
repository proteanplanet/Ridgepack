
clear
close all

cases={'v2.NARRM.historical_0101'}
labels={'NARRM 0101'}

fields={'timeDaily_avg_uVelocityGeo','timeDaily_avg_vVelocityGeo'};

startyear=1980;
endyear=2014;

for j=1:length(cases)

 casen=char(cases{j});

 cd(['/Users/afroberts/data/MODEL/E3SM/v2/',casen,'/data/ice/hist'])


% get velocity on native grid
for i=minmonth:maxmonth

 cd(['/Users/afroberts/data/MODEL/E3SM/v2/',casenow,'/data/ice/hist'])
 infile=['mpassi.hist.am.',casenow,'.state.mean.',...
                      num2str(i,'%2.2i'),'.',...
                      num2str(startyear,'%4.4i'),'_',...
                      num2str(endyear,'%4.4i'),'.nc'];

 if i==minmonth
  nc=ridgepack_clone(infile,...
                   {'timeMonthly_avg_uVelocityGeo',...
                    'timeMonthly_avg_vVelocityGeo',...
                    'timeMonthly_avg_iceAreaCell'});
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

% create masks to use only all area north of 40N (or south of 40S)
% mask at 10% concentration for calculating streamlines
maskc=find(hemisphere*nclatc.latCell.data>40*pi/180 & ...
            nc.timeMonthly_avg_iceAreaCell.data>0.01);
maskv=NaN*zeros([7 length(maskc)]);
for i=1:length(maskc)
  maskv(1:ncvonc.nEdgesOnCell.data(maskc(i)),i)=...
           ncvonc.verticesOnCell.data(...
           1:ncvonc.nEdgesOnCell.data(maskc(i)),maskc(i));
end
maskv=unique(sort(maskv(~isnan(maskv(:)))));

% project onto a polar stereographic projection
[X,Y]=ridgepack_geodetictoxy(nclatv.latVertex.data(maskv)*180/pi,...
                             nclonv.lonVertex.data(maskv)*180/pi,...
                             hemisphere);
[Xc,Yc]=ridgepack_geodetictoxy(nclatc.latCell.data(maskc)*180/pi,...
                               nclonc.lonCell.data(maskc)*180/pi,...
                               hemisphere);
Xi=Xc; Yi=Yc;

U=nc.timeMonthly_avg_uVelocityGeo.data(maskv);
V=nc.timeMonthly_avg_vVelocityGeo.data(maskv);
CONC=nc.timeMonthly_avg_iceAreaCell.data(maskc);
Speed=sqrt(nc.timeMonthly_avg_uVelocityGeo.data.^2+...
           nc.timeMonthly_avg_vVelocityGeo.data.^2);

% change mask to 15% concentration for speed 
maskc=find(hemisphere*nclatc.latCell.data>40*pi/180 & ...
           nc.timeMonthly_avg_iceAreaCell.data>0.15);

% generate speed on native grid, interpolating from vertices
% (N.B. this needs to be improved with intra-cell weightings)
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
  ridgepack_polarm('seaice','noland','grid')
 end
 textm(42,-10,[num2str(0.01*median(Speedc(~isnan(Speedc))),'%4.2f'),...
               '~m~s$^{-1}$'],'FontSize',7,'Color','b',...
               'HorizontalAlignment','right')
else
 if cols==1 & rows==1
  ridgepack_polarm('antarctic','grid','label','noland')
 else
  ridgepack_polarm('antarctic','noland','grid')
 end
 textm(-87,45,[num2str(0.01*median(Speedc(~isnan(Speedc))),'%4.2f'),...
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
[c,d] = mfwdtran(gcm,[nccoast.latitude.data],[nccoast.longitude.data]);
h1=line(c,d,'Color',0.5*[1 1 1]);
clear c d 

% plot streamlines on the map
for i=1:length(vertices)
 if ~isempty(vertices{i})
  xy=vertices{i};
  [lat,lon]=ridgepack_xytogeodetic(xy(:,1),xy(:,2),hemisphere);
  [c,d]=mfwdtran(gcm,lat,lon,gca,'surface');
  if lk==2
   plot(c,d,'Color',0*[1 1 1],'LineWidth',0.4)
  else
   plot(c,d,'Color',0*[1 1 1],'LineWidth',0.4)
  end
 end
end

% plot streamline arrow on the map
for i=1:length(arrowvertices)
 if ~isempty(arrowvertices{i})
  xy=arrowvertices{i};
  [lat,lon]=ridgepack_xytogeodetic(xy(:,1),xy(:,2),hemisphere);
  [c,d]=mfwdtran(gcm,lat,lon,gca,'surface');
  if lk==2
   plot(c,d,'Color',0*[1 1 1],'LineWidth',0.4)
  else
   plot(c,d,'Color',0*[1 1 1],'LineWidth',0.4)
  end
 end
end

% colorbar
if speedcolor & rows==1 & cols==1
 ridgepack_colorbar(cont,{'\times 10^{-2}','m~s^{-1}'});
 ridgepack_multicb(gca);
end

% titles
if rows==1 
 title(char(headers{cols}),'FontSize',10);
end
if cols==1
 ylabel(char(labels{rows}),'FontSize',10);
end

end
end

% line up the plots
ridgepack_multialign(gcf)

cd('/Users/afroberts/work')

if hemisphere==1
 nams='.north.streamline.';
else
 nams='.south.streamline.';
end

ridgepack_fprint('png',[ridgepack_cellcat(cases,'_'),nams,num2str(startyear,'%4.4i'),'_',num2str(endyear,'%4.4i'),'.png'],1,1)

end % lk





