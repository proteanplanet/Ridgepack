function ridgepack_e3smstream(ncu,varu,ncv,varv,ncc,varc,...
                              ncvert,nccell,hemisphere,density,...
                              cont,ref)

% switch on ice speed colors
speedcolor=true;
%speedcolor=false;

% set contour colormap, with grey between 10 and 15
[cmap]=ridgepack_colormap(cont,ref);

u=ncu.(varu).data;
v=ncv.(varv).data;
lat=nccell.latitude.data;
lon=nccell.longitude.data;

% create masks to use only all area north of 40N (or south of 40S)
% mask at 10% concentration for calculating streamlines
maskc=find(hemisphere*nccell.latitude.data>40*pi/180 & ...
           ncc.(varc).data>0.01);
maskv=NaN*zeros([ncvert.maxEdges.data(end) length(maskc)]);
for i=1:length(maskc)
  maskv(1:ncvert.nEdgesOnCell.data(maskc(i)),i)=...
           ncvert.verticesOnCell.data(...
           1:ncvert.nEdgesOnCell.data(maskc(i)),maskc(i));
end
maskv=unique(sort(maskv(~isnan(maskv(:)))));

% project onto a polar stereographic projection
[X,Y]=ridgepack_geodetictoxy(ncvert.latitude.data(maskv)*180/pi,...
                             ncvert.longitude.data(maskv)*180/pi,...
                             hemisphere);
[Xc,Yc]=ridgepack_geodetictoxy(nccell.latitude.data(maskc)*180/pi,...
                               nccell.longitude.data(maskc)*180/pi,...
                               hemisphere);
Xi=Xc; Yi=Yc;

U=ncu.(varu).data(maskv);
V=ncv.(varv).data(maskv);
CONC=ncc.(varc).data(maskc);
Speed=sqrt(ncu.(varu).data.^2+ncv.(varv).data.^2);

% change mask to 15% concentration for speed 
maskc=find(hemisphere*nccell.latitude.data>40*pi/180 & ...
           ncc.(varc).data>0.15);

% generate speed on native grid, interpolating from vertices
% (N.B. this needs to be improved with intra-cell weightings)
Speedc=NaN*zeros([length(maskc) 1]);
Latc=NaN*zeros([length(maskc) 8]);
Lonc=NaN*zeros([length(maskc) 8]);
for i=1:length(maskc)
 maxidx=ncvert.nEdgesOnCell.data(maskc(i));
 maxv=ncvert.verticesOnCell.data(1:maxidx,maskc(i));
 Speedc(i)=100*sum(Speed(maxv))./maxidx;
 Latc(i,1:maxidx)=ncvert.latitude.data(maxv)*180/pi;
 Lonc(i,1:maxidx)=ncvert.longitude.data(maxv)*180/pi;
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

if hemisphere==1

 ridgepack_polarm('seaice','grid','label','noland')

elseif hemisphere==-1

 ridgepack_polarm('antarctic','grid','label','noland')

else

 error('hemisphere needs to be either 1 (north) or -1 (south)')

end

xlabel(['Median: ',...
       num2str(0.01*median(Speedc(~isnan(Speedc))),'%4.2f'),...
       '~m~s$^{-1}$'],'FontSize',10,'Color',0.25*[1 1 1])

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
nccoast=ridgepack_psatcoaste3sm(ncvert,false);
[c,d] = mfwdtran(gcm,nccoast.latitude.data,nccoast.longitude.data);
h1=line(c,d,'Color',0.5*[1 1 1]);
clear c d 

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

% plot ice edge
%nc=ridgepack_pthresholde3sm(ncc,varc,0.15,ncvert);
%[c,d]=mfwdtran(gcm,nc.latitude.data,nc.longitude.data,gca,'surface');
%plot(c,d,'Color','m')

% colorbar
if speedcolor 
 ridgepack_colorbar(cont,{'\times 10^{-2}','m~s^{-1}'});
end



