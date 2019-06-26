function ridgepack_psatcole3sm(nc,var,ncvert,nccell,cont,ref,...
                               centlat,centlon,horizon,altitude)

% define colormap
[cmap]=ridgepack_colormap(cont,20);

% reduce the data use to the plotting area to speed things up
% and fine plotting edge limit of cells
maxth=deg2rad(horizon);
for i=1:length(nccell.latitude.data)

 maxidx=ncvert.nEdgesOnCell.data(i);

 la=ncvert.latitude.data(ncvert.verticesOnCell.data(1:maxidx,i));
 lo=ncvert.longitude.data(ncvert.verticesOnCell.data(1:maxidx,i));

 [x,y,z,ph,th]=ridgepack_satfwd(rad2deg(la),rad2deg(lo),...
                                centlat,centlon,horizon,altitude);

 % filter cells no in frame, and find cropping limit
 if all(isnan(x)) 
  nc.(var).data(i)=NaN;
 elseif any(isnan(x)) & ~all(isnan(x))
  [x,y,z,ph,th]=ridgepack_satfwd(rad2deg(la),rad2deg(lo),...
                     centlat,centlon,horizon,altitude,false);
  maxt=max(th(:));
  maxth=max(maxth,maxt);
 end

end

% now shade the regions
for j=1:length(cont)-1

 if j<length(cont)-1
  id=find(nc.(var).data>cont(j) & nc.(var).data<=cont(j+1));
 else
  id=find(nc.(var).data>cont(j));
 end

 [zindex,truecol]=ridgepack_colorindex(nc.(var).data(id),cont,20);

 if length(id)>0

  lat=zeros([1 8]);
  lon=zeros([1 8]);
  xl=zeros(length(id),8);
  yl=zeros(length(id),8);
  zl=zeros(length(id),8);
  phl=zeros(length(id),8);
  thl=zeros(length(id),8);
  cl=zeros(3,length(id));

  for i=1:1:length(id)

   maxidx=ncvert.nEdgesOnCell.data(id(i));

   la=ncvert.latitude.data(ncvert.verticesOnCell.data(1:maxidx,id(i)));
   lo=ncvert.longitude.data(ncvert.verticesOnCell.data(1:maxidx,id(i)));

   lat(1:maxidx)=la;
   lon(1:maxidx)=lo;

   lat(maxidx+1:8)=la(1);
   lon(maxidx+1:8)=lo(1);

   [xl(i,:),yl(i,:),zl(i,:),phl(i,:),thl(i,:)]=...
    ridgepack_satfwd(rad2deg(squeeze(lat(:))),...
                     rad2deg(squeeze(lon(:))),...
                     centlat,centlon,2*horizon,altitude);

  end

  patch(xl',yl',zl',truecol(1,:),'EdgeColor','none')
  hold on
  drawnow

 end

 clear zindex truecolor xl yl zl phl thl lon lat

end

% crop by overlaying a white ring
N = 100;
thetavec = linspace(deg2rad(horizon),maxth,N);
phivec = linspace(0,2*pi,N);
[th, ph] = meshgrid(thetavec,phivec);
R = ones(size(th)); % should be your R(theta,phi) surface in general
cx = R.*sin(th).*cos(ph);
cy = R.*sin(th).*sin(ph);
cz = 1.09*R.*cos(th);
c1 = ones(size(cx));
clear cc
cc(:,:,1)=c1;
cc(:,:,2)=c1;
cc(:,:,3)=c1;
surf(cx,cy,cz,cc,'EdgeColor','none');

% add black frame
ph=deg2rad([0:0.001:361]);
th=deg2rad(horizon*ones(size(ph)));
R = ones(size(th)); % should be your R(theta,phi) surface in general
cx = R.*sin(th).*cos(ph);
cy = R.*sin(th).*sin(ph);
cz = 1.095*R.*cos(th);
plot3(cx,cy,cz,'k')

% make axes equal and tight, set viewing angle
axis equal
view([0 0 0.4])
axis tight
axis off

% add lighting from infinite sources directly overhead
hl=light('Position',[0 0 10000],'Style','local')
material dull


