function ridgepack_e3smsatcol(nc,var,ncvert,cont,ref,...
                              centlat,centlon,horizon,altitude,...
                              lighting,reversed,logscale,colors)



% define colormap
[cmap]=ridgepack_colormap(cont,ref,colors,logscale,reversed);

% reduce the data use to the plotting area to speed things up
% and fine plotting edge limit of cells
maxth=deg2rad(horizon);
for i=1:length(ncvert.nCells.data)

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
mavar=max(nc.(var).data(:));

for j=1:length(cont)

 if j==1
  id=find(nc.(var).data<cont(2));
 elseif j==length(cont)
  id=find(nc.(var).data>=cont(j-1));
 else
  id=find(nc.(var).data>=cont(j) & nc.(var).data<cont(j+1));
 end

 [zindex,truecol]=ridgepack_colorindex(nc.(var).data(id),cont,ref);

 % find indices only within the plotting contour range 
 idx=id(~isnan(zindex));
 zind=find(~isnan(zindex));

 if length(idx)>0

  idmax=ncvert.maxEdges.data(end);

  lat=zeros([1 idmax]);
  lon=zeros([1 idmax]);
  xl=zeros(length(idx),idmax);
  yl=zeros(length(idx),idmax);
  zl=zeros(length(idx),idmax);
  phl=zeros(length(idx),idmax);
  thl=zeros(length(idx),idmax);
  cl=zeros(3,length(idx));

  for i=1:1:length(idx)

   maxidx=ncvert.nEdgesOnCell.data(idx(i));

   la=ncvert.latitude.data(ncvert.verticesOnCell.data(1:maxidx,idx(i)));
   lo=ncvert.longitude.data(ncvert.verticesOnCell.data(1:maxidx,idx(i)));

   lat(1:maxidx)=la;
   lon(1:maxidx)=lo;

   lat(maxidx+1:idmax)=la(1);
   lon(maxidx+1:idmax)=lo(1);

   [xl(i,:),yl(i,:),zl(i,:),phl(i,:),thl(i,:)]=...
    ridgepack_satfwd(rad2deg(squeeze(lat(:))),...
                     rad2deg(squeeze(lon(:))),...
                     centlat,centlon,2*horizon,altitude);

%   thealtitude=(1-0.05*(nc.(var).data(idx(i))./mavar))*altitude;
%
%   [xl(i,:),yl(i,:),zl(i,:),phl(i,:),thl(i,:)]=...
%    ridgepack_satfwd(rad2deg(squeeze(lat(:))),...
%                     rad2deg(squeeze(lon(:))),...
%                     centlat,centlon,2*horizon,thealtitude);
%                     centlat,centlon,2*horizon,altitude);

  end

  patch(xl',yl',zl',truecol(zind(1),:),'EdgeColor','none')
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
if lighting
 hl=light('Position',[0 0 10000],'Style','local')
 material dull
end

