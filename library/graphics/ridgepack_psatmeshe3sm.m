function [cells]=ridgepack_psatmeshe3sm(nc,var,ncvert,...
                      cont,ref,centlat,centlon,horizon,altitude)

gridheight=1.01; % height of grid superimposed on plot

% reduce the data use to the plotting area to speed things up
% and find plotting edge limit of cells
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

% allocate cell indices included in the horizon
if nargout==1
 disp('Outputing cells within the satellite view')
 cells=find(~isnan(nc.(var).data));
 return
end

% now draw the grid
for j=1:length(cont)-1

 if j<length(cont)-1
  id=find(nc.(var).data>cont(j) & nc.(var).data<=cont(j+1));
 else
  id=find(nc.(var).data>cont(j));
 end

 c=[];
 d=[];

 if length(id)>0

  for i=1:1:length(id)

   maxidx=ncvert.nEdgesOnCell.data(id(i));

   la=ncvert.latitude.data(ncvert.verticesOnCell.data(1:maxidx,id(i)));
   lo=ncvert.longitude.data(ncvert.verticesOnCell.data(1:maxidx,id(i)));

   la(end+1)=la(1);
   lo(end+1)=lo(1);

   la(end+1)=NaN;
   lo(end+1)=NaN;

   c=[c;la];
   d=[d;lo];

  end

  % find values within twice the horizon since we've already
  % weeded out cells that only fall within the horizon
  [dx,dy,dz,phi,theta]=ridgepack_satfwd(rad2deg(c),rad2deg(d),...
                                  centlat,centlon,2*horizon,1);

  % draw grid
  plot3(dx,dy,dz,'b','LineWidth',0.5)
  hold on

 end

end

% crop by overlaying a white ring
N = 100;
thetavec = linspace(deg2rad(horizon),maxth,N);
phivec = linspace(0,2*pi,N);
[th, ph] = meshgrid(thetavec,phivec);
R = ones(size(th)); % should be your R(theta,phi) surface in general
cx = R.*sin(th).*cos(ph);
cy = R.*sin(th).*sin(ph);
cz = gridheight*R.*cos(th);
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
cz = gridheight*R.*cos(th);
plot3(cx,cy,cz,'k')

% make axes equal and tight, set viewing angle
axis equal
view([0 0 0.4])
axis tight
axis off

% add lighting from infinite sources directly overhead
light('Position',[0 0 10000],'Style','local')
material dull


