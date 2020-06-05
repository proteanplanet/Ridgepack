function [cells]=ridgepack_e3smsatmeshv(nccell,centlat,centlon,...
                                        horizon,altitude,ncmask,var)

% apply universal mask if non supplied
if nargin<6
 ncmask=nccell;
 var='nVertices';
end

% height of grid superimposed on plot
gridheight=1.01; 

% reduce the data use to the plotting area to speed things up
% and find plotting edge limit of cells
maxth=deg2rad(horizon);
maxidx=nccell.vertexDegree.data(end);
for i=1:length(nccell.nVertices.data)

 idx=nccell.cellsOnVertex.data(1:maxidx,i);
 idx=idx(idx>0);

 la=nccell.latCell.data(idx);
 lo=nccell.lonCell.data(idx);

 [x,y,z,ph,th]=ridgepack_satfwd(rad2deg(la),rad2deg(lo),...
                                centlat,centlon,...
                                horizon,altitude);

 % filter cells not in frame, and find cropping limit
 if all(isnan(x)) 
  ncmask.(var).data(i)=NaN;
 elseif any(isnan(x)) & ~all(isnan(x))
  [x,y,z,ph,th]=ridgepack_satfwd(rad2deg(la),rad2deg(lo),...
                     centlat,centlon,horizon,altitude,false);
  maxt=max(th(:));
  maxth=max(maxth,maxt);
 end

end

% allocate cell indices included in the horizon
verts=find(~isnan(ncmask.(var).data));
if nargout==1
 disp('Outputing cells within the satellite view')
 return
end

% now draw the grid
if length(verts)>0

 % Mesh triangulation
 SHA=ridgepack_e3smeshv(nccell,verts);

 % Find values within twice the horizon since we've already weeded
 % out cells that only fall within the horizon so as to plot entire 
 % cells, not just the vertices within the horizon.
 [dx,dy,dz,phi,theta]=ridgepack_satfwd([SHA.Latitude],...
                                       [SHA.Longitude],...
                                       centlat,centlon,2*horizon,1);

 % draw grid
 plot3(dx,dy,dz,'Color',[0.7 0.2 0],'LineWidth',0.2)
 hold on

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

if nargout<1
 clear cells
end


