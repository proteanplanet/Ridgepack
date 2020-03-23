function [cells]=ridgepack_psatmeshe3sm(ncvert,centlat,centlon,...
                                        horizon,altitude,ncmask,var)

% apply universal mask if non supplied
if nargin<6
 ncmask=ncvert;
 var='nCells';
end

% height of grid superimposed on plot
gridheight=1.01; 

% reduce the data use to the plotting area to speed things up
% and find plotting edge limit of cells
maxth=deg2rad(horizon);
for i=1:length(ncvert.nCells.data)

 maxidx=ncvert.nEdgesOnCell.data(i);

 la=ncvert.latitude.data(ncvert.verticesOnCell.data(1:maxidx,i));
 lo=ncvert.longitude.data(ncvert.verticesOnCell.data(1:maxidx,i));

 [x,y,z,ph,th]=ridgepack_satfwd(rad2deg(la),rad2deg(lo),...
                                centlat,centlon,horizon,altitude);

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
cells=find(~isnan(ncmask.(var).data));
if nargout==1
 disp('Outputing cells within the satellite view')
 return
end

% now draw the grid
c=[];
d=[];

if length(cells)>0

 for i=1:1:length(cells)

  maxidx=ncvert.nEdgesOnCell.data(cells(i));

  la=ncvert.latitude.data(ncvert.verticesOnCell.data(1:maxidx,cells(i)));
  lo=ncvert.longitude.data(ncvert.verticesOnCell.data(1:maxidx,cells(i)));

  la(end+1)=la(1);
  lo(end+1)=lo(1);

  la(end+1)=NaN;
  lo(end+1)=NaN;

  c=[c;la];
  d=[d;lo];

 end

 % find values within twice the horizon since we've already wedded
 % out cells that only fall within the horizon so as to plot entire 
 % cells, not just the vertices within the horizon.
 [dx,dy,dz,phi,theta]=ridgepack_satfwd(rad2deg(c),rad2deg(d),...
                                  centlat,centlon,2*horizon,1);

 % draw grid
 plot3(dx,dy,dz,'b','LineWidth',0.5)
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


