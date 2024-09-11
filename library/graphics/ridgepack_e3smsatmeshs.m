function [cells]=ridgepack_e3smsatmeshs(ncvert,centlat,centlon,...
                                        horizon,altitude,ncmask,var)

% ridgepack_e3smsatmeshs - plot a tesselated E3SM MPAS mesh on the globe
%
% function [cells]=ridgepack_e3smsatmeshs(ncvert,centlat,centlon,...
%                                         horizon,altitude,ncmask,var)
%
% This function plots the tesselation of an MPAS E3SM mesh on the globe
% for a chosen horizon specified in degrees from tthe center latitude
% and longitude of the center of the plot.
% 
% INPUT:
%
% ncvert   - A netcdf structure containing all of the necessesary grid
%            information from MPAS
% centlat  - The central latitude of the plot
% centlon  - The central longitude of the plot
% horizon  - The reach of the plot, expressed as a polar angle from the
%            central latitude and longitude (degrees)
% altitude - Altitude of the points to be plotted on the sphere
% ncmask   - Cell indices of cells to be plotted [optional]
% var      - variable in ncmask used for masking [optional]
%
% OUTPUT:
%
% cells - indices plotted within the horizon of the figure
%
% Ridgepack Version 2.0
% Andrew Roberts, LANL, 2020 (afroberts@lanl.gov) 

% apply universal mask if non supplied
if nargin<6
 ncmask=ncvert;
 var='nCells';
end

% height of grid superimposed on plot
gridheight=1.01; 

% grid color
gridcolor=[0 0.301 0.435];

% reduce the data use to the plotting area to speed things up
% and find plotting edge limit of cells
maxth=deg2rad(horizon);
for i=1:length(ncvert.nCells.data)

 maxidx=ncvert.nEdgesOnCell.data(i);

 la=ncvert.latVertex.data(ncvert.verticesOnCell.data(1:maxidx,i));
 lo=ncvert.lonVertex.data(ncvert.verticesOnCell.data(1:maxidx,i));

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
if length(cells)>0

 % plot great circles (slower)
 greatcircles=false
 %greatcircles=true

 % Mesh triangulation
 SHA=ridgepack_e3smeshs(ncvert,cells);

 if greatcircles

  % split into 1NM great circle segments
  for i=2:length(SHA.Latitude)

   if SHA.Latitude(i-1)~=SHA.Latitude(i) & SHA.Longitude(i-1)~=SHA.Longitude(i)

    % find points on great circle of mesh
    [dist,angl,phi,tracklat,tracklon,tracklen]=...
          ridgepack_greatcircle(SHA.Latitude(i-1),SHA.Longitude(i-1),...
                                SHA.Latitude(i),SHA.Longitude(i));

    % Find values within twice the horizon since we've already weeded
    % out cells that only fall within the horizon so as to plot entire 
    % cells, not just the vertices within the horizon.
    [dx,dy,dz,phi,theta]=ridgepack_satfwd([tracklat],[tracklon],...
                                          centlat,centlon,2*horizon,1);

    % draw grid
    plot3(dx,dy,dz,'Color',gridcolor,'LineWidth',0.2)
    hold on

   end

  end

 else

  % Find values within twice the horizon since we've already weeded
  % out cells that only fall within the horizon so as to plot entire 
  % cells, not just the vertices within the horizon.
  [dx,dy,dz,phi,theta]=ridgepack_satfwd([SHA.Latitude],...
                                        [SHA.Longitude],...
                                        centlat,centlon,2*horizon,1);

  % draw grid
  plot3(dx,dy,dz,'Color',gridcolor,'LineWidth',0.2)
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
if horizon>30
 light('Position',[0 0 10000],'Style','local')
 material dull
end

if nargout<1
 clear cells
end


