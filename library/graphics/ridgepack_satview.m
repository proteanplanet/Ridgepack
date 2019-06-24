function [h]=ridgepack_satview(centlat,centlon,horizon,surface)

% input error checking
if nargin<3
 error('missing input variables')
elseif nargin==3
 surface=1;
elseif nargin==4
 if surface<=1
  disp('Adding an underlying surface at a given altitude')
 elseif surface==0
  disp('Not including underlying surface')
 elseif (surface>1 | surface<0)
  error('surface can take the values of 0 or 1 only')
 end
end

% set paramater space
surfh=0.999 % altiude of underlying surface
gridheight=1.01 % height of grid superimposed on plot

clf

% plot underlying surface
if surface>0 & surface<=1
 N = 500;
 thetavec = linspace(0,deg2rad(horizon),N);
 phivec = linspace(0,2*pi,N);
 [th, ph] = meshgrid(thetavec,phivec);
 R = surfh*ones(size(th)); % should be your R(theta,phi) surface in general
 cx = R.*sin(th).*cos(ph);
 cy = R.*sin(th).*sin(ph);
 cz = R.*cos(th);
 theta=abs(th);
 maxtheta=max(theta(:));
 c1=cos(theta/max(1,horizon/45)); % 3D shading (halved fade-off a 90 deg)
 cc(:,:,1)=c1;
 cc(:,:,2)=c1;
 cc(:,:,3)=c1;
 surf(cx,cy,cz,cc,'EdgeColor','none');
 hold on
end

% add parallels 
for lat=-0:10:80
 lons=[0:0.01:360];
 lats=lat.*ones(size(lons));
 [x,y,z]=ridgepack_satmap(lats,lons,centlat,centlon,horizon,1);
 plot3(x,y,gridheight*z,':','Color',0.5*[1 1 1])
 hold on
end

% add meridians
for lon=-180:30:180

 if lon==-90 | lon==0 | lon==90 | lon==180
  lats=[-90:0.01:90];
 else
  lats=[-80:0.01:80];
 end

 lons=lon.*ones(size(lats));
 [x,y,z,phi,theta]=ridgepack_satmap(lats,lons,centlat,centlon,horizon,1);
 plot3(x,y,gridheight*z,':','Color',0.5*[1 1 1])
 hold on

 if lon==-90 | lon==0 | lon==90 | lon==180

  if lon==180 | lon==0
   ending='$^{\circ}$E';
  elseif lon<0
   ending='$^{\circ}$W';
  elseif lon>0
   ending='$^{\circ}$E';
  end

  labellat=35;

  [x,y,z]=ridgepack_satmap(labellat,lon,centlat,centlon,horizon,1);

  rotation=lon-centlon;
  if rotation<91 & rotation>-91
   text(x,y,1.05*z,[num2str(abs(lon)),ending],...
         'HorizontalAlignment','center',...
         'VerticalAlignment','base',...
         'FontSize',8,...
	 'Color',0.5*[1 1 1],...
         'Rotation',rotation-90);
  end

 end

end

% plot frame
ph = deg2rad([0:0.01:360]); 
th = deg2rad(horizon*ones(size(ph)));
x = sin(th).*cos(ph);
y = sin(th).*sin(ph);
z = gridheight*cos(th);
plot3(x,y,z,'k-')
hold on

% fix axes
axis off
axis equal
view([0 0 0.4])
axis tight

