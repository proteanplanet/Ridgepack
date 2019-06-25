function [h]=ridgepack_satview(centlat,centlon,horizon,surface)

% first pass input error checking 
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

% second pass error checking
if centlat>90 | centlat<-90 
 error('latitude must be between -90 and 90')
elseif centlon>180 | centlon<-180
 error('latitude must be between -180 and 180')
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
 R = surfh*ones(size(th)); % should be R(theta,phi) surface in general
 cx = R.*sin(th).*cos(ph);
 cy = R.*sin(th).*sin(ph);
 cz = R.*cos(th);
 theta=abs(th);
 maxtheta=max(theta(:));
 %c1=cos(theta/max(1,horizon/45)); % 3D shading (halved fade-off a 90 deg)
 c1=cos(theta/2); % 3D shading (halved fade-off a 90 deg)
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

  if labellat>88 | labellat<-88
   error('labellat should be less than abs(89)')
  end

  % get coordinates of meredian label near frame edge
  lats=[-75:10:75];
  lons=lon*ones(size(lats));
  [x,y,z,ph,th]=ridgepack_satmap(lats,lons,centlat,centlon,horizon,1);
  idx=find(min(abs(rad2deg(th)-50))==abs(rad2deg(th)-50));

  if ~isempty(idx) & rad2deg(th(idx))<60

   rad2deg(th(idx))

   % find local angle of labels
   [xl,yl,zl]=ridgepack_satmap([lats(idx)-2 lats(idx)+2],[lon lon],...
                               centlat,centlon,horizon,1);

   % find angle of grid label and rotate to readable angle
   rotation=atan2d(diff(yl),diff(xl));
   if rotation>91 | rotation<-91
    rotation=mod(rotation+180,0)
   end
 
   % add grid label
   if ~isnan(rotation)
    text(x(idx),y(idx),1.02*z(idx),[num2str(abs(lon)),ending],...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom',...
               'FontSize',7,...
	       'Color',0.5*[1 1 1],...
               'Rotation',rotation);
   end

  end

 end
end

% plot frame
ph = deg2rad([0:0.001:361]); 
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

