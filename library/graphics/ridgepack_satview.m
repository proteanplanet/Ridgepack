function [h]=ridgepack_satview(centlat,centlon,horizon,surface,gridon)

% set paramater space
surfh=0.999; % altiude of underlying surface
gridheight=1.01; % height of grid superimposed on plot
gridcolor=0.45*[1 1 1]; % color of grid lines and labels
labelfontsize=6; % font size of labels

% first pass input error checking 
if nargin<3
 error('missing input variables')
elseif nargin==3
 surface=1;
 gridon=2;
elseif nargin==4
 if surface<=1
  disp(['Adding an underlying surface at altitude ',num2str(surfh)])
 elseif surface==0
  disp('Not including underlying surface')
 elseif (surface>1 | surface<0)
  error('surface can take the values of 0 or 1 only')
 end
 gridon=2;
end

% check for grid options
if nargin==5
 if gridon~=0 & gridon~=1 & gridon~=2
  error('gridon must be equal to zero, one or two')
 end
end

if gridon==0
 disp('No grid')
elseif gridon==1
 disp('Grid on but grid labeling off')
elseif gridon==2
 disp('Grid on and grid labeling on')
end

% second pass error checking
if centlat>90 | centlat<-90 
 error('latitude must be between -90 and 90')
elseif centlon>180 | centlon<-180
 error('latitude must be between -180 and 180')
end

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
if gridon>0
 for lat=-80:10:80
 
  % plot grid lines
  lons=[0:0.01:360];
  lats=lat.*ones(size(lons));
  [x,y,z]=ridgepack_satfwd(lats,lons,centlat,centlon,horizon,1);
  plot3(x,y,gridheight*z,':','Color',gridcolor)
  hold on

  % label parallels
  if mod(lat-10,20)==0 & gridon==2

   lons=centlon-mod(centlon-30,30)-15;
   lats=lat;

   [x,y,z,ph,th]=ridgepack_satfwd(lats,lons,centlat,centlon,horizon,1);
   diffhoriz=(0.8*horizon)-rad2deg(th);
   idx=find(diffhoriz>0);
 
   if ~isempty(idx)

    idx=idx(1);

    % find local angle of labels
    [xl,yl,zl,phl,thl]=ridgepack_satfwd([lats lats],...
                                       [lons(idx)-2 lons(idx)+2],...
                                       centlat,centlon,horizon,1);
    if max(rad2deg(thl))<horizon

     % find angle of grid label and rotate to readable angle
     rotation=atan2d(diff(yl),diff(xl));
     if rotation>91 | rotation<-91
      rotation=mod(rotation+180,0);
     end

     % add grid label
     if ~isnan(rotation)

      if lats>0
       ending='$^{\circ}$N';
      elseif lats<0
       ending='$^{\circ}$S';
      else
       ending='$^{\circ}$';
      end

      % get coordinates of meredian label near frame edge
      text(x(idx),y(idx),2*gridheight*z(idx),...
               [num2str(abs(lats)),ending],...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom',...
               'FontSize',labelfontsize,...
               'Color',gridcolor,...
               'Rotation',rotation);
     end

    end

   end

  end

 end

 % add meridians
 for lon=-180:30:180

  % plot grid lines
  if lon==-90 | lon==0 | lon==90 | lon==180
   lats=[-90:0.01:90];
  else
   lats=[-80:0.01:80];
  end
  lons=lon.*ones(size(lats));
  [x,y,z,phi,theta]=ridgepack_satfwd(lats,lons,centlat,centlon,horizon,1);
  plot3(x,y,gridheight*z,':','Color',gridcolor)
  hold on

  % label grid lines
  if (lon==-90 | lon==90 | lon==0  | lon==180) & gridon==2

   % get coordinates of meredian label near frame edge
   if centlat>0
    lats=[-75:10:75];
   else
    lats=[75:-10:-75];
   end
   lons=lon*ones(size(lats));

   [x,y,z,ph,th]=ridgepack_satfwd(lats,lons,centlat,centlon,horizon,1);
   diffhoriz=(0.8*horizon)-rad2deg(th);
   idx=find(diffhoriz>0);

   if ~isempty(idx) 
    idx=idx(1);
    % find local angle of labels
    [xl,yl,zl,phl,thl]=ridgepack_satfwd([lats(idx)-2 lats(idx)+2],...
                                       [lon lon],...
                                       centlat,centlon,horizon,1);
    if max(rad2deg(thl))<horizon

     % find angle of grid label and rotate to readable angle
     rotation=atan2d(diff(yl),diff(xl));
     if rotation>91 | rotation<-91
      rotation=mod(rotation+180,0);
     end
 
     % add grid label
     if ~isnan(rotation)

      if lon==180 | lon==0
       ending='$^{\circ}$E';
      elseif lon<0
       ending='$^{\circ}$W';
      elseif lon>0
       ending='$^{\circ}$E';
      end

      text(x(idx),y(idx),2*gridheight*z(idx),...
               [num2str(abs(lon)),ending],...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom',...
               'FontSize',labelfontsize,...
	       'Color',gridcolor,...
               'Rotation',rotation);
     end
    end
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

