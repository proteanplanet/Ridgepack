function [h]=ridgepack_satview(centlat,centlon,horizon,surface,gridon)


global debug;
%debug=true;
if debug; disp(['Entering ',mfilename,'...']); end


% set paramater space
surfh=0.999; % altiude of underlying surface
gridheight=1.01; % height of grid superimposed on plot
gridcolor=0.80*[1 1 1]; % color of grid lines and labels
labelcolor=0.80*[1 1 1]; % color of grid lines and labels
scalecolor=0.70*[1 1 1]; % color of 
labelfontsize=5; % font size of labels
scalefontsize=8; % font size of labels

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

 % insep indicates lat-lon labels within the figure
 insep=false;
 if horizon<2
  latspace=0.5;
  lonspace=1;
  latlabsp=1;
  lonlabsp=1;
 elseif horizon<5
  latspace=1;
  lonspace=2;
  latlabsp=1;
  lonlabsp=1;
 elseif horizon<10
  latspace=2;
  lonspace=5;
  latlabsp=1;
  lonlabsp=1;
 elseif horizon<30
  latspace=5;
  lonspace=10;
  latlabsp=1;
  lonlabsp=1;
 else
  insep=true;
  latspace=10;
  lonspace=30;
  latlabsp=2;
  lonlabsp=3;
 end
 
 % initialization of psl, which records location of lat labels
 latk=1;
 psl(latk)=270; 

 % now add latitude lines
 for lat=-90+latspace:latspace:90-latspace
 
  % plot grid lines
  lons=[0:0.01:360];
  lats=lat.*ones(size(lons));
  [xs,ys,zs,ps,ts]=ridgepack_satfwd(lats,lons,centlat,centlon,horizon,1);
  plot3(xs,ys,gridheight*zs,':','Color',gridcolor)
  hold on

  % label parallels
  if mod(lat-latspace,latlabsp*latspace)==0 & gridon==2

   lons=centlon-mod(centlon-lonspace,lonspace)-lonspace/2;
   lats=lat;

   [x,y,z,ph,th]=ridgepack_satfwd(lats,lons,centlat,centlon,horizon,1);
   if insep
    diffhoriz=(0.8*horizon)-rad2deg(th);
   else
    diffhoriz=horizon-rad2deg(th);
   end
   idx=find(diffhoriz>0);
 
   if ~isempty(idx)

    idx=idx(1);

    % find local angle of labels
    [xl,yl,zl,phl,thl]=ridgepack_satfwd([lats lats],...
                                       [lons(idx)-0.2 lons(idx)+0.2],...
                                       centlat,centlon,horizon,1);
    % add grid label
    if max(rad2deg(thl))<horizon

     if lats>0
       ending='$^{\circ}$N';
     elseif lats<0
       ending='$^{\circ}$S';
     else
       ending='$^{\circ}$';
     end

     % get coordinates of meredian label near frame edge
     if insep

      % find grid label angle & rotate to readable angle
      rotation=atan2d(diff(yl),diff(xl));
      if isnan(rotation)
        rotation=90;
      end

      if ~isnan(rotation)
       % plot labels inside frame
       rotation=wrapTo180(rotation);
       if rotation>91 | rotation<-91
        rotation=mod(rotation+180,0);
       end
       text(x(idx),y(idx),2*gridheight*z(idx),...
               [num2str(abs(lats)),ending],...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom',...
               'FontSize',labelfontsize,...
               'Color',labelcolor,...
               'Rotation',rotation);
      end

     else

      % find edge of frame and plot on edge
      for i=2:length(xs)
       if isnan(xs(i-1)) & ~isnan(xs(i))
         latk=latk+1;
         psl(latk)=ps(i); 
         rotation=90+ps(i)*180/pi;
         rotation=wrapTo180(rotation);
         if rotation>90 | rotation<-90
          rotation=mod(rotation+180,0);
         end
         if ys(i)<0 
          vposition='top';
         elseif ys(i)>=0
          vposition='bottom';
         end
         
         text(xs(i),ys(i),gridheight,...
             [num2str(abs(lats)),ending],...
             'HorizontalAlignment','center',...
             'VerticalAlignment',vposition,...
             'FontSize',labelfontsize,...
             'Color',labelcolor,...
             'Rotation',rotation);
 
       end
      end

     end

    end
   end
  end
 end

 % add meridians
 lonlabel=[-180+lonlabsp*lonspace:lonlabsp*lonspace:180];
 for lon=-180:lonspace:180

  % plot grid lines
  if lon==-90 | lon==0 | lon==90 | lon==180
   lats=[-90:0.01:90];
  else
   lats=[-80:0.01:80];
  end
  lons=lon.*ones(size(lats));
  [xs,ys,zs,ps,ts]=ridgepack_satfwd(lats,lons,...
                              centlat,centlon,horizon,1);
  plot3(xs,ys,gridheight*zs,':','Color',gridcolor)
  hold on

  if any(lonlabel==lon) & gridon==2

   % get coordinates of meredian label near frame edge
   if centlat>0
    lats=[-90+lonspace+lonspace/2:lonspace:90-lonspace-lonspace/2];
   else
    lats=[90-lonspace-lonspace/2:-lonspace:-90+lonspace+lonspace/2];
   end
   lons=lon*ones(size(lats));

   [x,y,z,ph,th]=ridgepack_satfwd(lats,lons,centlat,centlon,horizon,1);
   if insep
    diffhoriz=(0.8*horizon)-rad2deg(th);
   else
    diffhoriz=horizon-rad2deg(th);
   end
   idx=find(diffhoriz>0);

   if ~isempty(idx) 

    idx=idx(1);

    % find local angle of labels
    [xl,yl,zl,phl,thl]=ridgepack_satfwd([lats(idx)-0.2 lats(idx)+0.2],...
                                       [lon lon],...
                                       centlat,centlon,horizon,1);

    if max(rad2deg(thl))<horizon

     % add grid label
     if lon==180 | lon==0
       ending='$^{\circ}$E';
     elseif lon<0
       ending='$^{\circ}$W';
     elseif lon>0
       ending='$^{\circ}$E';
     end

     if insep

      % find angle of grid label and rotate to readable angle
      rotation=atan2d(diff(yl),diff(xl));
      if isnan(rotation)
       rotation=90;
      end

      if ~isnan(rotation)

       % plot labels inside frame
       rotation=wrapTo180(rotation);
       if rotation>90 | rotation<-90
        rotation=mod(rotation+180,0);
       end
       text(x(idx),y(idx),2*gridheight*z(idx),...
               [num2str(abs(lon)),ending],...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom',...
               'FontSize',labelfontsize,...
	       'Color',gridcolor,...
               'Rotation',rotation);

      end

     else

      % find edge of frame and plot on edge
      for i=2:length(xs)
       if isnan(xs(i-1)) & ~isnan(xs(i))
        rotation=90+ps(i)*180/pi;
        rotation=wrapTo180(rotation);
        if rotation>91 | rotation<-91
          rotation=mod(rotation+180,0);
        end

        % avoid overwriting latitude labels
        if (all(abs(ps(i)-psl)>5*pi/180))
          if ys(i)<=0 
           vposition='top';
          elseif ys(i)>0
           vposition='bottom';
          end

          text(xs(i),ys(i),gridheight,...
             [num2str(abs(lon)),ending],...
             'HorizontalAlignment','center',...
             'VerticalAlignment',vposition,...
             'FontSize',labelfontsize,...
             'Color',labelcolor,...
             'Rotation',rotation);
        end
       end
      end
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

% plot scale in km in lower right
if ~insep
 xrange=range(x);
 diamkm=horizon*60*1.852;
 mnxs=max(x)/2;
 mxxs=29*max(x)/30;
 mnys=29*min(y)/30;
 mxys=mnys;
 maxscale=mxxs-mnxs;
 scalelen=(diamkm*maxscale)/(xrange);
 newscale=10*floor(scalelen/10);
 xdist=newscale*(mxxs-mnxs)/scalelen;
 nmnxs=mnxs+(mxxs-mnxs)*xdist/2;
 nmxxs=mxxs-(mxxs-mnxs)*xdist/2;
 xscale=[nmnxs nmnxs NaN nmnxs nmxxs NaN nmxxs nmxxs];
 yscale=[0.99*mnys 1.01*mnys NaN mnys mxys NaN 0.99*mnys 1.01*mnys];
 plot3(xscale,yscale,gridheight*ones(size(xscale)),'Color',scalecolor)
 text((nmnxs+nmxxs)/2,mnys,gridheight,...
      [num2str(newscale),'~km'],...
      'HorizontalAlignment','center',...
      'VerticalAlignment','bottom',...
      'FontSize',scalefontsize,...
      'Color',scalecolor,...
      'Rotation',0);
end

% fix axes
axis off
axis equal
view([0 0 0.4])
axis tight

% clear drawing
drawnow

if debug; disp(['Leaving ',mfilename,'...']); end

