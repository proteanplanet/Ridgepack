function [h,nc]=ridgepack_e3smsatthreshold(ncvert,centlat,centlon,horizon,color,linestyle,label)

if nargin<2
 centlat=90;
elseif ~isnumeric(centlat)
 error('centlat should be a real number between -90 and 90')
elseif centlat>90 | centlat<-90
 error('centlat should be between -90 and 90')
end

if nargin<3
 centlon=0;
elseif ~isnumeric(centlon)
 error('centlon should be a real number between -90 and 90')
elseif centlat>180 | centlat<-180
 error('centlon should be between -180 and 180')
end

if nargin<4
 horizon=90;
elseif ~isnumeric(horizon)
 error('horizon must be a real number between 0 and 90')
elseif horizon<0 | horizon>90
 error('horizon must be between 0 and 90')
end

if nargin<5
 color=0.25*[0 0 1];
end

if nargin<6
 linestyle='-';
end

if nargin<7
 label=[];
end

% altitude
altitude=1.0001;

% generate coastline if required
if strcmp(ncvert.attributes.title,'E3SM-MPAS Edge Definition')
 nc=ncvert;
else
 nc=ridgepack_e3smcoastm(ncvert);
end

% plot coastline within the horizon
[x,y,z,phi,theta]=ridgepack_satfwd(nc.xlatitude.data,...
                                   nc.xlongitude.data,...
                                   centlat,centlon,...
                                   horizon,altitude);

h=plot3(x,y,z,'Color',color,'LineWidth',0.6-sin(deg2rad(horizon))*0.4,'LineStyle',linestyle);

% plot label if requested
if ~isempty(label)
 exit('no current way to plot labels')
end

% find lines that cross the horizon
k=0;
sidx=[];
for i=2:length(x)
 if (isnan(x(i)) & ~isnan(x(i-1)))  | (~isnan(x(i)) & isnan(x(i-1)))
  k=k+1;
  sidx(k)=[i-1];
 end
end

% truncate coastal edges that cross the boundary of the plot
for k=1:length(sidx)

 idx=[sidx(k):sidx(k)+1];

 [x,y,z,phi,theta]=ridgepack_satfwd(nc.xlatitude.data(idx(2)),...
                                    nc.xlongitude.data(idx(2)),...
                                    nc.xlatitude.data(idx(1)),...
                                    nc.xlongitude.data(idx(1)),...
                                    180,altitude);

 % divide offending line into 100ths, and draw to closest as 
 % possible to edge of plot. Phi is constant along line
 incrt=theta/100;
 thetax=[0:incrt:theta];
 phix=phi*ones(size(thetax));
 
 % now generate lats and lons along the line
 [lats,lons]=ridgepack_satinv(phix,thetax,...
                              nc.xlatitude.data(idx(1)),...
                              nc.xlongitude.data(idx(1)));

 % NaNs will be set for points outside of plot
 [x,y,z,phi,theta]=ridgepack_satfwd(lats,lons,...
                                    centlat,centlon,...
                                    horizon,altitude);

 plot3(x,y,z,'Color',color,...
             'LineWidth',0.6-sin(deg2rad(horizon))*0.4,...
             'LineStyle',linestyle)

end

% clear if there is not nargout
if nargout==0
 clear nc
 clear h
end


