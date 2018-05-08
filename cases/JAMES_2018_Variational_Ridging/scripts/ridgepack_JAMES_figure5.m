% ridgepack_JAMES_figure5 - Generates Figure 5 in JAMES Variation Ridging paper
% 
% This script generates Figure 5 from:
%
% Roberts, A.F., E.C. Hunke, S.M. Kamal, W.H. Lipscomb, C. Horvat, W. Maslowski (2018),
% Variational Method for Sea Ice Ridging in Earth System Models, Part I: Theory, 
% submitted to J. Adv. Model Earth Sy.
%
% Andrew Roberts, Naval Postgraduate School, April 2018 (afrobert@nps.edu)

clear

% Create figure
h=figure(1);
clf

% set up colors
cols=lines(10);

boxleftx=0.05;
boxrightx=1-boxleftx;

for frame=1:2 % figures (a) and (b)
%for frame=1:1 % figures (a) and (b)

boxcenterxt=(boxleftx+boxrightx)/2;
boxcenterxb=(boxleftx+boxrightx)/2;

figLk=boxrightx-boxleftx;
levelextent=0.0;
sealeft=-0.1
searight=1.1
aspectratio=1;
textoffset=0.012;
arrowhead=0.006;
envelope=0.025;
bluecol=[0 0.447058826684952 0.74117648601532];
gcol=cols(7,:);
greycol=0.0*[1 1 1];
icecol=0.75*[1 1 1];
fricol=11;

% parameter settings of scheme
rhoi=917.0; % density of ice (kg/m^3)
rhos=330.0; % density of snow (kg/m^3)
rhow=1026.0; % density of seawater (kg/m^3)

% level ice and snow thickness left side
hfi1=2.0; % level ice thickness
hfs1=0.3; % snow thickness

if frame==1
 % level ice and snow thickness right side
 hfi2=hfi1; % level ice thickness
 hfs2=hfs1; % snow thickness
 alpha2=22; % angle of ridge/keel repose
 alpha1=22; % angle of ridge/keel repose
elseif frame==2
 % level ice and snow thickness right side
 hfi2=0.5; % level ice thickness
 hfs2=0.0; % snow thickness
 alpha1=22; % angle of ridge/keel repose
 alpha2=35; % angle of ridge/keel repose
end

% deformed ice and snow thickness
if frame==1
 hdi1=3.0; % ridged ice thickness
 hds1=0.3; % ridged snow thickness
elseif frame==2
 hdi1=3.0; % ridged ice thickness
 hds1=0.3; % ridged snow thickness
end

porosity=0.2; % porosity of ridge and keel complex
%porosity=0.3; % porosity of ridge and keel complex

epsilon1=(hfi1-hdi1*(1-porosity))/(hdi1*(1-porosity))
epsilon2=(hfi2-hdi1*(1-porosity))/(hdi1*(1-porosity))

% calculate freeboard and draft of level ice
hfd1=(rhoi*hfi1+rhos*hfs1)/rhow; % level draft
hff1=(hfi1+hfs1)-hfd1; % level freeboard

% calculate freeboard and draft of level ice
hfd2=(rhoi*hfi2+rhos*hfs2)/rhow; % level draft
hff2=(hfi2+hfs2)-hfd2; % level freeboard

% calculate freeboard and draft of deformed ice 
hdd1=(rhoi*hdi1+rhos*hds1)/rhow; % ridged draft
hdf1=(hdi1+hds1)-hdd1; % ridged freeboard

% check for bounds of level ice
if hfd1<0
 hfd1=0;
 hff1=0;
elseif hfd1>hdd1/(1-porosity) | hff1>hdf1/(1-porosity)
 hdd1=hfd1;
 hdf1=hff1;
end

% calculate depth of keel relative to sea level
Hk=(2*hdd1/(1-porosity))-hfd1;

% calculate horizontal extent of keel structure 
Lk1=2*(Hk-hfd1)/tan(alpha1*pi/180);

% calculate horizontal extent of keel structure 
Lk2=2*(Hk-hfd2)/tan(alpha2*pi/180);

% calculate height of ridge
Hr=hff1+sqrt((tan(alpha1*pi/180))*(((hdf1*Lk1)/(1-porosity))-hff1*Lk1));

% calculate horizontal extent of ridge structure 
Lr1=2*(Hr-hff1)/tan(alpha1*pi/180);

% calculate horizontal extent of ridge structure 
Lr2=2*(Hr-hff2)/tan(alpha2*pi/180);

% use this to determine the scale factor for the diagram
scalefactorx=figLk/Lk1;
scalefactory=min(aspectratio*scalefactorx,1/(Hk+Hr));

% calculate level ice box and surface snow, centered slightly to left
sealevely=0.6; % rhoi/rhow
ylsnow1=sealevely+hff1*scalefactory;
ylbottom1=sealevely-hfd1*scalefactory;
xlleft=boxleftx+(boxrightx-boxleftx-Lk1*scalefactorx)/2;

ylsnow2=sealevely+hff2*scalefactory;
ylbottom2=sealevely-hfd2*scalefactory;
xlright=boxrightx-(boxrightx-boxleftx-Lk2*scalefactorx)/2;

% calculate ridge top and keel bottom
yrtop=ylsnow1+(Hr-hff1)*scalefactory;
ykbottom=ylbottom1-(Hk-hfd1)*scalefactory;

% start axis
if frame==1
 axes('Position',[0.025 0.25 0.95 0.95],...
      'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1],...
      'Xlim',[0 1],'Ylim',[0 1],'ActivePositionProperty','position',...
      'Box','off','Visible','off','Clipping','off','Xtick',[],'XTickLabel',[],...
      'Ytick',[],'YTickLabel',[])
else
 axes('Position',[0.025 -0.225 0.95 0.95],...
      'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1],...
      'Xlim',[0 1],'Ylim',[0 1],'ActivePositionProperty','position',...
      'Box','off','Visible','on','Clipping','off','Xtick',[],'XTickLabel',[],...
      'Ytick',[],'YTickLabel',[])
end
 
ridgepack_clearax

ymin=get(gca,'Ylim');

set(gca,'FontSize',12)
labelsize=get(gca,'FontSize');

% calculate area above sea level in plotting coordinates for idealized shape
freeboardarea1=0.5*scalefactorx*Lr1*(yrtop-ylsnow1);
draftarea1=0.5*scalefactorx*Lk1*(ylbottom1-ykbottom);

freeboardarea2=0.5*scalefactorx*Lr2*(yrtop-ylsnow2);
draftarea2=0.5*scalefactorx*Lk2*(ylbottom2-ykbottom);

% plot roughened ridge structure with these parameters
maxx=25; % grid points along half-keel
if frame==1
 jx=[2:5:24]; % randomly grab index to adjust
 delt=5;  % size of delta to draw keels
else
 jx=[2:5:24]; % randomly grab index to adjust
 delt=3;  % size of delta to draw keels
end
if max(jx)>maxx
 error('maxx is too small')
end


% calculate left half keel

xincrement=(boxcenterxb-xlleft)./maxx;
yincrement=(ylbottom1-ykbottom)./maxx;
yincrementfirst=yincrement;

x=zeros([1 maxx+1]);
y=zeros([1 maxx+1]);

while any(y(2:end)>=ylbottom1 | y(1:end-1)<=ykbottom)
 area=zeros([maxx 1]);
 for i=1:maxx+1
  x(i)=xlleft+(i-1)*xincrement;
  coeff=(i/(maxx+1))^4;
  if i==1
   y(i)=ylbottom1;
  elseif i==maxx+1
   y(i)=ykbottom;
   area(i)=(x(i)-x(i-1))*(ylbottom1-0.5*(y(i)+y(i-1)));
  else
   y(i)=y(i-1)-delt*yincrement*(rand-0.49);
   y(i)=max(ykbottom+3*yincrement,min(y(i),ylbottom1-3*yincrement));
   y(i)=(1-coeff)*y(i)+coeff*ykbottom;
   area(i)=(x(i)-x(i-1))*(ylbottom1-0.5*(y(i)+y(i-1)));
  end
 end
 totalareabelow=sum(area(:));

 for idx=jx
  newarea=sum(area(idx:idx+1))+(draftarea1/2-totalareabelow)/length(jx);
  y1=(ylbottom1-y(idx-1));
  y2=(ylbottom1-y(idx));
  y3=(ylbottom1-y(idx+1));
  x1=x(idx-1);
  x2=x(idx);
  x3=x(idx+1);
  y(idx)=ylbottom1-(2*newarea-y1*(x2-x1)-y3*(x3-x2))./(x3-x1);
 end
end

plot(x,y,'b')

hold on

xxb=x;
yyb=y;

% calculate right half keel

xincrement=(scalefactorx*Lk2)./(2*maxx);
yincrement=(ylbottom2-ykbottom)./maxx;

x=zeros([1 maxx+1]);
y=zeros([1 maxx+1]);

while any(y(1:end-1)>=ylbottom2 | y(2:end)<=ykbottom)
 area=zeros([maxx 1]);

 for i=1:maxx+1
  x(i)=boxcenterxb+(i-1)*xincrement;
  coeff=(i/(maxx+1))^4;
  if i==1
   y(i)=ykbottom;
  elseif i==maxx+1
   y(i)=ylbottom2;
  else
   y(i)=y(i-1)+delt*yincrement*(rand-0.49);
   y(i)=max(ykbottom+3*yincrement,min(y(i),ylbottom2-3*yincrement));
   y(i)=(1-coeff)*y(i)+coeff*ylbottom2;
   area(i)=(x(i)-x(i-1))*(ylbottom2-0.5*(y(i)+y(i-1)));
  end
 end
 totalareabelow=sum(area(:));

 for idx=jx
  newarea=sum(area(idx:idx+1))+(draftarea2/2-totalareabelow)/length(jx);
  y1=(ylbottom2-y(idx-1));
  y2=(ylbottom2-y(idx));
  y3=(ylbottom2-y(idx+1));
  x1=x(idx-1);
  x2=x(idx);
  x3=x(idx+1);
  y(idx)=ylbottom2-(2*newarea-y1*(x2-x1)-y3*(x3-x2))./(x3-x1);
 end
end

plot(x,y,'b')

xxb=[xxb x(2:end)];
yyb=[yyb y(2:end)];


% calculate fill keel/ridge in with rubble dots
ya=ykbottom:yincrementfirst:yrtop;
xk=zeros([length(xxb) length(ya)]);
yk=zeros([length(xxb) length(ya)]);
for i=1:length(xxb)
for j=1:length(ya)
 if ya(j)<=yyb(i)-yincrement | ya(j)>yrtop
  xk(i,j)=NaN;
  yk(i,j)=NaN;
% elseif ya(j)>=ylbottom1 & xxb(i)<(boxcenterxt-scalefactorx*Lr1/2)
%  xk(i,j)=NaN;
%  yk(i,j)=NaN;
% elseif ya(j)>=ylbottom2 & xxb(i)>(boxcenterxt+scalefactorx*Lr2/2)
%  xk(i,j)=NaN;
%  yk(i,j)=NaN;
 elseif mod(i-mod(j,2),2)==0
  xk(i,j)=NaN;
  yk(i,j)=NaN;
 else
  xk(i,j)=xxb(i);
  yk(i,j)=ya(j);
 end
end
end


% calculate left half ridge

xincrement=(boxcenterxt-(boxcenterxt-scalefactorx*Lr1/2))./maxx;
yincrement=(yrtop-ylsnow1)./maxx;

x=zeros([1 maxx+1]);
y=zeros([1 maxx+1]);

while any(y(2:end)<=ylsnow1 | y(1:end-1)>=yrtop)
 area=zeros([maxx 1]);
 for i=1:maxx+1
  x(i)=(boxcenterxt-scalefactorx*Lr1/2)+(i-1)*xincrement;
  coeff=(i/(maxx+1))^4;
  if i==1
   y(i)=ylsnow1;
  elseif i==maxx+1
   y(i)=yrtop;
  else
   y(i)=(y(i-1)+delt*yincrement*(rand-0.49));
   y(i)=max(ylsnow1+3*yincrement,min(y(i),yrtop-3*yincrement));
   y(i)=(1-coeff)*y(i)+coeff*yrtop;
   area(i)=(x(i)-x(i-1))*(0.5*(y(i)+y(i-1))-ylsnow1);
  end
 end
 totalareaabove=sum(area(:));

 for idx=jx
  newarea=sum(area(idx:idx+1))+((freeboardarea1/2)-totalareaabove)/length(jx);
  y1=(y(idx-1)-ylsnow1);
  y2=(y(idx)-ylsnow1);
  y3=(y(idx+1)-ylsnow1);
  x1=x(idx-1);
  x2=x(idx);
  x3=x(idx+1);
  y(idx)=ylsnow1+(2*newarea-y1*(x2-x1)-y3*(x3-x2))./(x3-x1);
 end
end

plot(x,y,'b')

xxt=[xlleft x];
yyt=[ylsnow1 y];


% calculate right half ridge

xincrement=(boxcenterxt-(boxcenterxt-scalefactorx*Lr2/2))./maxx;
yincrement=(yrtop-ylsnow1)./maxx;

x=zeros([1 maxx+1]);
y=zeros([1 maxx+1]);

while any(y(1:end-1)<=ylsnow2 | y(2:end)>=yrtop)
 area=zeros([maxx 1]);
 for i=1:maxx+1
  x(i)=boxcenterxt+(i-1)*xincrement;
  coeff=(i/(maxx+1))^4;
  if i==1
   y(i)=yrtop;
  elseif i==1
   y(i)=yrtop;
  else
   y(i)=(y(i-1)-delt*yincrement*(rand-0.49));
   y(i)=max(ylsnow2+3*yincrement,min(y(i),yrtop-3*yincrement));
   y(i)=(1-coeff)*y(i)+coeff*ylsnow2;
   area(i)=(x(i)-x(i-1))*(0.5*(y(i)+y(i-1))-ylsnow2);
  end
 end
 totalareaabove=sum(area(:));

 for idx=jx
  newarea=sum(area(idx:idx+1))+((freeboardarea2/2)-totalareaabove)/length(jx);
  y1=(y(idx-1)-ylsnow2);
  y2=(y(idx)-ylsnow2);
  y3=(y(idx+1)-ylsnow2);
  x1=x(idx-1);
  x2=x(idx);
  x3=x(idx+1);
  y(idx)=ylsnow2+(2*newarea-y1*(x2-x1)-y3*(x3-x2))./(x3-x1);
 end
end

plot(x,y,'b')

xxt=[xxt x xlright];
yyt=[yyt y ylsnow2];

x=[xxt fliplr(xxb)];
y=[yyt fliplr(yyb)];

patch(x,y,icecol,'EdgeColor',icecol)

hold on

% locate centroid
yb=sealevely-scalefactory*hdd1/(1-porosity);
yt=sealevely+scalefactory*(hdf1-hds1)/(1-porosity);
ys=sealevely+scalefactory*hdf1/(1-porosity);

yc=(0.5*(yt+yb)*rhoi + 0.5*(ys+yt)*rhos)./(rhoi+rhos);
xc=mean([xlleft xlright]);

% put in rubble descriptor, removing nearby stipples
ht=text(xc,ylbottom1-2*textoffset,{'rubble'},...
	   'Interpreter','Latex',...
           'VerticalAlignment','middle',...
           'HorizontalAlignment','center',...
           'Fontsize',labelsize+1,...
           'Color',0.99999*[1 1 1],...
           'FontWeight','bold',...
           'EdgeColor','none');
ext=get(ht,'Extent');
mask=(xk>ext(1)-textoffset & xk<=ext(1)+ext(3)+textoffset & ...
      yk>=ext(2)-0.0*ext(4) & yk<=(ext(2)+1.25*ext(4)));
xk(mask)=NaN;
yk(mask)=NaN;


% label centroid
ht=text(xc,sealevely-1.1*textoffset,'centroid',...
	   'Interpreter','Latex',...
           'VerticalAlignment','top',...
           'HorizontalAlignment','center',...
           'Fontsize',labelsize-1,...
           'Color',cols(4,:));
ext=get(ht,'Extent');
mask=(xk>ext(1)-1.0*textoffset & xk<=ext(1)+ext(3)+1.0*textoffset & ...
      yk>=ext(2)-0.5*textoffset & yk<=(ext(2)+ext(4)-0.5*textoffset));
xk(mask)=NaN;
yk(mask)=NaN;

% plot sea level text across figure
st=text(xc,sealevely-0.5*textoffset,...
           {'sea level'},...
	   'Interpreter','Latex',...
           'VerticalAlignment','bottom',...
           'HorizontalAlignment','center',...
           'Fontsize',labelsize,...
           'Color',bluecol,...
           'EdgeColor','none');
ext=get(st,'Extent');
mask=(xk>ext(1)-1.0*textoffset & xk<=ext(1)+ext(3)+1.0*textoffset & ...
      yk>=ext(2) & yk<=(ext(2)+ext(4)+1.0*textoffset));
xk(mask)=NaN;
yk(mask)=NaN;


% draw rubble zone stipples
plot(xk,yk,'.','MarkerEdgeColor',0.999999*[1 1 1],'MarkerSize',5)

% draw centroid
plot(xc,yc,'o','Color',cols(4,:),'MarkerSize',10,'LineWidth',0.75)
hp=plot(xc,yc,'x','Color',cols(4,:),'MarkerSize',18,'LineWidth',0.75)

drawnow



% plot feeding level ice plates left
hold on
x=[xlleft sealeft sealeft xlleft];
y=[ylsnow1 ylsnow1 ylbottom1 ylbottom1];
patch(x,y,icecol,'EdgeColor',icecol)
plot([sealeft sealeft],[ylsnow1 ylbottom1],':','Color',0.99999*[1 1 1],'LineWidth',3)

% plot feeding level ice plates right
hold on
x=[searight xlright xlright searight];
y=[ylsnow2 ylsnow2 ylbottom2 ylbottom2];
patch(x,y,icecol,'EdgeColor',icecol)
plot([searight searight],[ylsnow2 ylbottom2],':','Color',0.99999*[1 1 1],'LineWidth',3)

% plot feeding plate bottom line right
hold on
x=[xlright boxcenterxt+scalefactorx*Lr2/2];
y=[ylbottom2 ylbottom2];
plot(x,y,':','Color',greycol,'LineWidth',0.5)

% plot feeding plate bottom line left
hold on
x=[boxcenterxt-scalefactorx*Lr1/2 xlleft];
y=[ylbottom1 ylbottom1];
plot(x,y,':','Color',greycol,'LineWidth',0.5)

% add feeding plate left arrows
x=[(boxcenterxt-scalefactorx*Lk1/2)-2*envelope ...
   (boxcenterxt-scalefactorx*Lk1/2)-2*arrowhead ...
   (boxcenterxt-scalefactorx*Lk1/2)-2*arrowhead ...
   (boxcenterxt-scalefactorx*Lk1/2) ...
   (boxcenterxt-scalefactorx*Lk1/2)-2*arrowhead ...
   (boxcenterxt-scalefactorx*Lk1/2)-2*arrowhead];
y=[(ylbottom1+ylsnow1)/2 ...
   (ylbottom1+ylsnow1)/2 ...
   (ylbottom1+ylsnow1)/2+arrowhead ...
   (ylbottom1+ylsnow1)/2 ...
   (ylbottom1+ylsnow1)/2-arrowhead ...
   (ylbottom1+ylsnow1)/2];

plot(x,y,'k')
if frame==1
 text((boxcenterxt-scalefactorx*Lk1/2)-2*envelope-textoffset/2,y(1),...
      '$\hat{v}_{F_a}$',...
      'VerticalAlignment','middle',...
      'HorizontalAlignment','right',...
      'Interpreter','Latex',...
      'Fontsize',labelsize-1,...
      'Color','k');

 text((sealeft+xlleft)/2,ylbottom1,...
      {'Floe $F_a$'},...
      'Interpreter','Latex',...
      'VerticalAlignment','bottom',...
      'HorizontalAlignment','center',...
      'Fontsize',labelsize,...
      'Color',0.99*[1 1 1]);
else
 text((boxcenterxt-scalefactorx*Lk1/2)-2*envelope-textoffset/2,y(1),...
      '$\hat{v}_{F_a}$',...
      'VerticalAlignment','middle',...
      'HorizontalAlignment','right',...
      'Interpreter','Latex',...
      'Fontsize',labelsize-1,...
      'Color','k');

% text((sealeft+xlleft)/2,ylbottom1,...
%      {'\makebox[1in][c]{Thick Floe}'},...
%      'Interpreter','Latex',...
%      'VerticalAlignment','bottom',...
%      'HorizontalAlignment','center',...
%      'Fontsize',labelsize,...
%      'Color',0.5*[1 1 1]);
end

% add feeding plate right arrows
x=[(boxcenterxt+scalefactorx*Lk2/2)+2*envelope ...
   (boxcenterxt+scalefactorx*Lk2/2)+2*arrowhead ...
   (boxcenterxt+scalefactorx*Lk2/2)+2*arrowhead ...
   (boxcenterxt+scalefactorx*Lk2/2) ...
   (boxcenterxt+scalefactorx*Lk2/2)+2*arrowhead ...
   (boxcenterxt+scalefactorx*Lk2/2)+2*arrowhead];
y=[(ylbottom2+ylsnow2)/2 ...
   (ylbottom2+ylsnow2)/2 ...
   (ylbottom2+ylsnow2)/2+arrowhead ...
   (ylbottom2+ylsnow2)/2 ...
   (ylbottom2+ylsnow2)/2-arrowhead ...
   (ylbottom2+ylsnow2)/2];

plot(x,y,'k')

if frame==1
 text((boxcenterxt+scalefactorx*Lk2/2)+2*envelope+textoffset/2,y(1),...
      '$\hat{v}_{F_b}$',...
      'VerticalAlignment','middle',...
      'HorizontalAlignment','left',...
      'Interpreter','Latex',...
      'Fontsize',labelsize-1,...
      'Color','k');
 text((searight+xlright)/2,ylbottom2,...
      {'Floe $F_b$'},...
      'Interpreter','Latex',...
      'VerticalAlignment','bottom',...
      'HorizontalAlignment','center',...
      'Fontsize',labelsize,...
      'Color',0.99*[1 1 1]);
else
 text((boxcenterxt+scalefactorx*Lk2/2)+2*envelope+textoffset/2,y(1),...
      '$\hat{v}_{F_b}$',...
      'VerticalAlignment','middle',...
      'HorizontalAlignment','left',...
      'Interpreter','Latex',...
      'Fontsize',labelsize-1,...
      'Color','k');
% text((searight+xlright)/2,ylbottom2-textoffset/2,...
%      {'\makebox[1in][c]{Thin Floe}'},...
%      'Interpreter','Latex',...
%      'VerticalAlignment','top',...
%      'HorizontalAlignment','center',...
%      'Fontsize',labelsize,...
%      'Color',0.5*[1 1 1]);
end

% annotate angles of repose
radius=(scalefactorx*Lr1/2)*0.25
theta=-[0:1:atan2d(ylbottom1-ykbottom,(boxrightx-boxleftx)/2)];
x=boxleftx+radius*cosd(theta);
y=ylbottom1+radius*sind(theta);
line(x,y,'Color','k','LineWidth',0.5)
line([boxleftx boxleftx+1.25*radius],[ylbottom1 ylbottom1],...
     'LineWidth',0.5,...
     'Color',greycol) 

if frame==1
 alphatext1='$\alpha_{R_a}$';
 alphatext2='$\alpha_{R_b}$';
elseif frame==2
 alphatext1='$\alpha_{R_a}$';
 alphatext2='$\alpha_{R_b}$';
end

text(boxleftx+radius+0.5*textoffset,mean(y)-0.01*textoffset,...
           alphatext1,...
           'Color',greycol,...
           'VerticalAlignment','middle',...
	   'HorizontalAlignment','left',...
           'Interpreter','Latex',...
           'Fontsize',labelsize,...
           'EdgeColor','none');

theta=180+[0:1:atan2d(ylbottom2-ykbottom,scalefactorx*Lk2/2)];
x=boxcenterxt+scalefactorx*Lk2/2+radius*cosd(theta);
y=ylbottom2+radius*sind(theta);
line(x,y,'Color','k','LineWidth',0.5)
line([boxcenterxt+scalefactorx*Lk2/2 boxcenterxt+scalefactorx*Lk2/2-1.25*radius],...
     [ylbottom2 ylbottom2],...
     'LineWidth',0.5,...
     'Color',greycol)

text(boxcenterxt+scalefactorx*Lk2/2-radius-0.5*textoffset,mean(y)-0.00*textoffset,...
           alphatext2,...
           'Color',greycol,...
           'VerticalAlignment','middle',...
           'HorizontalAlignment','right',...
           'Interpreter','Latex',...
           'Fontsize',labelsize,...
           'EdgeColor','none');

theta=[0:1:atan2d(ylbottom1-ykbottom,(boxrightx-boxleftx)/2)];
x=boxcenterxt-scalefactorx*Lr1/2+radius*cosd(theta);
y=ylsnow1+radius*sind(theta);
line(x,y,'Color','k','LineWidth',0.5)
line([boxcenterxt-scalefactorx*Lr1/2 boxcenterxt-scalefactorx*Lr1/2+1.25*radius],...
     [ylsnow1 ylsnow1],...
     'LineWidth',0.5,...
     'Color',greycol) 

text(boxcenterxt-scalefactorx*Lr1/2+radius+0.5*textoffset,mean(y)+0.35*textoffset,...
           alphatext1,...
           'Color',greycol,...
           'VerticalAlignment','middle',...
           'HorizontalAlignment','left',...
           'Interpreter','Latex',...
           'Fontsize',labelsize,...
           'EdgeColor','none');

theta=180-[0:1:atan2d(ylbottom2-ykbottom,scalefactorx*Lk2/2)];
x=boxcenterxt+scalefactorx*Lr2/2+radius*cosd(theta);
y=ylsnow2+radius*sind(theta);
line(x,y,'Color','k','LineWidth',0.5)
line([boxcenterxt+scalefactorx*Lr2/2 boxcenterxt+scalefactorx*Lr2/2-1.25*radius],...
     [ylsnow2 ylsnow2],...
     'LineWidth',0.5,...
     'Color',greycol) 

text(boxcenterxt+scalefactorx*Lr2/2-radius-0.5*textoffset,mean(y)+0.49*textoffset,...
           alphatext2,...
           'Color',greycol,...
           'VerticalAlignment','middle',...
           'HorizontalAlignment','right',...
           'Interpreter','Latex',...
           'Fontsize',labelsize,...
           'EdgeColor','none');

% plot sea level line
x=[sealeft searight]; 
y=[sealevely sealevely];
line(x,y,'Color',bluecol,'LineStyle','-.');

% plot isostatic unit
x=[xlleft, xlleft, boxcenterxt-scalefactorx*Lr1/2, boxcenterxt,...
    boxcenterxt+scalefactorx*Lr2/2, xlright, xlright, boxcenterxb, xlleft];
y=[ylbottom1, ylsnow1, ylsnow1, yrtop, ylsnow2, ylsnow2, ylbottom2, ykbottom, ylbottom1];
plot(x,y,'-','Color',cols(2,:))

% These lines put a circle at the vertices
%plot([xlleft xlleft boxcenterxb xlright xlright boxcenterxt+scalefactorx*Lr2/2 boxcenterxt boxcenterxt-scalefactorx*Lr1/2],...
%     [ylsnow1 ylbottom1 ykbottom ylbottom2 ylsnow2 ylsnow2 yrtop ylsnow1],...
%     'o','Color',cols(2,:))

axis off
axis equal

hold on

% label isostatic unit
%labelquanta=true;
labelquanta=false;
if labelquanta
 if frame==1 
  text(xlleft,ylsnow1,{'isostatic quanta'},...
	   'Interpreter','Latex',...
           'VerticalAlignment','bottom',...
           'HorizontalAlignment','left',...
           'Fontsize',labelsize,...
           'Color',cols(2,:));
  text(boxcenterxt+scalefactorx*Lk2/2,ylsnow1,{'symmetric'},...
	   'Interpreter','Latex',...
           'VerticalAlignment','bottom',...
           'HorizontalAlignment','right',...
           'Fontsize',labelsize,...
           'Color',cols(2,:));
 elseif frame==2
  text(xlleft,ylsnow1,{'isostatic quanta'},...
	   'Interpreter','Latex',...
           'VerticalAlignment','bottom',...
           'HorizontalAlignment','left',...
           'Fontsize',labelsize,...
           'Color',cols(2,:));
  text(boxcenterxt+scalefactorx*Lk2/2,ylsnow2,{'asymmetric'},...
	   'Interpreter','Latex',...
           'VerticalAlignment','bottom',...
           'HorizontalAlignment','right',...
           'Fontsize',labelsize,...
           'Color',cols(2,:));
 end
end

% locate centroid of keel
%xk=mean([xlleft boxcenterxb xlright])
%yk=mean([ylbottom ykbottom ylbottom])
%plot(xk,yk,'o')

% plot information text in lower left corner
info=true;
if info & frame==1
 textbox={...
  ['$\phi_R=',num2str(porosity,'%8.1f'),'$'],... %  ['$\theta_R=0$'],...
  ['$\alpha_{R_a}=\alpha_{R_b}=',num2str(alpha1),'^{\circ}$'],...
  ['$h_{F_a}=h_{F_b}=',num2str(hfi1,'%8.1f'),'$ m'],...
  ['$h_{{Fs}_a}=h_{{Fs}_b}=',num2str(hfs1,'%8.1f'),'$ m'],...
  ['$\epsilon_{R_{I_a}}=\epsilon_{R_{I_b}}=',num2str((hfi1-hdi1)/(hdi1),'%8.2f'),'$'],...
  ['$H_K=',num2str(Hk,'%8.1f'),'$ m, $H_S=',num2str(Hr,'%8.1f'),'$ m'],...
  ['$L_{K_a}=L_{K_b}=',num2str(Lk1,'%8.1f'),'$ m']
  };
elseif info & frame==2
 textbox={...
  ['$\phi_R=',num2str(porosity,'%8.1f'),'$'],...  %  ['$\theta_R=0$'],...
  ['$\alpha_{R_a}=',num2str(alpha1),'^{\circ}, \alpha_{R_b}=',num2str(alpha2),'^{\circ}$'],...
  ['$h_{F_a}=',num2str(hfi1,'%8.1f'),'$ m, $h_{F_b}=',num2str(hfi2,'%8.1f'),'$ m'],...
  ['$h_{{Fs}_a}=',num2str(hfs1,'%8.1f'),'$ m, $h_{{Fs}_b}=',num2str(hfs2,'%8.1f'),'$ m'],...
  ['$\epsilon_{R_{I_a}}=',num2str((hfi1-hdi1)/(hdi1),'%8.2f'),'$, ',...
   '$\epsilon_{R_{I_b}}=',num2str((hfi2-hdi1)/(hdi1),'%8.2f'),'$'],...
  ['$H_K=',num2str(Hk,'%8.1f'),'$ m, $H_S=',num2str(Hr,'%8.1f'),'$ m'],...
  ['$L_{K_a}=',num2str(Lk1,'%8.1f'),'$ m, $L_{K_b}=',num2str(Lk2,'%8.1f'),'$ m']
  };
end

if info
 %text(sealeft+2*textoffset,ykbottom,textbox,...
 text(sealeft+2*textoffset,ykbottom-2*textoffset,textbox,...
	    'Interpreter','Latex',...
            'VerticalAlignment','bottom',...
            'HorizontalAlignment','left',...
            'Fontsize',labelsize-2,...
            'Color','k',...
            'EdgeColor','none');
end

if frame==2

 % provide horizontal and vertical scale bar in lower right corner
 x=[searight-2*scalefactorx searight]-2*textoffset;
 y=[ykbottom ykbottom]-2*textoffset;
 line(x,y,'Color',0.0*[1 1 1],'LineStyle','-');
 line([x(1) x(1)],[y(1) y(1)+arrowhead],...
          'Color',0.0*[1 1 1],'LineStyle','-');
 text(sum(x)/2,sum(y)/2-0.5*textoffset,...
           ['2 m'],...
	   'Interpreter','Latex',...
           'VerticalAlignment','bottom',...
           'Fontsize',labelsize-1,...
           'HorizontalAlignment','center',...
           'EdgeColor','none');

 x=[searight searight]-2*textoffset;
 y=[ykbottom ykbottom+2*scalefactory]-2*textoffset;
 line(x,y,'Color',0.0*[1 1 1],'LineStyle','-');
 line([x(1)-arrowhead x(1)],[y(2) y(2)],...
          'Color',0.0*[1 1 1],'LineStyle','-');
 text(sum(x)/2+0.5*textoffset,sum(y)/2,...
           ['2 m'],...
	   'Interpreter','Latex',...
           'VerticalAlignment','bottom',...
           'HorizontalAlignment','center',...
           'Fontsize',labelsize-1,...
           'Rotation',90,...
           'EdgeColor','none');

end

% place alphabetical indicator for frame
if frame==1
 text(sealeft,ylsnow1+textoffset,'a)',...
      'Interpreter','Tex',...
      'VerticalAlignment','bottom',...
      'Fontsize',labelsize,...
      'HorizontalAlignment','left',...
      'EdgeColor','none',...
      'FontName','Helvetica');
elseif frame==2
 text(sealeft,ylsnow1+textoffset,'b)',...
      'Interpreter','Tex',...
      'VerticalAlignment','bottom',...
      'Fontsize',labelsize,...
      'HorizontalAlignment','left',...
      'EdgeColor','none',...
      'FontName','Helvetica');
end

end

% print out file in Ridgepack/cases/JAMES_2018_Variational_Ridging/output
dir=fileparts(which(mfilename));
cd([dir(1:strfind(dir,'scripts')-1),'output']);
ridgepack_fprint('png',['figure5'],1,2)
ridgepack_fprint('epsc',['figure5'],1,2)

