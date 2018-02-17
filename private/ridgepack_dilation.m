%function [STRAIN,PHI,PE]=ridge_expected_path(strain,phi,pe)

clear
clf

% constants
[rhoi,rhos,rhow,delrho,g]=ridge_constants;

hfii=[2.0]; % thickness of initial ice
hfsi=[0.0]; % thickness of snow on the ice
spacing=0.01;

% create strain and split strain coordinates
stii=[-1+spacing:spacing:-spacing];
dstii=[-1+1.5*spacing:spacing:-1.5*spacing];

% create phi and split phi coordinates
phii=[0+spacing:spacing:1-spacing];
dphii=[0+1.5*spacing:spacing:1-1.5*spacing];

for k=1:length(stii);
for j=1:length(phii);

 strain(k,j)=stii(k);
 phi(k,j)=phii(j);

 [pe(k,j),LK(k,j),LS(k,j),HK(k,j),HS(k,j),ALPHA(k,j)]=...
           paper_ridge_energetics(hfii,hfsi,stii(k),phii(j));

end
end

% create strain and phi split coordinates
[dstrain,dphi]=meshgrid(dstii,dphii);

% calculate vector gradient
dVdst=(pe(2:end,:)-pe(1:end-1,:))./(strain(2:end,:)-strain(1:end-1,:));
dVdph=(pe(:,2:end)-pe(:,1:end-1))./(phi(:,2:end)-phi(:,1:end-1));
dVdstrain=(dVdst(:,1:end-1)+dVdst(:,2:end))./2;
dVdphi=(dVdph(1:end-1,:)+dVdph(2:end,:))./2;

% calculate d2/d1
d2ond1=flipud(dVdphi./dVdstrain);

% calculate streamline
STRAIN=stii(end:-1:2);
jdx=1;
PHI(1)=phii(jdx); % initial condition
PE(1)=pe(end,jdx);
for i=2:length(STRAIN)
 PHI(i)=PHI(i-1)-d2ond1(i,jdx)*spacing;
 jdx=find(min(abs(PHI(i)-phii(:)))==abs(PHI(i)-phii(:)));
 PE(i)=pe(i,jdx);
end

% magnitude of vectors
mag=sqrt(dVdstrain.^2 + dVdphi.^2);

% plot PE as color underlay, streamline field, and then path
AKS=10.^[2:1:11];
cmap=nccolormap(AKS,0,'bluered',true);
cmap=flipud(cmap);
colormap(cmap);
[zindex,truecol]=nccolorindex(pe,AKS,0);
surface(strain,phi,pe,truecol,'FaceAlpha',0.95,'EdgeColor','k')
shading flat
hold on
ncquiverref(dstrain,dphi,dVdstrain./mag,dVdphi./mag,'','median','w')
plot3(STRAIN,PHI,PE,'r')

ncfprint('png','test',1,1)


%dphidstrain=

%[STRAIN,PHI]=ode45(






