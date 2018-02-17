function [STRAIN,PHI,PE,strainsplit,phisplit,dVdstrain,dVdphi]=paper_ridging_path(strain,phi,pe)

% set debug flag
debug=false;

% create strain and split strain coordinates
stspacing=abs(strain(2,1)-strain(1,1));
stii=strain(:,1);
dstii=(stspacing./2)+stii(1:end-1);

% create phi and split phi coordinates
phspacing=abs(phi(1,2)-phi(1,1));
phii=phi(1,:);
dphii=(phspacing./2)+phii(1:end-1);

% calculate vector gradient
dVdst=(pe(2:end,:)-pe(1:end-1,:))./(strain(2:end,:)-strain(1:end-1,:));
dVdph=(pe(:,2:end)-pe(:,1:end-1))./(phi(:,2:end)-phi(:,1:end-1));
dVdstrain=(dVdst(:,1:end-1)+dVdst(:,2:end))./2;
dVdphi=(dVdph(1:end-1,:)+dVdph(2:end,:))./2;

% calculate streamline from initial conditions
STRAIN=stii(1:1:end-1);
jdx=1;
PHI(1)=phii(jdx); % initial condition
PE(1)=pe(1,jdx);
for i=2:length(STRAIN)
 PHI(i)=PHI(i-1)-(dVdphi(i,jdx)./dVdstrain(i,jdx)).*stspacing;
 jdx=find(min(abs(PHI(i)-phii(:)))==abs(PHI(i)-phii(:)));
 if jdx>size(dVdphi,2) | jdx>size(dVdstrain,2)
  error('Max phi exceeded')
 end
 PE(i)=pe(i,jdx);
end

% plot PE as color underlay, streamline field, and then path
if debug
 figure(2)
 clf
 AKS=10.^[2:1:11];
 cmap=nccolormap(AKS,0,'bluered',true);
 cmap=flipud(cmap);
 colormap(cmap);
 [zindex,truecol]=nccolorindex(pe,AKS,0);
 surface(strain,phi,pe,truecol,'FaceAlpha',0.95,'EdgeColor','k')
 shading flat
 hold on
 mag=sqrt(dVdstrain.^2 + dVdphi.^2);
 ncquiverref(strainsplit,phisplit,dVdstrain./mag,dVdphi./mag,'','median','w')
 plot3(STRAIN,PHI,PE,'r')
 ncfprint('png','path_test',2,1)
end

% prepare vector field for output
[strainsplit,phisplit]=meshgrid(dstii,dphii);
dVdstrain=dVdstrain';
dVdphi=dVdphi';

