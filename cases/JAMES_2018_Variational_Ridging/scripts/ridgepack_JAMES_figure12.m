function ridgepack_dilationplot(titleflag)

% function ridgepack_dilationplot
%
% INPUT: 
%
% titleflag - logical turning on title
%
% This function is part of Ridgepack Version 1.0.
% It generates a graphical example of the state space trajectory and 
% dilation field for a ridge with parent sheet thickness of 2m and no snow
% cover. This function is useful for testing the code in ridgepack_trajectory
% and that is why it was created.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

% plot on figure 2
clear
figure(2)
clf

if nargin<1
 titleflat=false;
end

% initial ice thickness (m)
hf=2.0;

% initial snow thickness (m)
hfs=0.0;

% retrieve constants
[rhoi,rhos,rhow,delrho,g]=ridgepack_constants;

% calculate trajectory for given thickness
[EPSILON,PHI,ALPHAHAT,VR,HK,HS,LK,LS,epmesh,phmesh,vr,epsplitmesh,phsplitmesh,d1,d2]=...
            ridgepack_trajectory(hf,hfs);

% set log color scale for potential energy density
contourss=10.^[floor(log10(min(vr(:)))):0.25:ceil(log10(max(vr(:))))];
cmap=colormap(parula(length(contourss)));
[zindex,truecol]=ridgepack_colorindex(vr,contourss,0);

% plot potential energy
surface(epmesh,phmesh,-0.1*ones(size(vr)),truecol,'EdgeColor','k')
shading flat
hold on

% plot the dilation field over the top as streamlines
hstr=streamslice(epsplitmesh,phsplitmesh,d1,d2);
set(hstr,'Color',[1 1 1])

% only plot zeta-hat trajector up to a min strain of -0.96
idx=find(EPSILON>=-0.96);
EPSILON=EPSILON(idx);
PHI=PHI(idx);
VR=VR(idx);

% plot ridging path
plot3(EPSILON,PHI,0.001*ones(size(VR)),'r','LineStyle','-','LineWidth',1.0)

% annote trajectory
text(EPSILON(end)+0.025,PHI(end),0.9,...
           ['$\hat{\zeta}$'],...
            'Color','r',...
            'Fontsize',10,...
            'Rotation',10,...
            'Margin',0.5,...
            'VerticalAlignment','bottom',...
            'HorizontalAlignment','center',...
            'Interpreter','Latex')

% arrow head of mean path
text(EPSILON(end)+0.008,PHI(end),0.9,...
                 ['$\bigtriangleup$'],...
                 'Color','r',...
                 'Fontsize',9,...
                 'Margin',0.0001,...
                 'Rotation',90,...
                 'VerticalAlignment','bottom',...
                 'HorizontalAlignment','center',...
                 'Interpreter','Latex')

% set the dynamic limits
xlim([-0.99 0])
ylim([0 0.99])
zlim([-1 1])

% plot the colorbar
ridgepack_colorbar(contourss,['J m^{-2}'],0.1,'vertical');

% fiddle to get the grid to plot
grid
grid off
grid on

% set log scale for z-axis 
set(gca,'Box','on','TickLabelInterpreter','Latex','Layer','top')

% set tick marks
set(gca,'Xtick',[-0.9:0.1:0])
set(gca,'Ytick',[0:0.1:0.9])

% axis labels
xlabel('Strain, $\epsilon_{R_I}$','Interpreter','Latex','fontsize',12)
ylabel('Porosity, $\phi_R$','Interpreter','Latex','fontsize',12)
zlabel('Potential Energy Density, $\mathcal{V}_R$ (J m$^{-2}$)','Interpreter','Latex','fontsize',12)

if titleflat

 % title
 title(['Dilation field and state space trajectory for $h_f=$',num2str(hf),'m, $h_{f_s}=$',num2str(hfs),'m'])
 
 % fix axes
 axis square

else

 % changes the position
 set(gca,'Position',[0.1 0.08 0.7 0.9])

end

% determine directory for read/write of zeta-hat plane data
writedir=[fileparts(which('ridgepack')),'/figures'];
cd(writedir)

% print figure
ridgepack_fprint('png','ridgepack_dilationplot',2,2)
ridgepack_fprint('epsc','ridgepack_dilationplot',2,2)

