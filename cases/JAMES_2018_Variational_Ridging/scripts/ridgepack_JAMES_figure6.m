% ridgepack_JAMES_figure6 - Generates Figure 6 in JAMES Variation Ridging paper
% 
% This script generates Figure 6 from:
%
% Roberts, A.F., E.C. Hunke, S.M. Kamal, W.H. Lipscomb, C. Horvat, W. Maslowski (2018),
% Variational Method for Sea Ice Ridging in Earth System Models, Part I: Theory, 
% submitted to J. Adv. Model Earth Sy.
%
% This function requires a png version of figure 16c from:
% 
% Davis, N. R., and P. Wadhams, 1995: A statistical analysis of Arctic pressure 
% ridge morphology. J. Geophys. Res., 100, 10915-10925, doi:10.1029/95JC00007. 
%
% The png file must be cropped to the exact axis limits of that figure. The 
% name of the file is set with 'fname', and its location is set with 'floc'.
% The file cannot be supplied with Ridgepack due to copyright restrictions.
%
% Andrew Roberts, Naval Postgraduate School, April 2018 (afrobert@nps.edu)

set(0,'DefaultTextInterpreter','Latex')

clear
close all

fonts=12;

notation=['ab'];

% set grey color tones
colormap('jet');
colmap=colormap(gray(2));
colmap(1,:)=0.65*[1 1 1];
colormap(colmap);

% read in image of Peter Wadham's paper
fname='DavisWadhamsFig16c.png';
floc='/Users/aroberts/science/publications/2015_Variational_Principle/Data/Wadhams_Data';
try
 cd(floc);
 im=imread(fname);
 im=im(1:end,1:end,1);
 im(end-50:end-5,5:50)=1;
catch
 disp('This function uses AGU copyrighted and cropped figure 16c from Davis and Wadhams')
 disp(['Unable to find ',fname,' at ',floc]);
 error('Figure 16c from Davis and Wadhams (1996) must be supplied with this function')
end

% switch off Wadhams tick marks (comment out this to check for alignment)
im(end-15:end,1:end)=1;
im(1:end,1:15)=1;

% set x and y coordinates
x=[1:size(im,1)]*15/size(im,1);
y=[size(im,2):-1:1]*90/size(im,2);

cls=ones([length(x) length(y) 3]);
for i=1:length(x)
for j=1:length(y)
 if abs(im(i,j))>0 & abs(im(i,j))<=255
  cls(i,j,:)=[1 1 1];
 else
  cls(i,j,:)=colmap(1,:);
 end
end
end

for ncol=1:2

 ridgepack_multiplot(1,2,1,ncol,[notation(ncol),')'])

 % plot wadhams data
 hi=image(x,y,cls);
 hold on

 hp=plot(1,100,'s','Color',colmap(1,:),'MarkerSize',10,'MarkerFaceColor',colmap(1,:));

 axis xy

 xmax=15;
 xmin=0;

 ylim([0 90])
 xlim([xmin xmax])

 set(gca,'Color','none');
 set(gca,'FontSize',fonts);

 cols(1,:)=[0 0 1];
 cols(2,:)=[0 0 1];
 cols(3,:)=[0 0 1];
 cols(4,:)=[0 0 1];
 cols(5,:)=[0 0 1];
 cols(6,:)=[0 0 1];
 cols(7,:)=[0 0 1];
 
 % calculate porosity along paths for hfi=2.0, hfs=0
 cont=[0:0.05:0.35];
 ridgepack_colormap(cont,0,'bluered');

 [ratio,alpha]=ridgeshape(0);
 [strainp,phip,alphap,VR,HK,HS,LK,LS]=ridgepack_trajectory(2.0,0.0);

 if ncol==1
  thetad=0;
 else
  thetad=55;
 end

 [strainp,phip,alphap,VR,HK,HS,LK,LS]=ridgepack_trajectory(2.0,0.0);

 [ratio1,alpha1]=ridgeshape(0);
 x1=ratio1;
 y1=squeeze(alpha1)';
 for i=1:length(alpha1)
  if alphap(end)>=alpha1(i)*cosd(thetad)
   idx=find(min(abs(alphap-alpha1(i)*cosd(thetad)))==abs(alphap-alpha1(i)*cosd(thetad)));
   phi1(i)=phip(idx);
  else
   phi1(i)=NaN;
  end
 end
 col1 = phi1';  % This is the color, vary with x in this case.

 [ratio2,alpha2]=ridgeshape(0.5);
 x2=ratio2;
 y2=squeeze(alpha2)';
 for i=1:length(alpha2)
  if alphap(end)>=alpha2(i)*cosd(thetad)
   idx=find(min(abs(alphap-alpha2(i)*cosd(thetad)))==abs(alphap-alpha2(i)*cosd(thetad)));
   phi2(i)=phip(idx);
  else
   phi2(i)=NaN;
  end
 end
 col2 = phi2';  % This is the color, vary with x in this case.

 z=zeros(size(x1));
 surface([x1 x2],[y1 y2],[z z],[col1 col2],'FaceAlpha',0.5);
 shading flat

 if ncol==2 
  ridgepack_colorbar(cont,'\phi_R')
 end

 % now plot the model
 [ratio,alpha]=ridgeshape(0);
 h1=plot(ratio(:),squeeze(alpha(:)),'color',cols(1,:),'LineWidth',1);

 [ratio,alpha]=ridgeshape(0.5);
 h2=plot(ratio(:),squeeze(alpha(:)),':','color',cols(1,:),'LineWidth',1);

 if ncol==2 
  legend([hp h1 h2],{'Sonar',...
   '$h_{Fd}=0$',...
   '$h_{Fd}=H_K/2$'},...
   'location','NorthEast',...
   'FontSize',fonts,...
   'Color','none',...
   'Interpreter','Latex')
  legend('boxoff')
 end

 xlabel('${L_K}/2({H_K}-{h_{Fd}})$','Interpreter','Latex','FontSize',fonts)
 text(2,83.5,['$\theta_R=',num2str(180-thetad),'^{\circ}$'],'FontSize',fonts)
 %title(['$\theta_R=',num2str(thetad),'^{\circ}$'],'FontSize',fonts)

 if ncol==1
  ylabel('$\alpha_K$ (degrees)','Interpreter','Latex')
 else
  set(gca,'YTickLabels',[])
  ylabel('$\;$','Interpreter','Latex','FontSize',2)
 end

 axis square

 set(gca,'BoxStyle','full','Layer','top')

end

ridgepack_multialign

% determine directory for read/write
dir=fileparts(which(mfilename));
cd([dir(1:strfind(dir,'scripts')-1),'output']);

% determine filename
x=strfind(mfilename,'_');
thisfilename=mfilename;
graphicsout=thisfilename(x(end)+1:end);

% output
disp(['Writing graphics output ',graphicsout,' to:',char(13),' ',pwd])

% print figure
ridgepack_fprint('epsc',graphicsout,1,2)


function [ratio,alpha]=ridgeshape(beta)

 ratio=[0.5:0.1:15]';

 for i=1:length(ratio)

  alpha(i)=acotd(ratio(i)/(1+beta));

 end

end
