% ridgepack_JAMES_figure16 - Generates Figure 16 in JAMES Variational Ridging paper
% 
% This script generates Figure 16 from:
%
% Roberts, A.F., E.C. Hunke, S.M. Kamal, W.H. Lipscomb, C. Horvat, 
% W. Maslowski (2019), Variational Method for Sea Ice Ridging in 
% Earth System Models, J. Adv. Model Earth Sy.
%
% The method used to integrate g(h,phi) is the "slow" method.  A faster
% computational method is available that uses convolusions (an FFT, in fact),
% and an example is available in the Mathematica script for a single ridge, 
% ridgepack_JAMES_figure13.nb
%
% VERSION/CASES: Ridgepack 1.0.1/JAMES_2019_VarRidging
%
% CONTACT: Andrew Roberts, afroberts@lanl.gov 
%
% FILE HISTORY:
% Author: Andrew Roberts, Naval Postgraduate School, April 2018 
% Update: Andrew Roberts, Los Alamos National Laboratory, December 2018

% version check
[v d] = version;
if str2num(d(end-3:end))<2018
 warning('This script designed for MATLAB 2018a or more recent version')
end

clear
close all

% set resolution of Epsilon (resolution 0.005 require for JAMES paper)
%resolution=0.001
resolution=0.005

% if reseting resolution, you must regenerate
%generate=false;
generate=true;

% determine directory for read/write of zeta-hat plane data
dir=fileparts(which(mfilename));
writedir=[dir(1:strfind(dir,'scripts')-1),'output'];
[status,msg]=mkdir(writedir);
cd(writedir);


% calculate or load zeta-hat plane
if generate
 disp(['Generating data in: ',char(13),' ',writedir])
 [HF,EPSILON,PHI,ALPHAHAT,VR,HK,HS,LK,LS]=ridgepack_zetahatplane(resolution);
 save('ridgepack_zetahatplane','HF','EPSILON','PHI','ALPHAHAT','VR','HK','HS','LK','LS')
else
 disp(['Reading from:',char(13),' ',writedir])
 try
  load ridgepack_zetahatplane
 catch
  error('Unable to load zeta-hat plane - switch to generate')
 end
end

% inititalize the thickness distribution, thickness grid, porosity grid and epsilon grid
[hincr,eincr,hgrid]=ridgepack_gridinit(resolution);

% use PHI on the zetahatplane as the porosity coordinate for ghphi
% This is done because it is numerically preferable in this circumstance
phigrid=PHI(1,:);
phincr(1)=diff(phigrid(1:2));
phincr(2:length(phigrid))=diff(phigrid);

% create 2-D arrays of increments for later use
[phincrs,hincrs]=meshgrid(phincr,hincr);

% initialize thickness distribution with zeros where ghphi is MxN 
% where M is the length of hgrid, and N is the length of phigrid
[ghinit] = 0*meshgrid(phigrid,hgrid);
[ghphi] = 0*meshgrid(phigrid,hgrid);

% print diagnostics for checking
disp(['Size of ghphi is ',num2str(size(ghphi))])

% set initial thickness distribution of sea ice field (this allows for extension
% to different initial shapes of thickness distirbutions if desired)
gshape='delta';
hinitial = 5.0;

% Case of the initial thickness distribution as a delta function on discrete grid
if strcmp(gshape,'delta')
 idx=find(min(abs(hgrid(:)-hinitial))==abs(hgrid(:)-hinitial));
 ghphi(idx,1,1)=1/(hincr(idx)*phincr(1));
else
 error('gshape specification not currently available')
end

clf

%integrate=false;
integrate=true;

bivariate=false;
%bivariate=true;

if integrate

 % intitial condition
 ghphihist(:,:,1)=ghphi;

 %for i=2:7
 for i=2:15

  disp(['Step: ',num2str(i)]) 

  [ghphi(:,:)]=ridgepack_redistribution(hgrid,hincr,phigrid,phincr,...
                     EPSILON,PHI,VR,HK,HS,LK,LS,ghphi);

  % history
  ghphihist(:,:,i)=ghphi;

  if bivariate

   clf
   surf(hgrid,phigrid,ghphi','FaceAlpha',0.75)
   shading flat
   set(gca,'Zscale','log')
   ylim([0 0.4])
   xlim([0 10])
   zlim([10^-8 10^4])
   view(135,30)

  else 

   % integrate through the porosity dimension
   % (equation 1 in Roberts et al. 2019)
   gh=squeeze(ghphihist(:,:,i)).*phincrs;
   gh=sum(gh/sum(gh(:)),2);

   % calculate thickness distribution
   semilogy(hgrid,gh)
   hold on

  end

  drawnow

 end

 save('ridgepack_ghphi','ghphihist')

else % read file to plot

 load ridgepack_ghphi

end

if bivariate

   % plot up 2D thickness distribution 
   clf

   ridgepack_multiplot(1,2,1,1,'a)')

   % draw initial redistribution
   i=2
   surf(hgrid,phigrid,squeeze(ghphihist(:,:,i))','FaceAlpha',0.75)
   shading flat
   hold on

   set(gca,'Zscale','log')
   ylim([0 0.4])
   xlim([0 10])
   zlim([10^-5 10^4])

   % draw initial delta function, including vertical arrow
   plot3([hinitial hinitial],[0 0],[10^-5 10^4],'b')
   text(hinitial,0,10^4,'$\uparrow$','Color','b',...
                        'HorizontalAlignment','center',...
                        'VerticalAlignment','cap')

   ylabel('$\phi$','Interpreter','Latex',...
          'HorizontalAlignment','center',...
          'Position',[13.5 0.25 10^-5],'fontsize',14,'Rotation',-08)
   zlabel('$g(h,\phi)$')
   set(gca,'XTick',[0:2:10])
   set(gca,'YTick',[0:0.10:0.4])
   set(gca,'ZTick',10.^[-4:2:4])
   grid on
   grid off
   grid on
   view(120,22)

   title('Initial Ridging over $A$')

   ridgepack_multiplot(1,2,1,2,'b)')

   i=4
   surf(hgrid,phigrid,squeeze(ghphihist(:,:,i))','FaceAlpha',0.75)
   shading flat
   hold on
   set(gca,'Zscale','log')
   ylim([0 0.4])
   xlim([0 10])
   zlim([10^-5 10^4])

   % draw initial delta function, including vertical arrow
   plot3([hinitial hinitial],[0 0],[10^-5 10^4],'b')
   text(hinitial,0,10^4,'$\uparrow$','Color','b',...
                        'HorizontalAlignment','center',...
                        'VerticalAlignment','cap')

   xlabel('$h$ (m)','Interpreter','Latex',...
          'HorizontalAlignment','center',...
          'Position',[7 0.5 10^-5],'fontsize',14,'Rotation',28)
   ylabel('$\phi$','Interpreter','Latex',...
          'HorizontalAlignment','center',...
          'Position',[13.5 0.25 10^-5],'fontsize',14,'Rotation',-08)

   %zlabel('$g(h,\phi)$')

   set(gca,'XTick',[0:2:10])
   set(gca,'YTick',[0:0.10:0.4])
   set(gca,'ZTick',10.^[-4:2:4])
   set(gca,'ZTickLabel',[])
   grid on
   grid off
   grid on
   view(120,22)

   title('Repeated Ridging over $A$')

   ridgepack_multialign(gcf)

   nameinfo='';

else 

   % plot up thickness distribution along thickness axis

   % integrate through the porosity dimension

   % apply mask
   idx=find(phigrid<0.35);

   clf

   ridgepack_multiplot(1,2,1,1,'a)')

   set(gca,'FontSize',12')

   k=0

   % pick of lines to use for JAMES paper
   if size(ghphihist,3)<9
    error('not enough to replicate JAMES paper')
   else
    % pick of lines to show
    yourpick=[1 4 6 7 9];

    % set colors of lines
    colines=colormap(lines(length(yourpick)));
    cols(1,:)=colines(1,:);
    cols(2,:)=colines(2,:);
    cols(3,:)=colines(end,:);
    cols(4,:)=colines(3,:);
    cols(5,:)=colines(end-1,:);
   end
   
   % (equation 1 in Roberts et al. 2019)
   for i=yourpick

    gh=squeeze(ghphihist(:,idx,i)).*phincrs(:,idx);
    gh=sum(gh/sum(gh(:)),2);

    % calculate thickness distribution
    k=k+1
    h(k)=plot(hgrid,gh,'Color',cols(k,:))
    hold on

   end

   xlim([0 10])

   ylim([0 2.5*10^-2])

   ylabel('$g(h)=\int\limits_{0}^{1} g(h,\phi_R)\;d\phi_R$')

   xlabel('h (m)')

   % fix up a funky order due to corrupted history file
   legend([h(1) h(end)],{'initial distribution','final distribution'},...
          'FontSize',12)

   legend('boxoff')

   ridgepack_multiplot(1,2,1,2,'b)')

   set(gca,'FontSize',12')

   k=0
   
   % (equation 1 in Roberts et al. 2019)
   for i=yourpick
   %for i=1:3:size(ghphihist,3)

    gh=squeeze(ghphihist(:,idx,i)).*phincrs(:,idx);
    gh=sum(gh/sum(gh(:)),2);

    % calculate thickness distribution
    k=k+1;
    plot(hgrid,gh,'Color',cols(k,:))

    hold on

   end

   set(gca,'YScale','log')

   xlim([0 10])

   ylim([10^-6 10^-1])

   xlabel('h (m)')

   ridgepack_multialign(gcf)

   nameinfo='';

end


% determine directory for read/write
dir=fileparts(which(mfilename));
outdir=[dir(1:strfind(dir,'scripts')-1),'output'];
[status,msg]=mkdir(outdir);
cd(outdir);

% determine filename
x=strfind(mfilename,'_');
thisfilename=mfilename;
graphicsout=[thisfilename(x(end)+1:end),nameinfo];

% output
disp(['Writing graphics output ',graphicsout,' to:',char(13),' ',pwd])

% print figure
ridgepack_fprint('epsc',graphicsout,1,2)
ridgepack_fprint('png',graphicsout,1,2)

