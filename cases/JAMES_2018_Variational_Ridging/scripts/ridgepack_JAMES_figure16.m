% ridgepack_JAMES_figure16 - Generates Figure 16 in JAMES Variational Ridging paper
% 
% This script generates Figure 16 from:
%
% Roberts, A.F., E.C. Hunke, S.M. Kamal, W.H. Lipscomb, C. Horvat, W. Maslowski (2018),
% Variational Method for Sea Ice Ridging in Earth System Models, Part I: Theory, 
% submitted to J. Adv. Model Earth Sy.
%
% Andrew Roberts, Naval Postgraduate School, April 2018 (afrobert@nps.edu)

% version check
[v d] = version;
if str2num(d(end-3:end))<2018
 warning('This script designed for MATLAB 2018a or more recent version')
end

clear
clf
%close all

% set resolution of Epsilon
resolution=0.001
%resolution=0.005
%resolution=0.01

% if reseting resolution, you must regenerate
generate=false;
%generate=true;

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
hinitial = 1.0;

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

 for i=2:7

  disp(['Step: ',num2str(i)]) 

  [ghphi(:,:)]=ridgepack_redistribution(hgrid,hincr,phigrid,phincr,...
                     EPSILON,PHI,VR,HK,HS,LK,LS,ghphi);

  % history
  ghphihist(:,:,i)=ghphi;
  size(ghphihist)

  if bivariate

   clf
   surf(hgrid,phigrid,ghphi','FaceAlpha',0.75)
   shading flat
   set(gca,'Zscale','log')
   ylim([0 0.4])
   xlim([0 10])
   zlim([10^-10 10^2])
   view(135,30)

  else 

   % integrate through the porosity dimension
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

return

if bivariate

   % plot up 2D thickness distribution 
   clf
   i=3
   surf(hgrid,phigrid,squeeze(ghphihist(:,:,i))','FaceAlpha',0.75)
   shading flat
   set(gca,'Zscale','log')
   ylim([0 0.4])
   xlim([0 10])
   zlim([10^-5 10^2])
   view(135,30)

else 

   % plot up thickness distribution along thickness axis

   % integrate through the porosity dimension
   
   for i=1:size(ghphihist,3)

    gh=squeeze(ghphihist(:,:,i)).*phincrs;
    gh=sum(gh/sum(gh(:)),2);

    % calculate thickness distribution
    semilogy(hgrid,gh)
    hold on

   end

   ylim([10^-7 10^-1])


   

  %legend({'initial distribution','final distribution'})

end

return


% determine directory for read/write
dir=fileparts(which(mfilename));
outdir=[dir(1:strfind(dir,'scripts')-1),'output'];
[status,msg]=mkdir(outdir);
cd(outdir);

% determine filename
x=strfind(mfilename,'_');
thisfilename=mfilename;
graphicsout=thisfilename(x(end)+1:end);

% output
disp(['Writing graphics output ',graphicsout,' to:',char(13),' ',pwd])

% print figure
ridgepack_fprint('epsc',graphicsout,1,2)
ridgepack_fprint('png',graphicsout,1,2)

