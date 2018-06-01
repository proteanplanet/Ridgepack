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

% determine directory for read/write of zeta-hat plane data
dir=fileparts(which(mfilename));
writedir=[dir(1:strfind(dir,'scripts')-1),'output'];
[status,msg]=mkdir(writedir);
cd(writedir);

generate=false;
%generate=true;

% calculate or load zetahat plane
if generate
 disp(['Generating data in: ',char(13),' ',writedir])
 [HF,EPSILON,PHI,ALPHAHAT,VR,HK,HS,LK,LS]=ridgepack_zetahatplane(resolution);
 save('ridgepack_zetahatplane','HF','EPSILON','PHI','ALPHAHAT','VR','HK','HS','LK','LS')
else
 disp(['Reading from:',char(13),' ',writedir])
 try
  load ridgepack_zetahatplane
 catch
  disp('Unable to load zeta-hat plane - switch to generate')
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

epsilondot=-10^-6;
dt=10;

clf

% intitial condition
ghphihist(:,:,1)=ghphi;

for i=2:10

 i

 [ghphi(:,:)]=ridgepack_redistribution(hgrid,hincr,phigrid,phincr,...
                     EPSILON,PHI,VR,HK,HS,LK,LS,ghphi,epsilondot,dt);

 % history
 ghphihist(:,:,i)=ghphi;

if i==0

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

return

clf
surf(hgrid,phigrid,ghphi','FaceAlpha',0.75)
shading flat
set(gca,'Zscale','log')
ylim([0 0.4])
xlim([0 10])
zlim([10^-7 10^5])
view(135,30)


return

% integrate through the porosity dimension
gh=ghphi.*phincrs;
gh=sum(gh/sum(gh(:)),2);

% calculate thickness distribution
clf
semilogy(hgrid,gh)
hold on
drawnow

semilogy(hgrid,sum(ghphi,2))
hold on
drawnow

legend({'initial distribution','final distribution'})

% plot up thickness distribution along thickness axis




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

