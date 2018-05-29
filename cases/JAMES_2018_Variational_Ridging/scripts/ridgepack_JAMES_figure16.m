% ridgepack_JAMES_figure5 - Generates Figure 16 in JAMES Variational Ridging paper
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
close all


% set resolution of Epsilon
resolution=0.01

% inititalize the thickness distribution, thickness grid, porosity grid and epsilon grid
[hincr,eincr,hgrid,epsilongrid,phigrid,epsilonsplit,phisplit,ghphi]=...
      ridgepack_gridinit(resolution);

% print diagnostics for checking
disp(['Size of ghphi is ',num2str(size(ghphi))])

% set initial thickness distribution of sea ice field (this allows for extension
% to different initial shapes of thickness distirbutions if desired)
gshape='delta';
hinitial = 2.0;

% Case of the initial thickness distribution as a delta function on discrete grid
if strcmp(gshape,'delta')
 idx=find(min(abs(hgrid(:)-hinitial))==abs(hgrid(:)-hinitial));
 ghphi(idx,1)=hinitial/hincr(idx);
else
 error('gshape specification not currently available')
end

epsilondot=-10^-6;
dt=10;

[ghphinew]=ridgepack_redistribution(ghphi,resolution,epsilondot,dt);


return

% calculate thickness distribution
for i=1:100
 [ghphi]=ridgepack_redistribution(ghphi,hgrid,phisplit);
 disp(num2str(i))
 if i==1
  clf
  semilogy(hgrid,sum(ghphi,2))
  hold on
  drawnow
 end
end

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

