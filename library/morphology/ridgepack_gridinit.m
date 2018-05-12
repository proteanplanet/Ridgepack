function [hincr,eincr,hgrid,epsilongrid,phigrid,epsilonsplit,phisplit,ghphi]=...
            ridgepack_gridinit


% epsilon  - ridge strain coordinate on alphahat plane (dimensionless)
% phi      - ridge macroporosity coordinate on alphahat plane (dimensionless)


global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% resolution of epsilon and phi in calculating zeta hat plane (dimensionless)
eincr = 0.001;
%eincr = 0.005;
%eincr = 0.01;

% minimim abs(strain) and porosity
minstrain = 0.01;

% maximim abs(strain) and porosity
maxstrain = 0.99-eincr-minstrain;

% minimum thickness on zeta-hat plane trajectory plane (m)
minthick = 0.01;

% maximum thickness on zeta-hat plane trajectory plane (m)
maxthick = 100;

% check resolution settings
if eincr>0.01
 error('eincr too large for convergence')
end

%logthickness=false;
logthickness=true;
if logthickness
 % set initial log thickness resolution on zeta-hat plane (m)
 hincr = (log10(10)-log10(minthick+0.01))/1000;

 % ice thickness grid
 hgrid = [minthick 10.^[log10(0.01):hincr:log10(maxthick)]];
else
 % set linear thickness resolution on zeta-hat plane (m)
 hincr = 0.1;

 hgrid = [minthick hincr:hincr:maxthick];
end


% determine final increment
hincr(1)=hgrid(1);
hincr(2:length(hgrid))=diff(hgrid);

% strain grid and splot grid
epsilongrid = [-minstrain:-eincr:-maxstrain];
epsilonsplit = [-minstrain-eincr/2:-eincr:-maxstrain+eincr/2];

% porosity grid
phigrid = [minstrain:eincr:maxstrain];
phisplit = [minstrain+eincr/2:eincr:maxstrain-eincr/2];

% initialize thickness distribution with zeros
[ghphi] = 0*meshgrid(phisplit,hgrid);

% debug information
if debug
 disp(['size of hgrid is ',num2str(size(hgrid))])
 disp(['size of epsilongrid is ',num2str(size(epsilongrid))])
 disp(['size of epsilonsplit is ',num2str(size(epsilonsplit))])
 disp(['size of phigrid is ',num2str(size(phigrid))])
 disp(['size of phisplit is ',num2str(size(phisplit))])
 disp(['size of ghphi is ',num2str(size(ghphi))])
end

if debug; disp(['...Leaving ',mfilename]); end

