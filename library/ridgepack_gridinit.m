function [hgrid,epsilongrid,phigrid,epsilonsplit,phisplit,ghphi]=ridgepack_gridinit


% epsilon  - ridge strain coordinate on alphahat plane (dimensionless)
% phi      - ridge macroporosity coordinate on alphahat plane (dimensionless)


global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% retrieve constants
[rhoi,rhos,rhow,delrho,g,eincr,hincr,minthick,maxthick]=ridgepack_constants;

% check resolution settings
if eincr>0.01
 error('eincr too large for convergence')
end

% Ice thickness grid
hgrid=10.^[log10(0.01):hincr:log10(maxthick)];

% strain grid and splot grid
epsilongrid=[-eincr:-eincr:-0.99];
epsilonsplit=[-eincr-eincr/2:-eincr:-0.99+eincr/2];

% Porosity grid
phigrid=[eincr:eincr:0.99];
phisplit=[eincr+eincr/2:eincr:0.99-eincr/2];

% initialize thickness distribution with zeros
[ghphi]=0*meshgrid(epsilonsplit,hgrid);

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


