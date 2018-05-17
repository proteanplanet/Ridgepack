function [HINCR,EINCR,HGRID,EPSILONGRID,PHIGRID,EPSILONSPLIT,PHISPLIT,GHPHI]=...
            ridgepack_gridinit

% RIDGEPACK_GRIDINIT - set up the initial thickness, strain and porosity grid 
%
% function [HINCR,EINCR,HGRID,EPSILONGRID,PHIGRID,EPSILONSPLIT,PHISPLIT,GHPHI]=...
%             ridgepack_gridinit
% 
% OUTPUT:
%
% HINCR        - increment in thickness grid (m)
% EINCR        - resolution of strain and porosity (dimensionless)
% HGRID        - thickness grid (m)
% EPSILONGRID  - strain grid (dimensionless)
% PHIGRID      - porosity grid (dimensionless)
% EPSILONSPLIT - split strain grid (dimensionless)
% PHISPLIT     - split porosity grid (dimensionless)
% GHPHI        - initialized thickness distribition g(h,phi)
%
% Ridgepack Version 1.0.
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% resolution of epsilon and phi in calculating zeta hat plane (dimensionless)
EINCR = 0.001;
%EINCR = 0.005;
%EINCR = 0.01;

% minimim abs(strain) and porosity
minstrain = 0.01;

% maximim abs(strain) and porosity
maxstrain = (1-EINCR)-EINCR-minstrain;

% minimum thickness on zeta-hat plane trajectory plane (m)
minthick = 0.01;

% maximum thickness on zeta-hat plane trajectory plane (m)
maxthick = 100;

% check resolution settings
if EINCR>0.01
 error('EINCR too large for convergence')
end

% flag for logarithmic thickness grid
%logthickness=false;
logthickness=true;

if logthickness
 % set initial log thickness resolution on zeta-hat plane (m)
 HINCR = (log10(10)-log10(minthick))/1000;

 % ice thickness grid
 HGRID = [10.^[log10(minthick):HINCR:log10(maxthick)]];
else
 % set linear thickness resolution on zeta-hat plane (m)
 HINCR = 0.1;

 HGRID = [minthick:HINCR:maxthick];
end

% determine final increment
HINCR(1)=HGRID(1);
HINCR(2:length(HGRID))=diff(HGRID);

% strain grid and splot grid
EPSILONGRID = [-minstrain:-EINCR:-maxstrain];
EPSILONSPLIT = [-minstrain-EINCR/2:-EINCR:-maxstrain+EINCR/2];

% porosity grid
PHIGRID = [minstrain:EINCR:maxstrain];
PHISPLIT = [minstrain+EINCR/2:EINCR:maxstrain-EINCR/2];

% initialize thickness distribution with zeros
[GHPHI] = 0*meshgrid(PHISPLIT,HGRID);

% debug information
if debug
 disp(['size of HGRID is ',num2str(size(HGRID))])
 disp(['size of EPSILONGRID is ',num2str(size(EPSILONGRID))])
 disp(['size of EPSILONSPLIT is ',num2str(size(EPSILONSPLIT))])
 disp(['size of PHIGRID is ',num2str(size(PHIGRID))])
 disp(['size of PHISPLIT is ',num2str(size(PHISPLIT))])
 disp(['size of GHPHI is ',num2str(size(GHPHI))])
end

if debug; disp(['...Leaving ',mfilename]); end

