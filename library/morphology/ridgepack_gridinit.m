function [HINCR,EINCR,HGRID,EPSILONGRID,PHIGRID,EPSILONSPLIT,PHISPLIT,GHPHI]=...
            ridgepack_gridinit(resolution,minstrain,maxstrain)

% ridgepack_gridinit - set up the initial thickness, strain and porosity grid 
%
% function [HINCR,EINCR,HGRID,EPSILONGRID,PHIGRID,EPSILONSPLIT,PHISPLIT,GHPHI]=...
%             ridgepack_gridinit(resolution,minstrain,maxstrain)
%
% Initializes the numerical grid for calculating the state-space trajectory
% of a ridge and for integrating thickness distribution changes. The thickness
% distribution GHPHI is also initialized.
%
% INTPUT:
%
% resolution - resolution of the strain and porosity grid with typical 
%              values between 0.01 and 0.001 (optional)
% minstrain  - minimum absolute strain in the epsilon grid (optional)
% maxstrain  - maximum absolute strain on the epsilon grid (optional)
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
% REFERENCE: 
%
% Roberts, A., E. Hunke, S. Kamal, W. Lipscomb, C. Horvat, W. Maslowski (2019),
% A Variational Method for Sea Ice Ridging in Earth System Models, 
% J. Adv. Model Earth Sy. 
% 
% VERSION/LIBRARY: Ridgepack 1.0.1/MORPHOLOGY
%
% CONTACT: Andrew Roberts, afroberts@lanl.gov 
%
% FILE HISTORY:
% Author: Andrew Roberts, Naval Postgraduate School, March 2018 
% Reviewed: by Samy Kamal, Naval Postgraduate School, May 2018
% Update: Andrew Roberts, Los Alamos National Laboratory, December 2018

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% resolution of epsilon and phi in calculating zeta hat plane (dimensionless)
if nargin==1
 EINCR = resolution;
else
 EINCR = 0.001;
 %EINCR = 0.005;
 %EINCR = 0.01;
end

if debug
 if EINCR<0.001
  disp('WARNING: This could be computationally expensive for little gain')
 elseif EINCR>0.01
  disp('WARNING: This may produce unconverged results')
 end
end

% minimim abs(strain) and porosity
if nargin<2
 minstrain = EINCR;
elseif isnumeric(minstrain)
 if maxstrain<0 | maxstrain>1
  error('minstrain is the absolute value and has limits 0 < minstrain < 1')
 else
  disp(['Setting minimum strain manually to ',num2str(-minstrain)])
 end
else
 error('minstrain is not numeric')
end

% maximim abs(strain) and porosity
if nargin<3
 maxstrain = (1-EINCR)-EINCR-minstrain;
elseif isnumeric(maxstrain)
 if maxstrain<0 | maxstrain>1
  error('maxstrain is the absolute value and has limits 0 < maxstrain < 1')
 elseif maxstrain<=minstrain
  error('maxstrain is less than minstrain')
 else
  disp(['Setting maximum strain manually to ',num2str(-maxstrain)])
 end
else
 error('maxstrain is not numeric')
end

% minimum thickness on zeta-hat plane trajectory plane (m)
minthick = 0.01;

% maximum thickness on zeta-hat plane trajectory plane (m)
maxthick = 100;

% flag for logarithmic thickness grid
logthickness=false;
%logthickness=true;

if logthickness
 % set initial log thickness resolution on zeta-hat plane (m)
 HINCR = (log10(10)-log10(minthick))/1000;

 % ice thickness grid
 HGRID = [10.^[log10(minthick):HINCR:log10(maxthick)]];
else
 % set linear thickness resolution on zeta-hat plane (m)
 HINCR = 0.001;

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
 disp(['...Leaving ',mfilename]); end
end

