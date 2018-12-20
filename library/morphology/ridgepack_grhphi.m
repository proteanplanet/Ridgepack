function [GRHPHI,hgrid]=ridgepack_grhphi(hF,hFs,epsilon,phi)

% ridgepack_grhphi - gR distribution for strain, porosity and initial conditions
%
% function [GRHPHI,hgrid]=ridgepack_grhphi(hF,hFs,epsilon,phi)
%
% This function calculates the thickness distribution step function GRHPHI on the 
% thickness axis given the initial sea ice thickness, hF, snow thickness on 
% top, hFs, as well as strain and porosity of the ice in the given ridge. If the 
% step function exceeds the thickness grid, it is returned as NaNs.
%
% INPUT:
% 
% hF      - parent sheet sea ice thickness (m) 
% hFs     - parent sheet snow thickness on sea ice (m)
% epsilon - ridge strain with value range (-1,0) (dimensionless)
% phi     - ridge macroporosity with value range (0,1)
%
%
% OUTPUT:
%
% GRHPHI - step thickness distribution on hgrid with porosity phi
% hgrid  - thickness grid (m)
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

% error check
if nargin~=4
 error('incorrect number of inputs')
elseif length(hF)>1 | length(hFs)>1 
 size(hF)
 size(hFs)
 error('length of hF or hFs is greater than 1')
elseif length(epsilon)>1 | length(phi)>1 
 error('length of epsilon or phi is greater than 1')
end

% initialize
[hincr,eincr,hgrid]=ridgepack_gridinit;
GRHPHI=zeros(size(hgrid));

% get dimensions of ridge based on initial conditions, epsilon and porosity
[VR,ALPHAHAT,HK,HS,LK,LS]=ridgepack_energetics(hF,hFs,epsilon,phi);

% Find indices for of steps on thickness distribution grid
h1idx=find(hgrid>=hF & hgrid<=hF+(LK-LS)*tand(ALPHAHAT)/2);
h2idx=find(hgrid>hF+(LK-LS)*tand(ALPHAHAT)/2 & hgrid<=hF+(LK+LS)*tand(ALPHAHAT)/2);

% check you are within range thickness, otherwise hand back no distribution
if ((hF+(LK+LS)*tand(ALPHAHAT)/2)>hgrid(end))

 disp(['Reached end of thickness grid, epsilon:',num2str(epsilon)])
 GRHPHI=NaN*GRHPHI;

else

 % Thickness distribution
 GRHPHI(h1idx) = 2/(LK*tand(ALPHAHAT));
 GRHPHI(h2idx) = 1/(LK*tand(ALPHAHAT));

 % Normalize to account for numerical inaccuracy of given ice thickness grid
 area=sum(hincr.*GRHPHI);
 GRHPHI(h1idx)=GRHPHI(h1idx)/area;
 GRHPHI(h2idx)=GRHPHI(h2idx)/area;

end

if debug
 disp(['Total area of gR is: ',num2str(sum(hincr.*GRHPHI))]);
 disp(['...Leaving ',mfilename]); 
end

