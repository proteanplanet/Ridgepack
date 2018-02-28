function [VR,ALPHAHAT,HK,HS,LK,LS]=ridgepack_energetics(hf,hfs,epsilon,phi)

% RIDGEPACK_ENERGETICS - calculate potential energy density of a ridge
%
% function [VR,ALPHAHAT,HK,HS,LK,LS]=ridgepack_energetics(hf,hfs,epsilon,phi)
% 
% This calculates the potential energy density, angle of repose, keel depth
% sail height, and cross-sectional keel and sail widths given the floe ice 
% and snow thickness, and ridge strain and porosity. 
%
% INPUT:
%
% hf      - floe ice thickness (m)
% hfs     - floe snow thickness (m)
% epsilon - ridge strain (dimensionless)
% phi     - ridge porosity (dimensionless)
%
% 
% OUTPUT:
%
% VR       - potential energy density of a ridge (J m^-2)
% ALPHAHAT - angle of repose of a ridge (degrees)
% HK       - draft of a keel (m)
% HS       - height of a sail (m)
% LK       - cross-sectional width of a keel (m)
% LS       - cross-sectional width of a sail (m)
%
% Ridgepack Version 1.0.
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

% check there are sufficient inputs
if nargin~=4
 error('incorrect number of inputs')
elseif any(size(epsilon)~=size(phi))
 error('size of epsilon not the same as phi')
elseif length(hf)>1 | length(hfs)>1
 error('length of hf or hfs is greater than 1')
end

% get constants
[rho,rhos,rhow,delrho,ghat]=ridge_constants;

% calculate thickness of deformed ice mass from strain
hd=hf./(1+epsilon);

% set snow thickness on ridge same as on level ice 
hds=hfs;

% calculate angle of repose of the ridge
ALPHAHAT=ridgepack_alphahat(epsilon,phi,hf,hd);

% get the morphological shape of the ridge
[hfd,hff,hdd,hdf,HK,HS,LK,LS]=ridgepack_morphology(hf,hfs,hd,hds,phi,ALPHAHAT);

% calculate potential energy density 
VR = delrho*ghat*(1-phi).*(hfd.*LK+0.25*(LK.^2).*tand(ALPHAHAT)) + ...
        rho*ghat*(1-phi).*(hff.*LK+0.25*(LS.^2).*tand(ALPHAHAT));


