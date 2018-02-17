function [VR,ALPHAHAT,HK,HS,LK,LS]=ridgepack_energetics(hf,hfs,stii,phii)

% function [VR,ALPHA,HK,HS,LK,LS]=ridgepack_energetics(hf,hfs,stii,phii)
% 
% This function is part of Ridgepack Version 1.0.
% It calculates the potential energy density, angle of repose, keel depth
% sail height, cross-sectional keel and sail widths given the floe ice 
% and snow thickness, and ridge strain and porosity. 
%
% INPUT:
%
% hf   - floe ice thickness (m)
% hfs  - floe snow thickness (m)
% stii - ridge strain (dimensionless)
% phii - ridge porosity (dimensionless)
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
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

% get constants
[rho,rhos,rhow,delrho,g]=ridge_constants;

% calculate thickness of deformed ice mass from strain
hdi=hfi./(1+stii);

% set snow thickness on ridge same as on level ice 
hds=hfs;

% calculate angle of repose of the ridge
ALPHA=ridgepack_alpha(stii,phii,hfi,hdi);

% get the morphological shape of the ridge
[hfd,hff,hdd,hdf,HK,HS,LK,LS]=ridge_morphology(hfi,hfs,hdi,hds,phii,ALPHA);

% calculate potential energy density 
VR = delrho*g*(1-phii)*(hfd*LK+0.25*(LK.^2)*tand(ALPHA)) + ...
        rho*g*(1-phii)*(hff*LK+0.25*(LS.^2)*tand(ALPHA));


