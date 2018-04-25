function [HF,EPSILON,PHI,ALPHAHAT,VR,HK,HS,LK,LS]=ridgepack_zetahatplane

% RIDGEPACK_ZETAHATPLANE - Calculate ridge state on zeta-hat plane
%
% function [HF,EPSILON,PHI,ALPHAHAT,VR,HK,HS,LK,LS]=ridgepack_zetahatplane
%
% This calculates strain, porosity and the angle of repose of ridges
% on the the zeta hat trajectory plane of potential energy density. 
% The entire ridge state space may be determined from these variables,
% and is output as part of the calculations.
%
% OUTPUT:
%
% HF       - thickness of initial ice (m, length N)
% EPSILON  - strain of ridge (dimensionless, length M)
% PHI      - porosity of a ridge (dimensionless, size MxN)
% ALPHAHAT - angle of repose of a ridge (degrees, size MxN)
% VR       - potential energy density of a ridge (J m^-2, size MxN)
% HK       - draft of a keel (m, size MxN)
% HS       - height of a sail (m, size MxN)
% LK       - cross-sectional width of a keel (m, size MxN)
% LS       - cross-sectional width of a sail (m, size MxN)
%
% Please note that this function is dependent on snow cover, but for the 
% purpose of the paper that Ridgepack Version 1.0 supports, this dependency
% has been removed.
%
% Ridgepack Version 1.0.
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

% initialze grid
[hincr,eincr,HF]=ridgepack_gridinit;

% remove snow for Version 1.0
HFs=zeros(size(HF));

% step through initial ice thicknesses
for i=1:length(HF) 

 % calculate zeta-hat trajectory for ice and snow thickness HF and HFs
 [EPSILON,PHI(i,:),ALPHAHAT(i,:),VR(i,:),HK(i,:),HS(i,:),LK(i,:),LS(i,:)]=...
      ridgepack_trajectory(HF(i),HFs(i));

end

