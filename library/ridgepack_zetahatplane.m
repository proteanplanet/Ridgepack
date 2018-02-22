function [HF,EPSILON,PHI,ALPHAHAT,VR,HK,HS,LK,LS]=ridgepack_zetahatplane

% function [HF,EPSILON,PHI,ALPHAHAT,VR,HK,HS,LK,LS]=ridgepack_zetahatplane
%
% This calculates strain, porosity and the angle of repose of ridges
% on the the zeta hat trajectory plane of potential energy density. 
% The entire ridge state space may be determined from these variables,
% and is output as part of the calculations.
%
% OUTPUT:
%
% HF       - thickness of initial ice (m)
% EPSILON  - strain of ridge (dimensionless)
% PHI      - porosity of a ridge (dimensionless)
% ALPHAHAT - angle of repose of a ridge (degrees)
% VR       - potential energy density of a ridge (J m^-2)
% HK       - draft of a keel (m)
% HS       - height of a sail (m)
% LK       - cross-sectional width of a keel (m)
% LS       - cross-sectional width of a sail (m)
%
% Ridgepack Version 1.0.
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

% retrieve constants
[rhoi,rhos,rhow,delrho,g,eincr,hincr,minthick,maxthick]=ridgepack_constants;

% set HF grid from 1cm to maxthick thick ice
HF=10.^[log10(0.01):hincr:log10(maxthick)];

% assume there is no snow in v1 (although the library allows it)
HFs=zeros(size(HF));

% check resolution settings
if eincr>0.01
 error('eincr too large for convergence')
end

% create strain and porosity grid
epsiloni=[-eincr:-eincr:-0.99];
phii=[eincr:eincr:0.50];
[epsilon,phi]=meshgrid(epsiloni,phii);

% step through initial ice thicknesses
for i=1:length(HF) 

 % calculate zeta-hat trajectory for ice and snow thickness HF and HFs
 [EPSILON,PHI(i,:),ALPHAHAT(i,:),VR(i,:),HK(i,:),HS(i,:),LK(i,:),LS(i,:)]=...
      ridgepack_trajectory(HF(i),HFs(i),epsilon,phi);

end

