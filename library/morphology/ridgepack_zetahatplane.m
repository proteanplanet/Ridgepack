function [HF,EPSILON,PHI,ALPHAHAT,VR,HK,HS,LK,LS]=ridgepack_zetahatplane(res)

% ridgepack_zetahatplane - Calculate ridge state on zeta-hat plane
%
% function [HF,EPSILON,PHI,ALPHAHAT,VR,HK,HS,LK,LS]=ridgepack_zetahatplane(res)
%
% This calculates strain, porosity and the angle of repose of ridges
% on the the zeta hat trajectory plane of potential energy density. 
% The entire ridge state space may be determined from these variables,
% and is output as part of the calculations.
%
% INPUT:
%
% hf  - parent ice thickness (m)
% hfs - thickness of snow on parent ice (m)
% res - resolution of the strain and porosity grid (optional)
%       with typical values between 0.01 and 0.001.
%
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
% Reviewed by Samy Kamal, Naval Postgraduate School, May 2018

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if nargin<1
 res=0.005;
elseif res>0.01
 error('Resolution is incorretly set')
end

% initialze grid
[hincr,eincr,HF]=ridgepack_gridinit(res);

% remove snow for Version 1.0 (this step can be removed in future versions)
HFs=zeros(size(HF));

% step through initial ice thicknesses
tic;
for i=1:length(HF) 

 % calculate zeta-hat trajectory for ice and snow thickness HF and HFs
 [EPSILON(i,:),PHI(i,:),ALPHAHAT(i,:),VR(i,:),HK(i,:),HS(i,:),LK(i,:),LS(i,:)]=...
      ridgepack_trajectory(HF(i),HFs(i),res);

end
elapsedtime=toc;

EPSILON=EPSILON(1,:);

disp(['Zetahat plane calculated at ',num2str(res),' resolution in ',...
      num2str(elapsedtime),' seconds'])

if debug; disp(['...Leaving ',mfilename]); end

