function [ALPHAHAT]=ridgepack_alphahat(epsilon,phi,hf,hd)

% RIDGEPACK_ALPHAHAT - Calculate angle of repose in the alpha-hat manifold
%  
% function [ALPHAHAT]=ridgepack_alphahat(epsilon,phi,hf,hd)
%
% This calculates the angle of repose for which the gain in potential 
% energy density is minimized for a ridge given Coulombic deformation
% within a ridge, based on strain, porosity, floe ice and deformed ice
% thickness.
% 
% INPUT:
% 
% epsilon - compressive ridge strain (dimensionless)
% phi     - ridge porosity (dimensionless)
% hf      - floe ice thickness (m)
% hd      - deformed ice thickness (m)
%
%
% OUTPUT:
%
% ALPHAHAT - angle of repose observing the Principle of Virtual Work (degrees)
%
% Ridgepack Version 1.0.
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

% check there are sufficient inputs
if nargin~=4
 error('incorrect number of inputs')
end

% calculate angle of repose 
if hd>hf

 gamma=(phi+phi.*epsilon-epsilon);

 ALPHAHAT=2*atan(sqrt(((5*(gamma.^2)+6*gamma+3)-...
            sqrt(((5*(gamma.^2)+6*gamma+3).^2)-4*(gamma.^4)))./(2*(gamma.^2))));
 
 ALPHAHAT=180*ALPHAHAT/pi;

elseif hd<hf

 error('Deformed ice mean thickness is less than floe ice mean hickness (hd<hf)')

else

 ALPHAHAT=0;

end

