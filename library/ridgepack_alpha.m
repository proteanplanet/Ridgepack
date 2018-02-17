function [alphahat]=paper_ridge_alpha(epsilon,phi,hf,hd)

% function [alphahat]=paper_ridge_alpha(epsilon,phi,hf,hd)
%
% This function is part of Ridgepack Version 1.0.
% It calculates the angle of repose for which the gain is potential 
% energy density is minimized for a ridge given Coulombic deformation
% within a ridge.
% 
% INPUT:
% 
% epsilon - ridge strain (dimensionless)
% phi     - ridge porosity (dimensionless)
% hf      - floe ice thickness (m)
% hd      - deformed ice thickness (m)
%
%
% OUTPUT:
%
% alphahat - angle of repose observing the Principle of Virtual Work (degrees)
%
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

% calculate angle of repose 
if hd>hf

 gamma=(phi+phi*epsilon-epsilon);

 alphahat=2*atan(sqrt(((5*(gamma.^2)+6*gamma+3)-...
              sqrt(((5*(gamma.^2)+6*gamma+3)^2)-4*(gamma.^4)))/(2*(gamma.^2))));
 
 alphahat=180*alpha/pi;

else

 alphahat=0;

end

