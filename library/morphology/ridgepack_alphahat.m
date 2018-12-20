function [ALPHAHAT]=ridgepack_alphahat(epsilon,phi,hf,hr)

% ridgepack_alphahat - Calculate angle of repose in the alpha-hat manifold
%  
% function [ALPHAHAT]=ridgepack_alphahat(epsilon,phi,hf,hr)
%
% This calculates the angle of repose for which the gain in potential 
% energy density is minimized for a ridge given Coulombic deformation
% within a ridge, based on strain, porosity, floe ice and deformed ice
% thickness. Calculations are from Appendix D, Roberts et al. (2019)
% 
% INPUT:
% 
% epsilon - ridge strain with value range (-1,0) (dimensionless)
% phi     - ridge porosity (dimensionless)
% hf      - floe ice thickness (m)
% hr      - ridged ice thickness (m)
%
% OUTPUT:
%
% ALPHAHAT - angle of repose observing the Principle of Virtual Work (deg)
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

% check there are sufficient inputs
if nargin~=4
 error('incorrect number of inputs')
end

% calculate angle of repose (Equation D.5, Roberts et al. 2019)
if hr>hf

 gamma=(phi+phi.*epsilon-epsilon);

 ALPHAHAT=2*atan(sqrt(((5*(gamma.^2)+6*gamma+3)-...
            sqrt(((5*(gamma.^2)+6*gamma+3).^2)-4*(gamma.^4)))./(2*(gamma.^2))));
 
 ALPHAHAT=180*ALPHAHAT/pi;

elseif hr<hf

 error('Deformed ice mean thickness is less than floe ice mean hickness (hr<hf)')

else

 ALPHAHAT=0;

end

if debug
  disp(['...Leaving ',mfilename]);
end


