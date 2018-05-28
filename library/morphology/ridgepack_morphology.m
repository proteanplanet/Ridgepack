function [HFD,HFF,HRD,HRF,HK,HS,LK,LS,L0,EPSILON]=...
                    ridgepack_morphology(hf,hfs,hr,hrs,phi,alpha)

% ridgepack_morphology - Calculate ridge dimension and strain
%
% function [HFD,HFF,HRD,HRF,HK,HS,LK,LS,L0,EPSILON]=...
%                    ridgepack_morphology(hf,hfs,hr,hrs,phi,alpha)
%
% This calculates the floe draft and freeboard, mean ridge draft and freeboard,
% keel depth, sail height, undeformed cross-sectional length of a floe,
% and strain based upon the floe ice and snow thickness, ridge mean
% ice and snow thickness, macroporosity and angle of repose.
% 
% INPUT:
% 
% hf    - floe ice thickness (m)
% hfs   - floe snow thickness (m)
% hr    - ridge mean ice thickness (m)
% hrs   - ridge mean snow thickness (m)
% phi   - macroporosity of ridge (dimensionless)
% alpha - angle of repose of a ridge (degrees)
%
%
% OUTPUT:
%
% HFD     - floe draft (m)
% HFF     - floe freeboard (m)
% HRD     - ridge mean draft (m)
% HRF     - ridge mean freeboard (m)
% HK      - keel depth (m)
% HS      - sail height (m)
% LK      - cross-sectional keel width (m)
% LS      - cross-sectional sail width (m)
% L0      - original cross-sectional width of undeformed parent ice (m)
% EPSILON - compressional strain or ridge
% 
% Ridgepack Version 1.0.
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% Reviewed by Samy Kamal, Naval Postgraduate School, May 2018

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% check there are sufficient inputs
if nargin~=6
 error('incorrect number of inputs')
end

% error checking
if any(hr<hf)
 error('ERROR: hr<hf')
elseif any(hrs<0)
 error('ERROR: hrs<0')
elseif any(hf<0)
 error('ERROR: hf<0')
elseif any(hr<0)
 error('ERROR: hr<0')
end

% get constants
hc=ridgepack_astroconstants;
rho=hc.rhoi.const;  % density of ice (kg/m^3)
rhos=hc.rhos.const; % density of snow (kg/m^3)
rhow=hc.rhow.const; % density of seawater (kg/m^3)
delrho=hc.rhow.const-hc.rhoi.const; % difference of water and ice densities
ghat=hc.ghat.const; % acceleration due to gravity

% calculate freeboard and draft of level ice
HFD=(rho*hf+rhos*hfs)/rhow; % level draft
HFF=(hf+hfs)-HFD; % level freeboard

% check that answers add up for floe ice
if any((HFD+HFF)~=(hf+hfs))
 error('HFD+HFF is not equal fo hf+hfs')
elseif HFF<hfs
 error('snow is submerged on parent ice!')
end

% calculate freeboard and draft of deformed ice 
HRD=(rho*hr+rhos*hrs)/rhow; % ridged draft
HRF=(hr+hrs)-HRD; % ridged freeboard

% check that answers add up for ridged ice
if any((HRD+HRF)~=(hr+hrs))
 error('HRD+HRF is not equal fo hr+hrs')
end

% calculate depth of keel relative to sea level
HK=(2*HRD./(1-phi))-HFD;

% calculate height of ridge
HS=HFF+2*sqrt(((HRD./(1-phi))-HFD).*(((HRF./(1-phi))-HFF)));

% convert alpha to radians
alpha=alpha*pi/180;

% calculate horizontal extent of keel structure 
if alpha==0
 LK=0;
else
 LK=2*(HK-HFD).*cot(alpha);
end

% calculate horizontal extent of ridge structure 
if alpha==0
 LS=0;
else
 LS=2*(HS-HFF).*cot(alpha);
end

% initial length of sea ice in ridging
L0=LK.*hr./hf;

% strain 
if alpha==0
 EPSILON=0;
else
 EPSILON=(L0-LK)./L0;
end

if debug; disp(['...Leaving ',mfilename]); end


