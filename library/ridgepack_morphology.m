function [HFD,HFF,HDD,HDF,HK,HS,LK,LS,L0,EPSILON]=...
                    ridgepack_morphology(hf,hfs,hd,hds,phi,alpha)

% function [HFD,HFF,HDD,HDF,HK,HS,LK,LS,L0,EPSILON]=...
%                    ridgepack_morphology(hf,hfs,hd,hds,phi,alpha)
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
% hd    - ridge mean ice thickness (m)
% hds   - ridge mean snow thickness (m)
% phi   - macroporosity of ridge (dimensionless)
% alpha - angle of repose of a ridge (degrees)
%
%
% OUTPUT:
%
% HFD     - floe draft (m)
% HFF     - floe freeboard (m)
% HDD     - ridge mean draft (m)
% HDF     - ridge mean freeboard (m)
% HK      - keel depth (m)
% HS      - sail height (m)
% LK      - cross-sectional keel width (m)
% LS      - cross-sectional sail width (m)
% L0      - original cross-sectional width of undeformed parent ice (m)
% EPSILON - compressional strain or ridge
% 
% Ridgepack Version 1.0.
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

% check there are sufficient inputs
if nargin~=6
 error('incorrect number of inputs')
end

% error checking
if any(hd<hf)
 error('ERROR: hd<hf')
elseif any(hds<0)
 error('ERROR: hds<0')
elseif any(hf<0)
 error('ERROR: hf<0')
elseif any(hd<0)
 error('ERROR: hd<0')
end

% get constants
[rhoi,rhos,rhow,delrho,g]=ridgepack_constants;

% calculate freeboard and draft of level ice
HFD=(rhoi*hf+rhos*hfs)/rhow; % level draft
HFF=(hf+hfs)-HFD; % level freeboard

% check that answers add up for floe ice
if any((HFD+HFF)~=(hf+hfs))
 error('HFD+HFF is not equal fo hf+hfs')
elseif HFF<hfs
 error('snow is submerged on parent ice!')
end

% calculate freeboard and draft of deformed ice 
HDD=(rhoi*hd+rhos*hds)/rhow; % ridged draft
HDF=(hd+hds)-HDD; % ridged freeboard

% check that answers add up for ridged ice
if any((HDD+HDF)~=(hd+hds))
 error('HDD+HDF is not equal fo hd+hds')
end

% calculate depth of keel relative to sea level
HK=(2*HDD./(1-phi))-HFD;

% calculate height of ridge
HS=HFF+2*sqrt(((HDD./(1-phi))-HFD).*(((HDF./(1-phi))-HFF)));

% convert alpha to radians
alpha=alpha*pi/180;

% calculate horizontal extent of keel structure 
LK=2*(HK-HFD).*cot(alpha);

% calculate horizontal extent of ridge structure 
LS=2*(HS-HFF).*cot(alpha);

% initial length of sea ice in ridging
L0=LK.*hd./hf;

% strain 
EPSILON=(L0-LK)./L0;

