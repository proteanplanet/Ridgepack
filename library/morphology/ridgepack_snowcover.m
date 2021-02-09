function [hfs,hrs,deltas]=ridgepack_snowcover(hf,hfs)

% ridgepack_snowcover - Calculate updated snow thickness 
%
% function [hfs]=ridgepack_snowcover(hf,hfs)
%
%  

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% error checking
if nargin~=2
 error('incorrect number of inputs')
end

% get constants
hc=ridgepack_astroconstants;
rho=hc.rhoi.const;  % density of ice (kg/m^3)
rhos=hc.rhos.const; % density of snow (kg/m^3)
rhow=hc.rhow.const; % density of seawater (kg/m^3)
delrho=hc.rhow.const-hc.rhoi.const; % difference of water and ice densities
ghat=hc.ghat.const; % acceleration due to gravity

% calculate freeboard and draft of level ice
% (Equations A.3 and A.4 in Roberts et al. 2019)
% N.B. This does not make white ice (snow ice) as would be 
% required in an Earth System Model.
HFD=(rho.*hf+rhos.*hfs)./rhow; % level draft
if HFD>hf
 hfs=delrho.*hf./rhos;
 HFD=hf;
end
HFF=(hf+hfs)-HFD; % level freeboard

% calculate change in snow cover
deltas=hfs-hfs;

% set snow thickness on ridge same as on level ice 
% (equation A.11, Roberts et al. 2019)
hrs=hfs;

if debug; disp(['...Leaving ',mfilename]); end

