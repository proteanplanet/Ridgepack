function [PE,LK,LS,HK,HS,ALPHA]=paper_ridge_energetics(hfi,hfs,strain,phi)

% constants
[rho,rhos,rhow,delrho,g]=ridge_constants;

% calculate thickness of deformed ice mass from strain
hdi=hfi/(1+strain);

% set snow thickness on ridge same as on level ice (as indicated in the paper)
hds=hfs;

ALPHA=paper_ridge_alpha(strain,phi,hfi,hdi);

% get morphological shape
[hfd,hff,hdd,hdf,HK,HS,LK,LS,L0,STRAIN]=ridge_morphology(hfi,hfs,hdi,hds,phi,ALPHA);

% Calculate potential energy
PE = delrho*g*(1-phi)*(0.5*hfd*LK+0.125*(LK.^2)*tan(ALPHA*pi/180)) + ...
        rho*g*(1-phi)*(0.5*hff*LK+0.125*(LS.^2)*tan(ALPHA*pi/180));







