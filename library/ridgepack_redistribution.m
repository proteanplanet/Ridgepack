function [g]=ridgepack_redistribution(g,epsilondot,dt)

% RIDGEPACK_REDISTRIBUTION - Redistribution function Psi
%
% function [g]=ridgepack_redistribution(g,epsilondot,dt)
%
% Thid function calculate Psi in the sea ice mass conservation 
% equation for a bivariate thickness distribution that is a function
% of both thickness and macroporosity.
%
% INPUT:
%
% gin        - bivariate sea ice thickness distribution at time ti
% epsilondot - strain rate in an area A (/second)
% deltat     - timestep deltat = tf - ti (seconds)
%
%
% OUTPUT:
%
% gout - bivariate sea ice thickness distribution at time tf
%
% Ridgepack Version 1.0.
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

% constants
[rhoi,rhos,rhow,delrho,g]=ridgepack_constants;

% Calculate the zeta-hat plane. Please note that this function is dependent 
% on snow cover, but for the purpose of the paper that Ridgepack Version 1.0 
% supports, this dependency has been removed.
[HF,EPSILON,PHI,ALPHAHAT,VR,HK,HS,LK,LS]=ridgepack_zetahatplane;

% work per ridge shape M x N, where gin(1,:) is the concentration 
% of ice with zero porosity, and is therefore unridged
energyratio=sum(gin(1,:).*LK(:,:).*VR(:,:))./(gin(1,:).*LK(:,:).*VR(:,:));

% probability of a ridge forming with strain and HF
probability=energyratio./sum(energyratio(:));

% calculate ar dependent on strain for the area change from non-porous ice
for j=1:length(HF)
 for i=1:length(EPSILON)
  ar(i,h<HF(j))=0;
  ar(i,h>=HF(j) & h<=HF(j)+LK(i,j)-LS(i,j))=...
                         probability(i)*2*(1+EPSILON(i))/(2*LK(i,j));
  ar(i,h>HF(j)+LK(i,j)-LS(i,j) & h<=HF(j)+LK(i,j)+LS(i,j))=...
                         probability(i)*1*(1+EPSILON(i))/(2*LK(i,j));
  ar(i,h>HF(j)+LK(i,j)+LS(i,j))=0;
 end
end










