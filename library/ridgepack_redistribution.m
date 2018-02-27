function [GHPHI]=ridgepack_redistribution(ghphi,hgrid,phigrid,epsilondot,dt)

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
% ghphi      - bivariate sea ice thickness distribution at time ti
% epsilondot - strain rate in an area A (/second)
% deltat     - timestep deltat = tf - ti (seconds)
%
%
% OUTPUT:
%
% GHPHI - bivariate sea ice thickness distribution at time tf
%
% Ridgepack Version 1.0.
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% check there are sufficient inputs
if nargin~=3
 error('incorrect number of inputs')
end

% initialize ar and GHPHI
ar=zeros(size(ghphi));
GHPHI=zeros(size(ghphi));

% constants
[rhoi,rhos,rhow,delrho,g]=ridgepack_constants;

% Calculate the zeta-hat plane. Please note that this function is dependent 
% on snow cover, but for the purpose of the paper that Ridgepack Version 1.0 
% supports, this dependency has been removed. In this instance, there is a 
% one-to-one mapping of EPSILON to PHI, so it can be assumed that the first
% PHI dimension corresponds to EPSILON.  This assumption cannot be made 
% if snow is included in the calculation, and a further mapping step is required
% that is not included here.
[HF,EPSILON,PHI,ALPHAHAT,VR,HK,HS,LK,LS]=ridgepack_zetahatplane;

% work per ridge shape M x N, where gin(1,:) is the concentration 
% of ice with zero porosity, and is therefore unridged
denominator=LK(:,:).*VR(:,:);
energyratio=sum(denominator(:))./denominator;

% probability of a ridge forming with strain and HF
numerator=ghphi(:,1).*energyratio;
probability=numerator./sum(numerator(:));

% debug information
if debug
 disp(['size of ghphi is ',num2str(size(ghphi))])
 disp(['size of GHPHI is ',num2str(size(GHPHI))])
 disp(['size of ar is ',num2str(size(ar))])
 disp(['size of LK is ',num2str(size(LK))])
 disp(['size of LS is ',num2str(size(LS))])
 disp(['size of VR is ',num2str(size(VR))])
 disp(['size of HF is ',num2str(size(HF))])
 disp(['size of EPSILON is ',num2str(size(EPSILON))])
 disp(['size of energyratio is ',num2str(size(energyratio))])
 disp(['size of probability is ',num2str(size(probability))])
end

% only do calculations where there is ice
hidx=find(ghphi(:,1)>0)

% Calculate ar dependent on strain for the area change from non-porous ice
% where i represents both the EPSILON and PHI indicy is no snow is included.
for i=hidx
 for j=1:length(EPSILON)

  % calculate position on hgrid where steps occur
  h1idx=find(hgrid>=HF(i) & hgrid<=HF(i)+LK(i,j)-LS(i,j));
  h2idx=find(hgrid>HF(i)+LK(i,j)-LS(i,j) & hgrid<=HF(i)+LK(i,j)+LS(i,j));

  % calculate ar function for an indidividual parent sheet thickness
  ar(h1idx,j) = probability(i,j)*(1+EPSILON(j))/(LK(i,j)) + ar(h1idx,j);
  ar(h2idx,j) = probability(i,j)*(1+EPSILON(j))/(2*LK(i,j)) + ar(h2idx,j);

 end
end

% check ar function 
clf

size(ar)
size(hgrid)

plot(hgrid',ar(:,20))

return


% finish determining by determining gout       
GHPHI(:,2:end)=ar(:,2:end)+ghphi(:,2:end);
GHPHI(:,1)=ar(:,1);

if debug; disp(['...Leaving ',mfilename]); end


