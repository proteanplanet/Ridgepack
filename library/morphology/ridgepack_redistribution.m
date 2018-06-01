function [GHPHI]=ridgepack_redistribution(hgrid,hincr,phigrid,phincr,...
                     EPSILON,PHI,VR,HK,HS,LK,LS,ghphi,epsilondot,dt)

% ridgepack_redistribution - Redistribution function Psi
%
% function [g]=ridgepack_redistribution(g,epsilondot,dt)
%
% This function calculates Psi in the sea ice mass conservation 
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
% Reviewed by Samy Kamal, Naval Postgraduate School, May 2018

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% check there are sufficient inputs
if nargin~=14
 error('incorrect number of inputs')
end

[phigrids,hgrids]=meshgrid(phigrid,hgrid);

% initialize ar and GHPHI
ar=zeros(size(ghphi));
GHPHI=zeros(size(ghphi));

% Rerrange to HF variable for legibility with paper
HF=hgrid;

% set snow cover to zero for version 1
HFs=zeros(size(hgrid));

% Only include thickness within the range of the chosen ridge (this
% can be adjusted for Earth System Model with variable thicknesses grid.)
% This is achieved numerically by zero-ing out the potential energy off-grid
% thicknesses
VR(HK+HS>max(hgrid))=0;

% work per ridge shape M x N, where gin(1,:) is the concentration 
% of ice with zero porosity, and is therefore unridged
denominator=LK(:,:).*VR(:,:);
energyratio=sum(denominator(:))./denominator;

% probability of a ridge forming with strain and parent ice thickness HF
% with discrete distribution ghphi with minimim porosity.
numerator=ghphi(:,1).*energyratio;
numerator(VR==0)=0;
probability=numerator./sum(numerator(:));

size(probability)

sum(probability(:))

% calculate a normalized version of ghphi
[phincrs,hincrs]=meshgrid(phincr,hincr);
scratch=ghphi.*phincrs.*hincrs;
ghphinormal=scratch./sum(scratch(:));

%ghphinormal=ghphi./sum(ghphi(:));

scratch=ghphi.*phincrs.*hincrs;
disp(['Integrated g(h) before: ',num2str(sum(scratch(:)))])

% only do calculations where there is ice
hidx=find(ghphi(:,1)>0);

% Integral transform:

% Calculate ar dependent on strain for the area change from non-porous ice
% where i represents the undeformed ice index, and j is the strain index.
for i=hidx'
 for j=1:length(phigrid)

  if VR(i,j)>0

   % calculate step function for an indidividual ridge of the given strain
   GRHPHI=ridgepack_grhphi(HF(i),HFs(i),EPSILON(j),PHI(i,j))/phincr(j);
   
   % now determine total area based on strain and probability of occurrence
   % spread between the two closest numerical categories
   ar(:,j)=probability(i,j)*(1+EPSILON(j))*GRHPHI(:) + ar(:,j);
  
  end

 end
end


% only do calculations where there is ice
hidx=find(ghphi(:,2)>0);

if ~isempty(hidx)

for i=hidx'
 for j=2:2 %length(EPSILON)

  % work per ridge shape M x N, where gin(1,:) is the concentration 
  % of ice with zero porosity, and is therefore unridged
  denominator=LK(i,j:end).*VR(i,j:end);
  energyratio=sum(denominator(:))./denominator;

  % probability of a ridge forming with strain and parent ice thickness HF
  % with discrete distribution ghphi with minimim porosity.
  numerator=ghphi(:,j).*energyratio;
  numerator(VR(i,j:end)==0)=0;
  probability=numerator./sum(numerator(:));
 
  size(probability)

  disp(['Sum over probability: ',num2str(sum(probability(:)))])


%  for k=j:length(EPSILON)
%  if VR(i,j)>0

%   % calculate step function for an indidividual ridge of the given strain
%   GRHPHI=ridgepack_grhphi(HF(i),HFs(i),EPSILON(j),PHI(i,j))/phincr(j);
   
%   % now determine total area based on strain and probability of occurrence
%   % spread between the two closest numerical categories
%   probability(i,k)
%   EPSILON(j)
%   EPSILON(k)
%   ar(:,k)
%   ar(:,k)=probability(i,k)*(1+EPSILON(j)-EPSILON(k))*GRHPHI(:) + ar(:,k);
%  end
 end 
end

end


%scratch=ar.*phincrs.*hincrs;
%disp(['Integrated g(h) after: ',num2str(sum(scratch(:)))])

%GHPHI=ar;

return

% Now calculate the transform of porous ice, summing over all 
for j=2:length(phigrid)-1

 % vectorize second loop for efficiency
 jdx=find(phigrid>=phigrid(j));

 weight(j)=sum(ghphinormal(:,j));

 %ar(:,j)=ar(:,j)-weight(j)*ghphi(:,j);

% ar(:,j)=ar(:,j)-ghphi(:,j);

% ghphinew(:,jdx)=weight(j)*ghphi(:,jdx)

 ar(:,jdx)=weight(j).*ghphi(:,jdx) + ar(:,jdx);

end

% finish determining by determining gout       
GHPHI=ar;

scratch=ar.*phincrs.*hincrs;
disp(['Integrated g(h) after: ',num2str(sum(scratch(:)))])

% normalize the distribution, which is the equivalent of advecting in ice
% with the same distribution

if debug; 
 disp(['size of ghphi is ',num2str(size(ghphi))])
 disp(['size of ar is ',num2str(size(ar))])
 disp(['size of LK is ',num2str(size(LK))])
 disp(['size of LS is ',num2str(size(LS))])
 disp(['size of VR is ',num2str(size(VR))])
 disp(['size of HF is ',num2str(size(HF))])
 disp(['size of EPSILON is ',num2str(size(EPSILON))])
 disp(['size of energyratio is ',num2str(size(energyratio))])
 disp(['size of probability is ',num2str(size(probability))])
 disp(['...Leaving ',mfilename]); end
end


