function [GHPHI]=ridgepack_redistribution(hgrid,phigrid,...
                     EPSILON,PHI,VR,HK,HS,LK,LS,ghinit,ghphi,epsilondot,dt)

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
if nargin~=13
 error('incorrect number of inputs')
end

% initialize ar and GHPHI, which each have size MxN
% where M is the length of hgrid, and N is the length of phigrid
ar=zeros(size(ghphi));
GRHPHI=zeros(size(ghphi));
size(GRHPHI)

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


% Integral transform:

% filter results
hidx=find(ghinit(:)>0);
if isempty(hidx)
 error('Nothing positive found in the thickness distribution')
end

% prepare step function building block in continuous (h,phi) space
for i=hidx'
 for j=1:length(phigrid)
  if VR(i,j)>0
   GRHPHI(:,j)=ridgepack_grhphi(hgrid(i),HFs(i),EPSILON(j),PHI(i,j)); 
  end
 end
end

% calculate normalized ghphi at first timestep
ghphinormal=GRHPHI./sum(GRHPHI(:));

% integrate out
for i=hidx'

 for j=1:length(phigrid)

  if VR(i,j)>0

   if ghphi(i,j)>0

    weight(j)=sum(ghphinormal(:,j));

    numerator=ghphi(:,j).*energyratio(:,j:end);
    numerator(VR(:,j:end)==0)=0;

    probability=weight(j)*numerator./sum(numerator(:));

    for k=j:length(phigrid)-1

     if VR(i,k)>0

      % now determine total area based on strain and probability of occurrence
      ar(:,k)=probability(i,k-j+1)*(1+EPSILON(k))*GRHPHI(:,k) + ar(:,k);

     end
 
    end

   end
  
  end
 
 end

end

GHPHI=ar;

if debug; disp(['...Leaving ',mfilename]); end


