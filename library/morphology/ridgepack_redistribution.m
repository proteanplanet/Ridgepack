function [GHPHI]=ridgepack_redistribution(hgrid,hincr,phigrid,phincr,...
                     EPSILON,PHI,VR,HK,HS,LK,LS,ghphi)

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
% OUTPUT:
%
% GHPHI - bivariate sea ice thickness distribution at time tf
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
if nargin~=12
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
maskin=(HK+HS<max(hgrid));
maskout=(HK+HS>=max(hgrid));

% work per ridge shape M x N, where gin(1,:) is the concentration 
% of ice with zero porosity, and is therefore unridged
denominator=LK(:,:).*VR(:,:);
energyratio=sum(denominator(:))./denominator;

% probability of a ridge forming with strain and parent ice thickness HF
% with discrete distribution ghphi with minimim porosity.
numerator=ghphi(:,1).*energyratio;
numerator(maskout)=0;
probability=numerator./sum(numerator(:));

% calculate a normalized version of ghphi
[phincrs,hincrs]=meshgrid(phincr,hincr);
scratch=ghphi.*phincrs.*hincrs;
ghphinormal=scratch./sum(scratch(:));

scratch=ghphi.*phincrs.*hincrs;

disp(['Integrated g(h) before redistribution: ',num2str(sum(scratch(:)))])

% only do calculations where there is ice in the initial category
hidx=find(ghphi(:,1)>0);

% Integral transform:
%
% Calculate ar dependent on strain for the area change from non-porous ice
% where i represents the undeformed ice index, and j is the strain index.

weight=sum(ghphinormal(:,1))

for i=hidx'
 for j=1:length(phigrid)

  if maskin(i,j)

   % calculate step function for an indidividual ridge of the given strain
   GRHPHI=ridgepack_grhphi(HF(i),HFs(i),EPSILON(j),PHI(i,j))./phincr(j);
   
   % now determine total area based on strain and probability of occurrence
   ar(:,j)=weight*probability(i,j)*(1+EPSILON(j))*GRHPHI(:) + ar(:,j);
  
  end

 end
end

for i=hidx'
 for j=2:length(EPSILON)-1

  if ghphi(i,j)>0

   % work per ridge shape M x N, where gin(1,:) is the concentration 
   % of ice with zero porosity, and is therefore unridged

   denominator=LK(i,j:end).*VR(i,j:end);
   energyratio=sum(denominator(:))./denominator;

   % probability of a ridge forming with strain and parent ice thickness HF
   % with discrete distribution ghphi with minimim porosity.
   numerator=ghphi(:,j).*energyratio;
   numerator(maskout(:,j:end))=0;
   probability=numerator./sum(numerator(:));

   weight(j)=sum(ghphinormal(:,j));

   %disp(['Sum over probability: ',num2str(sum(probability(:)))])

   for k=j:length(EPSILON)-1

    if maskin(i,k)

     % calculate step function for an indidividual ridge of the given strain
     GRHPHI=ridgepack_grhphi(HF(i),HFs(i),EPSILON(k),PHI(i,k))./phincr(k);
   
     % now determine total area based on strain and probability of occurrence
     ar(:,k)=weight(j)*probability(i,k-j+1)*(1+EPSILON(j)-EPSILON(k))*GRHPHI(:) + ar(:,k);

    end


   end 

  end

 end

end


scratch=ar.*phincrs.*hincrs;
disp(['Integrated g(h) after redistribution: ',num2str(sum(scratch(:)))])

% Advect in ice of the same distribution
scratch=ar.*phincrs.*hincrs;
GHPHI=scratch./sum(scratch(:));
GHPHI=GHPHI./(phincrs.*hincrs);

scratch=GHPHI.*phincrs.*hincrs;
disp(['Integrated g(h) after advection: ',num2str(sum(scratch(:)))])

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

