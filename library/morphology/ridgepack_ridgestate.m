function [EPSILON,PHI,ALPHAHAT,HR,HRS,HK,HS,LK,LS]=...
                       ridgepack_ridgestate(hf,hfs,epsilon)

% ridgepack_ridgestate - Calculate ridge morphology along trajectory
%
% function [EPSILON,PHI,ALPHAHAT,HR,HRS,HK,HS,LK,LS]=...
%                       ridgepack_ridgestate(hf,hfs,epsilon)
%
% This function calculate ridge metrics along the state-space trajectory,
% given an initial thickness and snow cover of a floe.  Optionally, a 
% selected strain-rate vector can be provided to calculate ridge statistics
% at those locations along the state-space trajectory.
%
% INPUT:
%
% hf      - parent ice sheet thickness (m)
% hfs     - thickness of snow on parent ice (m)
% epsilon - ridge strain (-1,0) as a single number or vector of values for 
%           which ridge statistics are required (dimensionless; optional)
%            
%
% OUTPUT:
%
% EPSILON     - ridge strain along zeta-hat trajectory (dimensionless)
% PHI         - ridge porosity along zeta-hat trajectory (dimensionless)
% ALPHAHAT    - ridge alpha-hat along zeta-hat trajectory (degrees)
% HR          - mean thickness of ice in the ridge (m)
% HRS         - mean thickness of snow on the ridge (m)
% HK          - draft of a keel (m)
% HS          - height of a sail (m)
% LK          - cross-sectional width of a keel (m)
% LS          - cross-sectional width of a sail (m)
%
% Ridgepack Version 1.0.
% Andrew Roberts, Naval Postgraduate School, May 2018 (afrobert@nps.edu)
% Reviewed by Samy Kamal, Naval Postgraduate School, May 2018

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% obtain trajectory for given floe thickness and snow cover
if nargin<2
 error('Insufficient arguments')
else
 [epsilonout,phiout,alphaout]=ridgepack_trajectory(hf,hfs);
end

% Get ice thickness
HR=hf./(1+epsilon);

% Set snow thickness to be identical to that of level ice
HRS=hfs*ones(size(epsilon));

% If calculating for desired strain rates, find nearest match on 
% the state-space trajectory.
if nargin==3
 
 for i=1:length(epsilon)

  jdx=find(abs(epsilonout-epsilon(i))==min(abs(epsilonout-epsilon(i))));

  EPSILON(i)=epsilon(i);

  PHI(i)=mean(phiout(jdx));
 
  ALPHAHAT(i)=mean(alphaout(jdx));

  [HFD,HFF,HRD,HRF,HK(i),HS(i),LK(i),LS(i)]=...
              ridgepack_morphology(hf,hfs,HR(i),HRS(i),PHI(i),ALPHAHAT(i));

 end

% Otherwise calculate ridge state all the way along the trajectory 
else

 EPSILON=epsilonout;

 PHI=phiout;
 
 ALPHAHAT=alphaout;

 for i=1:length(epsilonout)

  [HFD,HFF,HRD,HRF,HK(i),HS(i),LK(i),LS(i)]=...
              ridgepack_morphology(hf,hfs,HR,HRS,PHI(i),ALPHAHAT(i));

 end

end

if debug; disp(['...Leaving ',mfilename]); end
 

