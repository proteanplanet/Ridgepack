function [HF,EPSILON,PHI,ALPHA,VR,HK,HS,LK,LS]=ridgepack_zetahatplane

% function [HF,EPSILON,PHI,ALPHA,VR,HK,HS,LK,LS]=ridgepack_zetahatplane
%
% This function is part of Ridgepack Version 1.0.
% It calculates strain, porosity and the angle of repose of ridges
% on the the zeta hat trajectory plane of potential energy density. 
% The entire ridge state space may be determined from these variables,
% and is output as part of the calculations.
%
% OUTPUT:
%
% HF      - thickness of initial ice (m)
% EPSILON - strain of ridge (dimensionless)
% PHI     - porosity of a ridge (dimensionless)
% ALPHA   - angle of repose of a ridge (degrees)
% VR      - potential energy density of a ridge (J m^-2)
% HK      - draft of a keel (m)
% HS      - height of a sail (m)
% LK      - cross-sectional width of a keel (m)
% LS      - cross-sectional width of a sail (m)
%
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

% retrieve constants
[rhoi,rhos,rhow,delrho,g,epres,maxthick]=ridgepack_constants;

% set HF grid from 1cm to maxthick thick ice
incr=(log10(10)-log10(0.01))/1000;
HF=10.^[log10(0.01):incr:log10(maxthick)];

% assume there is no snow in v1 (although the library allows it)
HFs=zeros(size(HF));

% create strain and porosity grid
stii=[-epres:-epres:-0.99];
phii=[epres:epres:0.50];
[strain,porosity]=meshgrid(stii,phii);

% step through floe thickness
for i=1:length(HF) 

 % calculate potential energy field for a given strain
 vr=paper_ridge_energetics(HF(i),HFs(i),stii,phii);

 % calculate optimal path
 [EPSILON,PHI(i,:),VR(i,:)]=paper_ridging_path(strain,porosity,vr);

 for k=1:length(EPSILON)

  % calculate thickness of deformed ice mass from strain
  hdi(i,k)=HF(i)/(1+EPSILON(k));

  % set snow thickness on ridge same as on level ice (as indicated in the paper)
  hds(i)=HFs(i);

  % get angle ot repose
  ALPHA(i,k)=paper_ridge_ALPHA(EPSILON(k),PHI(i,k),HF(i),hdi(i,k));

  % get morphological shape
  [HFd,HFf,hdd,hdf,HK(i,k),HS(i,k),LK(i,k),LS(i,k)]=...
      ridge_morphology(HF(i),HFs(i),hdi(i,k),hds(i),PHI(i,k),ALPHA(i,k));

 end

end

