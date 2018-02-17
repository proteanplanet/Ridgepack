function [hf,epsilon,phi,alpha,vr,HK,HS,LK,LS]=ridgepack_zetahatplane

% This function calculates strain, porosity and angle of repose of ridges
% on the the zeta hat trajectory plane of potential energy density. 
%
% Output:
% -------
% hf 
% epsilon
% phi
% alpha
% vr
% HK
% HS
% LK
% LS

% Written by Andrew Roberts, March 2018

% retrieve constants
[rhoi,rhos,rhow,delrho,g,epres,maxthick]=ridgepack_constants;

% set hf grid from 1cm to maxthick thick ice
incr=(log10(10)-log10(0.01))/1000;
hf=10.^[log10(0.01):incr:log10(maxthick)];

% assume there is no snow in this case
hfs=zeros(size(hf));

% create strain and split strain coordinates
stii=[-epres:-epres:-0.99];
phii=[epres:epres:0.50];

% generate 2-D arrays from these vectors
[strain,phi]=meshgrid(stii,phii);

% step through floe thickness
for i=1:length(hf) 

 % calculate potential energy field
 for k=1:length(stii);
  for j=1:length(phii);
   [PE(k,j)]=paper_ridge_energetics(hf(i),hfs(i),stii(k),phii(j));
  end
 end

 % calculate optimal path
 [epsilon,phi(i,:),vr(i,:)]=paper_ridging_path(strain,phi,PE);

 for k=1:length(epsilon)

  % calculate thickness of deformed ice mass from strain
  hdi(i,k)=hf(i)/(1+epsilon(k));

  % set snow thickness on ridge same as on level ice (as indicated in the paper)
  hds(i)=hfs(i);

  % get angle ot repose
  alpha(i,k)=paper_ridge_alpha(epsilon(k),phi(i,k),hf(i),hdi(i,k));

  % get morphological shape
  [hfd,hff,hdd,hdf,HK(i,k),HS(i,k),LK(i,k),LS(i,k)]=...
      ridge_morphology(hf(i),hfs(i),hdi(i,k),hds(i),phi(i,k),alpha(i,k));

 end

end

