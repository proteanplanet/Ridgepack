function [EPSILON,PHI,ALPHAHAT,VR,HK,HS,LK,LS,epsilonsplit,phisplit,d1,d2]=...
           ridgpack_trajectory(epsilon,phi,alphahat,vr,hk,hs,lk,ls)

% function [EPSILON,PHI,ALPHAHAT,VR,HK,HS,LK,LS,epsilonsplit,phisplit,d1,d2]=...
%            ridgpack_trajectory(epsilon,phi,alphahat,vr,hk,hs,lk,ls)
%
% This function is part of Ridgepack Version 1.0.
% It calculates the zeta-hat trajectory on VR leafs of the alpha-hat
% plane of an individual ridge with potential energy density vr.
%
% INPUT:
%
% epsilon  - ridge strain coordinate on alphahat plane (dimensionless)
% phi      - ridge macroporosity coordinate on alphahat plane (dimensionless)
% alphahat - ridge alpha-hat angle of repose (degrees)
% vr       - potential energy density on alphahat plane (J m^-2)
% hk       - draft of a keel (m)
% hs       - height of a sail (m)
% lk       - cross-sectional width of a keel (m)
% ls       - cross-sectional width of a sail (m)
%
%
% OUTPUT:
%
% EPSILON      - ridge strain along zeta-hat trajectory (dimensionless)
% PHI          - ridge porosity along zeta-hat trajectory (dimensionless)
% ALPHAHAT     - ridge alpha-hat along zeta-hat trajectory (degrees)
% VR           - potential energy density along zeta-hat trajectory (J m^-2)
% HK           - draft of a keel (m)
% HS           - height of a sail (m)
% LK           - cross-sectional width of a keel (m)
% LS           - cross-sectional width of a sail (m)
% epsilonsplit - ridge strain on split epsilon grid (dimensionless)
% phisplit     - porosity on split epsilon grid (dimensionless)
% d1           - dilation field [epsilon-component] (J m^-2)
% d2           - dilation field [phi-component] (J m^-2)
%
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

% error check
if nargin~=8
 error('incorrect number of inputs')
elseif any(size(epsilon)~=size(phi))
 error('size of epsilon not the same as phi')
elseif any(size(epsilon)~=size(alphahat))
 error('size of alphahat not the same as epsilon')
elseif ndims(epsilon)>2
 error('epsilon and phi have a dimension greater than 2')
end

% create epsilon and split epsilon coordinates
try
 stspacing=abs(epsilon(2,1)-epsilon(1,1));
catch
 error('epsilon array is not the correct shape or size')
end
stii=epsilon(:,1);
dstii=(stspacing./2)+stii(1:end-1);

% create phi and split phi coordinates
try
 phspacing=abs(phi(1,2)-phi(1,1));
catch
 error('epsilon array is not the correct shape or size')
end
size(phi)
phii=phi(1,:);
dphii=(phspacing./2)+phii(1:end-1);

% calculate vector gradient
dVdep=(vr(2:end,:)-vr(1:end-1,:))./(epsilon(2:end,:)-epsilon(1:end-1,:));
dVdph=(vr(:,2:end)-vr(:,1:end-1))./(phi(:,2:end)-phi(:,1:end-1));
d1=(dVdep(:,1:end-1)+dVdep(:,2:end))/2;
d2=(dVdph(1:end-1,:)+dVdph(2:end,:))/2;

% initialize arrays
EPSILON=stii(1:1:end-1);
ALPHAHAT=zeros(size(EPSILON));
VR=zeros(size(EPSILON));
HK=zeros(size(EPSILON));
HS=zeros(size(EPSILON));
LK=zeros(size(EPSILON));
LS=zeros(size(EPSILON));

% calculate streamline from initial conditions
jdx=1;
PHI(1)=phii(jdx); % initial condition
ALPHAHAT(1)=alphahat(1,jdx);
size(ALPHAHAT)
size(alphahat)
VR(1)=vr(1,jdx);
HK(1)=hk(1,jdx);
HS(1)=hs(1,jdx);
LK(1)=lk(1,jdx);
LS(1)=ls(1,jdx);
for i=2:length(EPSILON)
 PHI(i)=PHI(i-1)-(d2(i,jdx)./d1(i,jdx)).*stspacing;
 jdx=find(min(abs(PHI(i)-phii(:)))==abs(PHI(i)-phii(:)));
 if jdx>size(d2,2) | jdx>size(d1,2)
  error('Max phi exceeded')
 end
 i
 jdx
 ALPHAHAT(i)=alphahat(i,jdx);
 VR(i)=vr(i,jdx);
 HK(i)=hk(i,jdx);
 HS(i)=hs(i,jdx);
 LK(i)=lk(i,jdx);
 LS(i)=ls(i,jdx);
end

% prepare vector field for output
[epsilonsplit,phisplit]=meshgrid(dstii,dphii);
d1=d1';
d2=d2';

