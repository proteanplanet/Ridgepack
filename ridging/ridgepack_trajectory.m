function [EPSILON,PHI,ALPHAHAT,VR,HK,HS,LK,LS,vr,epsilonsplit,phisplit,d1,d2]=...
            ridgepack_trajectory(hf,hfs,epsilon,phi)

% function [EPSILON,PHI,ALPHAHAT,VR,HK,HS,LK,LS,epsilonsplit,phisplit,d1,d2]=...
%            ridgepack_trajectory(hf,hfs,epsilon,phi)
%
% This calculates the zeta-hat trajectory on VR leafs of the alpha-hat
% plane of an individual ridge with potential energy density.
%
% INPUT:
%
% hf       - parent ice thickness (m)
% hfs      - thickness of snow on parent ice (m)
% epsilon  - ridge strain coordinate on alphahat plane (dimensionless)
% phi      - ridge macroporosity coordinate on alphahat plane (dimensionless)
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
% vr           - potential energy density on (epsilon,phi) grid (J m^-2)
% epsilonsplit - ridge strain on split epsilon grid (dimensionless)
% phisplit     - porosity on split epsilon grid (dimensionless)
% d1           - dilation field [epsilon-component] (J m^-2)
% d2           - dilation field [phi-component] (J m^-2)
%
% Ridgepack Version 1.0.
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

% error check
if nargin~=4
 error('incorrect number of inputs')
elseif any(size(epsilon)~=size(phi))
 error('size of epsilon not the same as phi')
elseif ndims(epsilon)>2
 error('epsilon and phi have a dimension greater than 2')
elseif length(hf)>1 | length(hfs)>1
 error('length of hf or hfs is greater than 1')
end

% create epsilon and split epsilon coordinates (could use ridgepack_constants here)
try
 stspacing=abs(epsilon(1,2)-epsilon(1,1));
catch
 error('epsilon array is not the correct shape or size')
end
stii=epsilon(1,:);
dstii=(stspacing./2)+stii(1:end-1);

% create phi and split phi coordinates (could use ridgepack_constants here)
try
 phspacing=abs(phi(1,2)-phi(1,1));
catch
 error('epsilon array is not the correct shape or size')
end
phii=phi(:,1);
dphii=(phspacing./2)+phii(1:end-1);

% get potential energy density for given ice thickness *
[vr]=ridgepack_energetics(hf,hfs,epsilon,phi);

% calculate dilation field
dVdep=(vr(:,2:end)-vr(:,1:end-1))./(epsilon(:,2:end)-epsilon(:,1:end-1));
dVdph=(vr(2:end,:)-vr(1:end-1,:))./(phi(2:end,:)-phi(1:end-1,:));
d1=(dVdep(1:end-1,:)+dVdep(2:end,:))/2;
d2=(dVdph(:,1:end-1)+dVdph(:,2:end))/2;

% initialize arrays
EPSILON=stii(1:1:end-1);
PHI=zeros(size(EPSILON));

% integrate streamline from initial conditions
jdx=1;
PHI(1)=phii(jdx); % initial condition
for i=2:length(EPSILON)
 PHI(i)=PHI(i-1)-(d2(jdx,i)./d1(jdx,i)).*stspacing;
 jdx=find(min(abs(PHI(i)-phii(:)))==abs(PHI(i)-phii(:)));
 if jdx>size(d2,1) | jdx>size(d1,1)
  error('Max phi exceeded')
 end
end

% determine entire state space (note this step can be added to * for efficiency)
[VR,ALPHAHAT,HK,HS,LK,LS]=ridgepack_energetics(hf,hfs,EPSILON,PHI);

% prepare dilation field for output
[epsilonsplit,phisplit]=meshgrid(dstii,dphii);
d1=d1;
d2=d2;

