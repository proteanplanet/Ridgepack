%function ridgepack


% Written by Andrew Roberts, March 2018

% inititalize thickness distribution  
[hgrid,epsilongrid,phigrid,epsilonsplit,phisplit,ghphi]=ridgepack_gridinit;


disp(['size of ghphi is ',num2str(size(ghphi))])

% set initial thickness distribution of sea ice field (meters)
gshape='delta';
if strcmp(gshape,'delta')
 hinitial = 2.0;
 idx=find(min(abs(hgrid(:)-hinitial))==abs(hgrid(:)-hinitial));
 ghphi(idx,1)=1;
else
 error('gshape is set incorrectly')
end

% calculate thickness distribution
[GHPHI]=ridgepack_redistribution(ghphi,hgrid,epsilonsplit);






