%function ridgepack


% Written by Andrew Roberts, March 2018

% inititalize thickness distribution  
[hincr,eincr,hgrid,epsilongrid,phigrid,epsilonsplit,phisplit,ghphi]=ridgepack_gridinit;

disp(['size of ghphi is ',num2str(size(ghphi))])

% set initial thickness distribution of sea ice field (meters)
gshape='delta';
% Case of the initial thickness distribution as a delta function on discrete grid
if strcmp(gshape,'delta')
 hinitial = 2.0;
 idx=find(min(abs(hgrid(:)-hinitial))==abs(hgrid(:)-hinitial));
 ghphi(idx,1)=hinitial/hincr;
else
 error('gshape is set incorrectly')
end

% calculate thickness distribution
[GHPHI]=ridgepack_redistribution(ghphi,hgrid,phisplit);


% plot up thickness distribution along thickness axis
clf
semilogy(hgrid,sum(GHPHI,2))






