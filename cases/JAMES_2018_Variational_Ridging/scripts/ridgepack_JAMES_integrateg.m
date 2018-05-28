% PLACEHOLDER

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
 ghphi(idx,1)=hinitial/hincr(idx);
else
 error('gshape is set incorrectly')
end

% calculate thickness distribution
for i=1:100
 [ghphi]=ridgepack_redistribution(ghphi,hgrid,phisplit);
 disp(num2str(i))
 if i==1
  clf
  semilogy(hgrid,sum(ghphi,2))
  hold on
  drawnow
 end
end

semilogy(hgrid,sum(ghphi,2))
hold on
drawnow

legend({'initial distribution','final distribution'})

% plot up thickness distribution along thickness axis

% determine directory for read/write
dir=fileparts(which(mfilename));
outdir=[dir(1:strfind(dir,'scripts')-1),'output'];
[status,msg]=mkdir(outdir);
cd(outdir);

% determine filename
x=strfind(mfilename,'_');
thisfilename=mfilename;
graphicsout=thisfilename(x(end)+1:end);

% output
disp(['Writing graphics output ',graphicsout,' to:',char(13),' ',pwd])

% print figure
ridgepack_fprint('epsc',graphicsout,1,2)
ridgepack_fprint('png',graphicsout,1,2)

