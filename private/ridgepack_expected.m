function [hfi,strainp,pep,phip,alphap,HK,HS,LK,LS]=paper_ridge_expected_plane(lowres,hfi,hfs)

% constants
[rhoi,rhos,rhow,delrho,g]=ridge_constants;

% calculate over an entire plane
if nargin==1
 % set grid thickness of deforming ice
 incr=(log10(10)-log10(0.01))/1000;
 hfi=10.^[log10(0.01):incr:log10(20)];
end

% assume there is no snow if not supplied
if nargin<3
 % assume now snow in this example
 hfs=zeros(size(hfi));
end

% Quick argument check
if nargin>3
 disp('There''s a problem with the inputs')
end

% create strain and split strain coordinates
if lowres | nargin==0
 spacing=0.001;
else
 spacing=0.0001;
end
stii=[-spacing:-spacing:-0.99];
phii=[spacing:spacing:0.50];

% generate 2-D arrays from these vectors
for k=1:length(stii);
for j=1:length(phii);
 strain(k,j)=stii(k);
 phi(k,j)=phii(j);
end
end

% step through floe thickness
for i=1:length(hfi) 

 % calculate potential energy field
 for k=1:length(stii);
  for j=1:length(phii);
   [PE(k,j)]=paper_ridge_energetics(hfi(i),hfs(i),stii(k),phii(j));
  end
 end

 % calculate optimal path
 [strainp,phip(i,:),pep(i,:)]=paper_ridging_path(strain,phi,PE);

 for k=1:length(strainp)

  % calculate thickness of deformed ice mass from strain
  hdi(i,k)=hfi(i)/(1+strainp(k));

  % set snow thickness on ridge same as on level ice (as indicated in the paper)
  hds(i)=hfs(i);

  % get angle ot repose
  alphap(i,k)=paper_ridge_alpha(strainp(k),phip(i,k),hfi(i),hdi(i,k));

  % get morphological shape
  [hfd,hff,hdd,hdf,HK(i,k),HS(i,k),LK(i,k),LS(i,k)]=...
      ridge_morphology(hfi(i),hfs(i),hdi(i,k),hds(i),phip(i,k),alphap(i,k));

 end

end

