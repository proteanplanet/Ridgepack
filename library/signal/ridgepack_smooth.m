function [nc]=ridgepack_smooth(nc,x,y,var,radius,mask)

% ridgepack_smooth - Field smoother using area averages
%
% function ridgepack_smooth(nc,x,y,var,radius,mask)
%
% This function smooths 2D fields using 'disk' averaging
% around a grid point, with the 'radius' specifying the 
% radius, in grid points, of the disk around each grid
% point that are averaged to form the value at each point
% of the new field. Points near the boundary use only as 
% many grid points as are available to provide an average,
% with the new field having the same size as the unsmoothed
% field. To avoid showing data with boundary effects, these
% data are set as NaNs. NaN's are handled by first 
% interpolating to a field with NaN's removed, soothing, 
% and overlaying the original NaN map.
%
% INPUT:
%
% nc     - nc structure
% x      - x dimension of variable to be smoothed
% y      - y dimension of variable to be smoothed
% var    - variable to be smoothed in nc with two dimensions
% radius - radius, in grid points, used to average each 
%          smoothed point as the output
% mask   - optional logical variable determining if the data
%          should be masked using a 'mask' in the nc structure.
%  
%
% OUTPUT:
%
% nc     - nc structure containing a new field called
%          nc.var_smoothed which is a smoothed version of the 
%          2D field.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% Naval Postgraduate School, Department of Oceanography
%

nc=ncshuffle(nc,{'time'});
Z=nc.(var).data;

if any(strcmpi(nc.(var).dimension,'time')) & ndims(Z)>3
 error([var,' has too many dimensions']);
elseif ~any(strcmpi(nc.(var).dimension,'time')) & ndims(Z)>2 
 error([var,' has too many dimensions']);
elseif any(strcmpi(nc.(var).dimension,'time')) & ndims(Z)<3
 error([var,' has too few dimensions']);
elseif ~any(strcmpi(nc.(var).dimension,'time')) & ndims(Z)<2
 error([var,' has too few dimensions']);
elseif any(strcmpi(nc.(var).dimension,'time'))
 ll=size(Z,3);
else
 ll=1;
end

if nargin<6; mask=0; end

for i=1:ll

 disp(['Smoothing record: ',num2str(i)]);

 if any(strcmpi(nc.(var).dimension,'time'))
  ZZ=squeeze(Z(:,:,i));
 else
  ZZ=Z;
 end

 if min(size(nc.(x).data))>1 & min(size(nc.(y).data))>1 & i==1
  xd=nc.(x).data(~isnan(ZZ));
  yd=nc.(y).data(~isnan(ZZ));
 elseif i==1
  [X,Y] = meshgrid(nc.(x).data,nc.(y).data);
  xd=X(~isnan(ZZ));
  yd=Y(~isnan(ZZ));
 end

 % remove NaNs
 zd=ZZ(~isnan(ZZ));

 % regrid data without NaNs
 zi=griddata(xd,yd,zd,X,Y,'nearest');

 % run disk filter over data
 zi=filter2(fspecial('disk',radius),zi,'same');

 % overlay the original NaN mask
 zi(isnan(ZZ))=NaN;

 % mask the variable if required
 if mask & isfield(nc,'mask')
  zi(~nc.mask.data)=NaN;
 end

 % fill in the boundaries
 zi(1:radius,:)=NaN;
 zi(end-radius+1:end,:)=NaN;
 zi(:,1:radius)=NaN;
 zi(:,end-radius+1:end)=NaN;

 if any(strcmpi(nc.(var).dimension,'time'))
  Z(:,:,i)=zi(:,:);
 else
  Z=zi;
 end

end

%dimensions=ncsetdiff(nc.(var).dimension,ncsetdiff(nc.(var).dimension,{x y}));

if ~isfield(nc.(var),'units'); nc.(var).units=''; end

nc=ncadd(nc,[var,'_smoothed'],Z,['Smoothed ',nc.(var).long_name],...
         nc.(var).dimension,nc.(var).units,nc.(var).type);


