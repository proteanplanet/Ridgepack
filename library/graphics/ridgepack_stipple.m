function ridgepack_stipple(nc,var,mat,cutoff,colors)

% ridgepack_stipple - Stipple above a cutoff value using the current axes 
%
% function ridgepack_stipple(nc,var,mat,val)
%
% This function stipples data above a cutoff value on the current axes, whether 
% the axes is a map or is in Cartesian coordinates. The spacing of the stipples
% is supplied via mat, so that one stipple represents the mean over mat x mat
% grid cells, and is centered above them. This is particular useful for indicating
% areas of statistical significance as calculate from ncttest2.
%
% INPUT:
%
% nc     - netcdf structure containing the variable var used to stipple 
%
% var    - Variable being stippled
%
% mat    - Square mat of grid points over which on stipple is to be supplied
%          to represent the cutoff above mean value of the mat x mat grid cells.
%          (optional, default is 8)
%
% cutoff - Cuttoff value above which stipples are applied.
%          (optional, default is 0.5).
% 
% colors  - color of stipples in RGB vector (e.g. white is [1 1 1])
%
% Note, if using output from ncttest2, one could use 'ridgepack_stipple(nc,'h',3,0.5)'
% which would stipple a square mat of 9 grid points that are statistically significant.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

if nargin<2
 error('incorrect number of inputs')
elseif ~isstruct(nc)
 error('nc is not a structure')
elseif ~ischar(var)
 error('var is not a character variable')
elseif nargin<3 
 mat=5;
 cutoff=0.5;
elseif ~isnumeric(mat)
 error('mat must be an integer value')
elseif nargin<4 
 cutoff=0.5;
elseif ~isnumeric(cutoff)
 error('cutoff must be a number')
end

if nargin<5
 colors=[0 0 0];
end

ha=get(gcf,'CurrentAxes');

% get x and y data for stippling
if isempty(ha)

 error('missing current axes')

elseif ismap(ha)

 disp('Stippling on current map axes')
 
 if ~isfield(nc,'x') | ~isfield(nc,'y')
   error('missing x and/or y data to stipple')
 elseif ~isfield(nc,'longitude') | ~isfield(nc,'latitude')
   error('missing latitude and/or longitude data to stipple')
 end

 % get x and y location on the map
 [x,y]=meshgrid(nc.x.data,nc.y.data);
 sz=size(x);
 mstruct=gcm;
 [x,y] = mfwdtran(mstruct,nc.latitude.data,nc.longitude.data,ha,'none');
 xz=size(x);
 if sz~=xz
  error('Change in size of x using mfwdtran. Try changing surface to none in the code')
 end

else

 if ~isfield(nc,'x') | ~isfield(nc,'y')
   error('missing x and/or y data to stipple')
 end

 [x,y]=meshgrid(nc.x.data,nc.y.data);

end

count=0;
xmean=[];
ymean=[];

for i=1:mat:size(x,1)-mat
for j=1:mat:size(x,2)-mat

 h=nc.h.data(i:i+mat-1,j:j+mat-1);

 if sum(isnan(h(:)))<(mat*mat)/2

  count=count+1;

  if nanmean(h(:))>cutoff
   xmean(count)=nanmean(reshape(x(i:i+mat-1,j:j+mat-1),[mat*mat 1]));
   ymean(count)=nanmean(reshape(y(i:i+mat-1,j:j+mat-1),[mat*mat 1]));
  end

 end

end
end

if ~isempty(xmean)
 plot(xmean,ymean,'.k','MarkerSize',1.0,'Color',colors)
end


 













 




