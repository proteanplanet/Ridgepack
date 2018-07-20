function [nc]=ridgepack_ttest2(nc1,var1,nc2,var2,conf)

% NCTTEST - Performs Welch's two-sided t-test between two datasets
% 
% function [nc]=ridgepack_ttest2(nc1,var1,nc2,var2,conf)
% 
% This function performs a two-sided t-test for two means with standard deviation
% and effective sample size, N, provided in each netcdf structure as output
% from the ncreduce function. This uses Welch's t-test, with effective 
% degrees of freedom given by df=N1+N2-2.  The difference between the two means
% is also provided in the output.
%
% INPUT:
%
%
% nc1      - nc structure containing the mean field, nc1.var1, the standard 
%            deviation, nc1.var1_std, and the effective sample size nc1.var1_equiv.
%
% var1     - The first variable mean being used in the t-test.
%
% nc2      - nc structure containing the mean field, nc2.var2, the standard 
%            deviation, nc2.var2_std, and the effective sample size nc2.var2_equiv.
%
% var2     - The second variable mean being used in the t-test.
%
% conf     - Confidence interval, provided as one of two options:
%            conf=1 using the 99% interval and conf=2 using the 95% interval 
%            for a two-sided test.  This argument is optional, and the 
%            default is conf=1 (99%) if it is omitted.
%
%
% OUTPUT:
%
%
% nc - netcdf structure containing a variable 'h' giving the results of the 
%      hypothesis test. h=1 indicates a significant difference between the 
%      two means, and h=0 indicates no significance at the desired confidence
%      interval. It also contains the field 'diff', which is the difference
%      of the two means.
%
% Notes:
%
% Output from this function can be used by ncstipple to stipple areas with 
% significant mean differences on a map. For example, one would use:
%
% ncstipple(nc,'h',3,0.5) 
% 
% which would stipple one square 3x3 grid point area with a mean h greater
% 0.5, which corresponds to more than half of the 3x3 square grid point mat
% being statistically significant.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

if nargin<4
 error('incorrect number of inputs')
elseif nargin==5 & (conf~=1 & conf~=2 & conf~=3 & conf~=4) 
 error('Confidence option is either 1, 2, 3 or 4 for 99%, 95%, 5% or 1%')
elseif nargin==4
 conf=1;
end

% extract standard deviation and equivalent degrees of freedom
% from mean dataset.
findx=0;
[variablenames1,numbervariables1]=ridgepack_name(nc1);
for i=1:numbervariables1
 name=char(variablenames1{i});
 if strcmpi(name,[var1,'_std'])
  std1=nc1.(name).data;
  findx=findx+1;
 elseif strcmpi(name,[var1,'_equiv'])
  equiv1=nc1.(name).data;
  findx=findx+1;
 elseif strcmpi(name,[var1])
  mu1=nc1.(name).data;
  findx=findx+1;
 end
end

if findx<3
 error('Unable to find mean, standard deviation or effective sample size in nc1')
end

findx=0;
[variablenames2,numbervariables2]=ridgepack_name(nc2);
for i=1:numbervariables2
 name=char(variablenames2{i});
 if strcmpi(name,[var2,'_std'])
  std2=nc2.(name).data;
  findx=findx+1;
 elseif strcmpi(name,[var2,'_equiv'])
  equiv2=nc2.(name).data;
  findx=findx+1;
 elseif strcmpi(name,[var2])
  mu2=nc2.(name).data;
  findx=findx+1;
 end
end

if findx<3
 error('Unable to find mean, standard deviation or effective sample size in nc2')
end

% calculate Welches t-statistic
t = (mu1-mu2)./sqrt(((std1.^2)./equiv1) + ((std2.^2)./equiv2));

% use degrees of freedom from Von Storch and Zwiers (1999)
df = (equiv1+equiv2-2);

% obtain critical t-value for desired two-sided confidence interval
if conf==1 % 99% confidence interval
 tcrit=tinv(0.995,df);
elseif conf==2 % 95% confidence interval
 tcrit=tinv(0.975,df);
elseif conf==3
 tcrit=tinv(0.025,df);
elseif conf==4
 tcrit=tinv(0.005,df);
else
 error('conf has the wrong value - internal error')
end

% test the hypothesis that the means are different
h=ones(size(t));
h(isnan(t))=0;
h(-tcrit<t & t<tcrit)=0;

% place the hypothesis test in an nc structure
nc=nc1;
nc.attributes.title=['Hypothesis test of difference between ',...
                     nc1.attributes.title,' and ',nc2.attributes.title];

if isfield(nc,'turn')
 nc=rmfield(nc,'turn');
end

if isfield(nc,'ANGLE')
 nc=rmfield(nc,'ANGLE');
end

if isfield(nc,[var1,'_samp'])
 nc=rmfield(nc,[var1,'_samp']);
end

nc.diff=nc.(var1);
nc.diff.data=mu1-mu2;
nc.diff.long_name=['DIFFERENCE OF MEANS ',nc.(var1).long_name];

nc.h=nc.(var1);
nc.h.units='';
nc.h.data=h;
nc.h.long_name=['HYPOTHESIS TEST ',nc.(var1).long_name];

nc=rmfield(nc,{[var1,'_std'],[var1,'_equiv'],var1});

nc=ridgepack_struct(nc);

