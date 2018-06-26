
% This script tests the rotary coherence script for coherence by adding white
% noise to two time series with the same signal.


nc.attributes.title='coherence test'

nc.time.data=[0:0.001:300];
nc.time.units='days';
nc.time.dimension={'time'};
nc.time.long_name='time';

nc.x.data=cos(0.5+10*nc.time.data)+sin(15*nc.time.data);
nc.y.data=sin(0.5+10*nc.time.data)+cos(15*nc.time.data);
nc.x.units='m/s';
nc.y.units='m/s';
nc.x.dimension={'time'};
nc.y.dimension={'time'};
nc.x.long_name='x';
nc.y.long_name='y';

nc2=nc;

nc.x.data=nc.x.data+random('poisson',0,1,size(nc.x.data));
nc.y.data=nc.y.data+random('poisson',0,1,size(nc.y.data));

nc2.x.data=nc.x.data+random('normal',0,1,size(nc.x.data));
nc2.y.data=nc.y.data+random('normal',0,1,size(nc.y.data));

ts=ridgepack_2vectseries(nc,'x',nc,'y','spline',[0 0 0 1 0 0.0]);
ts2=ridgepack_2vectseries(nc2,'x',nc2,'y','spline',[0 0 0 1 0 0.0]);

ridgepack_coherence(ts,ts2,0.125)


