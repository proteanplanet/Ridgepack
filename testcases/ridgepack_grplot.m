function ridgepack_grplot

hF=2;
hFs=0.3;
epsilon=-1/3;
phi=0;

hstep=0.01;
hmax=20;
hgrid=[hstep:hstep:hmax];
gr=zeros(size(hgrid));

[VR,ALPHAHAT,HK,HS,LK,LS]=ridgepack_energetics(hF,hFs,epsilon,phi);

% Thickness distribution on a grid
h1idx=find(hgrid>=hF & hgrid<=hF+(LK-LS)*tand(ALPHAHAT)/2);
h2idx=find(hgrid>hF+(LK-LS)*tand(ALPHAHAT)/2 & hgrid<=hF+(LK+LS)*tand(ALPHAHAT)/2);

gr(h1idx) = 2/(LK*tand(ALPHAHAT));
gr(h2idx) = 1/(LK*tand(ALPHAHAT));

area=sum(hstep*gr)

gr(h1idx)=gr(h1idx)/area;
gr(h2idx)=gr(h2idx)/area;

area=sum(hstep*gr)

clf
plot(hgrid,gr)
hold on

% Thickness distribution in continuous space
step(1)=hF;
g(1)=0;

step(2)=hF;
g(2)=2/(LK*tand(ALPHAHAT));

step(3)=hF+(LK-LS)*tand(ALPHAHAT)/2;
g(3)=2/(LK*tand(ALPHAHAT));

step(4)=hF+(LK-LS)*tand(ALPHAHAT)/2;
g(4)=1/(LK*tand(ALPHAHAT));

step(5)=hF+(LK+LS)*tand(ALPHAHAT)/2;
g(5)=1/(LK*tand(ALPHAHAT));

step(6)=hF+(LK+LS)*tand(ALPHAHAT)/2;
g(6)=0;

area=g(3)*(LK-LS)*tand(ALPHAHAT)/2+g(5)*(2*LS)*tand(ALPHAHAT)/2

plot(step,g)

legend({'Discrete','Continuous'})


















