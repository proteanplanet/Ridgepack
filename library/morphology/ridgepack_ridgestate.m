function [PHI,ALPHAHAT,HR,HRs]=ridgepack_ridgestate(hF,hFs,theta,epsilon)

HR=hF./(1-epsilon);

[EPSILON,PHI,ALPHAHAT]=ridgepack_trajectory(hF,hFs);

[HFD,HFF,HDD,HDF,HK,HS,LK,LS]=ridgepack_morphology(hf,hfs,hd,hds,phi,alpha);













