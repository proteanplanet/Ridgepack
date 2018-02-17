function [alpha]=paper_ridge_alpha(strain,phi,hfi,hdi)

[rhoi,rhos,rhow,delrho,g]=ridge_constants;

if hdi>hfi

 gamma=(phi+phi*strain-strain);

 alpha=2*atan(sqrt(((5*(gamma^2)+6*gamma+3)-...
              sqrt(((5*(gamma^2)+6*gamma+3)^2)-4*(gamma^4)))/(2*(gamma^2))));
 
 alpha=180*alpha/pi;

else

 alpha=0;

end

