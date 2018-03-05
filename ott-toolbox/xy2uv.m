function [u,v]=xy2uv(x,y,xi);
%Transforms (x,y) coordinates into elliptical coordinates (u,v) for a given
%ellipticity xi.
%
% Usage:
%
% [u,v] = xy2uv(x,y,xi);
%
% warning relative parameter sizes greater the 1e14 will prevent x from
% being resolved in calculation. This is due to the limitation of number
% representation of (1-eps) ~ 1 and cannot be resolved except with
% arbitrary precision.
%
% PACKAGE_INFO

rP=sqrt((x/xi+1).^2+(y/xi).^2);
rM=sqrt((x/xi-1).^2+(y/xi).^2);

plus=(rP+rM)/2;
minus=(rP-rM)/2;

u=acosh(plus);
v=acos(minus);

v(y<0)=-v(y<0);

%small xi breaks v... approximation is that v ~ theta
smallvar=sqrt(y(:).^2+x(:).^2)>1e14*xi;
v(smallvar)=atan2(y(smallvar),x(smallvar));

%big xi breaks u... approximation is that sinh(u) ~ u ~ abs(y)/xi
%this fails because of round off error (1-eps) ~ 1.
bigvar=(sqrt((y(:)).^2+(x(:)).^2))<1e-14*xi;
u(bigvar)=abs(y(bigvar))/xi;
