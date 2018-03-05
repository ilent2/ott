function [x,y]=uv2xy(u,v,xi);
%Transforms (u,v) coordinates into cartesian (x,y) for a given
%ellipticity xi.
%
% Usage:
%
% [x,y] = ux2xy(u,v,xi);
%
% when xi is really big the conversion from elliptical coordinates may
% prevent x having a non-zero value because (1-eps) ~ 1.
%
% PACKAGE_INFO

x=xi*cosh(u).*cos(v);
y=xi*sinh(u).*sin(v);
