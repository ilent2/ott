function [B,C,P] = vsh(n,m,theta,phi)
% VSH calculate vector spherical harmonics
%
% [B,C,P] = VSH(n,m,theta,phi) calculates vector spherical harmonics
% for the locations theta, phi.  Vector m allowed.  Scalar n for the moment.
%
% [B,C,P] = VSH(n,theta,phi) outputs for all possible m.
%
% If scalar m: B,C,P are arrays of size length(theta,phi) x 3
% If vector m: B,C,P are arrays of size length((theta,phi),m) x 3
% theta and phi can be vectors (of equal length) or scalar.
%
% The three components of each vector are [r,theta,phi]
%
% "Out of range" n and m result in return of [0 0 0]

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

import ott.utils.*
ott.warning('internal');

if length(n)>1
    ott.warning('external');
    error('n must be a scalar in this version')
end

if nargin<4
    phi=theta;
    theta=m;
    m=[-n:n];
end

% Convert a scalar theta or phi to a vector to match a vector
% partner
[theta,phi] = matchsize(theta,phi);

[Y,Ytheta,Yphi] = spharm(n,m,theta,phi);

%this makes the vectors go down in m for n. has no effect if old version
%code.

Z = zeros(size(Y));

B = [Z,Ytheta,Yphi];

C = [Z,Yphi,-Ytheta];

P = [Y,Z,Z];

ott.warning('external');
