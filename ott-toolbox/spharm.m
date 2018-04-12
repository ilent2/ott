function [Y,Ytheta,Yphi] = spharm(n,m,theta,phi)
% SPHARM scalar spherical harmonics and angular partial derivatives.
%
% Y = SPHARM(n,m,theta,phi) calculates scalar spherical harmonics.
%
% [Y,Ytheta,Yphi] = SPHARM(n,m,theta,phi) additionally, calculates
% the angular partial derivatives dY/dtheta and 1/sin(theta)*dY/dphi.
%
% SPHARM(n,theta,phi) as above but for all m.
%
% Scalar n for the moment.
% 
% If scalar m is used Y is a vector of length(theta,phi) and is
% completely compatible with previous versions of the toolbox. If vector m
% is present the output will be a matrix with rows of length(theta,phi) for
% m columns.
%
% "Out of range" n and m result in return of Y = 0
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

ott_warning('internal');

if length(n)>1
    ott_warning('external');
    error('n must be a scalar at present')
end

if nargin<4
    phi=theta;
    theta=m;
    m=[-n:n];
end

%this is a cop out meant for future versions.
if nargout>1
    mi=m;
    m=[-n:n];
end

m=m(abs(m)<=n);

[theta,phi] = matchsize(theta,phi);
input_length = length(theta);

% if abs(m) > n | n < 0
%    Y = zeros(input_length,1);
%    Ytheta = zeros(input_length,1);
%    Yphi = zeros(input_length,1);
%    return
% end

pnm = legendrerow(n,theta);
%pnm = pnm(abs(m)+1,:).';

% Why is this needed? Better do it, or m = 0 square integrals
% are equal to 1/2, not 1.
% This is either a bug in legendre or a mistake in the docs for it!
% Check this if MATLAB version changes! (Version 5.X)
%pnm(1,:) = pnm(1,:) * sqrt(2);

pnm = pnm(abs(m)+1,:); %pick the m's we potentially have.

[phiM,mv]=meshgrid(phi,m);

pnm = [(-1).^mv(m<0,:).*pnm(m<0,:);pnm(m>=0,:)];

expphi = exp(1i*mv.*phiM);

%N = sqrt((2*n+1)/(8*pi));

Y = pnm .* expphi;

% Do we want to calculate the derivatives?
if nargout <= 1
   Y=Y.';
   % Doesn't look like it
   ott_warning('external');
   return
end

% We use recursion relations to find the derivatives, choosing
% ones that don't involve division by sin or cos, so no need to
% special cases to avoid division by zero

% exp(i*phi),exp(-i*phi) are useful for both partial derivatives
expplus = exp(1i*phiM);
expminus = exp(-1i*phiM);

% theta derivative
% d/dtheta Y(n,m) = 1/2 exp(-i phi) sqrt((n-m)(n+m+1)) Y(n,m+1)
%                 - 1/2 exp(i phi) sqrt((n-m+1)(n+m)) Y(n,m-1)

ymplus=[Y(2:end,:);zeros(1,length(theta))];
ymminus=[zeros(1,length(theta));Y(1:end-1,:)];

Ytheta = sqrt((n-mv+1).*(n+mv))/2 .* expplus .* ymminus ...
         - sqrt((n-mv).*(n+mv+1))/2 .* expminus .* ymplus;

% phi derivative - actually 1/sin(theta) * d/dphi Y(n,m)
% Note that this is just i*m/sin(theta) * Y(n,m), but we use a
% recursion relation to avoid divide-by-zero trauma.
% 1/sin(theta) d/dphi Y(n,m) = 
% i/2 * [ exp(-i phi) sqrt((2n+1)(n+m+1)(n+m+2)/(2n+3)) Y(n+1,m+1)
%     + exp(i phi) sqrt((2n+1)(n-m+1)(n-m+2)/(2n+3)) Y(n+1,m-1) ]

Y2 = spharm(n+1,theta,phi).';

ymplus=Y2(3:end,:);
ymminus=Y2(1:end-2,:);

% size(ymplus)
% size(mv)
% size(expminus)

Yphi = 1i/2 * sqrt((2*n+1)/(2*n+3)) * ...
   ( sqrt((n+mv+1).*(n+mv+2)) .* expminus .* ymplus ...
   + sqrt((n-mv+1).*(n-mv+2)) .* expplus .* ymminus );

Y=Y(n+mi+1,:).';
Yphi=Yphi(n+mi+1,:).';
Ytheta=Ytheta(n+mi+1,:).';

ott_warning('external');

return
