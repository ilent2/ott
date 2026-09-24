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
% CONVENTION: these are the standard Condon-Shortley harmonics, i.e. the same
% sign convention as scipy.special.sph_harm, treams, and the usual tables:
%
%     Y_1^1 = -sqrt(3/(8*pi)) sin(theta) exp(i*phi)      (note the minus)
%
% Upstream OTT OMITS that phase: legendrerow implements a geodesy recursion
% in which it is conventionally absent, and nothing downstream restored it.
% The shipped output therefore differs from the harmonics above by a factor
% (-1)^m: measured ratio exactly -1 for odd m and +1 for even m, at both
% signs of m.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

import ott.utils.*

ott.warning('internal');

if length(n)>1
    ott.warning('external');
    error('n must be a scalar at present')
end

if nargin<4
  if nargin < 3
    error('Not enough input arguments');
  end
  
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

% NOTE: the (-1)^m below is NOT the Condon-Shortley phase.
% legendrerow only computes rows for |m|, so a negative-m request is handed 
% the same P_n^|m| as its positive; this factor supplies the negative-order 
% relation that turns one into the other, and applies to m < 0 ONLY.
pnm = [(-1).^mv(m<0,:).*pnm(m<0,:);pnm(m>=0,:)];

expphi = exp(1i*mv.*phiM);

%N = sqrt((2*n+1)/(8*pi));

Y = pnm .* expphi;

% ===================== FIX 1: CONDON-SHORTLEY PHASE =====================
% Up to this line Y is the PHASE-FREE harmonic: legendrerow never applies the
% (-1)^m that the Condon-Shortley convention carries inside P_l^m.  So Y here
% and the standard Y_n^m differ by a factor (-1)^m -- for (n,m) = (1,1) this
% line receives +sqrt(3/(8*pi)) sin(th) exp(i*ph) where the standard Y_1^1 is
% MINUS that.  Multiplying by (-1)^m supplies the missing phase and lands on
% the standard convention.
Y = Y .* (-1).^mv;

% Do we want to calculate the derivatives?
if nargout <= 1
   Y=Y.';
   % Doesn't look like it
   ott.warning('external');
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

% Ytheta and Yphi are built by recursions relating Y(n,m) to Y(n,m+-1).  Those
% neighbours now carry (-1)^(m+-1) = -(-1)^m, so each derivative acquires one
% extra minus sign relative to the patched Y.  Undo it.
Ytheta = -Ytheta;
Yphi   = -Yphi;

Y=Y(n+mi+1,:).';
Yphi=Yphi(n+mi+1,:).';
Ytheta=Ytheta(n+mi+1,:).';

ott.warning('external');
