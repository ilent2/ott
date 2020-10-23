function stBessel = sphGetXiPsi(nNmax, s, x, NB)
  %% sphGetXiPsi
% Calculates modified Bessel function products for spheroids
% 
% sphGetXiPsi(nNmax, s, x, NB) calculates the products
% xi(n,x)*psi(k,sx) and psi(n,x)*psi(k,sx) for n,k from 0 to nNmax+1. For
% xipsi this calculates only the part that does not exhibit cancellations,
% i.e. F^+_{nk}/x (see Eq. 46 of JQSRT 2013). NB>=N multipoles are used
% for the calculation (for improved precision).
%
% These are calculated for one s(wavelength), and multiple x, for
% (n,k)=0..nNmax+1. The two extremes (0 and nNmax+1) are needed to compute
% the products involving derivatives. Only n+k even are calculated and returned.
% The single regular Bessel (not products) are also returned.
%
% Input:
%         nNmax: [1 x 1] The maxmimum value N required for the integrals;
%                this returns a matrix up to N=nNmax+1, for calculating the
%                derivatives.
%         s: [1 x 1] The relative refractive index of the particle
%         x: [T x 1] The values of x to calculate this at.
%         NB: [1 x 1] The number of n that are used to
%                  calculate the Bessel products (NB>N is needed for large x
%                  to ensure accuracy)
%
% Output:
%         stBessel: A structure containing various Bessel functions and
%                   their products. The fields are
%             - xipsi: [N+2 x N+2 x T] The contributing part of xi_n psi_k
%             - psipsi: [N+2 x N+2 x T] The product psi_n psi_k
%             - chin: [N+2 x T] chi_n(x) (for n=0..N+1)
%             - psin: [N+2 x T] psi_n(x) (for n=0..N+1)
%             - chik: [N+2 x T] chi_k(sx) (for n=0..N+1)
%
% Dependency: 
% sphGetFpovx, vshRBpsi

import ott.tmatrix.smarties.*;

% This call takes care of the chipsi products (which are the ones causing
% problems for spheroids)
[chipsi,chin,psik] = sphGetFpovx(NB+1, s, x);
% chipsi is [NB+2 x NB+2 x X]
% chin, psik are [X x NB+2]

% The rest calculates the normal psipsi products using standard Matlab
% functions to compute the Bessel functions (only up to N+2)
stBessel.psipsi = zeros(nNmax+2, nNmax+2, length(x));

stBessel.psin = vshRBpsi(0:(nNmax+1), x); % [X x N+2]

% To reduce computations, only the products we need are calculated
% i.e. only n+k even (others are left as zero)
for nInd=1:(nNmax+2)
    kInds = (mod(nInd+1,2)+1):2:(nNmax+2);
    stBessel.psipsi(nInd, kInds, :) = transpose( ...
        bsxfun(@times, psik(:, kInds), stBessel.psin(:, nInd))); % 1 x K/2 x X
end

stBessel.xipsi  = stBessel.psipsi + 1i*chipsi(1:(nNmax+2),1:(nNmax+2),:);
stBessel.chin = chin(:,1:(nNmax+2));
stBessel.psik = psik(:,1:(nNmax+2));
end
