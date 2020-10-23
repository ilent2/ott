function [stXipsiAll, stPsipsiAll] = sphGetModifiedBesselProducts(nNmax, s, x, NB)
  %% sphGetModifiedBesselProducts
% Returns matrices of modified Bessel products
% 
% sphGetModifiedBesselProducts(nNmax, s, x, NB) Returns structures
% containing matrices of modified Bessel function products, suitable for
% calculating T-matrix integrals for spheroids. See JQSRT 123, 153 (2013)
% for further details: the modified Bessel products correspond to
% F^+_{nk}(s,x)/x (Eq. 46) and equivalents as defined in Sec. 4.3.
%
% The structure for psipsi also contains functions needed for the
% diagonals of P and Q.
%
% Input:
%       nNmax: [1 x 1] The maximum value of N required in the integrals
%       s:     [1 x 1] The relative refractive index of the particle
%       x:     [T x 1] The values of x at which the integrals are evaluated
%       NB:    [1 x 1] The number of N that should be used when
%                calculating the Bessel function products.
%                For large x, it may be necessary to use NB>N
%
%       stXiPsiAll: Structure containing the matrices for the product
%                   xipsi, as used in Q. The fields are (all [N x N x T])
%           - xipsi: The Bessel products xi_npsi_k
%           - xiprimepsi: xi'_npsi_k (for n+k odd, or n=k)
%           - xipsiprime: xi_npsi'_k (for n+k odd, or n=k)
%           - xiprimepsiprime: xi'_npsi'_k (for n+k even)
%           - xipsiOversxx: xi_npsi_k/(sx^2) (for n+k even)
%           - xiprimepsiprimePlusnnp1xipsiOversxx: (for n+k even)
%                       xi'_npsi'_k + n(n+1)xi_npsi_k/(sx^2)
%           - xiprimepsiprimePluskkp1xipsiOversxx: (for n+k even)
%                       xi'_npsi'_k + k(k+1)xi_npsi_k/(sx^2)
%           - fordiagLt1: [N x T] for Ltilde1 in diagonal
%           - fordiagLt2: [N x T] for Ltilde2 in diagonal
%           - fordiagLt3: [N x T] for Ltilde3 in diagonal
%       stPsiPsiAll: Same with all the xi replaced by psi (for P-matrix)
%                    Note that the names of the fields
%                    include xi to be consistent with the previous structure,
%                    but this is for psipsi.
%
% Dependency: 
% sphGetBesselProductsPrimes [private], sphGetXiPsi

import ott.tmatrix.smarties.*;

    % Start with F^+/x i.e. xipsi and psipsi
    stBessel= sphGetXiPsi(nNmax, s, x, NB);

    % Then deduce from them the equivalent with derivatives
    stXipsiAll = sphGetBesselProductsPrimes(stBessel.xipsi);
    stPsipsiAll = sphGetBesselProductsPrimes(stBessel.psipsi);

    % For diagonals of Q11 and P11 (for Ltilde 1, Eq. 25)
    % s xi'_n psi_n(sx) - xi_n psi'_n(sx)
    % = s xi_n psi_{n+1}(sx) - xi_{n+1} psi_n(sx) [Eq. 65]
    % Same for P
    psinp1psin = transpose(stBessel.psin(:,3:(nNmax+2)) .* stBessel.psik(:, 2:(nNmax+1))); % [N x X]
    psinpsinp1 = transpose(stBessel.psin(:,2:(nNmax+1)) .* stBessel.psik(:, 3:(nNmax+2))); % [N x X]
    xinp1psin = psinp1psin + 1i*transpose(stBessel.chin(:,3:(nNmax+2)) .* stBessel.psik(:, 2:(nNmax+1))); % [N x X]
    xinpsinp1 = psinpsinp1 + 1i*transpose(stBessel.chin(:,2:(nNmax+1)) .* stBessel.psik(:, 3:(nNmax+2))); % [N x X]
    stPsipsiAll.fordiagLt1 = s*psinpsinp1 -psinp1psin;
    stXipsiAll.fordiagLt1 = s*xinpsinp1 -xinp1psin;
    % For diagonals of Q22 and P22 (for Ltilde2 and Ltilde3, Eq. 26)
    psinpsin = transpose(stBessel.psin(:,2:(nNmax+1)) .* stBessel.psik(:, 2:(nNmax+1))); % [N x X]
    xinpsin = psinpsin + 1i*transpose(stBessel.chin(:,2:(nNmax+1)) .* stBessel.psik(:, 2:(nNmax+1))); % [N x X]
    nvec=(1:nNmax).';
    stPsipsiAll.fordiagLt2 = psinpsinp1 - s*psinp1psin + (s-1)*(s+1)/s * ((nvec+1) * (1./x.')).*psinpsin;
    stXipsiAll.fordiagLt2 = xinpsinp1 - s*xinp1psin + (s-1)*(s+1)/s * ((nvec+1) * (1./x.')).*xinpsin;
    stPsipsiAll.fordiagLt3 = bsxfun(@times, psinpsin, 1./(s*x.^2).');
    stXipsiAll.fordiagLt3 = bsxfun(@times, xinpsin, 1./(s*x.^2).');
end




function stBesselPrimes	= sphGetBesselProductsPrimes(prods)
% Calculates modified products involving Bessel function derivatives
% sphGetBesselProductsPrimes(xipsi) calculates the derivatives of the
% bessel function products needed in the integrals, as well as other terms
% calculated in a similar manner.
%
% Input:
%       prods: [N+2 x N+2 x X] Bessel function products (either
%       xi_n(x)psi_k(sx) or psi_n(x)psi_k(sx)), for n=0 to n=N+1, for n+k
%       even. These are obtained from sphGetXiPsi
%
% Output:
%      stBesselPrimes: A structure containing matrices
%                      for the derivatives product. Each matrix is [N x N x X] (for n,k=1..N)
%                      Field names use "xipsi" for convenience but may also correspond
%                       correspond to psipsi if prods=psipsi.
%              - xipsi: The input bessel products (without n=0,N+1) (for
%                       n+k even)
%              - xiprimepsi: xi'_npsi_k (for n+k odd)
%              - xipsiprime: xi_npsi'_k (for n+k odd)
%              - xiprimepsiprime: xi'_npsi'_k (for n+k even)
%              - xipsiOversxx: xi_npsi_k/(sx^2) (for n+k even)
%              - xiprimepsiprimePlusnnp1xipsiOversxx: (for n+k even)
%                       xi'_npsi'_k + n(n+1)xi_npsi_k/(sx^2)
%              - xiprimepsiprimePluskkp1xipsiOversxx: (for n+k even)
%                       xi'_npsi'_k + k(k+1)xi_npsi_k/(sx^2)

N=size(prods,1)-2;
X=size(prods,3);

stBesselPrimes.xiprimepsi = zeros(N,N,X); % [N x N]
stBesselPrimes.xipsiprime = zeros(N,N,X);
stBesselPrimes.xiprimepsiprime = zeros(N,N,X);

stBesselPrimes.xiprimepsiprimePlusnnp1xipsiOversxx = zeros(N,N,X);
stBesselPrimes.xiprimepsiprimePluskkp1xipsiOversxx = zeros(N,N,X);
stBesselPrimes.xipsiOversxx = zeros(N,N,X);

% For xipsi, we have the result already
% We don't want k=0, n=0, and we only want to go to nNmax, not nNmax+1
stBesselPrimes.xipsi = prods(2:end-1, 2:end-1, :);

% Now calculate the products involving derivatives
for nn = 1:N
    nInd=nn+1; % for the indexing of prods
    for kk=(2-mod(nn,2)):2:N % n+k even only
        kInd = kk+1; % for the indexing of prods
        kkp1 = kk*(kk+1);
        nnp1 = nn*(nn+1);

        % From Eq. 61
        stBesselPrimes.xiprimepsiprimePluskkp1xipsiOversxx(nn,kk, :) = (...
            (kk + (nn+1))*(kk+1)*prods(nInd-1, kInd-1, :)...
            +(kkp1 - kk*(nn+1))*prods(nInd-1, kInd+1, :)...
            +(kkp1 - (kk+1)*nn)*prods(nInd+1, kInd-1, :)...
            +(kkp1 + kk*nn)*prods(nInd+1, kInd+1, :))/((2*nn+1)*(2*kk+1));

        % From Eq. 62
        stBesselPrimes.xiprimepsiprimePlusnnp1xipsiOversxx(nn,kk, :) = (...
            (nnp1 + (kk+1)*(nn+1))*prods(nInd-1, kInd-1, :)...
            +(nn - kk)*(nn+1)*prods(nInd-1, kInd+1, :)...
            +((nn+1) - (kk+1))*nn*prods(nInd+1, kInd-1, :)...
            +(nnp1 + kk*nn)*prods(nInd+1, kInd+1, :))/((2*nn+1)*(2*kk+1));

    end
    for kk=(1+mod(nn,2)):2:N % n+k odd only
        kInd = kk+1; % for the indexing of prods
        % Eqs. 59-60
        stBesselPrimes.xiprimepsi(nn,kk, :) = (prods(nInd-1, kInd, :)*(nn+1) - nn*prods(nInd+1,kInd, :))/(2*nn+1);
        stBesselPrimes.xipsiprime(nn,kk, :) = (prods(nInd, kInd-1, :)*(kk+1) - kk*prods(nInd,kInd+1, :))/(2*kk+1);
    end
end

end
