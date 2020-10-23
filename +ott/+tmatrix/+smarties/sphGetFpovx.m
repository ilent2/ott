function [Fpovx, rbchi, rbpsi] = sphGetFpovx(nNmax, s, x)
  %% sphGetFpovx
% Calculate F^+/x (see Eq. 46 of JQSRT 2013)
% 
% sphGetFpovx(nNmax, s, x) calculates the matrix F^+/x where
%	Fpovx = P^+(x chi_n(x) psi_k(sx))/x for 0<=(n,k)<=nNmax, for one s (wavelength),
%	and x can be a vector of length [T x 1]
%
%	This function calculates F^+/x=F/x in regions where there are no
%	cancellations, and calls sphGetFpRow to calculate the last row
%	of F^+ for the region where there are cancellations. It then fills up
%	the matrix using the "westward" recursion scheme (solving for n, k-1,
%   see Fig. 3c). See Sec. 4.2 for further details
%
% Input:
%        nNmax: [1 x 1] The maximum value of n that is desired.
%        s:     [1 x 1] The relative refractive index of the particle
%        x:     [T x 1] The values of x at which the function is evaluated
%
% Output:
%        Fpovx: [N+1 x N+1 x T] The matrix F^+/x
%        rbchi: [T x N+1] The Riccati-Bessel function chi_n(x)
%        rbpsi: [T x N+1] The Riccati-Bessel functoin psi_k(sx)
%
% Dependency: 
% sphGetFpRow, vshRBchi, vshRBpsi

import ott.tmatrix.smarties.*;

numX = length(x);

rbpsi=vshRBpsi(0:(nNmax), s.*x); % [T x K+1]
rbchi=vshRBchi(0:(nNmax), x); % [T x N+1]

Fpovx = zeros(nNmax + 1, nNmax+1, numX); % [N+1 x K+1 x T]
% for xind=1:numX
%     % column-row product gives matrix [N+1 x K+1]
%    Fpovx(:,:,xind) = rbchi(xind,:).' * rbpsi(xind,:);
% end

% To initialize the Westward recurrence, we need the last row (n=N),
% and the subdiagonal terms.

% This first calculates the last row (n=N, all k) for Fpovx
FpRow=sphGetFpRow(nNmax, s, x);
Fpovx(nNmax+1, 1:(nNmax-4+1), :) = bsxfun(@times,FpRow,1./x.');

% do n, k-1 (West) recursion [Solving for F^+_{n,k-1} in Eq. 51 (dividing
% all terms by x]
% We only do n+k even
for kk=nNmax:-1:0
    kInd=kk+1; % kk is the k-value, kInd is the index
    % First fill in the matrix where there are no cancellations (n=0..k+2)
    % nMin=0 (if k even) or 1 (if k odd), nMax=min(k+2,N)
    indnDirect = (1+mod(kk,2)):2:(1+min(kk+2,nNmax)); % [1 x Nred]
    % Here F^+/x = F/x = chi_n(x) * psi_k(s*x)
    % x-values are placed on the third dimension
    Fpovx(indnDirect,kInd,:) = transpose(bsxfun(@times,rbchi(:,indnDirect),rbpsi(:,kInd))); % [Nred x X]
    % Now calculates all terms in column k-1 for all relevant problematic n
    if kk>0 % for all non zero k
        for nn=(kk+3):2:(nNmax-1) % n=k+2 to N-1 with only n+k odd
            nInd=nn+1; % nn is the n-value, nInd is the index
            Fpovx(nInd, kInd-1, :) = (Fpovx(nInd+1,kInd, :) + Fpovx(nInd-1,kInd, :))*(2*kk+1)/(2*nn+1)/s - Fpovx(nInd,kInd+1, :);
        end
    end
end

end
