function [Delta, err] = sphEstimateDelta(stGeometry, stParams, NQmax)
  %% sphEstimateDelta
% Estimates Delta from the convergence of T^{22,m=1}_{11} [see JQSRT2015]
% 
% sphEstimateDelta(stGeometry, stParams, NQmax, acc)
% finds an estimate for Delta by studying the convergence of the
% T^{22,m=1}_{11} matrix element. Delta+1 is the NQ for which this element
% has reached convergence.
% See JQSRT 160, 29 (2015) for further details.
% Also returns the estimated converged precision.
%
% Input:
%       stGeometry:   Structure containing geometric information, as from
%                   sphMakeGeometry
%       stParams:   Structure containing simulation parameters, as from
%                   tmsMakeParams or tmsMakeParamsLambda (only the first
%                   lambda is considered)
%       NQmax (optional - default is 80):
%                   The maximum number of multipoles
%
% Dependency: 
% rvhGetTRfromPQ, rvhTruncateMatrices, sphCalculatePQ

import ott.tmatrix.smarties.*;

%set default parameters
if nargin < 3
    NQmax = 80;
end

minAcc = 1e-4;

absmvec = 1; % only one m

% This works on only one wavelength, so we choose the largest k1 * s
% as representative of the worst case
% Find max and min relative refractive index
[~,ind] = max(abs(stParams.k1 .* stParams.s));

stParam1.s =stParams.s(ind);
stParam1.k1 =stParams.k1(ind);

% Number of points used for linear fit to check if error has reached a
% plateau
NforConv = 5;
Afit = [ones(NforConv,1), (1:NforConv).'];
% Store errors
T2211err = zeros(floor((NQmax+1)/2),1);
% Calculates P,Q
CstPQa = sphCalculatePQ(NQmax, absmvec, stGeometry, stParam1, NQmax);
T2211=0;

for N = 1:2:NQmax % Loop over truncation (only odd numbers)
    % Truncate P,Q to N and get corresponding T
    nCount=(N+1)/2;
    CstTRa = rvhGetTRfromPQ(rvhTruncateMatrices(CstPQa, N),false);
    T2211new = CstTRa{1}.st4MTeo.M22(1,1);

    % Relative error
    T2211err(nCount) = abs (T2211./T2211new-1);
%    fprintf('N=%d\t T2211=%.16g\t err=%.16g', N, T2211new,T2211err(nCount));

    T2211=T2211new;
    % Minimum requirement for convergence
    if nCount > NforConv && T2211err(nCount)<minAcc
        % Then test if the last NforConv errors have flattened out by a
        % linear regression to the last NforConv points
        coef=Afit \ log10(abs(T2211err((nCount-NforConv+1):nCount)));
        slope = coef(2);
        % if slope<0, then still converging
        if slope>-0.5 % this means less than one digit better over 2NforConv steps
            Delta = N-2*NforConv+2;
            err = mean(T2211err((nCount-NforConv+1):nCount));
            return;
        end
    end

end
%disp('Convergence for Delta was not found..')
Delta = NaN;
err=0;
