function [stCoa, CstTRa] = slvForT(stParams, stOptions, stGeometry)
  %% slvForT
% Calculates the T-matrix and orientation-averaged properties
% Input:
%       - stParams: struct
%              The following parameters should be defined in stParams:
%              - a: semi-axis along x,y
%              - c: semi-axis along z
%              - k1: wavevector in embedding medium (of refractive index nM) (k1=2*pi*nM/lambda)
%              - s: relative refractive index (s=n_Particle / nM)
%              - N: number of multipoles for T-matrix
%              - nNbTheta: number of thetas for quadratures
%       - stOptions: struct with optional parameters, see
%              slvGetOptionsFromStruct for details.
%       - stGeometry (optional): struct with geometry
%
% Output:
%       - stCoa: struct with orientation averaged cross-sections
%                Cext, Csca, Cabs
%       - CstTRa: cell {1 x M} of structs defining the T (and possibly R)
%                matrices
%
% Dependency: 
% rvhGetAverageCrossSections, rvhGetSymmetricMat, rvhGetTRfromPQ,
% rvhTruncateMatrices, slvGetOptionsFromStruct, sphCalculatePQ,
% sphEstimateDelta, sphEstimateNB, sphMakeGeometry

import ott.tmatrix.smarties.*;

c = stParams.c;
a = stParams.a;

s = stParams.s;
k1 = stParams.k1;

% For convenience, k1 and s are stored in a struct
stk1s.k1=k1;
stk1s.s=s;

N = stParams.N;
nNbTheta = stParams.nNbTheta;

[bGetR,Delta,NB,absmvec,bGetSymmetricT, bOutput] = slvGetOptionsFromStruct(stParams,stOptions);

stk1s.bOutput=bOutput;


% Make structure describing spheroidal geometry and quadrature points for
% numerical integrations
if nargin<3
    stGeometry = sphMakeGeometry(nNbTheta, a, c);
end

if Delta<0 % then need to estimate Delta
    [Delta, T2211err]= sphEstimateDelta(stGeometry, stk1s);
    if isnan(Delta)
        disp ('ERROR: Delta could not be found. Results are likely to be non-converged. Try choosing Delta manually instead.');
        return;
    end
    disp (['Delta estimated to \Delta=', int2str(Delta),' with relative error in T_{11}^{22,m=1} of ',  num2str(T2211err)]);
end

NQ = N+Delta;% NQ>=N: Maximum multipole order for computing P and Q matrices

% Estimating NB
if NB<=0
    NB=sphEstimateNB(NQ, stGeometry, stk1s);
end
if NB<NQ
    NB=NQ; % NB must be at least NQ
end

% Calculates P and Q
CstPQa = sphCalculatePQ(NQ, absmvec, stGeometry, stk1s, NB);

% Get T (and possibly R)
CstTRa = rvhGetTRfromPQ(CstPQa,bGetR);

% If needed, discard higher order multipoles
% (which are affected by the finite size of P and Q)
if NQ>N
    CstTRa = rvhTruncateMatrices(CstTRa, N);
 end
% T and R matrices now go up to N multipoles

% If required, symmetrize the T-matrix
if bGetSymmetricT
    CstTRa = rvhGetSymmetricMat(CstTRa, {'st4MT'});
end

% Calculate the (Ext, Abs, Sca) cross-sections for orientation-averaged excitation
stCoa = rvhGetAverageCrossSections(k1, CstTRa);

end
