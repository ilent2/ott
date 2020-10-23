function NB = sphEstimateNB(NQ, stGeometry, stParams, acc)
  %% sphEstimateNB
% Finds the number of n required for accurate modified Bessel products
% 
% sphEstimateNB(stGeometry, stParams, acc) finds an estimate for NB to ensure
% required accuracy in Bessel products
% Input:
%       NQ:         Minimum number of multipoles required
%       stGeometry: Structure containing geometric information, as from
%                   sphMakeGeometry
%       stParams:   Structure containing simulation parameters k1 and s
%       (possibly vectors here)
%       acc (optional): relative accuracy required, default is 1e-13
% Output: NB scalar with number of multipole required
%
% Dependency: 
% sphCheckBesselConvergence [private], sphGetFpovx 

import ott.tmatrix.smarties.*;

%set default accuracy if not specified
if nargin < 4
    acc = 1e-13;
end

s =  stParams.s;
k1 = stParams.k1;

% Find max size parameters
xmax=max(k1)*max(stGeometry.r);
% Find max and min relative refractive index
[~,ind]=min(abs(s));
smin=s(ind);
[~,ind]=max(abs(s));
smax=s(ind);

% find required NB for each of the 2 extreme cases
N1 = sphCheckBesselConvergence(NQ,smax,xmax,acc,NQ);
NB = sphCheckBesselConvergence(NQ,smin,xmax,acc,N1);

end

function NB = sphCheckBesselConvergence(Nreq, s, x, acc, Nmin)
% Determine required NB values for Bessel product
% sphCheckBesselConvergence(Nmin, s, x, acc) returns the required value
% of NB for calculations
% of the Bessel function product for arguments x to be valid up to a given
% acuracy acc. This is assessed by checking the relative accuracy of the
% last row of the F^+/x matrix.
%
%  Input:
%		Nreq	   The maximum value of N that is required to be sufficently
%					accurate
%		s			Relative refractive index of the particle
%		x			A vector containing the sizes to check, which can be
%					representative of the whole particle
%       acc         The required relative accuracy that the results are
%					meant to converge within
%       Nmin        Minimum NB (where the search starts)

import ott.tmatrix.smarties.*;

% initial guess of N is the smallest possible size
NBstart= max(Nmin,Nreq);
NB = NBstart;
prod = sphGetFpovx(NB, s, x);

ToContinue = true;

%The maximum value of N we are prepared to go to
maxN = 500;
NBstep=16;

while(ToContinue && NB < maxN)
    NBnext=NB+NBstep;
    prodNew = sphGetFpovx(NBnext, s, x);

        % Worst relative accuracy in all matrix up to n=Nreq+1 (+1 for
        % derivatives to also be accurate)
        relAccee  = max(max(abs(prod(1:2:(Nreq+1),1:2:(Nreq+1))./prodNew(1:2:(Nreq+1),1:2:(Nreq+1))-1)));
        relAccoo  = max(max(abs(prod(2:2:(Nreq+1),2:2:(Nreq+1))./prodNew(2:2:(Nreq+1),2:2:(Nreq+1))-1)));
        relAcc=max(relAccee,relAccoo);
    if relAcc < acc
        % we have found sufficiently high n, fine-tune NB by going
        % back and changing step
        ToContinue = false;
    else
        NB = NBnext;
        prod = prodNew;
    end
end

% Then repeat with small step to fine-tune NB
if NB>NBstart
    NB=NB-NBstep;
    NBstep=1;
    ToContinue=true;
    while(ToContinue && NB < maxN)
        NB=NB+NBstep;
        prod = sphGetFpovx(NB, s, x);

        % We now use prodNew from previous step as our reference

        % Worst relative accuracy in all matrix up to n=Nreq+1 (+1 for
        % derivatives to also be accurate)
        relAccee  = max(max(abs(prod(1:2:(Nreq+1),1:2:(Nreq+1))./prodNew(1:2:(Nreq+1),1:2:(Nreq+1))-1)));
        relAccoo  = max(max(abs(prod(2:2:(Nreq+1),2:2:(Nreq+1))./prodNew(2:2:(Nreq+1),2:2:(Nreq+1))-1)));
        relAcc=max(relAccee,relAccoo);

        if relAcc < acc
            % we have found sufficiently high n
            ToContinue = false;
        end
    end

end

if ToContinue
    error('Problem in sphEstimateNB: convergence was not achieved');
end

end
