function [S,lossPrecS] = sphGetFpRow(n, s, x)
  %% sphGetFpRow
% Calculates problematic Bessel products for one row
% 
% sphGetFpRow(n, s, x) Calculates the product F^+_nk = P^+(xchi_n(x)psi_k(sx))
% for one n and all k with cancellations, i.e. all 0<=k<=n-4 and n+k even,
% for one value of s (so one wavelength), and possibly many x's.
%
% This uses a series implementation to calculate the product, and should
% treat well both large and small (near 1) values of s, and large and small
% values of x. The result of this may be used along with a recursion scheme
% in order to calculate all the required terms for the EBCM integrals for
% spheroids.
% See appendix B of [JQSRT 2013] for further details
%
% Input:
%           n: [1 x 1] The value of n to use; typically we will do this for
%              the last row, i.e. n=N+1, where N is the
%              highest value of n required in calculations
%           s: [1 x 1] The relative refractive index of the particle.
%           x: [X x 1] The values of  x to calculate for.
%
% Output:
%           S: [K=n-3 x X] The values that are returned. This is for k=0 to
%           k=n-4, and entries for odd n+k are zero.
%           lossPrecS: [K=n-3 x X] Potential loss of precision of precision
%           estimated as abs(MaxTerm in series)/S. Problem if lossPrecS is
%           large, i.e. 10^4 means 4-digit loss.
%
% Dependency: 
% sphGetUforFp [private]

numX = length(x);

%qmax = max(3,n-1); % this is where the u branch finishes
xsquared = reshape(x.^2, 1, numX); % [1 x X] get x^2 as a row
u = sphGetUforFp(n); % [B+2 x B+1] where B=floor(n/2)

alphaBark= 1; % \bar{alpha}_k [Eq. B.23]

S = zeros(n-4 + 1, numX); % [K x X] This is the partial sum in the series
% in Eq. B.20 needed for F^+, we will add terms to it sequentially
maxTermS = zeros(n-4 + 1, numX); % [K x X] Keeps track of maximum term
% in the series for potential loss of precision
lossPrecS = zeros(n-4 + 1, numX); % [K x X] Keeps track of loss of precision

% Computing the beta_{jq} [Eq. B.14]
beta = zeros(n, n-1); % [j x q] 0<=j<=q, 1<=q<=qint-1
%fill the first two beta columns, that is q=1, q=2
beta(1, 1) = 1; % j=0, q=1
beta(2, 1) = beta(1,1) * (s-1)*(s+1); % j=1, q=1;
% (s-1)*(s+1) is used instead of s^2-1 to reduce problems when s is close
% to 1

beta(1, 2) = 1; % j=0, q=2
beta(2, 2) = beta(1,2) * 2 *(s-1)*(s+1); % j=1, q=2
beta(3, 2) = beta(2,2) * 1/2*(s-1)*(s+1); % j=2, q=2

numToTest = 3; % this is the number of consecutive tests required to consider
% convergence has been reached (a single term goes to zero by chance)

for kk = n-4:-2:0
    kInd = kk + 1;
    qmin = (n-kk)/2-1;
    qint = n-kk; % "intermediate" value of q, where we switch which method applies
    alphaBark = -alphaBark/(n-kk-2); % [Eq. B.24]

    %fill the last two beta columns here, as this is the first time that we need
    %them
    beta(1, qint-1) = 1; % j=0, q=qint-1
    for jj = 1:(qint-1)
        jInd = jj+1;
        beta(jInd, qint-1) = (qint -1- jj+1)/(jj)*beta(jInd-1,qint-1)*(s-1)*(s+1);
    end
    beta(1, qint) = 1; % j=0, q=qint
    for jj = 1:qint
        jInd = jj+1;
        beta(jInd, qint) = (qint - jj+1)/(jj)*beta(jInd-1,qint)*(s-1)*(s+1);
    end
    % more initializations
    alphaq = alphaBark*ones(1, numX); % [1 x X] alpha_{qnk}(x) [Eq. B.21] for q=qmin
    xNotConverged = true(1, numX);
    currentTerm = zeros(1, numX);
    counter = zeros(1, numX);
    test= false(1, numX);

    % Method 1 for qmin<=q<=n-k-1 (Section B.1)
    for qq = qmin:(qint-1)
        b = n-kk-qq-1;
        % Computes gamma_{qnk} from the sum in Eq. B.13
        gammaqnk=0;
        for jj=max(0,2*qq+kk-n+1):qq
            jInd = jj + 1;
            gammaqnk = gammaqnk + beta(jInd, qq) * u(qq-jj+1, b+1);
        end

        % For the finite series of Method 1, we only test convergence for
        % |s|>=2. If s is close to 1, the series may appear to converge but
        % the terms then increase again, so to be safe we take all the
        % terms

        if abs(s)<2
            % |s|<2, no convergence test
                currentTerm = alphaq * gammaqnk;
                S(kInd, :) = S(kInd, :) + currentTerm; % [1 x X]
                maxTermS (kInd,:) = max(abs(currentTerm),maxTermS (kInd,:));
%                disp (['k=', int2str(kk), ' q=', int2str(qq), ' S=', num2str(S(kInd,:)), ' term=', num2str(currentTerm)]);
        else
            % Test if convergence of the series S [Eq. B.20] has been obtained
            test(:)=false; % reset the test flag to Converged
            currentTerm(xNotConverged)  = alphaq(xNotConverged) * gammaqnk;
            if nnz(xNotConverged) > 0
                % converged if added term is less than epsilon times current
                % value of the series. test=false if this is the case
                % test=true for notConverged x-values
                test(xNotConverged) = abs(currentTerm(xNotConverged)) > eps() * abs(S(kInd, xNotConverged));
            end
            if nnz(test) > 0 % for the non-converged x values
                if ~isfinite(currentTerm(test))
                    disp('Problem (1) in sphGetFpovx...');
                end
                % Add the current term to the series
                S(kInd, test) = S(kInd, test) + currentTerm(test);
                maxTermS (kInd,test) = max(abs(currentTerm(test)),maxTermS (kInd,test));
            end

            % Full convergence is assumed if the convergence test is fulfilled
            % for numToTest (=3 here) consecutive terms
            counter  = counter + ~test; % count number of positive converged tests
            counter(test) = 0; % reset at zero if not-converged
            xNotConverged = counter<numToTest; % number of consecutive positive converged outcomes

            % Break if convergence of all terms has been obtained
            if nnz(xNotConverged) == 0
                break;
            end
        end
        % Update alphaq for next term (using Eq. B.22]
        alphaq = alphaq.*-xsquared /(2*(qq+1)); % [1 x X]

    end

    % We now move on to Method 2 for q>=qint [Sec. B.2]
    % Reinitializing counters - convergence is now tested in all cases
    % (otherwise the series would be infinite...)
    xNotConverged(:) = true;
    counter(:) = 0;

    % Initializing c_{qqnk} [Eq. B.19]
    cqqnk = 1/(2*kk+1); % this is for q=qint=n-k
    qq=qint;
    while true % upper limit of series will be decided by convergence test

        ciqnk=cqqnk; %i=q term
        gammaqnk = cqqnk; % initialize gammaqnk =cqqnk (this is i=q) and add the other
        % term according to Eq. B.3
        for ii = qq-1:-1:0 % sum from i=0 to i=qq-1
            % calculates ciqnk by reccurence using Eq. B.17
            ciqnk = ciqnk*s^2*(ii+1)*(2*ii+1-2*n)/(qq-ii)/(2*kk +2*qq-2*ii+1);
            gammaqnk = gammaqnk + ciqnk; % add term
        end

        % Convergence testing
        test(:)=false;
        currentTerm(xNotConverged)  =alphaq(xNotConverged) * gammaqnk;
        test(xNotConverged) = abs(currentTerm(xNotConverged)) > eps() * abs(S(kInd, xNotConverged));
        if nnz(test) > 0
            if ~isfinite(currentTerm(test))
                    disp('Problem (2) in sphGetFpovx...');
            end
            % Add the current term to the series
            S(kInd, test) = S(kInd, test) + currentTerm(test);
            maxTermS (kInd,test) = max(abs(currentTerm(test)),maxTermS (kInd,test));
%            disp (['k=', int2str(kk), ' q=', int2str(qq), ' S=', num2str(S(kInd,:)), ' term=', num2str(currentTerm)]);
        end


        counter  = counter + ~test;
        counter(test) = 0;
        xNotConverged = counter<numToTest;

        if nnz(xNotConverged) == 0
            break;
        end

        % Update alphaq for next term (using Eq. B.22]
        alphaq(xNotConverged) = alphaq(xNotConverged).*-xsquared(xNotConverged) /(2*(qq+1));
        % Update starting point for cqqnk using Eq. B.18
        cqqnk = cqqnk / (2*qq+1-2*n);
        qq = qq+1;
    end

    % estimated loss of precision
    lossPrecS(kInd, :) = abs(maxTermS(kInd,:)./abs(S(kInd,:)));
    % Finally, we scale everything with the missing factor (see Eq. B.20)
    S(kInd, :) = -s^(kk+1)*S(kInd,:);
end

end

function u = sphGetUforFp(n)
% Calculates the matrix u for series Bessel implementation
% sphGetUforFp(n) calculates the matrix, for one value of n, for the series
% implementation of the Bessel function product. This is valid for all k up
% to k=n.
%
% The matrix u is defined as (Eq. B.10 of JQSRT 2013]
%   u_{rb} = 2^b (d/dX)^b [X^(n-1/2)*(1-1/X)^r]|X=1
% which leads to the relations [Eq. B.16 of JQSRT 2013]
%   u_{0,0} = 0
%   u_{r,0} = 1 (r > 0)
%   u_{r,b+1} = (n-1/2-2r)u_{rb} - (n-1/2-r)u_{r+1,b}
%
% Input:
%        n: [1 x 1] The value of n required
%
% Output:
%        u: [B+2 x B+1] where B = floor(n/2), the matrix u as defined above

	bMax = floor(n/2);

	u = zeros(bMax+2, bMax+1); % [r x b], both from 0 to bMax
	u(1,1) = 1;

	for bInd = 1:(bMax)
		%bb = bInd - 1;
		%do r == 1 case first
		rInd = 1;
		rr = 0;
		u(rInd, bInd+1) = (2*n-1-4*rr)*u(rInd, bInd) - (2*n-1-2*rr)*u(rInd+1, bInd);
		% then do other r
		for rInd = 2:bInd+1
			rr = rInd - 1;
			u(rInd, bInd+1) = (2*n-1-4*rr)*u(rInd, bInd) - (2*n-1-2*rr)*u(rInd+1, bInd) + 2*rr*u(rInd-1, bInd);
		end
	end
end
