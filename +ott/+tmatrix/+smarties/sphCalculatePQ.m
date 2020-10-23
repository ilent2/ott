function CstPQa = sphCalculatePQ(nNmax, absmvec, stRtfunc, stParams, NB)
  %% sphCalculatePQ
% Calculates P,Q matrices for a spheroid using the algorithm of [JQSRT 123, 153 (2013)]
% 
% sphCalculatePQ(nNmax, absmvec, stRtfunc, params, NB) calculates the P
% and Q matrices for a spheroid for all m in absmvec.
%
% stRtfunc should contain the geometry, as obtained from sphMakeGeometry,
% and stParams should have fields k1 and s (both scalar). The theta range for
% the geometry should be [0;pi/2] or [0;pi] (either is acceptable)
%
% This makes use of modified Riccati-Bessel function products that are
% optimised for spheroids, i.e they do not contain terms that integrate to
% zero, thus removing the loss of precision which affects traditional
% codes. Details of the algorithm can be found in
% J. Quant. Spectrosc. Radiat. Transfer 123, 153-168 (2013).
%
%
% Input:
%       nNmax:      [1 x 1] The maximum multipole order to return
%       absmvec:    [1 x M] The values of m that the calculation should be carried
%                   out for. The order should be [0 1 2 ... mMax], omitting
%                   those entries not required
%       stRtfunc:   Structure containing geometric information, as from
%                   sphMakeGeometry
%       stParams:   Structure containing simulation parameters, as from
%                   tmsMakeParams or tmsMakeParamsLambda
%       NB:         (optional - default is nNmax)
%                   The number of multipoles that should be used to
%                   calculate the Bessel function products (NB >= N)
%
% Output:
%	This returns a cell {1 x M} with one structure for each m, each containing
%	structures for P (st4MPeo, st4MPoe) and for Q (st4MQeo, st4MQoe).
%	Each of these contains the blocks of the matrix, M11, M12, M21, M22, as
%	well as a field m which contains the value of m for those matrices, and
%	vectors ind1 and ind2 indicating which row and column indices are contained in
%	each block. Each cell also contains
%	a cell CsMatList, which is a list {'st4MQ', st4MP'}, showing that P and
%	Q are contained in it.
%
% Dependency: 
% isOctave, sphGetModifiedBesselProducts, vshPinmTaunm

import ott.tmatrix.smarties.*;

if nargin < 5 || NB < nNmax
    NB = nNmax;
end

% ignore divide-by-zero warnings in Octave
if isOctave()
    currentWarning = warning('off','Octave:divide-by-zero');
else
    currentWarning = warning();
end

M = length(absmvec); % number of m values

bOutput=true; % default
if isfield(stParams,'bOutput');
    bOutput = stParams.bOutput;
end

if bOutput
    disp(['sphCalculatePQ: Calculating P,Q for ', int2str(M), ' m-values with N_Q = ', ...
        int2str(nNmax), ', N_B = ', int2str(NB), ', N_Theta = ', ...
        int2str(stRtfunc.nNbTheta)]);
end

CstPQa = cell(1,M); % initializes cell to be returned

% the following are used for all m
k1 = stParams.k1;
s = stParams.s;
nNint = stRtfunc.nNbTheta;
T = nNint; % number of theta's
x = stRtfunc.r * k1; % [T x 1] x(theta)
xTheta = stRtfunc.drdt * k1; % [T x 1] derivative of x(theta)

% Angular functions, pi and tau
stPinmTaunm=vshPinmTaunm(nNmax, stRtfunc.theta);

sint = sin(stRtfunc.theta); % [T x 1]

dxdtwt = xTheta .* stRtfunc.wTheta; % [T x 1]

% An (for prefactors)
Anvec = sqrt((2*(1:nNmax)+1)./(2*(1:nNmax).*((1:nNmax)+1))); % [1 x N]

% AnAk matrix (for prefactors)
AnAk = Anvec.' * Anvec; % Matrix product gives [N x N] matrix

% The hard part is to get the modified (non-problematic) radial functions
% Note that these do not depend on m
% The following function uses the algorithm of [JQSRT 2013] to achieve this
[stXipsiAll, stPsipsiAll]=sphGetModifiedBesselProducts(nNmax, s, x, NB);
% These are [N x K=N x T] arrays

% We will use them as [N x T x K] for a given k, so we can speed up
% computations by reshaping the arrays now
xiprimepsi = zeros(nNmax,T,nNmax);
xipsiprime = zeros(nNmax,T,nNmax);
xipsi = zeros(nNmax,T,nNmax);
xiprimepsiprimePluskkp1xipsiOversxx = zeros(nNmax,T,nNmax);
xiprimepsiprimePlusnnp1xipsiOversxx = zeros(nNmax,T,nNmax);
for kk=1:nNmax
    xipsi(:,:,kk) = stXipsiAll.xipsi(:, kk, :);
    xiprimepsi(:,:,kk) = stXipsiAll.xiprimepsi(:, kk, :);
    xipsiprime(:,:,kk) = stXipsiAll.xipsiprime(:, kk, :);
    xiprimepsiprimePluskkp1xipsiOversxx(:,:,kk) = stXipsiAll.xiprimepsiprimePluskkp1xipsiOversxx(:, kk, :);
    xiprimepsiprimePlusnnp1xipsiOversxx(:,:,kk) = stXipsiAll.xiprimepsiprimePlusnnp1xipsiOversxx(:, kk, :);
    % We will also need those for the diagonals [N x T]
    %     xiprimepsiDiag(kk,:) = stXipsiAll.xiprimepsi(kk,kk,:);
    %     xipsiprimeDiag(kk,:) = stXipsiAll.xipsiprime(kk,kk,:);
    %     xipsisxxDiag(kk,:) =  stXipsiAll.xipsiOversxx(kk,kk,:);
end
forQdiagLt1 = stXipsiAll.fordiagLt1;
forQdiagLt2 = stXipsiAll.fordiagLt2;
forQdiagLt3 = stXipsiAll.fordiagLt3;
clear stXipsiAll
% Same for psipsi
psiprimepsi = zeros(nNmax,T,nNmax);
psipsiprime = zeros(nNmax,T,nNmax);
psipsi = zeros(nNmax,T,nNmax);
psiprimepsiprimePluskkp1psipsiOversxx = zeros(nNmax,T,nNmax);
psiprimepsiprimePlusnnp1psipsiOversxx = zeros(nNmax,T,nNmax);
% psiprimepsiDiag = zeros(nNmax,T);
% psipsiprimeDiag = zeros(nNmax,T);
% psipsisxxDiag = zeros(nNmax,T);
for kk=1:nNmax
    psipsi(:,:,kk) = stPsipsiAll.xipsi(:, kk, :);
    psiprimepsi(:,:,kk) = stPsipsiAll.xiprimepsi(:, kk, :);
    psipsiprime(:,:,kk) = stPsipsiAll.xipsiprime(:, kk, :);
    psiprimepsiprimePluskkp1psipsiOversxx(:,:,kk) = stPsipsiAll.xiprimepsiprimePluskkp1xipsiOversxx(:, kk, :);
    psiprimepsiprimePlusnnp1psipsiOversxx(:,:,kk) = stPsipsiAll.xiprimepsiprimePlusnnp1xipsiOversxx(:, kk, :);
    % We will also need those for the diagonals [N x T]
    %     psiprimepsiDiag(kk,:) = stPsipsiAll.xiprimepsi(kk,kk,:);
    %     psipsiprimeDiag(kk,:) = stPsipsiAll.xipsiprime(kk,kk,:);
    %     psipsisxxDiag(kk,:) =  stPsipsiAll.xipsiOversxx(kk,kk,:);
end
% This one is also needed for the diagonal of P11 [N x T] (Eq. 65)
forPdiagLt1 = stPsipsiAll.fordiagLt1;
forPdiagLt2 = stPsipsiAll.fordiagLt2;
forPdiagLt3 = stPsipsiAll.fordiagLt3;
clear stPsipsiAll

% The rest of the calculations are m-dependent, so we loop on m
for mInd = 1:M

    m = absmvec(mInd);
    nMin = max(m,1);
    Nm = nNmax - nMin + 1;% size of the matrices for a given m (since n,k>=m)
    nvec = (nMin:nNmax); % [1 x Nm] a list of n values to use
    pvec = nvec.*(nvec + 1) + m; % [1 x Nm] p-index for positive m

    % angular functions
    pinm = transpose(stPinmTaunm.pinm(:, pvec)); % [Nm x T]
    taunm = transpose(stPinmTaunm.taunm(:, pvec)); % [Nm x T]
    if m==0
        dn = transpose(stPinmTaunm.pn0(:, nvec)); % [Nm x T]
    else
        dn = bsxfun(@times,pinm,transpose(sint)/m); % [Nm x T]
    end

    dntimesnnp1= bsxfun(@times, dn, transpose(nvec.*(nvec+1))); % [Nm x T]

    K1 = zeros(Nm);
    K2 = zeros(Nm);
    L5 = zeros(Nm);
    L6= zeros(Nm);

    K1P = zeros(Nm);
    K2P = zeros(Nm);
    L5P = zeros(Nm);
    L6P = zeros(Nm);

    % Note that the following code does not take into account the zeros in half of the
    % terms in the matrices and hence calculates twice as many terms as
    % needed. However, because of the way Matlab handles matrices, it is as
    % fast as a more optimized code (the time gained in calculation is
    % offset by the extra array lookups, which are slow in Matlab).

    % The integrals are calculated as sums using Gaussian quadratures
    % This is carried out for a given index k for all n (i.e. for a given column)
    % by performing a matrix multiplication of a N x T matrix by a T x 1
    % vector.
    % We therefore loop on k
    for kk = nMin:nNmax
        kInd = kk -nMin + 1; % this is the index for k
        dk = transpose(dn(kInd, :)); % d_k(theta), [T x 1]
        tauk = transpose(taunm(kInd, :)); % tau_k(theta), [T x 1]

        % These [T x 1] vectors will be needed
        dxdttauksint	= dxdtwt .* tauk; % [T x 1]
        dxdtdksint		= dxdtwt .* dk; % [T x 1]

        % For matrices K1 and K2, we use Eqs. 11/53 and 12/54 of JQSRT2013
        pinxiprimepsi = pinm .* xiprimepsi(nMin:end, :, kk); % [Nm x T]
        pinxipsiprime = pinm .* xipsiprime(nMin:end, :, kk); % [Nm x T]
        % The integrals are carried out as sums over theta's by taking a matrix product
        % of a [Nm x T] matrix by a [T x 1] vector
        K1(:, kInd) = pinxipsiprime * dxdtdksint; % [Nm x 1] Eq. 11/53
        K2(:, kInd) = pinxiprimepsi * dxdtdksint; % [Nm x 1] Eq. 12/54
        % These are for Q, we also do the same for P. For efficiency, we use
        % the same variable names even if xi is in fact replaced by psi
        pinxiprimepsi = pinm .* psiprimepsi(nMin:end, :, kk); % [Nm x T]
        pinxipsiprime = pinm .* psipsiprime(nMin:end, :, kk); % [Nm x T]
        K1P(:, kInd) = pinxipsiprime * dxdtdksint; % [Nm x 1]
        K2P(:, kInd) = pinxiprimepsi * dxdtdksint; % [Nm x 1]

        % For L5, we use Eqs. 18/52
        dnxipsinnp1		= dntimesnnp1 .* xipsi(nMin:end, :, kk); % [Nm x T]
        taunxipsi		= taunm.*xipsi(nMin:end, :, kk);% [Nm x T]

        L5(:, kInd)		= dnxipsinnp1 * dxdttauksint - taunxipsi * dxdtdksint*kk*(kk+1);

        % Same for P
        dnxipsinnp1     = dntimesnnp1 .* psipsi(nMin:end, :,kk); % [Nm x T]
        taunxipsi		= taunm.* psipsi(nMin:end, :, kk); % [Nm x T]

        L5P(:, kInd)	= dnxipsinnp1 * dxdttauksint - taunxipsi * dxdtdksint*kk*(kk+1);

        % For L6, we use Eqs. 20-22/55-56
        dnxiprimepsiprimennp1plusxipsioversxxnnp1kkp1 =...
            dntimesnnp1 .* xiprimepsiprimePluskkp1xipsiOversxx(nMin:end, :, kk); % [Nm x T]
        taunxiprimepsiprimeplusxipsioversxxnnp1 = ...
            taunm.*xiprimepsiprimePlusnnp1xipsiOversxx(nMin:end, :, kk); % [Nm x T]

        L6(:, kInd) = dnxiprimepsiprimennp1plusxipsioversxxnnp1kkp1 * dxdttauksint - ...
                      taunxiprimepsiprimeplusxipsioversxxnnp1 * dxdtdksint*kk*(kk+1);

        % Same for P
        dnxiprimepsiprimennp1plusxipsioversxxnnp1kkp1 =...
            dntimesnnp1 .* psiprimepsiprimePluskkp1psipsiOversxx(nMin:end, :, kk); % [Nm x T]

        taunxiprimepsiprimeplusxipsioversxxnnp1 = ...
            taunm.*psiprimepsiprimePlusnnp1psipsiOversxx(nMin:end, :, kk); % [Nm x T]

        L6P(:, kInd) = dnxiprimepsiprimennp1plusxipsioversxxnnp1kkp1 * dxdttauksint - ...
                       taunxiprimepsiprimeplusxipsioversxxnnp1 * dxdtdksint*kk*(kk+1);

    end

    % Prefactors
    Prefactor1 = (s-1)*(s+1) / s * AnAk(nvec,nvec); % [Nm x Nm] (s^2-1)/s * An*Ak
    % Note that the result of the bsxfun below is [Nm x Nm] matrix of
    % n(n+1)-k(k+1)
    Prefactor2 = 1i * Prefactor1 ./ (bsxfun(@minus,transpose(nvec.*(nvec+1)),nvec.*(nvec+1)));

    % Get Q matrix with prefactors
    Q12 = Prefactor1 .* K1; % Eq. 15
    Q21 = - Prefactor1 .* K2; % Eq. 16
    Q11 = Prefactor2 .* L5; % Eq. 17
    Q22 = Prefactor2 .* L6; % Eq. 19
    % Get P matrix with prefactors
    P12 = Prefactor1 .* K1P; % Eq. 15
    P21 = - Prefactor1 .* K2P; % Eq. 16
    P11 = Prefactor2 .* L5P; % Eq. 17
    P22 = Prefactor2 .* L6P; % Eq. 19

    % To finish, we need to do the diagonal terms of Q11/P11 and Q22/P22
    % calculated differently to avoid problems

    % First get the diagonal terms indices (using linear indices rather
    % than matrix pairs)
    diagindices = (1:(Nm+1):(Nm^2));
    % this is a row [1 x Nm] whose indices correspond to the diagonal of
    % a [Nm x Nm] matrix ,
    % i.e. for A [Nm x Nm], A(diagindices) = diag(A).

    PrefactDiag1 = (-1i/s * (2*nvec+1)./ (2*nvec.*(nvec+1))); % [1 x Nm]
    PrefactDiag2 = ( -1i*(s-1)*(s+1)/s/2 *(2*nvec+1));
    pi2ptau2 = (pinm.^2 + taunm.^2); % [Nm x T]

    % Fill in the diagonals by calculating the integrals as matrix products
    % of [Nm x T] by [T x 1] as before

    % Diagonal of Q11, from Eqs 23, 25, 65 (all [1 x Nm])
    Q11(diagindices) = PrefactDiag1 .* ...
        transpose((pi2ptau2 .* forQdiagLt1(nvec,:)) * stRtfunc.wTheta);
    % Diagonal of Q22, from Eqs 24, 26, 27 (all [1 x Nm])
    Ltilde2 = transpose((pi2ptau2 .* forQdiagLt2(nvec,:)) * stRtfunc.wTheta);
    Ltilde3 = transpose((dn .* taunm .* forQdiagLt3(nvec,:))  * dxdtwt);
    Q22(diagindices) = PrefactDiag1 .* Ltilde2 + PrefactDiag2 .* Ltilde3;
    % Diagonal of P11, from Eqs 23, 25, 65 (all [1 x Nm])
    P11(diagindices) = PrefactDiag1 .* ...
        transpose((pi2ptau2 .* forPdiagLt1(nvec,:)) * stRtfunc.wTheta);
    % Diagonal of P22, from Eqs 24, 26, 27 (all [1 x Nm])
    Ltilde2 = transpose((pi2ptau2 .* forPdiagLt2(nvec,:)) * stRtfunc.wTheta);
    Ltilde3 = transpose((dn .* taunm .* forPdiagLt3(nvec,:))  * dxdtwt);
    P22(diagindices) = PrefactDiag1 .* Ltilde2 + PrefactDiag2 .* Ltilde3;


    % Store the results in the cell for return

    % The code below uses the symmetry and the even-odd formulation
    evenodd=mod(nMin,2);
    % if nMin is even then the first index (1) is even so need to swap even-odd
    inde=(1+evenodd):2:Nm;
    indo=(2-evenodd):2:Nm;

    CstPQa{1, mInd}.st4MQeo.M12 = Q12(inde,indo);
    CstPQa{1, mInd}.st4MQeo.M21 = Q21(indo,inde);
    CstPQa{1, mInd}.st4MQeo.M11 = Q11(inde,inde);
    CstPQa{1, mInd}.st4MQeo.M22 = Q22(indo,indo);
    CstPQa{1, mInd}.st4MQeo.m = m;
    CstPQa{1, mInd}.st4MQeo.ind1 = inde;
    CstPQa{1, mInd}.st4MQeo.ind2 = indo;

    CstPQa{1, mInd}.st4MQoe.M12 = Q12(indo,inde);
    CstPQa{1, mInd}.st4MQoe.M21 = Q21(inde,indo);
    CstPQa{1, mInd}.st4MQoe.M11 = Q11(indo,indo);
    CstPQa{1, mInd}.st4MQoe.M22 = Q22(inde,inde);
    CstPQa{1, mInd}.st4MQoe.m = m;
    CstPQa{1, mInd}.st4MQoe.ind1 = indo;
    CstPQa{1, mInd}.st4MQoe.ind2 = inde;

    CstPQa{1, mInd}.st4MPeo.M12 = P12(inde,indo);
    CstPQa{1, mInd}.st4MPeo.M21 = P21(indo,inde);
    CstPQa{1, mInd}.st4MPeo.M11 = P11(inde,inde);
    CstPQa{1, mInd}.st4MPeo.M22 = P22(indo,indo);
    CstPQa{1, mInd}.st4MPeo.m = m;
    CstPQa{1, mInd}.st4MPeo.ind1 = inde;
    CstPQa{1, mInd}.st4MPeo.ind2 = indo;

    CstPQa{1, mInd}.st4MPoe.M12 = P12(indo,inde);
    CstPQa{1, mInd}.st4MPoe.M21 = P21(inde,indo);
    CstPQa{1, mInd}.st4MPoe.M11 = P11(indo,indo);
    CstPQa{1, mInd}.st4MPoe.M22 = P22(inde,inde);
    CstPQa{1, mInd}.st4MPoe.m = m;
    CstPQa{1, mInd}.st4MPoe.ind1 = indo;
    CstPQa{1, mInd}.st4MPoe.ind2 = inde;
    CstPQa{1,mInd}.CsMatList={'st4MQ', 'st4MP'};
end

warning(currentWarning); % reactivates divide-by-zero warning

end

