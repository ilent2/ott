function stPinmTaunm = vshPinmTaunm(nMax, theta)
%% vshPinmTaunm
% Calculates angular functions pi and tau as defined in Mishchenko 2002
% 
%	vshPinmTaunm(nMax, theta) computes angle functions pi_nm(theta)
% and tau_nm(theta) for n=1..nMax, |m|<=n, using recurrence relations.
% Those functions are defined on p. 373 of Mishchenko 2002.
% It also returns the Legendre polynomials d_{n0}(theta)=P_n(cos(theta)),
% which are necessary for m=0.
%
% Input:
% - nMax: scalar integer [1 x 1]
% - theta: column vector [T x 1]
%          with theta (in radians)
%          all theta's must be between 0 and pi
%
% Output: structure with fields:
% - pinm: matrix [T x P] with pi_nm(theta)
% - taunm: matrix [T x P] with tau_nm(theta)
% - pn0: matrix [T x nMax] with d_n0(theta)=P_n(cos(theta))
% where P = nMax*(nMax + 2)
% The arrays pinm and taunm use the p-index which stores the possible values
% of(n,m) in a linear array using the following convention p = n*(n+1)/2.
%
% Dependency: 
% none

if size(theta,2)>1
    disp 'Warning: theta must be a column vector in vshPinTaunm...';
end
if ~isempty(find(theta<0 , 1))
    disp 'Warning: theta must be >=0 in vshPinTaunm...';
end

nrows=length(theta);
P = nMax*(nMax + 2); % maximum number of columns for pi_n and tau_n matrices
stPinmTaunm.pinm = zeros(nrows,P);
stPinmTaunm.taunm = zeros(nrows,P); % initialize both matrices to zero

% Initialize Am for m=0
Amsinmm1 = ones(nrows,1); % [T x 1]

muc=cos(theta); % [T x 1]
mus=sin(theta); % [T x 1]

% m=0 case is treated separately
% loop on m
for m=1:nMax

    % Am * sin(theta)^(m-1) is computed by recurrence on m
    if m>1
        Amsinmm1 = Amsinmm1*sqrt(((2*m-1)/(2*m))).*mus;
    else
        Amsinmm1 = Amsinmm1*sqrt(((2*m-1)/(2*m)));
    end

    % Initialize recurrence pi_{m-1,m}=0, pi_{m,m}=m*Am*sin(theta)^(m-1)
    ncols = nMax - m + 2;
    piaux = zeros(nrows,ncols);
    piaux(:,2) = m*Amsinmm1;
    % piaux contains pi_{m+j-2,m} j=1..(nMax-m+2),
    % i.e. pi_{m-1,m}..pi_{nMax,m}

    % Get pi_{m+1,m} to pi_nMax by recurrence
    for jj=3:ncols
        n = m + jj -2;
        % piaux(:,jj) is pi_{m+jj-2,m}
        piaux(:,jj) = (1/sqrt((n-m)*(n+m)))*((2*n-1)*muc.*piaux(:,jj-1) - sqrt((n-1-m)*(n-1+m)).*piaux(:,jj-2));
    end

    nvec = (m:nMax); % [1 x N2] with N2=ncols-1

    pvec = nvec.*(nvec + 1) + m;   % computes p for positive m
    pvecn = pvec - 2*m;  % computes p for negative m

    % fill in pi_nm matrix for positive or negative m
    stPinmTaunm.pinm(:,pvec) = piaux(:,2:(ncols));
    stPinmTaunm.pinm(:,pvecn)= (-1)^(m+1).*piaux(:,2:(ncols));

    % return tau_nm matrix for positive or negative m
     stPinmTaunm.taunm(:,pvec) = bsxfun(@times,piaux(:,1:(ncols-1)),-sqrt((nvec-m).*(nvec+m))/m) + (muc * (nvec/m)).*piaux(:,2:(ncols));
     stPinmTaunm.taunm(:,pvecn) = (-1)^m*stPinmTaunm.taunm(:,pvec);
end;

% Now do m=0 case
% Initialize recurrence p_0=1, p_1=muc, t_0=0, t_1=-mus
% pnm1 contains p_{n-1} n=1..nMax+1, same for tnm1
pnm1=ones(nrows,nMax+1); % [T x N+1]
pnm1(:,2)=muc; % [T x 1]
tnm1=zeros(nrows,nMax+1); % [T x N+1]
tnm1(:,2)=-mus; % [T x 1]

% Get p_2 to p_nMax and t_2 to t_nMax by recurrence
% p_n is pnm1(:,n+1), t_n is tnm1(:,n+1)
for n=2:(nMax)
    pnm1(:,n+1)=(2*n-1)/(n)*muc.*pnm1(:,n)-(n-1)/n*pnm1(:,n-1);
    tnm1(:,n+1)=muc.*tnm1(:,n)-n*mus.*pnm1(:,n);
end;

% return p_n matrix (except n=0)
stPinmTaunm.pn0=pnm1(:,2:(nMax+1));
% fill in taunm for m=0 (pinm=0 for m=0 is already set)
nvec=1:nMax;
pvec=nvec.*(nvec+1);
stPinmTaunm.taunm(:,pvec)=tnm1(:,2:(nMax+1));
