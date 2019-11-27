function [A,B,C] = translate_z(nmax,z, varargin)
% Calculates translation matrices for translation of VSWFs along z axis.
%
% Usage
%   [A,B] = translate_z(nmax,z) calculates translation matrices.
%   The matrices are use as:
%
%   .. math::
%      M' = A M + B N
%
%      N' = B M + A N
%
%   [A,B,C] = ott.utils.translate_z(nmax,z) additionally, calculates C,
%   the scalar SWF translation coefficients in 3d packed form.
%
% Parameters
%   - nmax (int) -- Determines the number of multipole terms to include
%     in the translation matrices (multipole order).  Can be a single
%     integer or two integers for the ``[row, column]`` nmax.
%     If the row, column indices don't match, A and B will not be square.
%   - z (numeric) -- Translation distance.
%
% A and B are sparse matrices, since only m' = m VSWFs couple
%
% If z is a vector/matrix only A's and B's will be outputted. A and B will
% be cells of matrices the length of the number of elements in z. To save
% time only use unique values of z.
%
% Time *may* be saved by taking the conjugate transpose instead of
% calculating translations in the positive or negative direction.
%
% Optional named parameters
%   - 'type' (enum)   --  Type of translation matrix to generate.
%   - 'method' (enum) --  Method to calculate translation matrix.
%
% Translation matrix types
%   - 'sbesselj'            regular to regular.  Default.  Used for
%     most particle or beam translations.
%   - 'sbesselh1'           outgoing to regular.  Can be useful for
%     doing multiple scattering calculations.
%   - 'sbesselh2'           incoming to regular
%   - 'sbesselh1farfield'   outgoing to regular far-field limit
%   - 'sbesselh2farfield'   incoming to regular far-field limit
%
% Methods
%   - 'gumerov'     --  Use Gumerov method (default, recommended)
%   - 'videen'      --  Use Videen method.  (not recommended for
%     large translations of the beam, unstable)
%
% Example usage
%   The following example calculates the A and B translation matrices
%   and applies them to a Gaussian beam.  A procedure similar to this
%   is done when calling ``Bsc.translateZ`` with a distance::
%
%     beam = ott.BscPmGauss('NA', 0.9, 'index_medium', 1.0, ...
%         'polarisation', [1, 0], 'wavelength0', 1064e-9);
%
%     z = 1.0e-6 ./ beam.wavelength;
%     [A, B] = translate_z([beam.Nmax, beam.Nmax], z);
%
%     % Apply translation matrices
%     new_beam = beam.translate(A, B);


% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Refs:
% N. A. Gumerov and R. Duraiswami, "Recursions for the computation of
% multipole translation and rotation coefficients for the 3-D Helmholtz
% equation", SIAM J. Sci. Comput. 25(4), 1344-1381 (2003)
%
% G. Videen, "Light scattering from a sphere near a plane interface",
% chapter 5 (pp. 81-96) % in F. Moreno and F. Gonzalez (eds), "Light
% Scattering from Microstructures", Lecture Notes in Physics 534,
% Springer-Verlag, Berlin, 2000

import ott.utils.*
ott.warning('internal');

p = inputParser;
p.addParameter('type', 'sbesselj');
p.addParameter('method', 'gumerov');
p.parse(varargin{:});

if numel(z)>1
    A=cell(numel(z),1);
    B=A;
    for ii=1:numel(z)
        [A{ii},B{ii}]=translate_z(nmax,z(ii), varargin{:});
    end
    C=0;
    ott.warning('external');
    return
end

% Calculate Nmax for each dimension
if numel(nmax) == 1
  nmax1 = nmax;
  nmax2 = nmax;
else
  nmax1 = nmax(1);
  nmax2 = nmax(2);
  nmax = max(nmax(1:2));
end

if z==0
    A=speye(nmax1^2+nmax1*2, nmax2^2+nmax2*2);
    B=sparse(nmax1^2+nmax1*2, nmax2^2+nmax2*2);
    C=A;
    ott.warning('external');
    return
end

% Calculate the scalar coefficients
switch p.Results.method
  case 'videen'
    C = translate_z_videen(nmax1, nmax2, nmax, abs(z), p);
  case 'gumerov'
    C = translate_z_gumerov(nmax1, nmax2, nmax, abs(z), p);
  otherwise
    error('OTT:UTILS:translate_z:unknown_method', 'Unknown method');
end

[A,B] = calculate_AB(C, nmax1, nmax2, nmax, z, p);

if nargout>2
    C=C(1:nmax+1,1:nmax+1,1:min(nmax1, nmax2)+1);
end

ott.warning('external');

end

function C = translate_z_videen(nmax1, nmax2, nmax, z, p)

import ott.utils.*

N = 3*nmax+5;
N3 = min(nmax1, nmax2) + 1;

C = zeros(N,N,N3);

% First calculate the scalar translation coeffs

% Starting values, for m=0 and k=any -> n=0
% Videen (38)
k = 0:(N-1);

switch p.Results.type
  case 'sbesselj'       % regular to regular
    if z < 0
      C(:,1,1) = sqrt(2*k+1) .* sbesselj(k,2*pi*abs(z)) .* (-1).^(k);
    else
      C(:,1,1) = sqrt(2*k+1) .* sbesselj(k,2*pi*z);
    end

  case 'sbesselh1'      % outgoing to regular
    if z < 0
      C(:,1,1) = sqrt(2*k+1) .* sbesselh1(k,2*pi*abs(z)) .* (-1).^(k);
    else
      C(:,1,1) = sqrt(2*k+1) .* sbesselh1(k,2*pi*z);
    end

  case 'sbesselh2'      % incoming to regular
    if z < 0
      C(:,1,1) = sqrt(2*k+1) .* sbesselh2(k,2*pi*abs(z)) .* (-1).^(k);
    else
      C(:,1,1) = sqrt(2*k+1) .* sbesselh2(k,2*pi*z);
    end

  case 'sbesselh1farfield'    % outgoing to regular

    if 2*pi*abs(z) <= (N-1).^2
      ott.warning('Farfield limit may not be satisfied');
    end

    h1limit = (-1i).^k./(1i*2*pi*abs(z)) .* exp(1i*2*pi*abs(z));

    if z < 0
      C(:,1,1) = sqrt(2*k+1) .* h1limit .* (-1).^(k);
    else
      C(:,1,1) = sqrt(2*k+1) .* h1limit;
    end

  case 'sbesselh2farfield'    % incoming to regular

    if 2*pi*abs(z) <= (N-1).^2
      ott.warning('Farfield limit may not be satisfied');
    end

    h2limit = (1i).^k./(-1i*2*pi*abs(z)) .* exp(-1i*2*pi*abs(z));

    if z < 0
      C(:,1,1) = sqrt(2*k+1) .* h2limit .* (-1).^(k);
    else
      C(:,1,1) = sqrt(2*k+1) .* h2limit;
    end

  otherwise
    error('OTT:UTILS:translate_z:type_error', 'Unknown translation type');
end

% Do n=1 as a special case (Videen (40) with n=0,n'=k)
kk = 1:(N-2);
kk = kk(:);
Cm = C(kk,1,1);
Cm = Cm(:);
Cp = C(kk+2,1,1);
Cp = Cp(:);
C(1,2,1) = -C(2,1,1);
C(kk+1,2,1) = sqrt(3./(2*kk+1)) .* ...
    ( kk.*sqrt(1./(2*kk-1)) .* Cm - (kk+1).*sqrt(1./(2*kk+3)) .* Cp );

% Now do the rest, up to n=N-1
% Videen (40), with n(Videen) = n-1, n' = k
% Note that only the k=0 term is needed for n=N-1
for n = 2:(N-2)
    kk = 1:(N-n-1);
    kk = kk(:);
    Cm = C(kk,n,1);
    Cm = Cm(:);
    Cp = C(kk+2,n,1);
    Cp = Cp(:);
    C0 = C(kk+1,n-1,1);
    C0 = C0(:);
    C(1,n+1,1) = (-1)^n * C(n+1,1,1);
    C(kk+1,n+1,1) = sqrt((2*n+1)./(2*kk+1))/n .* ...
        ( kk.*sqrt((2*n-1)./(2*kk-1)) .* Cm ...
        + (n-1)*sqrt((2*kk+1)/(2*n-3)) .* C0 ...
        - (kk+1).*sqrt((2*n-1)./(2*kk+3)) .* Cp );
end
n = N-1;
C(1,N,1) = sqrt(2*n+1)/n * ...
    ( (n-1)*sqrt(1/(2*n-3)) * C(1,n-1,1) - sqrt((2*n-1)/3) * C(2,n,1) );

% OK, now m other than m=0
% Only need to do positive m, since C(-m) = C(m)
% Videen (41)
for m = 1:min(nmax1, nmax2)
  nn = m:nmax1;
  kk = (m:N-2).';
  C0 = C(kk+1,nn+1,m);
  Cp = C(kk+2,nn+1,m);
  Cm = C(kk,nn+1,m);
  C(kk+1,nn+1,m+1) = sqrt(1./((2*kk+1)*((nn-m+1).*(nn+m)))) .* ...
      ( sqrt(((kk-m+1).*(kk+m).*(2*kk+1))).*C0 ...
      -2*pi*z*sqrt((((kk-m+2).*(kk-m+1)))./((2*kk+3))).*Cp ...
      -2*pi*z*sqrt((((kk+m).*(kk+m-1)))./((2*kk-1))).*Cm );
end

end % translate_z_videen

function [A,B] = calculate_AB(C, nmax1, nmax2, nmax, z, p)

import ott.utils.*

% OK, that's the scalar coefficients
% Time to find the vector coefficients - Videen (43) & (44)

nn = 1:nmax1;
kk = (1:nmax2).';

matrixm=sqrt(kk.*(kk+1)) ./ sqrt(nn.*(nn+1));

central_iterator1=[1:nmax1].*[2:nmax1+1];
central_iterator2=[1:nmax2].*[2:nmax2+1];

[ciy,cix]=meshgrid(central_iterator1,central_iterator2);

mmm=0;

C0 = C(2:(nmax2+1),2:(nmax1+1),mmm+1);
Cp = C(3:(nmax2+2),2:(nmax1+1),mmm+1);
Cm = C(1:nmax2,2:(nmax1+1),mmm+1);

t = matrixm.*(C0 - 2*pi*abs(z)./(kk+1) .* ...
    sqrt((kk-mmm+1).*(kk+mmm+1)./((2*kk+1).*(2*kk+3))) .* Cp - ...
    2*pi*abs(z)./kk.*sqrt((kk-mmm).*(kk+mmm)./((2*kk+1).*(2*kk-1))).*Cm);

toIndexy=(ciy(:));
toIndexx=(cix(:));
A=t(:);
B=zeros(size(A));

% Total size of A and B: sum((1:nmax).^2)*2 + nmax^2

for mmm=1:min(nmax1, nmax2)

    sz1 = mmm:nmax2;
    sz2 = mmm:nmax1;

    C0 = C((1+mmm):(nmax2+1),(1+mmm):(nmax1+1),mmm+1);
    Cp = C((2+mmm):(nmax2+2),(1+mmm):(nmax1+1),mmm+1);
    Cm = C((mmm):nmax2,(1+mmm):(nmax1+1),mmm+1);

    tt = matrixm(sz1, sz2).*(C0 - 2*pi*abs(z)./(kk(sz1)+1) .* ...
        sqrt((kk(sz1)-mmm+1).*(kk(sz1)+mmm+1) ...
        ./((2*kk(sz1)+1).*(2*kk(sz1)+3))) .* Cp - ...
        2*pi*abs(z)./kk(sz1).*sqrt((kk(sz1)-mmm) ...
        .*(kk(sz1)+mmm)./((2*kk(sz1)+1).*(2*kk(sz1)-1))).*Cm);

    ciys=ciy(mmm:end,mmm:end);
    cixs=cix(mmm:end,mmm:end);

    toIndexy=[toIndexy;(ciys(:)+mmm);(ciys(:)-mmm)];
    toIndexx=[toIndexx;(cixs(:)+mmm);(cixs(:)-mmm)];
    A=[A;tt(:);tt(:)];

    tt = mmm./(kk(sz1).*(kk(sz1)+1)).*matrixm(sz1, sz2) .* C0;
    B=[B;tt(:);-tt(:)];

end

% Keep B real until the end, makes things run faster
B = 1i*2*pi*abs(z)*B;

% This is faster than A = A + sparse(...) and A(sub2ind(...)) = [...]
if z < 0
  [n1, ~] = ott.utils.combined_index(toIndexy);
  [n2, ~] = ott.utils.combined_index(toIndexx);
  B=sparse(toIndexy,toIndexx,B.*(-1).^(n1-n2+1),nmax1*(nmax1+2),nmax2*(nmax2+2));
  A=sparse(toIndexy,toIndexx,A.*(-1).^(n1-n2),nmax1*(nmax1+2),nmax2*(nmax2+2));
else
  B=sparse(toIndexy,toIndexx,B,nmax1*(nmax1+2),nmax2*(nmax2+2));
  A=sparse(toIndexy,toIndexx,A,nmax1*(nmax1+2),nmax2*(nmax2+2));
end

end % calculate_AB

function C = translate_z_gumerov(nmax1, nmax2, nmax, r, p)
% Optimised implementation from Alex

import ott.utils.*

% Having pre-computed a_nm it's fast?
% nmax=3;
% r=0.1;

%some "setup" values:
mmax=min(nmax1,nmax2);
m=0;
fval=2*nmax+1;
nd=[m:fval];

kr=2*pi*r;

%compute seed functions:
switch p.Results.type
  case 'sbesselj'       % regular to regular
    C_nd00=[sqrt(2*nd+1).*sbesselj(nd,kr)];

  case 'sbesselh1'      % outgoing to regular
    C_nd00=[sqrt(2*nd+1).*sbesselh1(nd,kr)]./2;

  case 'sbesselh2'      % incoming to regular
    C_nd00=[sqrt(2*nd+1).*sbesselh2(nd,kr)]./2;

  otherwise
    error('OTT:UTILS:translate_z:type_error', 'Unknown translation type');
end

C_ndn0=zeros(length(nd)+1,length(nd)+1);
C_ndn0(1+[1:length(C_nd00)],2)=C_nd00;
C_ndn0(2,1+[1:length(C_nd00)])=((-1).^(nd).*C_nd00).';

%gumerov's zonal coefficients are m=0. Compute columns, limited by diagonal:
%compute lower diagonal first:
for jj=1:nmax
    ii=[jj:fval-jj].';
    C_ndn0(ii+2,ii(1)+2)=(anm_l(ii(1)-2,0).*C_ndn0(ii+2,ii(1))-anm_l(ii,0).*C_ndn0(ii+3,ii(1)+1)+anm_l(ii-1,0).*C_ndn0(ii+1,ii(1)+1))./anm_l(ii(1)-1,0);
    C_ndn0(ii(1)+2,ii+2)=((-1).^(jj+ii).*C_ndn0(ii+2,ii(1)+2)).';
end

%create "C":
C=zeros(nmax2+2,nmax1+1,mmax+1);
C(:,:,1)=C_ndn0(2:(nmax2+3),2:(nmax1+2));

%Having computed anm for m=0; cases we now can compute anm for all
%remaining cases:
ANM=anm_l([0:2*nmax+1].',[1:nmax]);
IANM=1./ANM;
for m=1:mmax

    %having computed the zonal coefficients we now compute the "diagonal ones"
    %(tesseral)
    %i.e. ones which generate m on the first column we then reproduce the same
    %commputation for the n nd recursion:

    nd=[m:fval-m].';
    C_nd1m=(bnm_l(nd,-m).*C_ndn0(nd+1,m+1)-bnm_l(nd+1,m-1).*C_ndn0(nd+3,m+1))./bnm_l(m,(-m));

    %having computed the first seed column we now recur the elements:
    C_ndn1=zeros(size(C_ndn0)); %make zero as we re-use
    C_ndn1([1:length(C_nd1m)]+m+1,m+2)=C_nd1m;
    C_ndn1(m+2,[1:length(C_nd1m)]+m+1)=((-1).^(nd+m).*C_nd1m).';

    for jj=m+1:nmax
        ii=[jj:fval-jj].';
%         C_ndn1(ii+2,ii(1)+2)=(anm(ii(1)-2,m).*C_ndn1(ii+2,ii(1))-anm(ii,m).*C_ndn1(ii+3,ii(1)+1)+anm(ii-1,m).*C_ndn1(ii+1,ii(1)+1))./anm(ii(1)-1,m);
        C_ndn1(ii+2,ii(1)+2)=(ANM(ii(1)-1,m).*C_ndn1(ii+2,ii(1))-ANM(ii+1,m).*C_ndn1(ii+3,ii(1)+1)+ANM(ii,m).*C_ndn1(ii+1,ii(1)+1)).*IANM(ii(1),m);
        C_ndn1(ii(1)+2,ii+2)=((-1).^(jj+ii).*C_ndn1(ii+2,ii(1)+2)).';
    end
    C_ndn0=C_ndn1;

    C(:,:,m+1)=C_ndn0(2:(nmax2+3),2:(nmax1+2));

end

end % translate_z_gumerov

function a_nm = anm_l(n,m);
% For translate_z_gumerov
fn=1./(2*n+1)./(2*n+3);
a_nm=sqrt((n+abs(m)+1).*(n-abs(m)+1).*fn);
a_nm(n<0)=0;
a_nm(abs(m)>n)=0;
end

function b_nm = bnm_l(n,m);
% For translate_z_gumerov
b_nm=(2*(m<0)-1).*sqrt((n-m-1).*(n-m)./(2*n-1)./(2*n+1));
b_nm(abs(m)>n)=0;
end

