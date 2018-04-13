function [A,B,C] = translate_z(nmax,z)
% TRANSLATE_Z calculates translation matricies for translation of
% VSWFs along z axis.
%
% [A,B] = translate_z(nmax,z) calculates translation matricies.
% The matricies are use as: M' = A M + B N; N' = B M + A N.
%
% [A,B,C] = translate_z(nmax,z) additionally, calculates C,
% the scalar SWF translation coefficients in 3d packed form.
%
% A and B are sparse matrices, since only m' = m VSWFs couple
%
% If z is a vector/matrix only A's and B's will be outputted. A and B will
% be cells of matricies the length of the number of elements in z. To save
% time only use unique values of z.
%
% Time *may* be saved by taking the conjugate transpose instead of
% calculating translations in the positive or negative direction.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

ott_warning('ott:translate_z:move', ...
    'This file will move to ott.utils.translate_z');

% Refs:
% N. A. Gumerov and R. Duraiswami, "Recursions for the computation of
% multipole translation and rotation coefficients for the 3-D Helmholtz
% equation", SIAM J. Sci. Comput. 25(4), 1344-1381 (2003)
%
% G. Videen, "Light scattering from a sphere near a plane interface",
% chapter 5 (pp. 81-96) % in F. Moreno and F. Gonzalez (eds), "Light
% Scattering from Microstructures", Lecture Notes in Physics 534,
% Springer-Verlag, Berlin, 2000

ott_warning('internal');

if numel(z)>1
    A=cell(numel(z),1);
    B=A;
    for ii=1:numel(z)
        [A{ii},B{ii}]=translate_z(nmax,z(ii));
    end
    C=0;
    ott_warning('external');
    return
end

if z==0
    A=sparse(1:(nmax^2+nmax*2),1:(nmax^2+nmax*2),1);
    B=sparse((nmax^2+nmax*2),(nmax^2+nmax*2));
    C=A;
    ott_warning('external');
    return
end

N = 3*nmax+5;

C = zeros(N,N,N);

% First calculate the scalar translation coeffs

% Starting values, for m=0 and k=any -> n=0
% Videen (38)
k = 0:(N-1);
C(:,1,1) = sqrt(2*k+1) .* sbesselj(k,2*pi*z);
if z < 0
    C(:,1,1) = sqrt(2*k+1) .* sbesselj(k,2*pi*abs(z)) .* (-1).^(k);
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
for m = 1:(N-2)
    kk = m:(N-2);
    kk = kk(:);
    nn = m:(N-2);
    row1 = ones(size(nn));
    C0 = C(kk+1,nn+1,m);
    Cp = C(kk+2,nn+1,m);
    Cm = C(kk,nn+1,m);
    %     C(kk+1,nn+1,m+1) = sqrt(1./((2*kk+1)*((nn-m+1).*(nn+m)))) .* ...
    %         ( sqrt(((kk-m+1).*(kk+m))*(2*nn+1)) .* C0 ...
    %         - 2*pi*z * sqrt(((kk-m+2).*(kk-m+1))*row1) .* Cp ...
    %         - 2*pi*z * sqrt(((kk+m).*(kk+m-1))*row1) .* Cm );
    C(kk+1,nn+1,m+1) = sqrt(1./((2*kk+1)*((nn-m+1).*(nn+m)))) .* ...
        ( sqrt(((kk-m+1).*(kk+m).*(2*kk+1)))*row1 .* C0 ...
        -2*pi*z*sqrt((((kk-m+2).*(kk-m+1)))./((2*kk+3)))*row1.*Cp ...
        -2*pi*z*sqrt((((kk+m).*(kk+m-1)))./((2*kk-1)))*row1.*Cm );
end

% OK, that's the scalar coefficients
% Time to find the vector coefficients - Videen (43) & (44)

[nn,kk]=meshgrid([1:nmax],[1:nmax]);

matrixm=sqrt(kk.*(kk+1)) ./ sqrt(nn.*(nn+1));

central_iterator=[1:nmax].*[2:nmax+1];

[ciy,cix]=meshgrid(central_iterator,central_iterator);

mmm=0;

C0 = C(2:(nmax+1),2:(nmax+1),mmm+1);
Cp = C(3:(nmax+2),2:(nmax+1),mmm+1);
Cm = C(1:nmax,2:(nmax+1),mmm+1);

t = matrixm.*(C0 - 2*pi*z./(kk+1) .* ...
    sqrt((kk-mmm+1).*(kk+mmm+1)./((2*kk+1).*(2*kk+3))) .* Cp - ...
    2*pi*z./kk.*sqrt((kk-mmm).*(kk+mmm)./((2*kk+1).*(2*kk-1))).*Cm);

toIndexy=(ciy(:));
toIndexx=(cix(:));
A=t(:);
B=zeros(size(A));

for mmm=1:nmax
    C0 = C(2:(nmax+1),2:(nmax+1),mmm+1);
    Cp = C(3:(nmax+2),2:(nmax+1),mmm+1);
    Cm = C(1:nmax,2:(nmax+1),mmm+1);
    
    t = matrixm.*(C0 - 2*pi*z./(kk+1) .* ...
        sqrt((kk-mmm+1).*(kk+mmm+1)./((2*kk+1).*(2*kk+3))) .* Cp - ...
        2*pi*z./kk.*sqrt((kk-mmm).*(kk+mmm)./((2*kk+1).*(2*kk-1))).*Cm);

    tt=t(mmm:end,mmm:end);
    ciys=ciy(mmm:end,mmm:end);
    cixs=cix(mmm:end,mmm:end);
    
    toIndexy=[toIndexy;(ciys(:)+mmm);(ciys(:)-mmm)];
    toIndexx=[toIndexx;(cixs(:)+mmm);(cixs(:)-mmm)];
    A=[A;tt(:);tt(:)];

    t = 1i*2*pi*z*mmm./(kk.*(kk+1)).*matrixm .* C0;
    tt=t(mmm:end,mmm:end);
    B=[B;tt(:);-tt(:)];

end

B=sparse(toIndexy,toIndexx,B,nmax*(nmax+2),nmax*(nmax+2));
A=sparse(toIndexy,toIndexx,A,nmax*(nmax+2),nmax*(nmax+2));

if nargout>2
    C=C(1:nmax+1,1:nmax+1,1:nmax+1);
end

ott_warning('external');

return
