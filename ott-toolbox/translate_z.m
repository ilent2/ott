function [A,B,C] = translate_z(nmax,r)
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

if numel(r)>1
    A=cell(numel(r),1);
    B=A;
    for ii=1:numel(r)
        [A{ii},B{ii}]=translate_z(nmax,r(ii));
    end
    C=0;
    ott_warning('external');
    return
end

if r==0
    A=sparse(1:(nmax^2+nmax*2),1:(nmax^2+nmax*2),1);
    B=sparse((nmax^2+nmax*2),(nmax^2+nmax*2));
    C=A;
    ott_warning('external');
    return
end

% Having pre-computed a_nm it's fast?
% nmax=3;
% r=0.1;

%some "setup" values:s
m=0;
fval=2*nmax+1;
nd=[m:fval];

kr=2*pi*r;

%compute seed functions:
C_nd00=[sqrt(2*nd+1).*sbesselj(nd,kr)];

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
C=zeros(nmax+2,nmax+1,nmax+1);
C(:,:,1)=C_ndn0(2:(nmax+3),2:(nmax+2));

%Having computed anm for m=0; cases we now can compute anm for all
%remaining cases:
ANM=anm_l([0:2*nmax+1].',[1:nmax]);
IANM=1./ANM;
for m=1:nmax
    
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
    
    C(:,:,m+1)=C_ndn0(2:(nmax+3),2:(nmax+2));
    
end

% OK, that's the scalar coefficients
% Time to find the vector coefficients - Videen (43) & (44)

[nn,kk]=meshgrid([1:nmax],[1:nmax]);

matrixm=sqrt(kk.*(kk+1)) ./ sqrt(nn.*(nn+1));

central_iterator=[1:nmax].*[2:nmax+1];

[ciy,cix]=meshgrid(central_iterator,central_iterator);

mmm=0;

C0 = C(2:(nmax+1),2:end,mmm+1);
Cp = C(3:end,2:end,mmm+1);
Cm = C(1:nmax,2:end,mmm+1);

t = matrixm.*(C0 - kr./(kk+1) .* ...
    sqrt((kk-mmm+1).*(kk+mmm+1)./((2*kk+1).*(2*kk+3))) .* Cp - ...
    kr./kk.*sqrt((kk-mmm).*(kk+mmm)./((2*kk+1).*(2*kk-1))).*Cm);

toIndexy=(ciy(:));
toIndexx=(cix(:));
A=t(:);
B=zeros(size(A));

for mmm=1:nmax
    
    C0 = C(2:(nmax+1),2:end,mmm+1);
    Cp = C(3:end,2:end,mmm+1);
    Cm = C(1:nmax,2:end,mmm+1);

    t = matrixm.*(C0 - kr./(kk+1) .* ...
        reshape(ANM(kk+1,mmm),size(kk)) .* Cp - ...
        kr./kk.*reshape(ANM(kk,mmm),size(kk)).*Cm);

    tt=t(mmm:end,mmm:end);
    ciys=ciy(mmm:end,mmm:end);
    cixs=cix(mmm:end,mmm:end);
    
    toIndexy=[toIndexy;(ciys(:)+mmm);(ciys(:)-mmm)];
    toIndexx=[toIndexx;(cixs(:)+mmm);(cixs(:)-mmm)];
    A=[A;tt(:);tt(:)];

    t = 1i*kr*mmm./(kk.*(kk+1)).*matrixm .* C0;
    tt=t(mmm:end,mmm:end);
    B=[B;tt(:);-tt(:)];

end
B=sparse(toIndexy,toIndexx,B,nmax*(nmax+2),nmax*(nmax+2));
A=sparse(toIndexy,toIndexx,A,nmax*(nmax+2),nmax*(nmax+2));

ott_warning('external');

end

function a_nm = anm_l(n,m);
fn=1./(2*n+1)./(2*n+3);
a_nm=sqrt((n+abs(m)+1).*(n-abs(m)+1).*fn);
a_nm(n<0)=0;
a_nm(abs(m)>n)=0;
end

function b_nm = bnm_l(n,m);
b_nm=(2*(m<0)-1).*sqrt((n-m-1).*(n-m)./(2*n-1)./(2*n+1));
b_nm(abs(m)>n)=0;
end
