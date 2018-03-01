function [A,B,C] = translate_z(nmax,z)
% translate_z.m - translation of VSWFs along z axis.
%
% Usage:
% [A,B] = translate_z(nmax,z);
% [A,B,C] = translate_z(nmax,z);
% where
% M' = A M + B N; N' = B M + A N
%
% C are the scalar SWF translation coefficients
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
% PACKAGE INFO

% Refs:
% N. A. Gumerov and R. Duraiswami, "Recursions for the computation of
% multipole translation and rotation coefficients for the 3-D Helmholtz
% equation", SIAM J. Sci. Comput. 25(4), 1344-1381 (2003)
%
% G. Videen, "Light scattering from a sphere near a plane interface",
% chapter 5 (pp. 81-96) % in F. Moreno and F. Gonzalez (eds), "Light
% Scattering from Microstructures", Lecture Notes in Physics 534,
% Springer-Verlag, Berlin, 2000

if numel(z)>1
    A=cell(numel(z),1);
    B=A;
    for ii=1:numel(z)
        [A{ii},B{ii}]=translate_z(nmax,z(ii));
    end
    C=0;
    return
end

if z==0
    A=sparse(1:(nmax^2+nmax*2),1:(nmax^2+nmax*2),1);
    B=0;
    C=A;
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
else
    C(:,1,1) = sqrt(2*k+1) .* sbesselj(k,2*pi*z);
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
        ( sqrt(((kk-m+1).*(kk+m).*(2*kk+1))*row1) .* C0 ...
        -2*pi*z*sqrt((((kk-m+2).*(kk-m+1))*row1)./((2*kk+3)*row1)).*Cp ...
        -2*pi*z*sqrt((((kk+m).*(kk+m-1))*row1)./((2*kk-1)*row1)).*Cm );
end

% OK, that's the scalar coefficients
% Time to find the vector coefficients - Videen (43) & (44)

kkk = ones(nmax,nmax,nmax+1);
nnn = ones(nmax,nmax,nmax+1);
for k = 1:nmax
    kkk(k,:,:) = k;
    nnn(:,k,:) = k;
end
mmm = ones(nmax,nmax,nmax+1);
for m = 0:nmax
   mmm(:,:,m+1) = m;
end
C0 = C(2:(nmax+1),2:(nmax+1),1:(nmax+1));
Cp = C(3:(nmax+2),2:(nmax+1),1:(nmax+1));
Cm = C(1:nmax,2:(nmax+1),1:(nmax+1));

A = C0 - 2*pi*z./(kkk+1) .* ...
    sqrt((kkk-mmm+1).*(kkk+mmm+1)./((2*kkk+1).*(2*kkk+3))) .* Cp - ...
    2*pi*z./kkk.*sqrt((kkk-mmm).*(kkk+mmm)./((2*kkk+1).*(2*kkk-1))).*Cm;

B = 1i*2*pi*z*mmm./(kkk.*(kkk+1)) .* C0;

% Videen uses a different normalisation, so adjust to our VSWFs
A = A .* sqrt(kkk.*(kkk+1)) ./ sqrt(nnn.*(nnn+1));
B = B .* sqrt(kkk.*(kkk+1)) ./ sqrt(nnn.*(nnn+1));

% Dump the extra values from C
C = C(1:(nmax+1),1:(nmax+1),1:(nmax+1));


% Testing, testing, 1, 2, 3 ...
%C
%sum(C.*C,2)
%VV = C(:,1,1);
%sum(VV.*VV)

% Now stick all these values in the N_orders x N_orders translations
% matrices. These will be sparse since they're diagonal in m.

A_elements = [];
B_elements = [];
C_elements = [];
v_cols = [];
v_rows = [];
s_cols = [];
s_rows = [];
for m = 0:nmax
    mm = max(1,m);
    k = mm:nmax;
    n = mm:nmax;
    k0 = m:nmax;
    n0 = m:nmax;
    v_col = combined_index(k,m);
    v_row = combined_index(n,m);
    s_col = combined_index(k0,m) + 1;
    s_row = combined_index(n0,m) + 1;
    v_col_matrix = v_col' * ones(size(v_row));
    v_row_matrix = ones(size(v_col')) * v_row;
    s_col_matrix = s_col' * ones(size(s_row));
    s_row_matrix = ones(size(s_col')) * s_row;
    A_elements_matrix = A(k,n,m+1);
    B_elements_matrix = B(k,n,m+1);
    C_elements_matrix = C(k0+1,n0+1,m+1);
    A_elements = [ A_elements; A_elements_matrix(:) ];
    B_elements = [ B_elements; B_elements_matrix(:) ];
    C_elements = [ C_elements; C_elements_matrix(:) ];
    v_cols = [ v_cols; v_col_matrix(:) ];
    v_rows = [ v_rows; v_row_matrix(:) ];
    s_cols = [ s_cols; s_col_matrix(:) ];
    s_rows = [ s_rows; s_row_matrix(:) ];
    if m > 0
        v_col = combined_index(k,-m);
        v_row = combined_index(n,-m);
        s_col = combined_index(k0,-m) + 1;
        s_row = combined_index(n0,-m) + 1;
        v_col_matrix = v_col' * ones(size(v_row));
        v_row_matrix = ones(size(v_col')) * v_row;
        s_col_matrix = s_col' * ones(size(s_row));
        s_row_matrix = ones(size(s_col')) * s_row;
        A_elements = [ A_elements; A_elements_matrix(:) ];
        B_elements = [ B_elements; -B_elements_matrix(:) ];
        C_elements = [ C_elements; C_elements_matrix(:) ];
        v_cols = [ v_cols; v_col_matrix(:) ];
        v_rows = [ v_rows; v_row_matrix(:) ];
        s_cols = [ s_cols; s_col_matrix(:) ];
        s_rows = [ s_rows; s_row_matrix(:) ];
    end
end

A = sparse(v_rows,v_cols,A_elements);
B = sparse(v_rows,v_cols,B_elements);
C = sparse(s_rows,s_cols,C_elements);

return
