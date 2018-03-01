function [D,D2] = wigner_rotation_matrix( nmax, R )
% wigner_rotation_matrix.m
% Rotation matrix for rotation of spherical harmonics or T-matrices
%
% Usage:
% D = wigner_rotation_matrix(nmax,R)
% [D,D2] = wigner_rotation_matrix(nmax,R)
% where
% R = cartesian coordinate rotation matrix
% D = wigner D matrix, a' = D a
% D2 = double matrix, D2 = [ D 0; 0 D], so [a';b'] = D [a;b]
%
% D (and D2) are sparse
%
% This method from Choi et al., J. Chem. Phys. 111: 8825-8831 (1999)
% Note change in notation - here, use x' = Rx (x is row vector),
% a' = Da (a is column vector) etc.
%
% PACKAGE INFO

% Transform cartesian rotation matrix to spinor(?) rotation matrix
%
% Note that the order n = 1 scalar spherical harmonics are equal to
% orthogonal spinor vectors
% e+ = -exp(i*phi) sin(theta) / sqrt(2) = ( - x - iy )/sqrt(2)/r
% e- = exp(-i*phi) sin(theta) / sqrt(2) = ( x - iy )/sqrt(2)/r
% e0 = cos(theta) = z/r

% So, to transform cartesian coords to spinor coords,
% s = x C / r

C = [  1/sqrt(2) 0 -1/sqrt(2);
      -i/sqrt(2) 0 -i/sqrt(2);
       0         1  0 ];
invC = [ 1/sqrt(2) i/sqrt(2) 0;
         0         0         1;
        -1/sqrt(2) i/sqrt(2) 0 ];

% Since x' = x R, s' -> x' C = s invC R C -> s' = s (invC R C)

D = invC * R * C;

% This is also the rotation sub-matrix for n = 1
% Although, since we store spherical harmonic coeffs as column vectors,
D1 = D.';
DD = D1;

maxci = combined_index(nmax,nmax);
D = sparse(maxci,maxci);

D(1:3,1:3) = D1;

% Calculate for n = 2:nmax by recursion
for n = 2:nmax
    
    DDD = zeros(2*n+1);
    
    % Fill in whole block except for top and bottom row

    m0 = ones(2*n-1,1) * (-n:n);
    m1 = ((-n+1):(n-1)).' * ones(1,2*n+1);
    a = sqrt( (n+m0) .* (n-m0) ./ ( (n+m1) .* (n-m1) ) );
    b = sqrt( (n+m0) .* (n+m0-1) ./ ( 2*(n+m1) .* (n-m1) ) );
    
    DDD(2:end-1,2:end-1) = D1(2,2) * a(:,2:end-1) .* DD;
    DDD(2:end-1,3:end) = DDD(2:end-1,3:end) + D1(2,3) * b(:,3:end) .* DD;
    DDD(2:end-1,1:end-2) = DDD(2:end-1,1:end-2) + D1(2,1) * fliplr(b(:,3:end)) .* DD;
    
    % Top row
    
    m0 = (-n:n);
    c = sqrt( (n+m0).*(n-m0)/(n*(2*n-1)) );
    d = sqrt( (n+m0).*(n+m0-1)/(2*n*(2*n-1)) );
    
    DDD(1,2:end-1) = D1(1,2) * c(2:end-1) .* DD(1,:);
    DDD(1,3:end) = DDD(1,3:end) + D1(1,3) * d(3:end) .* DD(1,:);
    DDD(1,1:end-2) = DDD(1,1:end-2) + D1(1,1) * fliplr(d(3:end)) .* DD(1,:);
    
    % Bottom row
    
    DDD(end,2:end-1) = D1(3,2) * c(2:end-1) .* DD(end,:);
    DDD(end,3:end) = DDD(end,3:end) + D1(3,3) * d(3:end) .* DD(end,:);
    DDD(end,1:end-2) = DDD(end,1:end-2) + D1(3,1) * fliplr(d(3:end)) .* DD(end,:);

    % Dump into final matrix
    
    minci = combined_index(n,-n);
    maxci = combined_index(n,n);
    
    D(minci:maxci,minci:maxci) = DDD;
    DD = DDD;
    
end

return
