function D = wigner_rotation_matrix( nmax, R )
% WIGNER_ROTATION_MATRIX rotation matrix for rotation of spherical
% harmonics or T-matrices.
%
% D = WIGNER_ROTATION_MATRIX(nmax,R) calculates the rotation matrix
% for the VSH given a 3x3 coordinate rotation matrix R.  Usage: a' = D a.
%
% This method from Choi et al., J. Chem. Phys. 111: 8825-8831 (1999)
% Note change in notation - here, use x' = Rx (x is row vector),
% a' = Da (a is column vector) etc.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Check inputs
assert(isnumeric(nmax) && isscalar(nmax), ...
    'nmax must be a numeric scalar')
assert(isnumeric(R) && isequal(size(R), [3,3]), ...
    'R must be a 3x3 rotation matrix')

ott.warning('internal');

% Transform cartesian rotation matrix to spinor(?) rotation matrix
%
% Note that the order n = 1 scalar spherical harmonics are equal to
% orthogonal spinor vectors
% e+ = -exp(i*phi) sin(theta) / sqrt(2) = ( - x - iy )/sqrt(2)/r
% e- = exp(-i*phi) sin(theta) / sqrt(2) = ( x - iy )/sqrt(2)/r
% e0 = cos(theta) = z/r

% So, to transform cartesian coords to spinor coords,
% s = x C / r

% C = [  1/sqrt(2) 0 -1/sqrt(2);
%       -i/sqrt(2) 0 -i/sqrt(2);
%        0         1  0 ];
% invC = [ 1/sqrt(2) i/sqrt(2) 0;
%          0         0         1;
%         -1/sqrt(2) i/sqrt(2) 0 ];

C = [  1/sqrt(2) 0 -1/sqrt(2);
      1i/sqrt(2) 0 1i/sqrt(2);
       0         1  0 ];
invC = [ 1/sqrt(2) -1i/sqrt(2) 0;
         0         0         1;
         -1/sqrt(2) -1i/sqrt(2) 0 ];

% Since x' = x R, s' -> x' C = s invC R C -> s' = s (invC R C)

D = invC * R * C;

% This is also the rotation sub-matrix for n = 1
% Although, since we store spherical harmonic coeffs as column vectors,
D1 = D.';
DD = D1;

X = {};
X{1} = sparse(D1);

for n = 2:nmax

    DDD = complex(zeros(2*n+1));

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

    X{end+1} = DDD;

    DD = DDD;

end

D = blkdiag(X{:});

ott.warning('external');

