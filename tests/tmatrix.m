
% Ensure the ott package is in our path
addpath('../');

% Create a T-matrix for a sphere
T = ott.Tmatrix.simple('sphere', 1.0, 'wavelength0', 1.0, ...
    'index_medium', 1.0, 'index_particle', 1.2);
tol = 1.0e-6;

Ttotal = T.toTotal;
Tscat = T.toScattered;

%% Check that the total field T-matrix conserves power
assert(all(abs(1.0 - sum(abs(T.toTotal.data).^2, 2)) < tol));

%% Check we can convert to and from scattered and total field
assert(strcmpi(Ttotal.type, 'total'), 'Conversion to total');
assert(strcmpi(Tscat.type, 'scattered'), ...
    'Conversion to scattered');

% Check conversions back
assert(strcmpi(Tscat.toTotal.type, 'total'), ...
    'Conversion back to total');
assert(strcmpi(Ttotal.toScattered.type, 'scattered'), ...
    'Conversion back to scattered');

%% Check resizing T-matrix works
Tnew1 = T;
Tnew1.Nmax = Tnew1.Nmax + 5;
assert(all(Tnew1.Nmax == T.Nmax + 5), ...
    'Faild to increase Nmax with vector size');
assert(all(size(Tnew1.data) > size(T.data)), ...
    'Tmatrix size not actually increased (vector input)');

Tnew2 = Tnew1;
Tnew2.Nmax = T.Nmax;
assert(all(Tnew2.Nmax == T.Nmax), ...
    'Failed to decrease Nmax with vector size');
assert(all(size(Tnew2.data) == size(T.data)), ...
    'Tmatrix size not actually decreased (vector input)');

Tnew1 = T;
Tnew1.Nmax = T.Nmax(1) + 5;
assert(all(Tnew1.Nmax == T.Nmax + 5), ...
    'Faild to increase Nmax (scalar input)');
assert(all(size(Tnew1.data) > size(T.data)), ...
    'Tmatrix size not actually increased (scalar input)');

Tnew1 = T;
Tnew1.Nmax = [T.Nmax(1) + 5, T.Nmax(2)];
assert(all(Tnew1.Nmax == [T.Nmax(1) + 5, T.Nmax(2)]), ...
    'Faild to increase Nmax (uneven input)');
assert(size(Tnew1.data, 1) > size(T.data, 1) ...
    && size(Tnew1.data, 2) == size(T.data, 2), ...
    'Tmatrix size not increased correctly (uneven input)');

Tnew1 = T;
Tnew1.Nmax(1) = T.Nmax(1) + 5;
assert(all(Tnew1.Nmax == [T.Nmax(1) + 5, T.Nmax(2)]), ...
    'Faild to increase Nmax (index input)');
assert(size(Tnew1.data, 1) > size(T.data, 1) ...
    && size(Tnew1.data, 2) == size(T.data, 2), ...
    'Tmatrix size not increased correctly (index input)');

%% Check resizing total conserves power
Ttotal.Nmax = Ttotal.Nmax + 5;
assert(all(abs(1.0 - sum(abs(Ttotal.data).^2, 2)) < tol));

%% Check resizing scattered conserves power
Tscat.Nmax = Tscat.Nmax + 5;
assert(all(abs(1.0 - sum(abs(Tscat.toTotal.data).^2, 2)) < tol));
