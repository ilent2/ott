
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

%% Check resizing total conserves power
Ttotal.Nmax = Ttotal.Nmax + 5;
assert(all(abs(1.0 - sum(abs(Ttotal.data).^2, 2)) < tol));

%% Check resizing scattered conserves power
Tscat.Nmax = Tscat.Nmax + 5;
assert(all(abs(1.0 - sum(abs(Tscat.toTotal.data).^2, 2)) < tol));
