% Example file for finding spring constants
%
% PACKAGE INFO

% Specify refractive indices
n_medium = 1.34;
n_particle = 1.59;
n_relative = n_particle/n_medium;

% If you want to give all measurements in wavelengths in the surrounding
% medium, then:
wavelength = 1;
% wavelength = wavelength0 / n_medium;
% else you can give it in any units you want. Only k times lengths matters
k = 2*pi/wavelength;

radius = 1.5;
Nmax = ka2nmax(k*radius);

if Nmax < 12
    Nmax = 12;
end

% a Gaussian beam: w0 = 2/(k*tan(theta))
w0 = 0.2671; % Convergence half-angle of 50 degrees
% Use lg_mode_w0 to find w0 from the angle

% Polarisation. [ 1 0 ] is plane-polarised along the x-axis, [ 0 1 ] is
% y-polarised, and [ 1 -i ] and [ 1 i ] are circularly polarised.
polarisation = [ 1 0 ];

% Location of the focal point relative to the particle. These are the
% [ x y z ] coordinates.
beam_offset = [ 0 0 0];

[n,m,a0,b0] = bsc_pointmatch_farfield(Nmax,1,[ 0 0 w0 1 polarisation 90 beam_offset ]);
[a,b] = make_beam_vector(a0,b0,n,m);

% If you're going to do a range of particles, then the T-matrix has to
% calculated inside the loop.

% To search for a refractive index, I recommend a bisection search. For an
% example of bisection search, see find_axial_equilibrium.m

T = tmatrix_mie(Nmax,k,k*n_relative,radius);

[z,k] = axial_equilibrium(T,a,b)

