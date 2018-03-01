function [nn,mm,a,b] = bsc_pointmatch_focalplane( nmax, beam_type, parameters )
% bsc_pointmatch_focalplane.m
% Uses an overdetermined point-matching method to find
% spherical harmonic expansion of a laser beam.
%
% Usage:
% [n,m,a,b] = bsc_pointmatch_focalplane( nmax, beam_type, parameters );
%
% Currently available types of beams:
% 0 Gauss-Hermite beam
%   parameters: [ n m w0 P xcomponent ycomponent ]
% 1 Laguerre-Gauss beam
%   parameters: [ p l w0 P xcomponent ycomponent ]
%
% PACKAGE INFO

zero_rejection_level = 1e-8;

speed_of_light = 3.00e8;
%medium_refractive_index = 1.33;
medium_refractive_index = 1;
%beam_wavelength0 = 1064e-9;
beam_wavelength0 = 1;
beam_power = 1e-3;
speed_in_medium = speed_of_light/medium_refractive_index;
epsilon0 = 8.854187817e-12;
kappa = medium_refractive_index^2;

w0 = parameters(3);
xcomponent = parameters(5);
ycomponent = parameters(6);

k = 2*pi * medium_refractive_index / beam_wavelength0;

% For focal plane point-matching, note that M(odd n) is zero,
% and N(even n) is zero. So we can leave these out.
%
% For beams with rotational symmetry of irradiance and wavevector
% about the beam axis, we only need to consider m = +/- 1.
% However, we will program the general case of m = -n:n

% Number of unknowns is nmax * 2 (2 amn/bmn coefficients per n)
% and number of equations is number_of_points * 3
% (3 vector components at each point).
%
% So, if we use a number of points = 2*nmax, we obtain an
% overdetermined system.
%
% Why not just take points at phi = 0 and 90 deg, at varying r?

% degree and order of all modes
total_modes = nmax^2 + 2*nmax;
[nn,mm] = combined_index((1:total_modes)');
% Following 3 lines for m = +/- case
%total_modes = 2 * nmax;
%nn = sort([ 1:nmax 1:nmax ]');
%mm = ((-1).^(1:(2*nmax)))';

% radius,azimuthal angle grid
nkr = 2*(nmax + 1);
nphi = 2*(nmax + 1);
% For m = +/- 1 case
%nphi = 3;
kr = 1:nkr;
phi = ((1:nphi)-1)/nphi * 2*pi;
kr = kr.' * ones(1,nphi);
phi = ones(nkr,1) * phi;
kr = kr(:);
phi = phi(:);

% Maybe nice to choose the radial position of the points to
% be at the maxima of the spherical Bessel functions
% which are at approximately at kr =
%jnmaxkr = [ 2.08 3.34 4.51 5.65 6.76 7.85 8.93 10.01 11.08 12.14 ];

np = length(kr);

theta = ones(np,1)*pi/2;

r = kr/k;
[x,y,z] = rtp2xyz(r,theta,phi);

coefficient_matrix = zeros(3*np,2*nmax);
e_field = zeros(3*np,1);

% Find electric field at all points
rayleigh_range = k * w0^2 / 2;
q0 = - i * rayleigh_range;
qz = q0 + z;
central_irradiance = 2*beam_power / (pi*w0^2);
central_amplitude = sqrt(2*central_irradiance / ...
   (speed_in_medium*kappa));
central_amplitude = 1; % FOR TESTING

% Gaussian
beam_envelope = sqrt(2/pi) * q0./(w0+qz) .* exp( i*k * ( z + (r.^2)./(2*qz) ) );

% Bi-Gaussian
%beam_envelope = sqrt(2/pi) * q0./(w0+qz) .* exp( i*k * ( z + ((x/2).^2+(y/1).^2)./(2*qz) ) );

%radial_mode = 0;
%azimuthal_mode = 3;
%beam_envelope = (sqrt(2)*r/w0).^azimuthal_mode .* exp(-r.^2/w0^2 + i*azimuthal_mode*phi);

Ex = xcomponent * beam_envelope * central_amplitude;
Ey = ycomponent * beam_envelope * central_amplitude;
%Ey = i*Ex;
Ez = zeros(size(Ex));

e_field = [ Ex(:); Ey(:); Ez(:) ];

for n = 1:length(nn)

   % Now find RgM, RgN as appropriate for each mode
   [M,N] = vswfcart(nn(n),mm(n),kr,theta,phi,3);
   if rem(nn(n),2) == 0
      % Even n
      MN = [ M(:,1); M(:,2); M(:,3) ];
   else
      % Odd n
      MN = [ N(:,1); N(:,2); N(:,3) ];
   end
   coefficient_matrix(:,n) = MN;

end

%tic
%fprintf(1,'Beginning solution of linear system ... ');
expansion_coefficients = coefficient_matrix \ e_field;
%fprintf(1,'done!\n');
%toc

non_zero_elements = find( abs(expansion_coefficients) > max(abs(expansion_coefficients)) * zero_rejection_level );
expansion_coefficients = expansion_coefficients(non_zero_elements);
nn = nn(non_zero_elements); mm = mm(non_zero_elements);

a = zeros(size(nn));
b = zeros(size(nn));

for n = 1:length(nn)

   if rem(nn(n),2) == 0
      a(n) = expansion_coefficients(n);
      b(n) = expansion_coefficients(n) * sign(mm(n));
   else
      a(n) = expansion_coefficients(n) * sign(mm(n));
      b(n) = expansion_coefficients(n);
   end
   
end

return


