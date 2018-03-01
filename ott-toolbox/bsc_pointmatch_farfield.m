function [nn,mm,a,b] = bsc_pointmatch_farfield( nmax, beam_type, parameters )
% bsc_pointmatch_farfield.m
% Uses an overdetermined point-matching method to find
% spherical harmonic expansion of a laser beam.
%
% Usage:
% [n,m,a,b] = bsc_pointmatch_farfield( nmax, beam_type, parameters );
%
% Currently available types of beams:
% 0 Gauss-Hermite beam
%   parameters: [ n m w0 P xcomponent ycomponent ]
% 1 Laguerre-Gauss beam
%   parameters: [ p l w0 P xcomponent ycomponent truncation_angle xoffset yoffset zoffset ]
%
% PACKAGE INFO

axisymmetry = 1;
%axisymmetry = 0;

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
k = 2*pi * medium_refractive_index / beam_wavelength0;
xcomponent = parameters(5);
ycomponent = parameters(6);
radial_mode = parameters(1);
azimuthal_mode = parameters(2);

% Truncation angle
if length(parameters) > 6
   truncation_angle = parameters(7);
else
   truncation_angle = 90;
end

% Offset of focal point from coordinate origin
if length(parameters) > 9
   offset = parameters(8:10);
else
   offset = [];
end

% degree and order of all modes
total_modes = nmax^2 + 2*nmax;
[nn,mm] = combined_index((1:total_modes)');

if axisymmetry
    nn = sort([ 1:nmax 1:nmax ]');
    mm = ((-1).^(1:(2*nmax)))' + azimuthal_mode;
    
    removeels=find(nn<abs(mm));
    nn(removeels)=[];
    mm(removeels)=[];
    
end

% Grid of points over sphere
ntheta = (nmax + 1);
nphi = 2*(nmax + 1);
if axisymmetry
    ntheta = 2*(nmax+1);
    nphi = 3;
end

[theta,phi] = angulargrid(ntheta,nphi);

np = length(theta);

e_field = zeros(2*np,1);

% Find electric field at all points
% In the far-field, we have:
% w = 2|z|/(k w0)     (cylindrical coords)
% r/w = kr w0 / 2 |z| (cylindrical coords)
% r = z tan(theta)    (cylindrical -> spherical conversion)
% r/w = k w0 |tan(theta)|/2 (spherical)

%central_irradiance = 2*beam_power / (pi*w0^2);
%central_amplitude = sqrt(2*central_irradiance / ...
%   (speed_in_medium*kappa));

central_amplitude = 1;

rw = k^2 * w0^2 * tan(theta).^2 / 2;

% Bi-Gaussian beam
%rw = rw .* abs(sin(phi)) + 1/9 * rw .* abs(cos(phi));
%rw = 1/9 * rw .* abs(sin(phi)) + rw .* abs(cos(phi));

%beam_envelope = rw.^(azimuthal_mode/2) .* L(p,l,rw) .* ...
%       exp(-rw/2  + i*azimuthal_mode*phi);
L = laguerre(radial_mode,abs(azimuthal_mode),rw);
beam_envelope = rw.^abs(azimuthal_mode/2) .* L .* exp(-rw/2 + i*azimuthal_mode*phi);
%beam_envelope = exp(-rw/2);

% Phase shift due to offset?
% Phase shift is exp(-i*k * offset.rhat)
if ~isempty(offset)
    rhat = rtpv2xyzv( ones(size(theta)), zeros(size(theta)), zeros(size(theta)), ones(size(theta)), theta, phi );
    [offset,rhat] = matchsize(offset,rhat);
    phase_shift = exp( -i * k * dot(offset,rhat,2) );
    beam_envelope = beam_envelope .* phase_shift;
end


outbeam = find(theta<pi*(180-truncation_angle)/180);
beam_envelope(outbeam) = 0;

% FOR TESTING
% figure;polar(theta,beam_envelope)

Ex = xcomponent * beam_envelope * central_amplitude;
Ey = ycomponent * beam_envelope * central_amplitude;

Etheta = - Ex .* cos(phi) - Ey .* sin(phi);
Ephi = - Ex .* sin(phi) + Ey .* cos(phi);

e_field = [ Etheta(:); Ephi(:) ];
coefficient_matrix = zeros(2*np,2*length(nn));

for n = 1:length(nn)

   % Now find E as appropriate for each mode
   [B,C,P] = vsh(nn(n),mm(n),theta,phi);
   coefficient_matrix(:,n) = [ C(:,2); C(:,3) ] * i^(nn(n)+1)/sqrt(nn(n)*(nn(n)+1));
   coefficient_matrix(:,n+length(nn)) = [ B(:,2); B(:,3) ] * i^nn(n)/sqrt(nn(n)*(nn(n)+1));

end

%tic
%fprintf(1,'Beginning solution of linear system ... ');
expansion_coefficients = coefficient_matrix \ e_field;
%fprintf(1,'done!\n');
%toc

a = expansion_coefficients(1:end/2);
b = expansion_coefficients(1+end/2:end);

non_zero_elements = find( abs(a) + abs(b) > max(abs(a)+abs(b)) * zero_rejection_level );
a = a(non_zero_elements);
b = b(non_zero_elements);
nn = nn(non_zero_elements);
mm = mm(non_zero_elements);

return

