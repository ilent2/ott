function [nn,mm,a,b] = bsc_bessel_farfield( nmax, beam_type, parameters )
% bsc_bessel_farfield.m
% Uses an overdetermined point-matching method to find
% spherical harmonic expansion of a laser beam. It uses a half sine like
% function to simulate a beam with annular k. This can lead to a bessel
% like beam in the focal plane. 
%
% The higher the truncation angle, the less Bessel-like the beam
% seems.
%
% At present, there is no systematic way to decide what the beam shape is.
% However, this being said the beam has the correct Poynting vector field.
%
% Usage:
% [n,m,a,b] = bsc_bessel_farfield( nmax, beam_type, parameters );
%
% Currently available types of beams:
% 0 Bessel 0
%   parameters: [ crossing l angl P xcomponent ycomponent field_truncation]
%
% crossing is the sine like curve you want. Because the code originally had
% a Bessel beam in the far-field we've kept this feature to generate the
% sine-like curves in the annulus. E.g. Crossing of 1 gives the zeroth peak of
% a Bessel function which leads to a Gaussian like function. Crossing
% numbers > 3 are recommended. You have to experiement to get the profile
% you want.
%
% l is the angular momentum mode.
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

trangle = parameters(3);
k = 2*pi * medium_refractive_index / beam_wavelength0;
xcomponent = parameters(5);
ycomponent = parameters(6);
crossing = parameters(1);
azimuthal_mode = parameters(2);

% Truncation angle
if length(parameters) > 6
    truncation_angle = parameters(7);
else
    truncation_angle = trangle;
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
    total_modes = 2 * nmax;
    nn = sort([ 1:nmax 1:nmax ]');
    mm = ((-1).^(1:(2*nmax)))' + azimuthal_mode;
    
    removeels=find(nn<abs(mm));
    nn(removeels)=[];
    mm(removeels)=[];
    
    total_modes=total_modes-length(removeels);
end

%low order polyfit to zero bessel. nu=0 for 32 crossings.
pl=[3.063230681485826e-13,-5.409370459846392e-11, 4.140250813983338e-09,-1.800381039704490e-07, 4.899736360978252e-06,-8.668392800340810e-05, 1.003283964223102e-03,-7.460513408785807e-03, 3.408757318457949e-02, 3.053582162968942e+00,-6.765173623967824e-01];
%high order polyfit to zero bessel. nu=0 for 33--3183.
ph=[9.849331224056770e-35,-1.677767643312824e-30, 1.233230135756182e-26,-5.120391079825781e-23, 1.320345035667242e-19,-2.190091523570770e-16, 2.341773420953136e-13,-1.574740179965436e-10, 6.302974829391556e-08, 3.141579061259407e+00,-7.840391685281115e-01];

if crossing<33
    pl=[2.404825504706033e+00,
        5.520078107770384e+00,
        8.653727905625018e+00,
        1.179153441945821e+01,
        1.493091769885970e+01,
        1.807106395526346e+01,
        2.121163662949842e+01,
        2.435247152470611e+01,
        2.749347912686714e+01,
        3.063460646148027e+01,
        3.377582021212756e+01,
        3.691709834862685e+01,
        4.005842576261244e+01,
        4.319979171267694e+01,
        4.634118836705095e+01,
        4.948260989430253e+01,
        5.262405183899220e+01,
        5.576551075352520e+01,
        5.890698392499183e+01,
        6.204846918939744e+01,
        6.518996479953800e+01,
        6.833146932928049e+01,
        7.147298160306086e+01,
        7.461450064317542e+01,
        7.775602562983886e+01,
        8.089755587054134e+01,
        8.403909077627365e+01,
        8.718062984288920e+01,
        9.032217263635329e+01,
        9.346371878096545e+01,
        9.660526794987861e+01,
        9.974681985740902e+01];
    kr=pl(crossing);
    
    if crossing>1
        krm1=pl(crossing-1);
    else
        krm1=0;
    end
else
    kr=polyval(ph,crossing);
    krm1=polyval(ph,crossing-1);
end

% Grid of points over sphere
ntheta = 4*(nmax + 1);
nphi = 2*(nmax + 1);
%ntheta = ceil( 1.5 * sqrt(total_modes) );
%nphi = ceil( sqrt(total_modes) );
if axisymmetry
    ntheta = 16*(nmax+1);
    nphi = 3;
end

[theta,phi] = angulargrid(ntheta,nphi);

np = length(theta);

coefficient_matrix = zeros(2*np,2*nmax);
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

unk=kr./tan(trangle*pi/180);

thetatrunc=atan2(krm1,unk);

rw =  unk * tan(theta);

beam_envelope = besselj(0,rw).*exp(1i*azimuthal_mode*phi);

outbeam = find(theta>(pi-thetatrunc));

beam_envelope(outbeam)=0;

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

for n = 1:length(nn)
    
    % Now find E as appropriate for each mode
    [B,C,P] = vsh(nn(n),mm(n),theta,phi);
    coefficient_matrix(:,n) = [ C(:,2); C(:,3) ] * i^(nn(n)+1)/sqrt(nn(n)*(nn(n)+1));
    coefficient_matrix(:,n+total_modes) = [ B(:,2); B(:,3) ] * i^nn(n)/sqrt(nn(n)*(nn(n)+1));
    
end

%tic
%fprintf(1,'Beginning solution of linear system ... ');
expansion_coefficients = coefficient_matrix \ e_field;
%fprintf(1,'done!\n');
%toc

a = expansion_coefficients(1:total_modes);
b = expansion_coefficients((1+total_modes):end);

non_zero_elements = find( abs(a) + abs(b) > max(abs(a)+abs(b)) * zero_rejection_level );
a = a(non_zero_elements);
b = b(non_zero_elements);
nn = nn(non_zero_elements);
mm = mm(non_zero_elements);

return

