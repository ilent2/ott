% Example of calculation of force in a Gaussian beam trap
%
% Figure 3 in Nieminen et al., Optical tweezers computational toolbox,
% Journal of Optics A 9, S196-S203 (2007)
%
% How long should this take?
% With the original settings, on a 3GHz PC, it took us about 1.5 hours
% It's slow because the force is calculated for particles a long way from
% the focus, forcing a large Nmax to be used. For a particle of radius=0.5
% wavelengths, and only going out 2 wavelengths away from the focus cuts
% the time to 1.5 seconds.
%
% This file is part of the package Optical tweezers toolbox 1.0.1
% Copyright 2006-2007 The University of Queensland.
% See README.txt or README.m for license and details.
%
% http://www.physics.uq.edu.au/people/nieminen/software.html

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

%radius = 2.5;
radius = 1.0;
Nmax = ka2nmax(k*radius);

diam_microns = radius * 1.064 * 2 / n_medium

if Nmax < 12
    Nmax = 12;
end

% Specify the beam width. We can either start with the numerical
% aperture (NA) or the beam convergence angle. Either way, we convert
% to the equivalent paraxial beam waist, which is the w0 we put into the
% paraxial beam to obtain the desired (non-paraxial) far field.
% For a Gaussian beam: w0 = 2/(k*tan(theta))
NA = 1.25;
beam_angle = asin(NA/n_medium)*180/pi
w0 = lg_mode_w0( [ 0 0 ], beam_angle );
             
% Polarisation. [ 1 0 ] is plane-polarised along the x-axis, [ 0 1 ] is
% y-polarised, and [ 1 -i ] and [ 1 i ] are circularly polarised.
polarisation = [ 1 i ];

% Location of the focal point relative to the particle. These are the
% [ x y z ] coordinates.
beam_offset = [ 0 0 0];

[n,m,a0,b0] = bsc_pointmatch_farfield(Nmax,1,[ 0 0 w0 1 polarisation 90 beam_offset ]);
[a,b] = make_beam_vector(a0,b0,n,m);

T = tmatrix_mie(Nmax,k,k*n_relative,radius);

%z = linspace(-8,8,80);
%r = linspace(-4,4,80);
z = linspace(-5,5,80);
r = linspace(-3,3,80);

fz = zeros(size(z));
fr = zeros(size(r));

for nz = 1:length(z)
    
    [Nbeam,dummy] = combined_index(length(a));
    [Nparticle,dummy] = combined_index(max(size(T))/2);

    equiv_ka = nmax2ka(Nbeam);
    Nbeam2 = ka2nmax( equiv_ka + 2*pi*abs(z(nz)) );

    Nmax = max(Nparticle,Nbeam2);
    total_orders = combined_index(Nmax,Nmax);
    [n,m] = combined_index((1:total_orders)');

    T = change_nmax(T,Nmax);
    a0 = change_nmax(a,Nmax);
    b0 = change_nmax(b,Nmax);

    [A,B] = translate_z(Nmax,z(nz));
    a2 = ( A*a0 + B*b0 );
    b2 = ( A*b0 + B*a0 );

    pq = T * [ a2; b2 ];
    p = pq(1:total_orders);
    q = pq((total_orders+1):end);
    power = sum( abs(a2).^2 + abs(b2).^2 );
    fz(nz) = force_z(n,m,a2,b2,p,q)/power;
    
end

% Find equilibrium point
nequi = min(find(fz<0));
zeq = interp1(fz((nequi-1):nequi),z((nequi-1):nequi),0)

[Nbeam,dummy] = combined_index(length(a));
[Nparticle,dummy] = combined_index(max(size(T))/2);

equiv_ka = nmax2ka(Nbeam);
Nbeam = ka2nmax( equiv_ka + 2*pi*abs(zeq) );

Nmax = max(Nparticle,Nbeam);
total_orders = combined_index(Nmax,Nmax);
[n,m] = combined_index((1:total_orders)');

T = change_nmax(T,Nmax);
a0 = change_nmax(a,Nmax);
b0 = change_nmax(b,Nmax);

[A,B] = translate_z(Nmax,zeq);
a1 = ( A*a0 + B*b0 );
b1 = ( A*b0 + B*a0 );

for nr = 1:length(r)
    
    [Nbeam,dummy] = combined_index(length(a));
    [Nparticle,dummy] = combined_index(max(size(T))/2);

    equiv_ka = nmax2ka(Nbeam);
    Nbeam2 = ka2nmax( equiv_ka + 2*pi*abs(r(nr)) );

    Nmax = max(Nparticle,Nbeam2);
    total_orders = combined_index(Nmax,Nmax);
    [n,m] = combined_index((1:total_orders)');

    T = change_nmax(T,Nmax);
    a0 = change_nmax(a1,Nmax);
    b0 = change_nmax(b1,Nmax);
    
    Rx = [ 0 0 -1; 0 1 0; 1 0 0 ];
    Dx = wigner_rotation_matrix(Nmax,Rx);
    [A,B] = translate_z(Nmax,r(nr));
    a2 = ( A * Dx*a0 + B * Dx*b0 );
    b2 = ( A * Dx*b0 + B * Dx*a0 );

    pq = T * [ a2; b2 ];
    p = pq(1:total_orders);
    q = pq((total_orders+1):end);
    power = sum( abs(a2).^2 + abs(b2).^2 );
    fr(nr) = force_z(n,m,a2,b2,p,q)/power;
    
end

figure; plot(z,fz);
xlabel('{\it z} (x\lambda)');
ylabel('{\it Q_z}');
aa = axis;
hold on
line(aa(1:2),[ 0 0 ]);
line([0 0],aa(3:4));

figure; plot(r,fr);
xlabel('{\it r} (x\lambda)');
ylabel('{\it Q_r}');
aa = axis;
hold on
line(aa(1:2),[ 0 0 ]);
line([0 0],aa(3:4));


