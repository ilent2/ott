function z = find_axial_equilibrium(T,a,b,z)
% find_axial_equilibrium.m - find equilibrium position along beam axis
%
% Usage:
% z = find_axial_equilibrium(T,a,b);
% z = find_axial_equilibrium(T,a,b,initial_guess);
% where T = T-matrix, a,b = multipole expansion of beam.
%
% This file is part of the package Optical tweezers toolbox 1.0
% Copyright 2006 The University of Queensland.
% See README.txt or README.m for license and details.
%
% http://www.physics.uq.edu.au/people/nieminen/software.html

newton = true;

z_precision = 1e-4;
short_distance = 1e-5;
zpoints = 45;

if nargin < 4
    z = 0;
end

z_old = z + 10*z_precision;

[Nbeam,dummy] = combined_index(length(a));
[Nparticle,dummy] = combined_index(max(size(T))/2);

equiv_ka = nmax2ka(Nbeam);

if newton
    
% Newton's method
while abs(z-z_old) > z_precision

    Nbeam2 = ka2nmax( equiv_ka + 2*pi*z );

    Nmax = max(Nparticle,Nbeam2);
    total_orders = combined_index(Nmax,Nmax);
    [n,m] = combined_index((1:total_orders)');

    T2 = change_nmax(T,Nmax);
    a0 = change_nmax(a,Nmax);
    b0 = change_nmax(b,Nmax);

    [A,B] = translate_z(Nmax,z);
    a2 = ( A*a0 + B*b0 );
    b2 = ( A*b0 + B*a0 );

    [Ashort,Bshort] = translate_z(Nmax,short_distance);

    pq = T2 * [ a2; b2 ];
    p = pq(1:total_orders);
    q = pq((total_orders+1):end);
    power = sum( abs(a2).^2 + abs(b2).^2 );
    f1 = force_z(n,m,a2,b2,p,q)/power;

    a3 = ( Ashort*a2 + Bshort*b2 );
    b3 = ( Ashort*b2 + Bshort*a2 );

    pq = T2 * [ a3; b3 ];
    p = pq(1:total_orders);
    q = pq((total_orders+1):end);
    f2 = force_z(n,m,a3,b3,p,q)/power;

    dz = short_distance * f1/(f1-f2);
    
    z_old = z;
    
    z = z + dz;
    
end

else % end of if newton

% Bisection method

% Need initial guess

% Guess the radius
radius = nmax2ka(Nparticle)/(2*pi);

z = linspace(-radius,3*radius,zpoints);

Nbeam2 = ka2nmax( equiv_ka + 2*pi*max(z) );
Nmax = max(Nparticle,Nbeam2);
total_orders = combined_index(Nmax,Nmax);
[n,m] = combined_index((1:total_orders)');

T2 = change_nmax(T,Nmax);
a0 = change_nmax(a,Nmax);
b0 = change_nmax(b,Nmax);

fz = zeros(size(z));

for nz = 1:zpoints
    
    [A,B] = translate_z(Nmax,z(nz));
    a2 = ( A*a0 + B*b0 );
    b2 = ( A*b0 + B*a0 );

    pq = T2 * [ a2; b2 ];
    p = pq(1:total_orders);
    q = pq((total_orders+1):end);
    power = sum( abs(a2).^2 + abs(b2).^2 );
    fz(nz) = force_z(n,m,a2,b2,p,q)/power;
    
    if fz(nz) < 0
        z1 = z(nz-1);
        z2 = z(nz);
        f1 = fz(nz-1);
        f2 = fz(nz);
        break
    end

end

if f1 == 0
    z = z1;
    return
end

if nz == zpoints
    error('Failed to find starting points for bisection search');
end
    
% Now the actual bisection search

while z2 - z1 > z_precision
    
    z3 = (z1+z2)/2;

    [A,B] = translate_z(Nmax,z3);
    a2 = ( A*a0 + B*b0 );
    b2 = ( A*b0 + B*a0 );

    pq = T2 * [ a2; b2 ];
    p = pq(1:total_orders);
    q = pq((total_orders+1):end);
    power = sum( abs(a2).^2 + abs(b2).^2 );
    f3 = force_z(n,m,a2,b2,p,q)/power;
   
    if f3 == 0
        z = z3;
        break
    end
    if f1*f3 < 0
        z1 = z3;
        f1 = f3;
    else
        z2 = z3;
        f2 = f3;
    end
    
end
    
z = z3;

end % end of if bisection (ie not newton)

return
