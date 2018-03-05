function [out1,out2,out3,out4] = spring_constant(T,a,b,z)
% spring_constant.m
%
% Finds (x,y,z) spring constant for particle in trap.
%
% Usage:
% k = spring_constant(T,a,b);
% [k,z] = spring_constant(T,a,b);
% [kx,ky,kz,z] = spring_constant(T,a,b);
% k = spring_constant(T,a,b,z);
%
% The beam coefficients a,b need to be vectors with the standard packing
% If you the n,m,a,b output by bsc_*, try:
% [a,b] = make_beam_vector(a0,b0,n,m);
%
% PACKAGE INFO

short_distance = 1e-5;
equilibrium_force_level = 1e-10;

if nargin < 4
    z = find_axial_equilibrium(T,a,b);
end

[Nbeam,dummy] = combined_index(length(a));
[Nparticle,dummy] = combined_index(max(size(T))/2);

equiv_ka = nmax2ka(Nbeam);
Nbeam2 = ka2nmax( equiv_ka + 2*pi*z );

Nmax = max(Nparticle,Nbeam2);
total_orders = combined_index(Nmax,Nmax);
[n,m] = combined_index((1:total_orders)');

T = change_nmax(T,Nmax);
a0 = change_nmax(a,Nmax);
b0 = change_nmax(b,Nmax);

[A,B] = translate_z(Nmax,z);
a = ( A*a0 + B*b0 );
b = ( A*b0 + B*a0 );

Rx = [ 0 0 -1; 0 1 0; 1 0 0 ];
Ry = [ 1 0 0; 0 0 -1; 0 1 0 ];
[A,B] = translate_z(Nmax,short_distance);
Dx = wigner_rotation_matrix(Nmax,Rx);
Dy = wigner_rotation_matrix(Nmax,Ry);

a2 = ( A * Dx*a + B * Dx*b );
b2 = ( A * Dx*b + B * Dx*a );

pq = T * [ a2; b2 ];
p = pq(1:total_orders);
q = pq((total_orders+1):end);
power = sum( abs(a2).^2 + abs(b2).^2 );
fx = force_z(n,m,a2,b2,p,q)/power;

a2 = ( A * Dy*a + B * Dy*b );
b2 = ( A * Dy*b + B * Dy*a );

pq = T * [ a2; b2 ];
p = pq(1:total_orders);
q = pq((total_orders+1):end);
power = sum( abs(a2).^2 + abs(b2).^2 );
fy = force_z(n,m,a2,b2,p,q)/power;

a2 = ( A*a + B*b );
b2 = ( A*b + B*a );

pq = T * [ a2; b2 ];
p = pq(1:total_orders);
q = pq((total_orders+1):end);
power = sum( abs(a2).^2 + abs(b2).^2 );
fz = force_z(n,m,a2,b2,p,q)/power;

% And now, just in case the particle is not actually in an equilibrium
% position, find the force when undisplaced
pq = T * [ a; b ];
p = pq(1:total_orders);
q = pq((total_orders+1):end);
power = sum( abs(a).^2 + abs(b).^2 );
[f,t] = forcetorque(n,m,a,b,p,q);
force0 = f/power;

if sum(abs(force0)) > equilibrium_force_level
    warning('Particle doesn''t seem to be in equilibrium position');
end

k = ( force0 - [ fx fy fz ] ) / short_distance;

if nargout > 2
    out1 = k(1);
    out2 = k(2);
    out3 = k(3);
    out4 = z;
else
    out1 = k;
    out2 = z;
    out3 = 0;
    out4 = 0;
end

return