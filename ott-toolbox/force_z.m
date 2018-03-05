function force = force_z(n,m,a,b,p,q)
% force_z.m
%
% Finds z component of optical force
%
% Usage:
% Fz = forcez(n,m,a,b,p,q)
%
% What units are you using for a,b,p,q?
% If you have simple units like (using incoming/outgoing):
%     power = sum( abs(a).^2 ... )
% then divide by c to get newtons, divide by omega to get N.m
% If you have
%    sum( abs(a).^2 + abs(b).^2 ) = 1
% then the force and torque are in units of the momentum per photon
% and hbar per photon.
%
% WARNING: This code will be set up to expect a,b,p,q
% in either the incident-scattered or incoming-outgoing
% formulations! Check that it matches the one you use!
%
% This file is part of the package Optical tweezers toolbox 1.0.1
% Copyright 2006-2007 The University of Queensland.
% See README.txt or README.m for license and details.
%
% http://www.physics.uq.edu.au/people/nieminen/software.html

% Uncomment one of the following:
incidentscattered = 1; % YES, I AM USING INCIDENT-SCATTERED FORMULATION
% incidentscattered = 0; % NO, I AM NOT, I USE INCOMING-OUTGOING FORMULATION

% The force/torque calculations are easiest in the incoming-outgoing
% formulation, so convert to it if necessary
if incidentscattered
    p = 2*p + a;
    q = 2*q + b;
end

force = forcez(n,m,a,b) - forcez(n,m,p,q);

return


% Find z-component of force
% Magic formula from Crichton
function fz = forcez(n,m,a,b)

ci = combined_index(n,m);

aa = zeros(max(ci),1);
bb = zeros(max(ci),1);
aap = zeros(max(ci),1);
bbp = zeros(max(ci),1);

[nn,mm] = combined_index((1:max(ci))');

aa(ci) = a;
bb(ci) = i*b;

n1 = find( n>1 & n>abs(m));

ci1 = combined_index(n(n1)-1,m(n1));

aap(ci1) = a(n1);
bbp(ci1) = i*b(n1);

fz = 2 * mm ./ nn ./ (nn+1) .* imag( conj(aa) .* bb ) ...
    - 2 ./ (nn+1) .* sqrt( nn .* (nn+2) .*  (nn-mm+1) .* (nn+mm+1) ./ (2*nn+1) ./ (2*nn+3) ) ...
    .* imag( aa.*conj(aap) + bb.*conj(bbp) );

fz = sum(fz);

return

