function [force,torque,spin] = forcetorque(n,m,a,b,p,q)
% FORCETORQUE finds optical force and torque
%
% [force,torque,spin] = FORCETORQUE(n,m,a,b,p,q) and
% [force,torque,spin] = FORCETORQUE(n,m,ab,pq) calculate force, torque
% and spin between incident ab beam and scattered pq beam.
% Beams can be specified as separate column vectors or combined vectors
% ab = [a;b], pq = [p;q], where a,b,p,q are the beam coefficients.
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
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

ott_warning('ott:forcetorque:depreciated', ...
    ['forcetorque.m will be replaced with ' ...
    'force_torque_farsund.m in ott1.4.']);

ott_warning('internal');

% Uncomment one of the following:
incidentscattered = 1; % YES, I AM USING INCIDENT-SCATTERED FORMULATION
% incidentscattered = 0; % NO, I AM NOT, I USE INCOMING-OUTGOING FORMULATION

if nargin==4
    if length(a)==length(b)
        labpq=length(a)/2;
        p=b(1:labpq);
        q=b(labpq+1:end);
        b=a(labpq+1:end);
        a=a(1:labpq);
    else
        ott_warning('external');
        error('[a;b] must be the same size as [p;q]!')
    end
end

% The force/torque calculations are easiest in the incoming-outgoing
% formulation, so convert to it if necessary
if incidentscattered
    p = 2*p + a;
    q = 2*q + b;
end

force = [ 0 0 0 ];
torque = [ 0 0 0 ];
spin = [ 0 0 0 ];

% z-component is easiest
force(3) = forcez(n,m,a,b) - forcez(n,m,p,q);
torque(3) = sum( m.*( abs(a).^2 + abs(b).^2 - abs(p).^2 - abs(q).^2 ) );
spin(3) = spinz(n,m,a,b) - spinz(n,m,p,q);

% Now find x,y components by rotating by 90 degrees and re-using the
% z-component formulae

ci = combined_index(n,m);
[n2,m2] = combined_index(1:combined_index(max(n),max(n)));
n2 = n2(:);
m2 = m2(:);

% First, rotate x axis onto z axis

R = calc_rotation_matrix([0 -pi/2 0]);
D = wigner_rotation_matrix(max(n),R);
D = D(:,ci);

a2 = D*a;
b2 = D*b;
p2 = D*p;
q2 = D*q;

force(1) = forcez(n2,m2,a2,b2) - forcez(n2,m2,p2,q2);
torque(1) = sum( m2.*( abs(a2).^2 + abs(b2).^2 - abs(p2).^2 - abs(q2).^2 ) );
spin(1) = spinz(n2,m2,a2,b2) - spinz(n2,m2,p2,q2);

% Finally, rotate (original) y axis onto z axis

R = calc_rotation_matrix([-pi/2 0 0]);
D = wigner_rotation_matrix(max(n),R);
D = D(:,ci);

a2 = D*a;
b2 = D*b;
p2 = D*p;
q2 = D*q;

force(2) = forcez(n2,m2,a2,b2) - forcez(n2,m2,p2,q2);
torque(2) = sum( m2.*( abs(a2).^2 + abs(b2).^2 - abs(p2).^2 - abs(q2).^2 ) );
spin(2) = spinz(n2,m2,a2,b2) - spinz(n2,m2,p2,q2);

ott_warning('external');

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
bb(ci) = 1i*b;

n1 = find( n>1 & n>abs(m));

ci1 = combined_index(n(n1)-1,m(n1));

aap(ci1) = a(n1);
bbp(ci1) = 1i*b(n1);

fz = 2 * mm ./ nn ./ (nn+1) .* imag( conj(aa) .* bb ) ...
    - 2 ./ (nn+1) .* sqrt( nn .* (nn+2) .*  (nn-mm+1) .* (nn+mm+1) ./ (2*nn+1) ./ (2*nn+3) ) ...
    .* imag( aa.*conj(aap) + bb.*conj(bbp) );

fz = sum(fz);

return

% Also magic formula for z-component of spin
function sz = spinz(n,m,a,b)

ci = combined_index(n,m);

aa = zeros(max(ci),1);
bb = zeros(max(ci),1);
aap = zeros(max(ci),1);
bbp = zeros(max(ci),1);

[nn,mm] = combined_index((1:max(ci))');

aa(ci) = a;
bb(ci) = 1i*b;

n1 = find( n>1 & n>abs(m));

ci1 = combined_index(n(n1)-1,m(n1));

aap(ci1) = a(n1);
bbp(ci1) = 1i*b(n1);

sz = mm ./ nn ./ (nn+1) .* ( abs(aa).^2 + abs(bb).^2 ) ...
    - 2 ./ (nn+1) .* sqrt( nn .* (nn+2) .*  (nn-mm+1) .* (nn+mm+1) ./ (2*nn+1) ./ (2*nn+3) ) ...
    .* real( aa.*conj(bbp) - bb.*conj(aap) );

sz = sum(sz);

return

