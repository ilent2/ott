function alph = polarizability_FCD(d,m,kvec,E0)
% Filtered coupled dipole polarizability
%
% alpha = polarizability_FCD(d, m, kvec, E0) calculates the polarizability
% using the filtered coupled dipole approach.
%
% m : N length vector containing relative refractive indices
%                                (isotropic version)   
% d : lattice spacing
% kvec : wave vector [kx ky kz]     e.g. [0 0 1] z-direction
% E0 : E-field polarization [Ex Ey Ez]   [1 0 0] x-polarized
%                                        [1 i 0] left-handed circ pol.  
%
% Author: Vincent Loke
% Affiliation: Physics Dept, School of Physical Sciences
%              The University of Queensland
% Version: Pre-release (2007)

k = 2*pi;
N = length(m); % number of dipoles

alpha_CM = 3*d^3/(4*pi)*(m.^2 - 1)./(m.^2 + 2); % Clausius-Mossotti
alpha_FCD = alpha_CM./(1 + (alpha_CM/d^3).*(4/3*(k*d)^2 + 2/3*(i + log((pi-k*d)/(pi+k*d))/pi)*k^3*d^3));

alph = zeros(3*N,1);
% assuming same polarizability in x, y & z directions
for j = 1:N
  alph(3*(j-1) + 1) = alpha_FCD(j);
  alph(3*(j-1) + 2) = alpha_FCD(j);
  alph(3*(j-1) + 3) = alpha_FCD(j);
end
