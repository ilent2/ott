function alph = polarizability_LDR(d,m,kvec,E0)
% Lattice dispersion relation polarizablity
%
% Polarizability calculation based on Draine & Goodman, 
% Beyond Clausius-Mossoti: wave propagation on a polarizable point lattice
% and the discrete dipole approximation,
% The Astrophysical Journal, 405:685-697, 1993 March 10
%
% alpha = polarizability_LDR(d, m, kvec, E0)
%
% m : N length vector containing relative refractive indices
%                                (isotropic version)   
% d : lattice spacing
% kvec : wave vector [kx ky kz]     e.g. [0 0 1] z-direction
%
% E0 : E-field polarization [Ex Ey Ez]   [1 0 0] x-polarized
%                                        [1 i 0] left-handed circ pol.  
%
% Author: Vincent Loke
% Affiliation: Physics Dept, School of Physical Sciences
%              The University of Queensland
% Version: Pre-release (2007)

k0 = 2*pi;
N = length(m); % number of dipoles
b1 = -1.8915316;
b2 = 0.1648469;
b3 = -1.7700004;
msqr = m.^2;
dcube = d^3;

if nargin > 3  % we have polarization info   
  a_hat = kvec/norm(kvec);
  e_hat = E0/norm(E0);
  S = 0;
  for j = 1:3
    S = S + (a_hat(j)*e_hat(j))^2;
  end
else           % use randomly-oriented value; also for non-plane wave
  S = .2;
end    

alpha_CM = 3*dcube/(4*pi)*(msqr - 1)./(msqr + 2); % Clausius-Mossotti
alpha_LDR = alpha_CM./(1 + (alpha_CM/dcube).*((b1+msqr*b2+msqr*b3*S)*(k0*d)^2-2/3*i*k0^3*dcube));

alph = zeros(3*N,1);
% assuming same polarizability in x, y & z directions
for j = 1:N
  alph(3*(j-1) + 1) = alpha_LDR(j);
  alph(3*(j-1) + 2) = alpha_LDR(j);
  alph(3*(j-1) + 3) = alpha_LDR(j);
end
