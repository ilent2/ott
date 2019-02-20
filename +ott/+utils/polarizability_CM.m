% Clausius-Mossoti Polarizability 

% Author: Vincent Loke
% Affiliation: Physics Dept, School of Physical Sciences
%              The University of Queensland
% Version: Pre-release (2007)

function alph = polarizability_CM(d, m)

% m : N length vector containing relative refractive indices
% d : lattice spacing

N = length(m); % number of dipoles
msqr = m.^2;
dcube = d^3;

alpha_CM = 3*dcube/(4*pi)*(msqr - 1)./(msqr + 2); % Clausius-Mossotti

alph = zeros(3*N,1);

% assuming same polarizability in x, y & z directions
for j = 1:N
  alph(3*(j-1) + 1) = alpha_CM(j);
  alph(3*(j-1) + 2) = alpha_CM(j);
  alph(3*(j-1) + 3) = alpha_CM(j);
end
