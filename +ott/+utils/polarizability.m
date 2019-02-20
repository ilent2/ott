% Polarizability calculation based on Draine & Goodman, 
% Beyond Clausius-Mossoti: wave propagation on a polarizable point lattice
% and the discrete dipole approximation,
% The Astrophysical Journal, 405:685-697, 1993 March 10

% Author: Vincent Loke
% Affiliation: Physics Dept, School of Physical Sciences
%              The University of Queensland
% Version: Pre-release (2007)

function alph = polarizability(d,m,k)

% d : lattice spacing
% m : N length vector containing relative refractive indices
% k : wave number

N = length(m); % number of dipoles
S = 0; % temporary
b1 = -1.8915316;
b2 = 0.1648469;
b3 = -1.7700004;
msqr = m.^2;
dcube = d^3;

alpha_CM = 3*dcube/(4*pi)*(msqr - 1)./(msqr + 2); % Clausius-Mossotti

%Draine & Goodman
alpha_LDR = alpha_CM./(1 + (alpha_CM/dcube).*((b1+msqr*b2+msqr*b3*S)*(k*d)^2-2/3*i*k^3*dcube));

%Draine, The Astrophysical Journal, 333:848-872, 1988 Oct 15
%alpha1 = alpha_CM./(1 - 2/3*i*k0^3*alpha_CM);

% b1 = -sqrt(4*pi/3);
% alpha1 = alpha0./(1 + (alpha0/dcube).*b1*(k*d)^2);
alph = zeros(3*N,1);

% assuming same polarizability in x, y & z directions
%alph = [alpha1; alpha1; alpha1];
for j = 1:N
%   alph(3*(j-1) + 1) = alpha1(j);
%   alph(3*(j-1) + 2) = alpha1(j);
%   alph(3*(j-1) + 3) = alpha1(j);
  alph(3*(j-1) + 1) = alpha_LDR(j);
  alph(3*(j-1) + 2) = alpha_LDR(j);
  alph(3*(j-1) + 3) = alpha_LDR(j);
end
