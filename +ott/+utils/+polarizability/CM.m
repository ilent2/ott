function alpha_CM = polarizability_CM(d, index)
% Clausius-Mossoti Polarizability 
%
% Usage
%   alpha = polarizability_LDR(spacing, index)
%   Calculates a Nx1 element vector containing the isotropic
%   polarisabilities for N dipoles.
%
%   alpha = polarizability_LDR(spacing, index, kvec, E0)
%   As above but specifies the polarisability information for use
%   with plane wave illumination.
%
% Parameters
%   - spacing (numeric scalar) -- lattice spacing parameter
%   - index (Nx1 numeric) -- Relative refractive indices for N dipoles.
%   - kvec (1x3 numeric) -- Wave vector [kx, ky, kz]
%   - E0 (1x3 numeric) -- E-field polarisation [Ex, Ey, Ez]

% Based on the script by Vincent Loke.
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

msqr = index(:).^2;
dcube = d^3;

alpha_CM = 3*dcube/(4*pi)*(msqr - 1)./(msqr + 2); % Clausius-Mossotti
