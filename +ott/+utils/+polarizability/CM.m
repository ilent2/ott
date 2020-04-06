function alpha_CM = CM(spacing, index)
% Clausius-Mossoti Polarizability 
%
% Usage
%   alpha = CM(spacing, index)
%   Calculates a Nx1 element vector containing the isotropic
%   polarisabilities for N dipoles.
%
% Parameters
%   - spacing (numeric scalar) -- lattice spacing parameter
%   - index (Nx1 numeric) -- Relative refractive indices for N dipoles.

% Based on the script by Vincent Loke.
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

msqr = index(:).^2;
dcube = spacing^3;

alpha_CM = 3*dcube/(4*pi)*(msqr - 1)./(msqr + 2); % Clausius-Mossotti

