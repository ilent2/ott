function alpha_FCD = FCD(spacing, index, varargin)
% Filtered coupled dipole polarizability
%
% Usage
%   alpha = FCD(spacing, index)
%   Calculates a Nx1 element vector containing the isotropic
%   polarizabilities for N dipoles.
%
% Parameters
%   - spacing (numeric scalar) -- lattice spacing parameter
%   - index (Nx1 numeric) -- Relative refractive indices for N dipoles.
%
% Optional named arguments
%   - k0 (numeric) -- Wavenumber to scale spacing by.  Default: ``2*pi``.

% Based on the script by Vincent Loke.
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

p = inputParser;
p.addParameter('k0', 2*pi);
p.parse(varargin{:});

k = p.Results.k0;
msqr = index(:).^2;

alpha_CM = 3*spacing^3/(4*pi)*(msqr - 1)./(msqr + 2); % Clausius-Mossotti
alpha_FCD = alpha_CM./(1 + (alpha_CM/spacing^3).*(4/3*(k*spacing)^2 + 2/3*(1i + log((pi-k*spacing)/(pi+k*spacing))/pi)*k^3*spacing^3));
