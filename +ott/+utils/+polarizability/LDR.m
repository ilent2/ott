function alpha_LDR = LDR(spacing,index,varargin)
% Lattice dispersion relation polarizablity
%
% Polarizability calculation based on
%
%   Draine & Goodman, Beyond Clausius-Mossoti: wave propagation
%   on a polarizable point lattice and the discrete dipole approximation,
%   The Astrophysical Journal, 405:685-697, 1993 March 10
%
% Usage
%   alpha = LDR(spacing, index, ...)
%   Calculates a Nx1 element vector containing the isotropic
%   polarisabilities for N dipoles.
%
%   alpha = LDR(spacing, index, kvec, E0, ...)
%   As above but specifies the polarisability information for use
%   with plane wave illumination.
%
% Parameters
%   - spacing (numeric scalar) -- lattice spacing parameter
%   - index (Nx1 numeric) -- Relative refractive indices for N dipoles.
%   - kvec (1x3 numeric) -- Wave vector [kx, ky, kz]
%   - E0 (1x3 numeric) -- E-field polarisation [Ex, Ey, Ez]
%
% Optional named arguments
%   - k0 (numeric) -- Wavenumber to scale spacing by.  Default: ``2*pi``.

% Based on the script by Vincent Loke.
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

p = inputParser;
p.addOptional('kvec', [], @isnumeric);
p.addOptional('E0', [], @isnumeric);
p.addParameter('k0', 2*pi);
p.parse(varargin{:});

k0 = p.Results.k0;
b1 = -1.8915316;
b2 = 0.1648469;
b3 = -1.7700004;
msqr = index(:).^2;
dcube = spacing^3;

if nargin >= 4  && ~isempty(p.Results.kvec) % we have polarization info   
  
  kvec = p.Results.kvec;
  E0 = p.Results.E0;
  
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
alpha_LDR = alpha_CM./(1 + (alpha_CM/dcube).*((b1+msqr*b2+msqr*b3*S)*(k0*spacing)^2-2/3*1i*k0^3*dcube));

