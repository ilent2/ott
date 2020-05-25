classdef DdaDipoleForce
% Evaluate the force from individual dipoles
%
% TODO: Should this all be part of a beam class?
%   This isn't really scattering, its more force calculation
%
% Evaluates
%
%   F_rad = F_{inc} + F_{sca} = \sum_i F_{i,inc} + F_{i, sca}
%
% where `sca` and `inc` denote the incident and scattered force
% contributions which are given by
%
%   F_{i, inc} = \frac{1}{2} \Re i k (p^* \cdot E) \exp(i k \cdot r)
%
% and
%
%   F_{i, sca} = \sum_{i\neq j} \frac{1}{2}\Re F_{i,j}
%
% and :math:`F_{i,j}` is a rather long term involving the radiation
% from every other dipole.
%
% This method is efficient for very few dipoles and is useful when
% the force on each dipole is of interest.  The method id described
% in
%
%   Radiation forces in the discrete-dipole approximation
%   A. G. Hoekstra, M. Frijlink, L. B. F. M. Waters, and P. M. A. Sloot
%   JOSA A, Vol. 18, Issue 8, pp. 1944-1953 (2001)
%   https://doi.org/10.1364/JOSAA.18.001944

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function F = calculateDipoleInc()

      F = 0.5 * real(1i.*k.*sum(conj(p) .* E, 1) .* exp(1i*k.*r));

    end

    function F = calculateDipoleSca()

      pis = ...;   % dipole
      pjs = ...;   % Target

      r_vec = dipole_xyz - targets_xyz;
      r_ij = vecnorm(r_vec);
      n_ij = r_vec ./ r_ij;
      k = 2*pi;


      Fij = exp(1i*k*r_ij) .* ((...
          (conj(pis) .* pjs).*n + conj(pis).*(n .* pjs) ...
          + (conj(pis) .* n).*pjs - 5.*(conj(pis).*n).*n.*(n .* pjs)) ...
          .* (-k^2/r_ij.^2 - (3i*k)./r_ij.^3 + 3./r_ij.^4) ...
          + ((conj(pis).*pjs).*n - (conj(pis).*n).*n.*(n.*pjs)) ...
          .* ((1i.*k^3)./r_ij - k.^2./r_ij.^2));

    end
  end
end
