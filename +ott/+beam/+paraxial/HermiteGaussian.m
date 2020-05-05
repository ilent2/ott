classdef HermiteGaussian < ott.beam.paraxial.Paraxial ...
    & ott.beam.abstract.HermiteGaussian ...
    & ott.beam.utils.Vector2Polarisation
% Hermite-Gaussian beam in the paraxial approximation
%
% Solution to paraxial Helmholtz equation in Cartesian coordinates::
%
%   \phi_0 = u_m(x, z) u_n(y, z) \exp(-ikz)
%
% Properties
%   - waist      -- Beam waist at focus
%   - mmode      -- Hermite mode order
%   - nmode      -- Hermite mode order

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Hidden)
    function u = uInternal(beam, x, z, mode)
      % Calculate the separable part of the solution

      % From the wikipedia page
      % https://en.wikipedia.org/wiki/Gaussian_beam

      zR = pi .* beam.waist.^2 .* beam.index_medium ./ beam.wavelength;

      q0 = 1i .* zR;
      qz = z + q0;

      waistz = beam.waist .* sqrt(1 + (z./zR).^2);

      u = sqrt(1./(2.^mode .* factorial(mode))) ...
          .* sqrt(q0./qz) .* (-conj(qz)./qz).^(mode./2) ...
          .* hermiteH(mode, sqrt(2).*x./waistz) ...
          .* exp(-i .* beam.wavenumber .* x.^2 ./ 2 ./ qz);

    end

    function E = efieldInternal(beam, xyz)
      % Calculate the electric field

      um = beam.uInternal(xyz(1, :), xyz(3, :), beam.mmode);
      un = beam.uInternal(xyz(2, :), xyz(3, :), beam.nmode);

      E = um .* un .* exp(-1i*beam.wavenumber.*xyz(3, :));

      E0 = sqrt(4 .* beam.power .* beam.speed0 ...
          ./ (pi .* beam.index_medium .* beam.waist.^2));

      % Add in power and polarisation
      E = E .* E0 .* [beam.polarisation(:); 0];

      % Package output
      E = ott.utils.FieldVector(xyz, E, 'Cartesian');

    end
  end

  methods
    function beam = HermiteGaussian(waist, mmode, nmode, varargin)
      % Construct a new Hermite-Gaussian beam
      %
      % Usage
      %   beam = HermiteGaussian(waist, mmode, nmode, ...)
      %
      % Parameters
      %   - waist (numeric) -- Beam waist
      %   - mmode (integer) -- Mode number
      %   - nmode (integer) -- Mode number
      %
      % Optional named arguments
      %   - polarisation (2 numeric) -- Polarisation.  Default: [1,0].
      %   - power (numeric) -- Beam power.  Default: 1.0.
      %
      % For optional parameters, see :class:`Properties`.

      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('polarisation', [1;0]);
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Call base constructor
      beam = beam@ott.beam.utils.Vector2Polarisation(...
          p.Results.polarisation);
      beam = beam@ott.beam.abstract.HermiteGaussian(waist, ...
          mmode, nmode, unmatched{:});
    end
  end
end

