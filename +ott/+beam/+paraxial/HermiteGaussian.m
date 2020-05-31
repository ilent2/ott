classdef HermiteGaussian < ott.beam.paraxial.Paraxial ...
    & ott.beam.properties.HermiteGaussian ...
    & ott.beam.properties.VariablePower
% Hermite-Gaussian beam in the paraxial approximation.
% Inherits from :class:`Paraxial`
% and :class:`ott.beam.properties.HermiteGaussian`.
%
% Solution to paraxial Helmholtz equation in Cartesian coordinates::
%
%   \phi_0 = u_m(x, z) u_n(y, z) \exp(-ikz)
%
% Properties
%   - waist      -- Beam waist at focus
%   - mmode      -- Hermite mode order
%   - nmode      -- Hermite mode order
%   - power      -- Beam power
%
% See :class:`ott.beam.Beam` for additional methods/properties.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Static)
    function args = likeProperties(other, args)
      % Construct an array of like-properties
      args = ott.beam.properties.VariablePower.likeProperties(other, args);
      args = ott.beam.properties.HermiteGaussian.likeProperties(other, args);
    end

    function beam = like(other, varargin)
      % Create a beam like another beam
      %
      % Usage
      %   beam = HermiteGaussian.like(other, ...)
      %
      % See constructor for arguments.

      args = ott.beam.paraxial.HermiteGaussian.likeProperties(...
          other, varargin);
      beam = ott.beam.paraxial.HermiteGaussian(args{:});
    end
  end

  methods
    function beam = HermiteGaussian(varargin)
      % Construct a new Hermite-Gaussian beam
      %
      % Usage
      %   beam = HermiteGaussian(waist, mmode, nmode, ...)
      %   Parameters can also be passed as named arguments.
      %
      % Parameters
      %   - waist (numeric) -- Beam waist
      %   - mmode (integer) -- Mode number
      %   - nmode (integer) -- Mode number
      %
      % Optional named arguments
      %   - power (numeric) -- Beam power.  Default: 1.0.
      %
      % For properties see :class:`ott.beam.properties.HermiteGaussian`.

      args = ott.utils.addDefaultParameter('power', 1.0, varargin);
      beam = beam@ott.beam.properties.HermiteGaussian(args{:});
    end
  end

  methods (Hidden)
    function u = uInternal(beam, x, z, mode)
      % Calculate the separable part of the solution

      % From the wikipedia page
      % https://en.wikipedia.org/wiki/Gaussian_beam

      zR = pi .* beam.waist.^2 .* beam.medium.index ./ beam.wavelength;

      q0 = 1i .* zR;
      qz = z + q0;

      waistz = beam.waist .* sqrt(1 + (z./zR).^2);

      u = sqrt(1./(2.^mode .* factorial(mode))) ...
          .* sqrt(q0./qz) .* (-conj(qz)./qz).^(mode./2) ...
          .* hermiteH(mode, sqrt(2).*x./waistz) ...
          .* exp(-1i .* beam.wavenumber .* x.^2 ./ 2 ./ qz);

    end

    function E = efieldInternal(beam, xyz)
      % Calculate the electric field

      um = beam.uInternal(xyz(1, :), xyz(3, :), beam.mmode);
      un = beam.uInternal(xyz(2, :), xyz(3, :), beam.nmode);

      E = um .* un .* exp(-1i*beam.wavenumber.*xyz(3, :));

      E0 = sqrt(4 .* beam.power .* beam.vacuum.speed ...
          ./ (pi .* beam.medium.index .* beam.waist.^2));

      % Add in power and polarisation
      E = E .* E0 .* [beam.polarisation(:); 0];

      % Package output
      E = ott.utils.FieldVector(xyz, E, 'Cartesian');

    end
  end
end
