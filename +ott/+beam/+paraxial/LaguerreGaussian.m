classdef LaguerreGaussian < ott.beam.paraxial.Paraxial ...
    & ott.beam.properties.LaguerreGaussian ...
    & ott.beam.utils.VariablePower
% Laguerre-Gaussian beam in the paraxial approximation.
% Inherits from :class:`ott.beam.paraxial.Paraxial` and
% :class:`ott.beam.properties.LaguerreGaussian`.
%
% Solution to paraxial Helmholtz equation in cylindrical coordinates::
%
%   \phi_0 = u_l(r, z) u_p(\phi, z) \exp(-ikz)
%
% where :math:`u_l` and :math:`u_p` are the separable parts of the
% solution.
%
% Properties
%   - waist      -- Beam waist at focus
%   - lmode      -- Azimuthal Laguerre mode order
%   - pmode      -- Radial Laguerre mode order
%   - power      -- Beam power
%
% Inherited properties
%   - permittivity  -- Relative permittivity of medium (default: 1.0)
%   - wavelength    -- Wavelength of beam in medium (default: 1.0)
%   - speed0        -- Speed of light in vacuum (default: 1.0)
%   - omega         -- Optical frequency of light
%   - index_medium  -- Refractive index in medium
%   - wavenumber    -- Wave-number of beam in medium
%   - speed         -- Speed of light in medium
%   - wavelength0   -- Vacuum wavelength of beam

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Static)
    function args = likeProperties(other, args)
      % Construct an array of like-properties
      args = ott.beam.utils.VariablePower.likeProperties(other, args);
      args = ott.beam.properties.LaguerreGaussian.likeProperties(other, args);
    end

    function beam = like(other, varargin)
      % Create a beam like another beam
      %
      % Usage
      %   beam = LaguerreGaussian.like(other, ...)
      %
      % See constructor for arguments.

      args = ott.beam.paraxial.LaguerreGaussian.likeProperties(...
          other, varargin);
      beam = ott.beam.paraxial.LaguerreGaussian(args{:});
    end
  end

  methods
    function beam = LaguerreGaussian(varargin)
      % Construct a new Laguerre-Gaussian beam
      %
      % Usage
      %   beam = LaguerreGaussian(waist, lmode, pmode, ...)
      %   Parameters can also be passed as named arguments.
      %
      % Parameters
      %   - waist (numeric) -- Beam waist
      %   - lmode (integer)     -- Azimuthal LG mode
      %   - pmode (integer > 0) -- Radial LG mode
      %
      % Optional named arguments
      %   - power (numeric) -- Beam power.  Default: 1.0.
      %
      % For properties see :class:`ott.beam.properties.Gaussian`.

      args = ott.utils.addDefaultParameter('power', 1.0, varargin);
      beam = beam@ott.beam.properties.LaguerreGaussian(args{:});
    end
  end

  methods (Hidden)
    function u = uRadial(beam, r, z)
      % Radial part of separable solution

      % Based on Wikipedia page
      % https://en.wikipedia.org/wiki/Gaussian_beam

      zR = pi .* beam.waist.^2 .* beam.medium.index ./ beam.wavelength;
      waistz = beam.waist .* sqrt(1 + (z./zR).^2);

      invRz = z ./ (z.^2 + zR.^2);
      psiz = (abs(beam.lmode) + 2*beam.pmode + 1) .* atan2(z, zR);

      u = beam.waist ./ waistz .* (r .* sqrt(2)./ waistz).^abs(beam.lmode) ...
          .* exp(-r.^2./waistz.^2) ...
          .* laguerreL(beam.pmode, abs(beam.lmode), 2*r.^2./waistz.^2) ...
          .* exp(-1i.*beam.wavenumber.*r.^2./2.*invRz) ...
          .* exp(1i*psiz);
    end

    function u = uAzimuthal(beam, phi, z)
      % Azimuthal part of separable solution

      % Based on Wikipedia page
      % https://en.wikipedia.org/wiki/Gaussian_beam

      u = exp(-1i.*beam.lmode .*phi);
    end

    function E = efieldInternal(beam, xyz)
      % Calculate the electric field

      rho = sqrt(xyz(1, :).^2 + xyz(2, :).^2);
      phi = atan2(xyz(2, :), xyz(1, :));

      % Separable parts of solution
      ul = beam.uRadial(rho, xyz(3, :));
      up = beam.uAzimuthal(phi, xyz(3, :));

      % Normalisation factor
      Clp = sqrt(2*factorial(beam.pmode) ...
          ./(pi*factorial(beam.pmode + abs(beam.lmode))));

      E0 = sqrt(2 .* beam.power .* beam.vacuum.speed ...
          ./ (beam.medium.index .* beam.waist.^2));

      E = Clp .* E0 .* ul .* up .* exp(-1i*beam.wavenumber.*xyz(3, :));

      % Add in power and polarisation
      E = E .* [beam.polarisation(:); 0];

      % Package output
      E = ott.utils.FieldVector(xyz, E, 'Cartesian');

    end
  end
end
