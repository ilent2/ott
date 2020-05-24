classdef GaussianDavis5 < ott.beam.Beam ...
    & ott.beam.properties.Gaussian ...
    & ott.beam.utils.VariablePower
% Fifth-order Davis approximation of a Gaussian beam
% Inherits from :class:`Beam` and :class:`Gaussian`.
%
% A Davis beam is a approximate solution to the vector Helmholtz equation
% with a power series expansion of the beam width parameter.
% This implementation uses the derivation in
%
%   J. P. Barton and D. R. Alexander,
%   Journal of Applied Physics 66, 2800 (1989)
%   https://doi.org/10.1063/1.344207
%
% Units of the fields depend on units used for the properties.
%
% Methods
%   - efield        -- Calculate the electric field
%   - hfield        -- Calculate the magnetic field
%   - ehfield       -- Calculate the electric and magnetic fields
%
% Properties
%   - waist         -- Beam waist radius
%   - permittivity  -- Relative permittivity of medium (default: 1.0)
%   - wavelength    -- Wavelength of beam in medium (default: 1.0)
%   - speed0        -- Speed of light in vacuum (default: 1.0)
%   - power         -- Beam power (default: 1.0)
%   - polarisation  -- Beam polarisation (default: [1; 0])
%
% Properties (Dependent)
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
      args = ott.beam.properties.Gaussian.likeProperties(other, args);
    end

    function beam = like(other, varargin)
      % Create a beam like another beam
      %
      % Usage
      %   beam = GaussianDavis5.like(other, ...)
      %
      % See constructor for arguments.

      args = ott.beam.abstract.Gaussian.likeProperties(other, varargin);
      beam = ott.beam.GaussianDavis5(args{:});
    end
  end

  methods
    function beam = GaussianDavis5(varargin)
      % Construct a new fifth-order Davis Gaussian beam approximation
      %
      % Usage
      %   beam = GaussianDavis(waist, ...)
      %   Parameters can also be passed as named arguments.
      %
      % Optional named arguments
      %   - power (numeric) -- Beam power, default 1.0.
      %
      % For properties see :class:`ott.beam.properties.Gaussian`.

      args = ott.utils.addDefaultParameter('power', 1.0, varargin);
      beam = beam@ott.beam.properties.Gaussian(args{:});
    end
  end

  methods (Hidden)
    function [E, H] = ehfieldInternal(beam, xyz)
      % Calculate the E and H fields

      assert(isnumeric(xyz) && size(xyz, 1) == 3, ...
          'xyz must be 3xN numeric');

      % Calculate normalised beam waist (Eq 1)
      s = 1 ./ (beam.wavenumber .* beam.waist);

      % Calculate normalised radial directions
      x = xyz(1, :) ./ beam.waist;
      y = xyz(2, :) ./ beam.waist;
      rho2 = x.^2 + y.^2;

      % Calculate normalised axial direction
      z = xyz(3, :) ./ (beam.wavenumber .* beam.waist.^2);

      % Calculate the paraxial expansion term (Eq 17)
      Q = 1 ./ (1i + 2.0.*z);
      psi0 = 1i .* Q .* exp(-1i .* rho2 .* Q);

      polarisation = [beam.polarisation(:); 0];
      A0 = psi0 .* exp(-1i .* z ./ s.^2);

      E0 = sqrt(4 .* beam.power .* beam.vacuum.speed ./ (pi ...
          .* beam.medium.index .* beam.waist.^2 .* (1 + s.^2 + 1.5.*s.^4)));

      Ax = @(x) (1 ...
        + s.^2 .* (-rho2.*Q.^2 + 1i.*rho2.^2.*Q.^3 - 2.*Q.^2.*x.^2) ...
        + s.^4 .* ( 2.*rho2.^2.*Q.^4 - 3i.*rho2.^3.*Q.^5 ...
        - 0.5.*rho2.^4.*Q.^6 ...
        + (8.*rho2.*Q.^4 - 2i.*rho2.^2.*Q.^5).*x.^2)) .* E0 .* A0;
      Ey = (s.^2.*(-2.*Q.^2.*x.*y) ...
        + s.^4.*(8.*rho2.*Q.^4 - 2i.*rho2.^2.*Q.^5).*x.*y) .* E0.*A0;
      Az = @(x) (s.*(-2.*Q.*x) ...
        + s.^3.*(6.*rho2.*Q.^3 - 2i.*rho2.^2.*Q.^4).*x ...
        + s.^5.*(-20.*rho2.^2.*Q.^5 + 10i.*rho2.^3.*Q.^6 ...
        + rho2.^4.*Q.^7).*x) .* E0.*A0;

      Ex = Ax(x);
      Ez = Az(x);

      Hx = beam.medium.index .* Ey;
      Hy = beam.medium.index .* Ax(y);
      Hz = beam.medium.index .* Az(y);

      E = [Ex; Ey; Ez];
      H = [Hx; Hy; Hz];

      % Package output
      E = ott.utils.FieldVector(xyz, E, 'cartesian');
      H = ott.utils.FieldVector(xyz, H, 'cartesian');
    end

    function E = efieldInternal(beam, xyz)
      % Calculate only E-field
      [E, ~] = beam.ehfieldInternal(xyz);
    end

    function H = hfieldInternal(beam, xyz)
      % Calculate only H-field
      [~, H] = beam.ehfieldInternal(xyz);
    end

    function E = efarfieldInternal(beam)
      error('Not yet implemented');
    end
    function H = hfarfieldInternal(beam)
      error('Not yet implemented');
    end
  end
end
