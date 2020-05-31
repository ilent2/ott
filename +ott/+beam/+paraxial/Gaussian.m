classdef Gaussian < ott.beam.paraxial.Paraxial ...
    & ott.beam.properties.Gaussian ...
    & ott.beam.properties.VariablePower
% Paraxial approximation of a Gaussian beam
% Inherits from :class:`Paraxial` and :class:`ott.beam.properties.Gaussian`.
%
% The paraxial beam is described by the potential::
%
%   \phi_0 = i Q \exp(-i\rho Q)
%
% where :math:`Q` and :math:`\rho` are related to the axial and radial
% positions normalised by the beam waist and diffraction length.  For
% further details see
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
      args = ott.beam.properties.VariablePower.likeProperties(other, args);
      args = ott.beam.properties.Gaussian.likeProperties(other, args);
    end

    function beam = like(other, varargin)
      % Create a beam like another beam
      %
      % Usage
      %   beam = Gaussian.like(other, ...)
      %
      % See constructor for arguments.

      args = ott.beam.paraxial.Gaussian.likeProperties(other, varargin);
      beam = ott.beam.paraxial.Gaussian(args{:});
    end
  end

  methods
    function beam = Gaussian(varargin)
      % Construct a new Gaussian paraxial beam representation
      %
      % Usage
      %   beam = paraxial.Gaussian(waist, ...)
      %   Parameters can also be passed as named arguments.
      %
      % Parameters
      %   - waist (numeric) -- Beam waist [L]
      %
      % Optional named arguments
      %   - power (numeric) -- Beam power.  Default: ``1.0``.
      %
      % For properties see :class:`ott.beam.properties.Gaussian`.

      args = ott.utils.addDefaultParameter('power', 1.0, varargin);
      beam = beam@ott.beam.properties.Gaussian(args{:});
    end
  end

  methods (Hidden)
    function E = efieldInternal(beam, xyz)
      % Calculate the E (and H) fields

      % Calculate normalised beam waist (Eq 1)
      s = 1 ./ (beam.wavenumber .* beam.waist);

      % Calculate normalised radial directions
      rho2 = sum((xyz(1:2, :)./beam.waist).^2, 1);

      % Calculate normalised axial direction
      z = xyz(3, :) ./ (beam.wavenumber .* beam.waist.^2);

      % Calculate the paraxial expansion term (Eq 17)
      Q = 1 ./ (1i + 2.0.*z);
      psi0 = 1i .* Q .* exp(-1i .* rho2 .* Q);

      % Calculate the scalar potential (Eq 11)
      polarisation = [beam.polarisation(:); 0];
      A = polarisation .* psi0 .* exp(-1i .* z ./ s.^2);

      E0 = sqrt(4 .* beam.power .* beam.vacuum.speed ...
          ./ (pi .* beam.medium.index .* beam.waist.^2));

      E = E0 .* A;

      % Package output
      E = ott.utils.FieldVector(xyz, E, 'cartesian');
    end
  end
end
