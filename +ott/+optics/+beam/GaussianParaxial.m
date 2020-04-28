classdef GaussianParaxial < ott.optics.beam.Beam & ott.optics.beam.Gaussian
% Paraxial approximation of a Gaussian beam
% Inherits from :class:`Beam` and :class:`Gaussian`.
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

  properties
    polarisation      % x and y polarisation
  end

  properties (Hidden)
    internalPower     % Actual beam power value
  end

  methods (Hidden)
    function [E, H] = ehfieldInternal(beam, xyz)
      % Calculate the E and H fields

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

      E0 = sqrt(16 .* beam.power ./ (beam.index_medium ...
          .* beam.speed0 .* beam.waist.^2 .* (1 + s.^2 + 1.5.*s.^4)));

      E = E0 .* A;
      H = E0 .* beam.index_medium .* A;

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

    function val = getBeamPower(beam)
      % Get the internal power value
      val = beam.internalPower;
    end
    function beam = setBeamPower(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'power must be numeric scalar');
      beam.internalPower = val;
    end
  end

  methods
    function beam = GaussianParaxial(waist, varargin)
      % Construct a new Gaussian paraxial beam representation
      %
      % Usage
      %   beam = GaussianParaxial(waist, ...)
      %
      % Parameters
      %   - waist (numeric) -- Beam waist [L]
      %
      % Optional named arguments
      %   - permittivity (numeric) -- Relative permittivity of medium
      %   - wavelength (numeric) -- Wavelength in medium [L]
      %   - speed0 (numeric) -- Speed of light in vacuum [L/T]
      %   - power (numeric) -- Beam power.  Default: ``1.0``.
      %   - polarisation (2x1 numeric) -- Polarisation.  Default: ``[1;0]``.
      %
      %   - omega (numeric) -- Optical frequency [2*pi/T]
      %   - index_medium (numeric) -- Refractive index in medium
      %   - wavenumber (numeric) -- Wave-number in medium [2*pi/L]
      %   - speed (numeric) -- Speed of light in medium [L/T]
      %   - wavelength0 (numeric) -- Wavelength in medium [L]

      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('power', 1.0);
      p.addParameter('polarisation', [1;0]);
      p.parse(varargin{:});

      % Call base constructor
      unmatched = [fieldnames(p.Unmatched).'; struct2cell(p.Unmatched).'];
      beam = beam@ott.optics.beam.Gaussian(waist, varargin{:});

      beam.power = p.Results.power;
      beam.polarisation = p.Results.polarisation;
    end
  end

  methods
    function beam = set.polarisation(beam, val)
      assert(isnumeric(val) && all(size(val) == [2, 1]), ...
        'polarisation must be 2x1 vector');
      beam.polarisation = val;
    end
  end
end

