classdef Vacuum
% Describes properties of the vacuum (i.e., the units for the optical medium).
% The vacuum can have any units but must have SI dimensions.
%
% Properties
%   - permittivity      -- Permittivity of the vacuum
%   - permeability      -- Permeability of the vacuum
%   - speed             -- Speed in vacuum (dependent)
%
% Constant properties
%   - BaseSi      -- Vacuum with base SI units
%   - Unitary     -- Vacuum with unitary properties
%
% Supported casts
%   - Material    -- Constructs a dielectric material (index = 1.0)
%   - Medium      -- Uses Material

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  properties (SetAccess=protected)
    permittivity
    permeability
  end

  properties (Dependent)
    speed
  end

  properties (Constant)

    % Unitary vacuum (permittivity = permeability = speed = 1)
    Unitary = ott.beam.medium.Vacuum(1.0, 1.0);

    % Vacuum with base SI units
    %   * Permittivity = 8.854e-12 F/m
    %   * Permeability = 1.257e-6
    %   * Speed approximately 3.0e8 m/s
    BaseSi = ott.beam.medium.Vacuum(8.854e-12, 1.257e-6);

  end

  methods
    function vac = Vacuum(varargin)
      % Construct a vacuum representation
      %
      % Usage
      %   vacuum = Vacuum(permittivity, permeability)
      %
      % Parameters
      %   - permittivity (numeric) -- Permittivity of vacuum.  Default: 1.0.
      %
      %   - permeability (numeric) -- Permeability of vacuum.  Default: 1.0.

      p = inputParser;
      p.addOptional('permittivity', 1.0, @isnumeric);
      p.addOptional('permeability', 1.0, @isnumeric);
      p.parse(varargin{:});

      vac.permittivity = p.Results.permittivity;
      vac.permeability = p.Results.permeability;
    end

    function mat = ott.beam.medium.Material(~)
      % Construct a vacuum material (dielectric with index = 1.0)
      mat = ott.beam.medium.Dielectric(1.0);
    end

    function med = ott.beam.medium.Medium(vac)
      % Construct a vacuum medium (casts to Medium via Material)
      mat = ott.beam.medium.Material(vac);
      med = ott.beam.medium.Medium(mat);
    end
  end

  methods % Getters/setters
    function val = get.speed(medium)
      val = 1./sqrt(medium.permittivity .* medium.permeability);
    end

    function vac = set.permittivity(vac, val)
      assert(isnumeric(val) && isscalar(val), ...
        'permittivity must be numeric scalar');
      vac.permittivity = val;
    end

    function vac = set.permeability(vac, val)
      assert(isnumeric(val) && isscalar(val), ...
        'permeability must be numeric scalar');
      vac.permeability = val;
    end
  end
end
