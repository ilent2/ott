classdef Vacuum < ott.beam.medium.Medium
% Vacuum medium base class.
% Inherits from :class:`Medium`.
%
% Properties
%   - permittivity      -- Permittivity of the vacuum
%   - permeability      -- Permeability of the vacuum
%   - speed             -- Speed in vacuum (dependent)
%   - index             -- Refractive index (constant, 1.0)
%
% Static methods
%   - SiUnits     -- Vacuum with SI units
%   - Unitary     -- Vacuum with unitary properties

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  properties
    permittivity
    permeability
  end

  properties (Dependent)
    speed
    index
  end

  methods (Static)
    function vac = SiUnits()
      % Construct a vacuum with SI units
      %
      % Parameters are approximate
      %   * Permittivity = 8.854e-12 F/m
      %   * Permeability = 1.257e-6
      %   * Speed approximately 3.0e8 m/s
      %
      % Usage
      %   vacuum = Vacuum.SiUnits()

      vac = ott.beam.medium.Vacuum('permittivity', 8.854e-12, ...
          'permeability', 1.257e-6);
    end

    function vac = Unitary()
      % Construct a vacuum with unitary units
      %
      % Parameters are approximate
      %   * Permittivity = 1.0
      %   * Permeability = 1.0
      %   * Speed = 1.0
      %
      % Usage
      %   vacuum = Vacuum.Unitary()

      vac = ott.beam.medium.Vacuum('permittivity', 1.0, ...
          'permeability', 1.0);
    end

    function vac = Simple(varargin)
      % Construct a new vacuum medium from a list of properties
      %
      % Usage
      %   vacuum = Vacuum.Simple(...)
      %
      % Named parameters
      %   - permittivity0 -- Vacuum permittivity
      %   - permeability0 -- Vacuum permeability
      %   - speed0 -- Wave speed in vacuum.
      %
      %   - like (ott.beam.medium.Medium) -- Medium to use for default
      %     property values.  Default: ``Unitary``.
      %
      % Supports partial name matching.
      %
      % Maximum two of permittivity, permeability and speed can be supplied.

      p = inputParser;
      p.addParameter('permittivity0', []);
      p.addParameter('permeability0', []);
      p.addParameter('speed0', []);
      p.addParameter('like', ott.beam.medium.Vacuum.Unitary);
      p.parse(varargin{:});

      non_args = isempty(p.Results.permittivity0) ...
          + isempty(p.Results.permeability0) ...
          + isempty(p.Results.speed0);
      assert(non_args ~= 0, 'Maximum two property arguments supported');

      if non_args == 1
        % No defaults
        if isempty(p.Results.permittivity0)
          vac = ott.beam.medium.Vacuum('speed', p.Results.speed0, ...
              'permeability', p.Results.permeability0);
        elseif isempty(p.Results.permeability0)
          vac = ott.beam.medium.Vacuum('speed', p.Results.speed0, ...
              'permittivity', p.Results.permittivity0);
        else
          vac = ott.beam.medium.Vacuum(...
              'permeability', p.Results.permeability0, ...
              'permittivity', p.Results.permittivity0);
        end
      else
        % One or two defaults
        if ~isempty(p.Results.permittivity0)
          vac = ott.beam.medium.Vacuum('speed', p.Results.like.speed0, ...
              'permittivity', p.Results.permittivity0);
        elseif ~isempty(p.Results.permeability0)
          vac = ott.beam.medium.Vacuum('speed', p.Results.like.speed0, ...
              'permeability', p.Results.permeability0);
        elseif ~isempty(p.Results.speed0)
          vac = ott.beam.medium.Vacuum('speed', p.Results.speed0, ...
              'permittivity', p.Results.like.permittivity0);
        else
          % Two defaults
          vac = ott.beam.medium.Vacuum(...
              'permeability', p.Results.like.permeability, ...
              'permittivity', p.Results.like.permittivity);
        end
      end
    end
  end

  methods
    function vac = Vacuum(varargin)
      % Construct a vacuum representation
      %
      % Usage
      %   vacuum = Vacuum(...)
      %
      % Named parameters
      %   - permittivity (numeric) -- Permittivity of vacuum.
      %     Default: ``[]``.
      %
      %   - permeability (numeric) -- Permeability of vacuum.
      %     Default: ``[]``.
      %
      %   - speed (numeric) -- Speed of wave in vacuum.
      %     Default: ``[]`` (computed from permittivity/permeability).
      %
      % Two parameters must be supplied.

      p = inputParser;
      p.addParameter('permittivity', []);
      p.addParameter('permeability', []);
      p.addParameter('speed', []);
      p.parse(varargin{:});

      assert(numel(p.UsingDefaults) == 1, ...
        'ott:beam:medium:Vacuum:wrong_arg_count', ...
        'Must supply two parameters');

      if isempty(p.Results.speed)
        vac.permittivity = p.Results.permittivity;
        vac.permeability = p.Results.permeability;
      elseif isempty(p.Results.permittivity)
        vac.permeability = p.Rseults.permeability;
        vac.permittivity = 1 ./ p.Results.speed.^2 ./ vac.permeability;
      elseif isempty(p.Results.permeability)
        vac.permittivity = p.Rseults.permittivity;
        vac.permeability = 1 ./ p.Results.speed.^2 ./ vac.permittivity;
      end
    end
  end

  methods % Getters/setters
    function val = get.speed(medium)
      val = 1./sqrt(medium.permittivity .* medium.permeability);
    end

    function val = get.index(medium)
      val = 1.0;
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
