classdef (Abstract) Properties
% A base class for Beam and abstract.Beam representations.
%
% This class defines the common properties and methods to these
% two classes.
%
% Properties
%   - power         -- The power of the beam (may be infinite)
%   - omega         -- Beam optical frequency
%   - medium        -- Medium where beam is propagating
%
% Dependent properties
%   - vacuum        -- Vacuum material linked to the Medium
%   - wavelength    -- Wave-length of the beam in medium
%   - wavenumber    -- Wave-number of beam in medium
%   - wavelength0   -- Wave-length of the beam in vacuum
%   - wavenumber0   -- Wave-number of beam in vacuum
%
% Abstract methods
%   - getBeamPower      -- get method called by dependent property power

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    omega            % Beam optical frequency
    medium           % Medium where beam is propagating
  end

  properties (Dependent)
    power           % The power of the beam (may be infinite)
    vacuum          % Vacuum material liked to the medium

    wavenumber      % Wave-number of beam in medium
    wavenumber0     % Wave-number of beam in vacuum
    wavelength      % Wavelength of beam
    wavelength0     % Vacuum wavelength of beam
  end

  methods (Abstract)
    getBeamPower    % get method called by dependent property power
  end

  methods
    function beam = Properties(varargin)
      % Initialize beam properties
      %
      % Usage
      %   beam = beam@ott.beam.Properties(...)
      %
      % Named arguments
      %   - wavelength|wavelength0 (numeric) -- Beam wavelength
      %   - wavenumber|wavenumber0 (numeric) -- Beam wave-number
      %   - omega (numeric) -- Optical frequency of beam.
      %     Default: ``2*pi``
      %
      %   - vacuum (ott.beam.medium.Medium) -- Vacuum medium.
      %     Default: ``ott.beam.medium.Vacuum.Unitary``.
      %     If set, `medium` must not be Material.
      %
      %   - medium (ott.beam.medium.Medium|Material) -- Medium or
      %     material describing optical properties of medium.
      %     Default: ``ott.beam.medium.Vacuum.Unitary``.
      %
      %   - like (ott.beam.Properties) -- Uses another beam
      %     group for default parameters.
      %
      % The medium is constructed from
      % :meth:`ott.beam.medium.Material.Simple`.  All unmatched parameters,
      % along with `vacuum` and `medium` are passed to this function.
      %
      % Only one of omega, wavelength or wave-number should be set.

      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('omega', []);
      p.addParameter('wavelength', []);
      p.addParameter('wavelength0', []);
      p.addParameter('wavenumber', []);
      p.addParameter('wavenumber0', []);
      p.addParameter('vacuum', []);
      p.addParameter('medium', []);
      p.addParameter('like', []);
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Get default values from like
      default_omega = 2*pi;
      default_vacuum = ott.beam.medium.Vacuum.Unitary();
      default_medium = default_vacuum;
      if ~isempty(p.Results.like)
        default_omega = p.Results.like.omega;
        default_vacuum = p.Results.like.vacuum;
        default_medium = p.Results.like.medium;
      end

      % Check number of omega related arguments
      num_omega = isempty(p.Results.omega) ...
          + isempty(p.Results.wavelength) ...
          + isempty(p.Results.wavenumber) ...
          + isempty(p.Results.wavelength0) ...
          + isempty(p.Results.wavenumber0);
      assert(num_omega >= 4, ...
        'Must only provide one of omega|wavelength|wavenumber');

      % Get default values for medium construction
      medium = p.Results.medium;
      vacuum = p.Results.vacuum;
      if ~isempty(vacuum) && ~isempty(medium)
        assert(~isa(medium, 'ott.beam.medium.Material'), ...
            'medium must not be a material if vacuum is specified');
      elseif isempty(vacuum) && isempty(medium)
        medium = default_medium;
        vacuum = default_vacuum;
      elseif isempty(medium)
        medium = default_medium;
      elseif isempty(vacuum)
        if isa(medium, 'ott.beam.medium.Material')
          vacuum = medium.vacuum;
        else
          vacuum = default_vacuum;
        end
      end

      % Construct medium
      beam.medium = ott.beam.medium.Material.Simple(...
          'vacuum', vacuum, 'like', medium, unmatched{:});

      % Calculate optical frequency
      if ~isempty(p.Results.omega)
        beam.omega = p.Results.omega;
      elseif ~isempty(p.Results.wavelength)
        beam.setFrequency('wavelength', p.Results.wavelength);
      elseif ~isempty(p.Results.wavelength0)
        beam.setFrequency('wavelength0', p.Results.wavelength0);
      elseif ~isempty(p.Results.wavenumber)
        beam.setFrequency('wavenumber', p.Results.wavenumber);
      elseif ~isempty(p.Results.wavenumber0)
        beam.setFrequency('wavenumber0', p.Results.wavenumber0);
      else
        beam.omega = 2*pi;
      end
    end

    function beam = setFrequency(beam, name, val)
      % Set the optical frequency from another named parameter
      %
      % Usage
      %   beam = beam.setFrequency(name, value)
      %
      % Parameters
      %   - name (enum) -- One of the named properties supported by
      %     the constructor.
      %
      %   - value (numeric) -- Value to use in calculation

      switch name
        case 'wavelength'
          beam.omega = 2*pi*beam.medium.speed./val;
        case 'wavelength0'
          beam.omega = 2*pi*beam.vacuum.speed./val;
        case 'wavenumber'
          beam.omega = beam.medium.speed.*val;
        case 'wavenumber0'
          beam.omega = beam.vacuum.speed.*val;
        otherwise
          error('Unknown parameter name');
      end
    end

    function data = arrayApply(beam, data, func)
      % Apply function to each array in the beam array output.
      %
      % This function is overloaded by Array types in order to
      % implement incoherent combination.

      data = func(data);
    end
  end

  methods (Hidden)
    function beam = setBeamPower(beam, val)
      % Function to set the beam power (if supported)
      % Override this function if your beam supports this feature
      error('Setting beam power not supported');
    end
  end

  methods % Getters/setters
    function val = get.power(beam)
      val = beam.getBeamPower();
    end
    function beam = set.power(beam, val)
      beam = beam.setBeamPower(val);
    end

    function val = get.vacuum(beam)
      val = beam.medium.vacuum;
    end

    function beam = set.omega(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
          'omega must be numeric scalar');
      beam.omega = val;
    end

    function val = get.wavenumber(beam)
      val = beam.omega ./ beam.medium.speed;
    end
    function val = get.wavenumber0(beam)
      val = beam.omega ./ beam.vecuum.speed;
    end

    function val = get.wavelength(beam)
      val = 2*pi.*beam.medium.speed ./ beam.omega;
    end
    function val = get.wavelength0(beam)
      val = 2*pi.*beam.vacuum.speed ./ beam.omega;
    end
  end
end
