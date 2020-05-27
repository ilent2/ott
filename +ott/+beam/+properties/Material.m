classdef Material < ott.beam.properties.Beam
% Defines material properties for beams.
% Inherits from :class:`ott.beam.properties.Beam`.
%
% Properties
%   - omega       -- Optical frequency of beam.
%   - medium      -- Properties of optical medium.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    omega        % Optical frequency of beam.
    medium       % Properties of optical medium.
  end

  methods (Static)
    function args = likeProperties(other, args)
      if isa(other, 'ott.beam.properties.Material')
        % TODO: Need to check for other arguments supported by constructor
        %   Probably want to clean up constructor so it doesn't have other
        %   arguments.  Not sure how to do that.
        args = ott.utils.addDefaultParameter('omega', other.omega, args);
        args = ott.utils.addDefaultParameter('medium', other.medium, args);
      end
      args = ott.beam.properties.Beam.likeProperties(other, args);
    end
  end

  methods
    function beam = Material(varargin)
      % Construct material properties
      %
      % Usage
      %   beam = beam@Material(...)

      p = inputParser;
      p.addParameter('omega', []);
      p.addParameter('wavelength', []);
      p.addParameter('wavelength0', []);
      p.addParameter('wavenumber', []);
      p.addParameter('wavenumber0', []);
      p.addParameter('vacuum', []);
      p.addParameter('medium', []);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Construct base
      beam = beam@ott.beam.properties.Beam(unmatched{:});

      % Get default values from like
      % TODO: Use 'UsingDefaults' instead of isempty
      default_omega = 2*pi;
      default_vacuum = ott.beam.medium.Vacuum.Unitary();
      default_medium = default_vacuum;

      % Check number of omega related arguments
      % TODO: We often want to supply two of these (to specify medium)
      num_omega = isempty(p.Results.omega) ...
          + isempty(p.Results.wavelength) ...
          + isempty(p.Results.wavenumber) ...
          + isempty(p.Results.wavelength0) ...
          + isempty(p.Results.wavenumber0);
      assert(num_omega >= 4, ...
        'Must only provide one or two of omega|wavelength|wavenumber');

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
      % TODO: We should remove like from Materials too?
      beam.medium = ott.beam.medium.Material.Simple(...
          'vacuum', vacuum, 'like', medium);

      % Calculate optical frequency
      if ~isempty(p.Results.omega)
        beam.omega = p.Results.omega;
      elseif ~isempty(p.Results.wavelength)
        beam = beam.setFrequency('wavelength', p.Results.wavelength);
      elseif ~isempty(p.Results.wavelength0)
        beam = beam.setFrequency('wavelength0', p.Results.wavelength0);
      elseif ~isempty(p.Results.wavenumber)
        beam = beam.setFrequency('wavenumber', p.Results.wavenumber);
      elseif ~isempty(p.Results.wavenumber0)
        beam = beam.setFrequency('wavenumber0', p.Results.wavenumber0);
      else
        beam.omega = default_omega;
      end
    end
  end

  methods % Getters/setters
    function beam = set.omega(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
          'omega must be numeric scalar');
      beam.omega = val;
    end

    function beam = set.medium(beam, val)
      assert(isa(val, 'ott.beam.medium.Medium'), ...
          'medium must be a medium.Medium');
      beam.medium = val;
    end
  end
end
