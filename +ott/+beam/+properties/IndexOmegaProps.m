classdef IndexOmegaProps
% Declares index_medium and omega properties for a beam.
%
% Properties
%   - index_medium    -- Medium refractive index
%   - omega           -- Optical speed in medium

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    index_medium   % Refractive index of the medium
    omega          % Optical angular frequency of light [1/s]
  end

  methods (Static)
    function [omega, index, unmatched] = parseArgs(varargin)
      % Parse inputs to get omega/index and unmatched
      %
      % Usage
      %   [omega, index, unmatched] = parseArgs(...)
      %
      % Optional named arguments
      %   - index_medium (numeric) -- Refractive index of the medium.
      %     Default: ``1.0``.
      %
      %   - omega (numeric) -- Optical angular frequency [Hz].
      %     Default: ``3e8/1064e-9*2*pi`` (i.e., default vacuum
      %     wavelength is 1064 nm).  Only one of omega and wavelength0
      %     should be set.
      %
      %   - wavelength0 (numeric) -- Wavelength in vacuum [m].
      %     Default: ``1064e-9``.  Only one of omega and wavelength0
      %     should be set.

      p = inputParser;
      p.addParameter('index_medium', 1.0);
      p.addParameter('omega', []);
      p.addParameter('scale', 1.0);
      p.addParameter('wavelength0', []);
      p.addParameter('wavelength', [], ...
        @(~) error('wavelength not supported, did you mean wavelength0?'));
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      speed0 = ott.beam.Beam.speed0;

      index = p.Results.index_medium;

      % Set omega from omega or wavelength0
      if isempty(p.Results.omega) && isempty(p.Results.wavelength0)
        wavelength0 = 1064e-9;
        omega = 2*pi*speed0/wavelength0;
      elseif ~isempty(p.Results.omega)
        omega = p.Results.omega;
      elseif ~isempty(p.Results.wavelength0)
        omega = 2*pi*speed0/p.Results.wavelength0;
      else
        error('Only one of omega and wavelength0 should be set');
      end
    end
  end

  methods
    function [beam, unmatched] = IndexOmegaProps(omega, index_medium)
      % Construct beam index_medium/omega props
      %
      % Usage
      %   beam = Beam(omega, index_medium)

      beam.omega = omega;
      beam.index_medium = index_medium;
    end
  end

  methods % Getters/setters
    function beam = set.index_medium(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
          'index_medium must be numeric scalar');
      beam.index_medium = val;
    end

    function beam = set.omega(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
          'omega must be numeric scalar');
      beam.omega = val;
    end
  end
end
