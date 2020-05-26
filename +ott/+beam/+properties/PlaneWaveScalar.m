classdef (Abstract) PlaneWaveScalar < ott.beam.properties.PlaneWave
% Properties for a scalar plane wave.
% Inherits from :class:`PlaneWave`.
%
% Properties
%   - field         -- Field parallel and perpendicular to polarisation
%   - position      -- (Inherited) Beam position
%   - rotation      -- (Inherited) Beam rotation
%
% Dependent properties
%   - direction     -- Beam direction
%   - origin        -- Position used to calculate beam phase offset
%   - polarisation1 -- Primary polarisation direction
%   - polarisation2 -- Secondary polarisation direction
%   - intensity     -- Intensity of the beam
%   - wavevector    -- (Inherited) Beam wave-vector
%
% Methods
%   - setData     -- Set field
%   - rotate*     -- Rotate the particle around the X,Y,Z axis
%   - translate*  -- Apply a phase offset (translation) to the beam

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (SetAccess=protected)
    field             % Field parallel and perpendicular to polarisation
  end

  properties (Dependent, SetAccess=protected)
    origin            % Origin for phase calculation
    direction         % Direction of the wave
    polarisation1     % Primary Polarisation direction
    polarisation2     % Secondary Polarisation direction
    intensity         % Intensity of the beam
  end

  methods (Static)
    function args = likeProperties(other, args)
      % Add like-properties to argument list
      if isa(other, 'ott.beam.PlaneWaveScalar')
        args = ott.utils.addDefaultParameter('field', other.field, args);
      end
      args = ott.beam.properties.PlaneWave.likeProperties(other, args);
    end
  end

  methods
    function beam = PlaneWaveScalar(varargin)
      % Construct plane wave properties
      %
      % Usage
      %   beam = PlaneWave(...)
      %
      % Optional named arguments
      %   - field (2 numeric) -- Field parallel and perpendicular to
      %     plane wave polarisation direction.
      %     Default: ``[1, 1i]``.

      p = inputParser;
      p.addParameter('field', [1, 1i]);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.properties.PlaneWave(unmatched{:});
      beam.field = p.Results.field;
    end
    
    function beam = setData(field)
      % Set field data
      %
      % Usage
      %   beam = beam.setData(field)
      %
      % Parameters
      %   - field (2 numeric) -- Field parallel and perpendicular to
      %     plane wave polarisation direction.
      %     Default: ``[1, 1i]``.
      
      beam.field = field;
    end
  end

  methods % Getters/setters
    % Properties
    %   - field         -- Field parallel and perpendicular to polarisation
    %   - direction     -- Beam direction
    %   - origin        -- Position used to calculate beam phase offset
    %   - polarisation1 -- Primary polarisation direction
    %   - polarisation2 -- Secondary polarisation direction
    %   - intensity     -- Intensity of the beam

    function direction = get.direction(beam)
      direction = beam.rotation(:, 3);
    end
    function direction = get.polarisation1(beam)
      direction = beam.rotation(:, 1);
    end
    function direction = get.polarisation2(beam)
      direction = beam.rotation(:, 2);
    end

    function position = get.origin(beam)
      position = beam.position;
    end

    function intensity = get.intensity(beam)
      intensity = vecnorm(beam.field);
    end

    function beam = set.field(beam, val)
      assert(isnumeric(val) && isvector(val) && numel(val) == 2, ....
          'field must be 3 element numeric vector');;
      beam.field = val(:);
    end
  end
end
