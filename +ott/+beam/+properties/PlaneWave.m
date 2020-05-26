classdef (Abstract) PlaneWave < ott.beam.properties.MaterialBeam
% Properties of scattered beams.
% Inherits from :class:`MaterialBeam`.
%
% Properties
%   - position      -- Position used to calculate beam phase offset
%   - rotation      -- Describes the direction of the beam.
%   - field         -- Field parallel and perpendicular to polarisation
%   - polarisation  -- Primary polarisation direction
%
% Dependent properties
%   - direction     -- Beam direction
%   - wavevector    -- Beam wave-vector
%   - intensity     -- Intensity of the beam
%
% Methods
%   - rotate*     -- Rotate the particle around the X,Y,Z axis
%   - translate*  -- Apply a phase offset (translation) to the beam

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    field             % Field parallel and perpendicular to polarisation
    polarisation      % Polarisation direction
  end

  properties (Dependent)
    wavevector        % Wave-vectors of plane wave components
    intensity         % Intensity of plane wave components
  end

  methods % Getters/setters
    % Properties
    %   - field         -- Field parallel and perpendicular to polarisation
    %   - polarisation  -- Primary polarisation direction

    function wv = get.wavevector(beam)
      % Get the plane wave wave-vector
      wv = beam.direction .* beam.wavenumber ./ vecnorm(beam.direction);
    end

    function intensity = get.intensity(beam)
      intensity = sum(abs(beam.field), 1);
    end

    function beam = set.field(beam, val)

      % Check type and rows
      assert(isnumeric(val) && any(size(val, 1) == [1, 2]), ...
          'field must be numeric 1xN or 2xN matrix');

      % Check length
      assert(any(size(val, 2) == [1, size(beam.direction, 2)]), ...
          'field must have length 1 or same length as direction');
      if size(val, 2) == 1
        val = repmat(val, 1, size(beam.direction, 2));
      end

      beam.field = val;
    end

    function beam = set.polarisation(beam, val)

      % Check type and rows
      assert(isnumeric(val) && size(val, 1) == 3, ...
          'polarisation must be numeric 3xN matrix');

      % Check length
      assert(any(size(val, 2) == [1, size(beam.direction, 2)]), ...
          'polarisation must have length 1 or same length as direction');
      if size(val, 2) == 1
        val = repmat(val, 1, size(beam.direction, 2));
      end

      beam.polarisation = val;
    end
  end
end
