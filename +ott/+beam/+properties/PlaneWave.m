classdef (Abstract) PlaneWave < ott.beam.properties.Material
% Properties of scattered beams.
% Inherits from :class:`Material`.
%
% Properties
%   - position      -- (Inherited) Beam position
%   - rotation      -- (Inherited) Beam rotation
%
% Abstract properties
%   - field         -- Field parallel and perpendicular to polarisation
%   - direction     -- Beam direction
%   - origin        -- Position used to calculate beam phase offset
%   - polarisation1 -- Primary polarisation direction
%   - polarisation2 -- Secondary polarisation direction
%   - intensity     -- Intensity of the beam
%
% Dependent properties
%   - wavevector    -- Beam wave-vector
%
% Methods
%   - rotate*     -- Rotate the particle around the X,Y,Z axis
%   - translate*  -- Apply a phase offset (translation) to the beam

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Abstract, SetAccess=protected)
    field             % Field parallel and perpendicular to polarisation
    direction         % Beam direction
    origin            % Position used to calculate beam phase offset
    polarisation1     % Primary polarisation direction
    polarisation2     % Secondary polarisation direction
    intensity         % Intensity of plane wave components
  end

  properties (Dependent)
    wavevector        % Wave-vectors of plane wave components
  end

  methods % Getters/setters
    % Properties
    %   - wavevector        % Wave-vectors of plane wave components

    function wv = get.wavevector(beam)
      % Get the plane wave wave-vector
      wv = beam.direction .* beam.wavenumber ./ vecnorm(beam.direction);
    end
  end
end
