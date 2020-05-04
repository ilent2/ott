classdef PlaneWave < ott.beam.abstract.Abstract & ott.utils.Vector
% Abstract representation of a plane wave beam
%
% Properties
%   - origin        -- Ray origins, 3xN array (default [0;0;0])
%   - direction     -- Direction of propagation (3xN Cartesian)
%   - field         -- Field parallel and perpendicular to polarisation
%   - polarisation  -- Primary polarisation direction
%
% Methods
%   - rotate      -- Rotate the direction and polarisation
%   - rotate*     -- Rotate the particle around the X,Y,Z axis
%   - size        -- Get size of beam array

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    field             % Field parallel and perpendicular to polarisation
    polarisation      % Polarisation direction
  end

  properties (Dependent)
    wavevector
  end

  methods
    function beam = PlaneWave(varargin)
      % Construct a new abstract plane wave representation
      %
      % Usage
      %   beam = PlaneWave(...)
      %
      % Optional named arguments
      %   - direction (3xN numeric) -- direction vectors (Cartesian)
      %     Default: ``[0;0;1]``.
      %
      %   - polarisation (3xN numeric) -- polarisation vectors (Cartesian)
      %     Default: ``[1;0;0]``.
      %
      %   - field (1xN|2xN numeric) -- Field vectors parallel and
      %     (optionally) perpendicular to the polarisation direction.
      %     Allows for 0 intensity with finite polarisation direction.
      %     Default: ``1``.
      %
      %   - origin (3xN numeric) -- Origin of plane waves.
      %     Default: ``[0;0;0]``.

      % Parse parameters
      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('direction', [0;0;1]);
      p.addParameter('polarisation', [0;0;1]);
      p.addParameter('origin', [0;0;0]);
      p.addParameter('field', 1.0);
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Get Vector to store most
      beam = beam@ott.utils.Vector(p.Results.origin, p.Results.direction);
      beam = beam@ott.beam.abstract.Abstract(unmatched{:});

      % Store remaining parameters
      beam.field = p.Results.field;
      beam.polarisation = p.Results.polarisation;
    end

    function beam = rotate(beam, varargin)
      % Rotate the beam and the polarisation vector
      %
      % Usage
      %   rbeam = beam.rotate(R, ...)
      %
      % Parameters
      %   - R (3x3 numeric) -- rotation matrix
      %
      % Optional named arguments
      %   - origin (logical) -- If true, the origin is rotated too.
      %     Default: ``false``.

      % Rotate the location (and origin)
      beam = rotate@ott.utils.Vector(beam, varargin{:});

      % Rotate the polarisation
      beam.polarisation = R * beam.polarisation;

    end

    function sz = size(vec, varargin)
      % Get the number of beams contained in this object
      %
      % The leading dimension is always 1.  May change in future.

      sz = size(vec.data(1, :), varargin{:});
    end
  end

  methods (Hidden)
    function p = getBeamPower(beam)
      % Returns infinite plane wave power
      p = Inf;
    end
  end

  methods % Getters/setters
    % Properties
    %   - field         -- Field parallel and perpendicular to polarisation
    %   - polarisation  -- Primary polarisation direction

    function wv = get.wavevector(beam)
      % Get the plane wave wave-vector
      wv = beam.direction .* beam.wavenumber;
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

