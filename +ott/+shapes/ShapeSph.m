classdef (Abstract) ShapeSph < ott.shapes.Shape
% Abstract class for shapes best described in a Spherical coordinate system.
% Inherits from :class:`Shape`.
%
% Abstract methods
%   - insideRtp     -- Determines if a point is inside the shape
%
% Methods
%   - insideXyz     -- Cartesian coordinate inputs, calls insideRtp

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  methods
    function b = insideXyz(shape, varargin)
      % Determine if point is inside the shape (Cartesian coordinates)
      %
      % Usage
      %   b = shape.insideXyz(x, y, z, ...) determine if the
      %   point in Cartesian coordinates is inside the shape.
      %
      %   b = shape.insideXyz(xyz, ...) as above, but uses a 3xN input.
      %
      % Optional arguments
      %   - origin (enum) -- Coordinate system origin.  Either 'world'
      %     or 'shape' for world coordinates or shape coordinates.

      % Parse inputs
      xyz = shape.insideXyzParseArgs(shape.position, varargin{:});

      % For planes, it's easier to work in Cartesian coordinates
      rtp = ott.utils.xyz2rtp(xyz);

      % Call Cartesian method
      b = shape.insideRtp(rtp, 'origin', 'shape');
    end
  end
end

