classdef (Abstract) ShapeCart < ott.shapes.Shape
% Abstract class for shapes best described in a Cartesian coordinate system.
% Inherits from :class:`Shape`.
%
% Abstract methods
%   - insideXyz     -- Determines if a point is inside the shape
%
% Methods
%   - insideRtp     -- Spherical coordinate inputs, calls insideXyz

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  methods
    function b = insideRtp(shape, varargin)
      % Determine if point is inside the shape (Spherical coordinates)
      %
      % Usage
      %   b = shape.insideRtp(radius, theta, phi, ...) determine if the
      %   point described by radius, theta (polar), phi (azimuthal)
      %   is inside the shape.
      %
      %   b = shape.insideRtp(rtp, ...) as above, but uses a 3xN input.
      %
      % Optional arguments
      %   - origin (enum) -- Coordinate system origin.  Either 'world'
      %     or 'shape' for world coordinates or shape coordinates.

      % Parse inputs
      rtp = shape.insideRtpParseArgs(shape.position, varargin{:});

      % For planes, it's easier to work in Cartesian coordinates
      xyz = ott.utils.rtp2xyz(rtp.').';

      % Call Cartesian method
      b = shape.insideXyz(xyz, 'origin', 'shape');
    end
  end
end

