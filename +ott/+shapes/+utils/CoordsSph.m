classdef (Abstract) CoordsSph
% A mixin class for shapes best defined in Spherical coordinates.
%
% Abstract methods
%   - insideRtpInternal     -- Determines if a point is inside the shape
%   - normalsRtpInternal    -- Calculate normals at a surface location
%
% Methods
%   - insideXyzInternal     -- Cartesian coordinate inputs, calls insideRtp
%   - normalsXyzInternal    -- Calculate normals, calls normalsRtp

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  methods (Hidden)
    function b = insideXyzInternal(shape, xyz, varargin)
      % Determine if point is inside the shape (Cartesian coordinates)

      % For planes, it's easier to work in Cartesian coordinates
      rtp = ott.utils.xyz2rtp(xyz);

      % Call Cartesian method
      b = shape.insideRtpInternal(rtp, 'origin', 'shape');
    end

    function nxyz = normalsXyzInternal(shape, xyz, varargin)
      % Calculate normals at the specified surface locations

      % Transform to spherical
      rtp = ott.utils.xyz2rtp(xyz);

      % Call Spherical method
      nxyz = shape.normalsRtpInternal(rtp, varargin{:});
    end
  end
end

