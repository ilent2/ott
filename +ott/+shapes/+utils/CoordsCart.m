classdef (Abstract) CoordsCart
% A mixin class for shapes best defined in Cartesian coordinates.
%
% Abstract methods
%   - insideXyzInternal     -- Determines if a point is inside the shape
%   - normalsXyzInternal    -- Calculate normals at a surface location
%
% Methods
%   - insideRtpInternal     -- Spherical coordinate inputs, calls insideXyz
%   - normalsRtpInternal    -- Calculate normals, calls normalsXyz

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  methods (Hidden)
    function b = insideRtpInternal(shape, rtp, varargin)
      % Determine if point is inside the shape (Spherical coordinates)

      % For planes, it's easier to work in Cartesian coordinates
      xyz = ott.utils.rtp2xyz(rtp);

      % Call Cartesian method
      b = shape.insideXyzInternal(xyz, varargin{:});
    end

    function nxyz = normalsRtpInternal(shape, rtp, varargin)
      % Calculate normals at the specified surface locations

      % Transform to Cartesian coordinates
      xyz = ott.utils.rtp2xyz(rtp);

      % Call Spherical method
      nxyz = shape.normalsXyzInternal(xyz, varargin{:});
    end
  end
end

