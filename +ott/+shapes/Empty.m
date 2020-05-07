classdef Empty < ott.shapes.Shape
% An empty shape (with no geometry)
%
% This is the default element in an empty Shape array.
%
% Properties
%   volume        % Shape has no volume
%   maxRadius     % Shape max radius is zero

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  methods
    function shape = Empty(varargin)
      % Construct a new empty shape
      %
      % Usage
      %   shape = Empty(...)
      %
      % Optional parameters
      %   - position (3xN numeric) -- Position of the shape.
      %     Default: ``[0;0;0]``.
      %
      %   - rotation (3x3N numeric) -- Orientation of the shape.
      %     Default: ``eye(3)``.

      shape = shape@ott.shapes.Shape(varargin{:});
    end
  end

  methods (Hidden)
    function b = insideRtpInternal(shape, rtp)
      % No point are inside empty shapes
      b = false(size(rtp(1, :)));
    end

    function b = insideXyzInternal(shape, xyz)
      % No point are inside empty shapes
      b = false(size(xyz(1, :)));
    end

    function nxyz = normalsXyzInternal(shape, xyz)
      % Normals are nan (no normals for empty shapes)
      nxyz = nan(size(xyz));
    end

    function nxyz = normalsRtpInternal(shape, rtp)
      % Normals are nan (no normals for empty shapes)
      nxyz = nan(size(rtp));
    end

    function r = get_maxRadius(shape)
      % Shape has no radius
      r = 0.0;
    end

    function v = get_maxVolume(shape)
      % Shape has no volume
      v = 0.0;
    end
  end
end
