classdef Empty < ott.shapes.Shape
% An empty shape (with no geometry)
%
% The element still inherits from Shape and has position and rotation
% properties.  The object can be visualised with surf, which draws a
% cube with transparent faces.
%
% This is the default element in an empty Shape array.
%
% Properties
%   - volume        -- Shape has no volume
%   - maxRadius     -- Shape max radius is zero
%
% Methods
%   - surf          -- Draws a hollow cube
%   - surfPoints    -- Returns an empty array

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  properties (Dependent)
    volume        % Constant: 0
    maxRadius     % Constant: 0
    boundingBox   % Constant: zeros
    starShaped    % Constant: true
    xySymmetry    % Constant: true
    zRotSymmetry  % Constant: 0
  end

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

    function [xyz, nxyz, dA] = surfPoints(varargin)
      % Generate empty arrays

      xyz = zeros(3, 0);
      nxyz = zeros(3, 0);
      dA = [];
    end
  end

  methods (Hidden)
    function b = insideRtpInternal(~, rtp)
      % No point are inside empty shapes
      b = false(size(rtp(1, :)));
    end

    function b = insideXyzInternal(~, xyz)
      % No point are inside empty shapes
      b = false(size(xyz(1, :)));
    end

    function nxyz = normalsXyzInternal(~, xyz)
      % Normals are nan (no normals for empty shapes)
      nxyz = nan(size(xyz));
    end

    function nxyz = normalsRtpInternal(~, rtp)
      % Normals are nan (no normals for empty shapes)
      nxyz = nan(size(rtp));
    end

    function [locs, norms, dist] = intersectAllInternal(~, vecs)
      % All nans
      locs = nan(3, numel(vecs));
      norms = nan(3, numel(vecs));
      dist = nan(1, numel(vecs));
    end

    function [locs, norms] = intersectInternal(~, vecs)
      % All nans
      locs = nan(3, numel(vecs));
      norms = nan(3, numel(vecs));
    end

    function S = surfInternal(~, varargin)
      % Constructs a cube for visualisation.  Sets face alpha to 0.
      %
      % Usage
      %   p = shape.surfInternal(...)
      %   Returns a handle to the patch.
      %
      % Optional named parameters
      %   - scale (numeric) -- Scale of the cube
      %
      % Unmatched arguments are passed to `Cube.surfInternal`.

      p = inputParser;
      p.addParameter('scale', 1.0);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      cube = ott.shapes.Cube(p.Results.scale);
      S = cube.surfInternal(unmatched{:});

      % Make transparent
      S.FaceAlpha = 0;
    end
  end

  methods % Getters/setters
    function r = get.maxRadius(~)
      % Shape has no radius
      r = 0.0;
    end

    function v = get.volume(~)
      % Shape has no volume
      v = 0.0;
    end

    function bb = get.boundingBox(~)
      bb = zeros(3, 2);
    end

    function b = get.starShaped(~)
      b = true;
    end
    function b = get.xySymmetry(~)
      b = true;
    end
    function q = get.zRotSymmetry(~)
      q = 0;
    end
  end
end
