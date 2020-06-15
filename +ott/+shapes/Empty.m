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

    function varargout = surf(shape, varargin)
      % Constructs a cube for visualisation, turns off the surfaces
      %
      % Usage
      %   p = shape.surf(...)
      %   Returns a handle to the patch.
      %
      % Optional named parameters
      %   - scale (numeric) -- Scale of the cube
      %
      % Unmatched arguments are passed to `cube.surf`.

      p = inputParser;
      p.addParameter('scale', 1.0);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      cube = ott.shapes.Cube(p.Results.scale, ...
          'position', shape.position, 'rotation', shape.rotation);
      sp = cube.surf(unmatched{:});

      % Make transparent
      sp.FaceAlpha = 0;

      % Assign outputs
      if nargout == 1
        varargout{1} = sp;
      end
    end

    function [xyz, nxyz, dA] = surfPoints(varargin)
      % Generate empty arrays

      xyz = zeros(3, 0);
      nxyz = zeros(3, 0);
      dA = [];
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
  end

  methods % Getters/setters
    function r = get.maxRadius(shape)
      % Shape has no radius
      r = 0.0;
    end

    function v = get.volume(shape)
      % Shape has no volume
      v = 0.0;
    end

    function bb = get.boundingBox(shape)
      bb = zeros(3, 2);
    end

    function b = get.starShaped(shape)
      b = true;
    end
    function b = get.xySymmetry(shape)
      b = true;
    end
    function q = get.zRotSymmetry(shape)
      q = 0;
    end
  end
end
