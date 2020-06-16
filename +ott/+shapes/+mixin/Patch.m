classdef Patch < ott.shapes.mixin.IntersectTriMesh
% Shapes describes by lists of vertices and faces.
% Inherits from :class:`IntersectTriMesh`.
%
% Abstract properties
%   - verts     -- 3xN list of N vertices
%   - faces     -- mxN list of patches
%
% Methods
%   - surf        -- Generate a surface using the patch function
%   - surfPoints  -- Cast to TriangularMesh and call surfPoints
%   - intersectInternal -- Casts to TriangularMesh
%   - intersectAllInternal -- Casts to TriangularMesh
%
% Supported casts
%   - TriangularMesh    -- Inherited (uses PatchMesh)
%   - PatchMesh

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (Abstract)
    verts
    faces
    position
    rotation
  end

  methods
    function varargout = surf(shape, varargin)
      % Generate a visualisation of the shape
      %
      % Usage
      %   p = shape.surf()
      %   Display visualsation in current figure.  Optionally, returns the
      %   patch handle.
      %
      % Optional named parameters
      %   - axes (handle) -- axis to place surface in (default: gca)
      %
      %   - surfOptions (cell) -- options to be passed to surf (default: {})
      %
      %   - showNormals (logical) -- Show surface normals (default: false)
      %
      %   - origin (enum) -- Coordinate origin for drawing.
      %     Can be 'global' or 'local'  Default: 'global'.

      p = inputParser;
      p.addParameter('surfOptions', {});
      p.addParameter('axes', []);
      p.addParameter('showNormals', false);
      p.addParameter('origin', 'global');
      p.parse(varargin{:});

      xyz = shape.verts;

      % Translate to world coordinates
      switch p.Results.origin
        case 'global'
          xyz = shape.local2global(xyz);
        case 'local'
          % Nothing to do
        otherwise
          error('Unknown origin specified');
      end

      % Place the surface in the specified axes
      our_axes = p.Results.axes;
      if isempty(our_axes)
        our_axes = gca();
      end

      % Patch doesn't watch for hold, so clear it ourselves
      isholdon = ishold();
      if ~isholdon
        cla(our_axes);
      end

      sp = patch(our_axes, 'Faces', shape.faces.', ...
          'Vertices', xyz.', 'FaceColor', 'flat', ...
          'FaceVertexCData', [0.9290 0.6940 0.1250], ...
          p.Results.surfOptions{:});

      % Clear the orientation/aspect if hold isn't on
      if ~isholdon
        view(our_axes, [60, 30]);
        daspect(our_axes, [1, 1, 1]);
      end

      % Compute and add surface normals at centre of faces
      if p.Results.showNormals

        % Find mean of each patch
        mXyz = zeros(3, size(shape.faces, 2));
        for ii = 1:size(shape.faces, 1)
          mXyz = mXyz + shape.verts(:, shape.faces(ii, :));
        end
        mXyz = mXyz ./ size(shape.faces, 1);

        nxyz = shape.normalsXyz(mXyz);

        % Preserve hold-no state
        if ~isholdon
          hold('on');
        end

        % Generate plot of surface normals
        if isfinite(shape.maxRadius)
          s = 0.1*shape.maxRadius;
        else
          s = 1.0;
        end
        quiver3(mXyz(1, :), mXyz(2, :), mXyz(3, :), ...
            s.*nxyz(1, :), s.*nxyz(2, :), s.*nxyz(3, :), 0);

        if ~isholdon
          hold('off');
        end
      end

      if nargout ~= 0
        varargout{1} = sp;
      end
    end

    function [xyz, nxyz, dA] = surfPoints(shape, varargin)
      % Cast to TriangularMesh and call surfPoints
      %
      % Usage
      %   [xyz, nxyz, dA] = shape.surfPoints(...)

      shape = ott.shapes.TriangularMesh(shape);
      [xyz, nxyz, dA] = shape.surfPoints(varargin{:});
    end

    function shape = ott.shapes.PatchMesh(shape)
      % Convert the shape to a PatchMesh

      shape = ott.shapes.PatchMesh(shape.verts, shape.faces, ...
          'position', shape.position, 'rotation', shape.rotation);
    end
  end
end
