classdef (Abstract) Shape < ott.utils.RotationPositionProp ...
    & matlab.mixin.Heterogeneous
% Shape abstract class for optical tweezers toolbox shapes.
% Inherits from :class:`ott.utils.RotationPositionProp` and
% `matlab.mixin.Hetrogeneous`.
%
% Properties disregard the rotation/position properties and describe
% the object in the local coordinates.
%
% Properties
%   - position   -- Location of shape ``[x, y, z]``
%   - rotation   -- Orientation of the particle (3x3 rotation matrix)
%
% Abstract properties
%   - maxRadius         -- Maximum particle radius
%   - volume            -- Particle volume
%   - boundingBox       -- Cartesian coordinate bounding box
%   - starShaped        -- True if the shape is star-shaped
%   - xySymmetry        -- True if shape is xy-mirror symmetric
%   - zRotSymmetry      -- z-axis rotational symmetry of particle
%
% Methods
%   - voxels          -- Generate array of voxels or voxel visualisation
%   - insideRtp       -- Determine if Spherical point is inside shape
%   - insideXyz       -- Determine if Cartesian point is inside shape
%   - normalsRtp      -- Calculate normals at surface location
%   - normalsXyz      -- Calculate normals at surface location
%   - writeWavefrontObj -- write shape to Wavefront OBJ file
%   - intersect       -- Calculate intersection between vectors and surface
%   - intersectAll    -- Calculate intersection between vectors and surface
%   - intersectBoundingBox -- Calculate intersection with bounding box
%   - getBoundingBox  -- Get the bounding box with transformations applied
%   - getBoundingBoxShape -- Get a shape representing the bounding box
%   - rotate*         -- Functions for rotating the entity
%   - translate*      -- Functions for translating the entity
%   - operator|       -- Union operator: creates a new set
%   - operator&       -- Intersection operator: creates a new set
%   - operator~       -- Inverse operator: creates a new :class:`Inverse`.
%
% Abstract methods
%   - surf                -- Generate surface visualisation
%   - surfPoints          -- Calculate points for surface integration
%   - intersectInternal   -- Method called by intersect
%   - intersectAllInternal -- Method called by intersectAll
%   - insideRtpInternal   -- Determine if Spherical point is inside shape
%   - insideXyzInternal   -- Determine if Cartesian point is inside shape
%   - normalsRtpInternal  -- Calculate normals at surface location
%   - normalsXyzInternal  -- Calculate normals at surface location
%
% Supported casts
%   - TriangularMesh -- Requires cast for PatchMesh

% Copyright 2018-2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (Abstract)
    maxRadius          % Maximum particle radius
    volume             % Particle volume
    boundingBox        % Cartesian coordinate bounding box (no rot/pos)
    starShaped         % True if the particle is star-shaped
    xySymmetry         % True if the particle is xy-plane mirror symmetric
    zRotSymmetry       % z-axis rotational symmetry of particle
  end

  methods (Abstract)
    surf(obj)            % Generate surface visualisation
    surfPoints(obj)      % Calculate points for surface integration
  end

  methods (Abstract, Hidden)
    insideRtpInternal(obj)    % Determine if point is inside shape (Spherical)
    insideXyzInternal(obj)    % Determine if point is inside shape (Cartesian)
    normalsRtpInternal(obj)   % Calculate normals (Spherical)
    normalsXyzInternal(obj)   % Calculate normals (Cartesian)
    intersectInternal(obj)    % Method called by intersect
    intersectAllInternal(obj)    % Method called by intersectAll
  end

  methods
    function shape = Shape(varargin)
      % Construct a new shape instance.
      %
      % This class cannot be instanced directly, use one of the other
      % shape descriptions to create a new shape.
      %
      % Usage
      %   shape = shape@ott.shapes.Shape(...)
      %
      % Optional named arguments
      %   - position (3 numeric) -- Position of the shape.
      %     Default: ``[0;0;0]``.
      %
      %   - rotation (3x3 numeric) -- Orientation of the shape.
      %     Default: ``eye(3)``.

      p = inputParser;
      p.addParameter('position', [0;0;0]);
      p.addParameter('rotation', eye(3));
      p.parse(varargin{:});

      Nposition = size(p.Results.position, 2);
      Nrotation = size(p.Results.rotation, 2);
      assert(mod(Nrotation, 3) == 0, ...
        'rotation must be a 3x3N matrix');
      Nrotation = Nrotation ./ 3;

      % Create array if required
      if Nposition ~= 1 || Nrotation ~= 1
        assert(Nposition == 1 || Nrotation == 1 || Nposition == Nrotation, ...
            'length of position and rotation must be same length');
        shape = repelem(shape, 1, max([Nposition, Nrotation]));
      end

      % Store initial position
      if Nposition > 1
        for ii = 1:Nposition
          shape(ii).position = p.Results.position(:, ii);
        end
      else
        [shape.position] = deal(p.Results.position);
      end

      % Store initial rotation
      if Nrotation > 1
        for ii = 1:Nrotation
          shape(ii).rotation = p.Results.rotation(:, (1:3) + (ii-1)*3);
        end
      else
        [shape.rotation] = deal(p.Results.rotation);
      end
    end

    function shape = ott.shapes.TriangularMesh(shape, varargin)
      % Cast to TriangularMesh via PatchMesh
      shape = ott.shapes.PatchMesh(shape, varargin{:});
      shape = ott.shapes.TriangularMesh(shape);
    end

    function shape = not(shape)
      % Take the inverse of the shape
      %
      % Usage
      %   shape = ~shape

      shape = ott.shapes.Inverse(shape);
    end

    function shape = or(a, b)
      % Create a union of two shapes
      %
      % Usage
      %   shape = shape1 | shape2

      shape = ott.shapes.Union([a, b]);
    end

    function shape = and(a, b)
      % Create a intersection of two shapes
      %
      % Usage
      %   shape = shape1 & shape2

      shape = ott.shapes.Intersection([a, b]);
    end

    function bb = getBoundingBox(shape, varargin)
      % Get bounding box after applying transformations
      %
      % Usage
      %   bb = shape.getBoundingBox(...)
      %
      % Optional named arguments
      %   - origin (enum) -- Coordinate origin.  'local' or 'global'.

      p = inputParser;
      p.addParameter('origin', 'global');
      p.parse(varargin{:});

      bb = shape.boundingBox;

      % Translate to world coordinates
      switch p.Results.origin
        case 'global'
          xyz = [bb(1, [1, 1, 1, 1, 2, 2, 2, 2]); ...
                 bb(2, [1, 1, 2, 2, 1, 1, 2, 2]); ...
                 bb(3, [1, 2, 1, 2, 1, 2, 1, 2])];
          xyz = shape.local2global(xyz);
          bb = [min(xyz, [], 2), max(xyz, [], 2)];
        case 'local'
          % Nothing to do
        otherwise
          error('Unknown origin specified');
      end
    end

    function bbshape = getBoundingBoxShape(shape, varargin)
      % Get shape representing the bounding box.
      %
      % Usage
      %   bb = shape.getBoundingBoxShape(...)
      %
      % Optional named arguments
      %   - origin (enum) -- Coordinate origin.  'local' or 'global'.

      p = inputParser;
      p.addParameter('origin', 'global');
      p.parse(varargin{:});

      bb = shape.boundingBox;
      bbsize = diff(bb, [], 2);
      bbshape = ott.shapes.RectangularPrism(bbsize);

      % Translate to world coordinates
      switch p.Results.origin
        case 'global'
          bbshape.position = shape.position + bb(:, 1) + bbsize./2;
          bbshape.rotation = shape.rotation;
        case 'local'
          bbshape.position = bb(:, 1) + bbsize./2;
        otherwise
          error('Unknown origin specified');
      end
    end

    function writeWavefrontObj(shape, filename)
      % Write representation of shape to Wavefront OBJ file
      %
      % Convert the shape to a TriangularMesh and write it to a file.
      %
      % Usage
      %   shape.writeWavefrontObj(filename)

      shape = ott.shapes.TriangularMesh(shape);
      shape.writeWavefrontObj(filename);
    end

    function [locs, norms] = intersectBoundingBox(shape, vecs, varargin)
      % Calculate intersections with bounding box.
      %
      % Internally, constructs a bounding box shape and calculates
      % all intersections with :meth:`intersectAll`.
      %
      % Usage
      %   [locs, norms] = shape.intersectBoundingBox(shape, vecs, ...)
      %
      % All unmatched arguments are passed to getBoundingBoxShape.

      bbshape = shape.getBoundingBoxShape(varargin{:});
      [locs, norms] = bbshape.intersectAll(vecs);
    end

    function [locs, norms, dist] = intersectAll(shape, vecs, varargin)
      % Calculate intersection with all faces in a direction
      %
      % Any rays/faces that don't intersect result in nans being returned.
      %
      % Usage
      %   [locs, norms, dist] = shape.intersectAll(shape, vecs, ...)
      %
      % Parameters
      %   - vecs (3xM Vector) -- vectors to intersect.  Must be castable
      %     to a ott.utils.Vector object.
      %
      % Returns
      %   - locs (3xNxM numeric) -- intersections with N locations
      %   - norms (3xNxM numeric) -- surface normals at N intersections
      %   - dist (1xNxM numeric) -- distance from vector origin
      %
      % Optional named arguments
      %   - origin (enum) -- Coordinate origin.  'local' or 'global'.
      %
      % Additional arguments passed to intersectAllInternal.

      p = inputParser;
      p.addParameter('origin', 'global');
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      if ~isa(vecs, 'ott.utils.Vector')
        vecs = ott.utils.Vector(vecs);
      end

      % Translate to local coordinates
      switch p.Results.origin
        case 'global'
          vecs.origin = shape.global2local(vecs.origin);
        case 'local'
          % Nothing to do
        otherwise
          error('Unknown origin specified');
      end

      % Defer to shape-specific implementation
      [locs, norms, dist] = shape.intersectAllInternal(vecs, unmatched{:});

      % Translate to global coordinatse
      switch p.Results.origin
        case 'global'
          sz = size(locs);
          locs = reshape(shape.local2global(reshape(locs, 3, [])), sz);
          norms = reshape(shape.rotation * reshape(norms, 3, []), sz);
        case 'local'
          % Nothing to do
        otherwise
          error('Unknown origin specified');
      end
    end

    function [locs, norms] = intersect(shape, vecs, varargin)
      % Calculate intersection locations and normals.
      %
      % Any rays that don't intersect shape result in nans.
      %
      % Usage
      %   [locs, norms, dist] = shape.intersectAll(shape, vecs, ...)
      %
      % Parameters
      %   - vecs (3xM Vector) -- vectors to intersect.  Must be castable
      %     to a ott.utils.Vector object.
      %
      % Returns
      %   - locs (3xN numeric) -- intersections with N locations
      %   - norms (3xN numeric) -- surface normals at N intersections
      %
      % Optional named arguments
      %   - origin (enum) -- Coordinate origin.  'local' or 'global'.
      %
      % Additional arguments passed to intersectInternal.

      p = inputParser;
      p.addParameter('origin', 'global');
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      if ~isa(vecs, 'ott.utils.Vector')
        vecs = ott.utils.Vector(vecs);
      end

      % Translate to local coordinates
      switch p.Results.origin
        case 'global'
          vecs.origin = shape.global2local(vecs.origin);
        case 'local'
          % Nothing to do
        otherwise
          error('Unknown origin specified');
      end

      % Defer to shape-specific implementation
      [locs, norms] = shape.intersectInternal(vecs, unmatched{:});

      % Translate to global coordinatse
      switch p.Results.origin
        case 'global'
          locs = shape.local2global(locs);
          norms = shape.rotation * norms;
        case 'local'
          % Nothing to do
        otherwise
          error('Unknown origin specified');
      end
    end

    function varargout = voxels(shape, varargin)
      % Generate an array of xyz coordinates for voxels inside the shape
      %
      % Usage
      %   voxels(spacing) shows a visualisation of the shape with
      %   circles placed at locations on a Cartesian grid.
      %
      %   xyz = voxels(spacing) returns the voxel locations.
      %
      % Optional named arguments
      %   - 'plotoptions'   Options to pass to the plot3 function
      %
      %   - 'visualise'     Show the visualisation (default: nargout == 0)
      %
      %   - origin (enum) -- Coordinate system origin.  Either 'global'
      %     or 'local' for world coordinates or shape coordinates.
      %
      %   - axis (handle) -- Axes hanlde to place plot in.
      %     Default: ``[]``, uses gca() when available.

      p = inputParser;
      p.addOptional('spacing', shape.maxRadius/10);
      p.addParameter('plotoptions', []);
      p.addParameter('visualise', nargout == 0);
      p.addParameter('scale', 1.0);
      p.addParameter('axes', []);
      p.addParameter('origin', 'local');
      p.addParameter('even_range', false);
      p.parse(varargin{:});

      plotoptions = p.Results.plotoptions;
      if isempty(plotoptions)
        plotoptions = {...
          'MarkerFaceColor', 'w', ...
          'MarkerEdgeColor', [.5 .5 .5], ...
          'MarkerSize', 20*p.Results.spacing/shape.maxRadius/p.Results.scale};
      end

      % Calculate range of dipoles
      numr = ceil(shape.maxRadius * p.Results.scale / p.Results.spacing);

      % Add an extra point so we don't have a point around zero
      if p.Results.even_range
        numr = numr + 0.5;
      end

      rrange = (-numr:numr)*p.Results.spacing;

      % Generate the voxel grid
      [xx, yy, zz] = meshgrid(rrange, rrange, rrange);
      xyz = [xx(:), yy(:), zz(:)].';

      % Remove points outside shape
      mask = shape.insideXyz(xyz ./ p.Results.scale, 'origin', 'local');
      xyz(:, ~mask) = [];

      % Translate to world coordinates
      switch p.Results.origin
        case 'global'
          xyz = shape.local2global(xyz);
        case 'local'
          % Nothing to do
        otherwise
          error('Unknown origin specified');
      end

      % Visualise the result
      if p.Results.visualise

        % Get the axes to use
        our_axes = p.Results.axes;
        if isempty(our_axes)
          our_axes = gca();
        end

        plot3(our_axes, xyz(1,:), xyz(2,:), xyz(3,:), 'o', plotoptions{:});
        axis(our_axes, 'equal');
        title(our_axes, ['spacing = ' num2str(p.Results.spacing) ...
            ', N = ' int2str(sum(mask))])
      end

      % Assign output
      if nargout ~= 0
        varargout = {xyz};
      end
    end

    function varargout = isosurface(shape, varargin)
      % Generate an isosurfrace for the shape from the voxel data
      %
      % This method is more intensive than the surf method and often
      % less accurate but provides a useful alternative when the shape
      % doesn't directly support a surf method.
      %
      % Usage
      %   FV = shape.isosurface(...)
      %
      % Optional named parameters
      %   - samples (3 numeric) -- Number of samples per Cartesian axes.
      %
      %   - visualise (logical) -- Show the visualisation.
      %     Default: ``nargout == 0``
      %
      %   - origin (enum) -- Coordinate system origin.  Either 'global'
      %     or 'local' for world coordinates or shape coordinates.
      %
      %   - axis (handle) -- Axes hanlde to place plot in.
      %     Default: ``[]``, uses gca() when available.

      p = inputParser;
      p.addParameter('samples', [50, 50, 50]);
      p.addParameter('visualise', nargout == 0);
      p.addParameter('axes', []);
      p.addParameter('origin', 'local');
      p.addParameter('padding', 0.1);
      p.parse(varargin{:});

      % Get shape bounding box
      bb = shape.boundingBox;

      % Compute padding
      pad = p.Results.padding .* diff(bb, [], 2);

      % Generate grid of points
      x = linspace(bb(1, 1)-pad(1), bb(1, 2)+pad(1), p.Results.samples(1));
      y = linspace(bb(2, 1)-pad(2), bb(2, 2)+pad(2), p.Results.samples(2));
      z = linspace(bb(3, 1)-pad(3), bb(3, 2)+pad(3), p.Results.samples(3));
      [X, Y, Z] = meshgrid(x, y, z);
      xyz = [X(:), Y(:), Z(:)].';

      % Generate voxels
      voxels = reshape(shape.insideXyz(xyz, 'origin', 'local'), size(X));

      % Smooth the data
      voxels = smooth3(voxels);

      % Translate to world coordinates
      switch p.Results.origin
        case 'global'
          xyz = shape.local2global(xyz);
          X = reshape(xyz(1, :), size(X));
          Y = reshape(xyz(2, :), size(Y));
          Z = reshape(xyz(3, :), size(Z));
        case 'local'
          % Nothing to do
        otherwise
          error('Unknown origin specified');
      end

      % Generate isosurface
      FV = isosurface(X, Y, Z, voxels, 0.5);

      % Generate visualisation
      if p.Results.visualise

        % Get the axes to use
        our_axes = p.Results.axes;
        if isempty(our_axes)
          our_axes = gca();
        end

        % Patch doesn't watch for hold, so clear it ourselves
        isholdon = ishold();
        if ~isholdon
          cla(our_axes);
        end

        sp = patch(our_axes, FV, 'FaceColor', 'flat', ...
            'FaceVertexCData', [0.9290 0.6940 0.1250]);
        isonormals(X,Y,Z, voxels, sp);
        sp.EdgeColor = 'none';
        camlight; lighting phong;

        % Clear the orientation/aspect if hold isn't on
        if ~isholdon
          view(our_axes, [60, 30]);
          daspect(our_axes, [1, 1, 1]);
        end
      end

      % Assign output
      if nargout ~= 0
        varargout{1} = FV;
      end
    end

    function b = insideRtp(shape, varargin)
      % Determine if point is inside the shape (Spherical coordinates)
      %
      % Usage
      %   b = shape.insideRtp(rtp, ...)
      %
      % Parameters
      %   - rtp (3xN numeric) -- Spherical coordinates.
      %
      % Optional arguments
      %   - origin (enum) -- Coordinate system origin.  Either 'global'
      %     or 'local' for world coordinates or shape coordinates.

      % Parse inputs
      rtp = shape.insideRtpParseArgs(varargin{:});

      b = shape.insideRtpInternal(rtp);
    end

    function nxyz = normalsRtp(shape, varargin)
      % Calculate normals at the specified surface locations
      %
      % Usage
      %   nxyz = shape.normalsRtp(rtp, ...)
      %   Result is in Cartesian coordinates.
      %
      % Parameters
      %   - rtp (3xN numeric) -- Spherical coordinates.
      %
      % Optional arguments
      %   - origin (enum) -- Coordinate system origin.  Either 'global'
      %     or 'local' for world coordinates or shape coordinates.

      % Parse inputs
      rtp = shape.insideRtpParseArgs(varargin{:});

      nxyz = shape.normalsRtpInternal(rtp);
    end

    function b = insideXyz(shape, varargin)
      % Determine if point is inside the shape (Cartesian coordinates)
      %
      % Usage
      %   b = shape.insideXyz(xyz, ...)
      %
      % Parameters
      %   - xyz (3xN numeric) -- Cartesian coordinates.
      %
      % Optional arguments
      %   - origin (enum) -- Coordinate system origin.  Either 'global'
      %     or 'local' for world coordinates or shape coordinates.

      % Get xyz coordinates from inputs and translated to origin
      xyz = shape.insideXyzParseArgs(varargin{:});

      b = shape.insideXyzInternal(xyz);
    end

    function nxyz = normalsXyz(shape, varargin)
      % Calculate normals at the specified surface locations
      %
      % Usage
      %   nxyz = shape.normalsXyz(xyz, ...)
      %   Result is in Cartesian coordinates.
      %
      % Parameters
      %   - xyz (3xN numeric) -- Cartesian coordinates.
      %
      % Optional arguments
      %   - origin (enum) -- Coordinate system origin.  Either 'global'
      %     or 'local' for world coordinates or shape coordinates.

      % Get xyz coordinates from inputs and translated to origin
      xyz = shape.insideXyzParseArgs(varargin{:});

      nxyz = shape.normalsXyzInternal(xyz);
    end
  end

  methods (Hidden)
    function xyz = insideXyzParseArgs(shape, xyz, varargin)
      % Helper for parsing arguments for insideXyz functions

      % Collect arguments
      p = inputParser;
      p.addParameter('origin', 'global');
      p.parse(varargin{:});

      % Translate to local coordinates
      switch p.Results.origin
        case 'global'
          % Translate to local coordinates, the Internal functions use local
          xyz = shape.global2local(xyz);
        case 'local'
          % Nothing to do
        otherwise
          error('Unknown origin specified');
      end
    end

    function rtp = insideRtpParseArgs(shape, rtp, varargin)
      % Argument parser helper for insideRtp function

      % Collect arguments
      p = inputParser;
      p.addParameter('origin', 'global');
      p.parse(varargin{:});

      assert(isnumeric(rtp) && size(rtp, 1) == 3, ...
          'rtp must be 3xN numeric matrix');
      assert(all(rtp(1, :) >= 0), 'Radii must be positive');

      % Translate to local coordinates
      switch p.Results.origin
        case 'global'
          % Translate to local coordinates, the Internal functions use local
          % Only do work if we need to
          if vecnorm(shape.position) ~= 0
            xyz = ott.utils.rtp2xyz(rtp);
            xyz = shape.global2local(xyz);
            rtp = ott.utils.xyz2rtp(xyz);
          end
        case 'local'
          % Nothing to do
        otherwise
          error('Unknown origin specified');
      end
    end
  end

  methods (Static, Sealed, Access = protected)
    function default_object = getDefaultScalarElement
      % Default object for a Shape array is an empty shape
      default_object = ott.shapes.Empty;
    end
  end
end
