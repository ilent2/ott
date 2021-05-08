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
%   - surf            -- Generate surface visualisation
%   - scale           -- Scale the geometry (uses `scaleInternal`)
%   - voxels          -- Generate array of voxels or voxel visualisation
%   - insideRtp       -- Determine if Spherical point is inside shape
%   - insideXyz       -- Determine if Cartesian point is inside shape
%   - normalsRtp      -- Calculate normals at surface location
%   - normalsXyz      -- Calculate normals at surface location
%   - writeWavefrontObj -- write shape to Wavefront OBJ file
%   - intersect       -- Calculate intersection between vectors and surface
%   - starRadii       -- Calculate radii of star shaped particle
%   - intersectAll    -- Calculate intersection between vectors and surface
%   - intersectBoundingBox -- Calculate intersection with bounding box
%   - getBoundingBox  -- Get the bounding box with transformations applied
%   - getBoundingBoxShape -- Get a shape representing the bounding box
%   - rotate*         -- Functions for rotating the entity
%   - translate*      -- Functions for translating the entity
%   - operator/       -- Scale the object (uses scale)
%   - operator*       -- Scale the object (uses scale)
%   - operator|       -- Union operator: creates a new set
%   - operator&       -- Intersection operator: creates a new set
%   - operator~       -- Inverse operator: creates a new :class:`Inverse`.
%
% Abstract methods
%   - scaleInternal       -- Scale the geometry (called by times/rdivide)
%   - surfInternal        -- Generate surface visualisation
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

% TODO: Whats the difference between voxels.scale and voxels.spacing?
%   Should voxels.scale exist?

  properties (Abstract)
    maxRadius          % Maximum particle radius
    volume             % Particle volume
    boundingBox        % Cartesian coordinate bounding box (no rot/pos)
    starShaped         % True if the particle is star-shaped
    xySymmetry         % True if the particle is xy-plane mirror symmetric
    zRotSymmetry       % z-axis rotational symmetry of particle
  end

  methods (Abstract)
    surfPoints(obj)      % Calculate points for surface integration
    scaleInternal(obj)           % Scale the object
  end

  methods (Abstract, Hidden)
    surfInternal(obj)         % Generate surface visualisation
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
      %   shape = shape@ott.shape.Shape(...)
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
      Nrotation = floor(Nrotation ./ 3);

      % Create array if required
      if Nposition ~= 1 || Nrotation ~= 1
        assert(Nposition == 1 || Nrotation == 1 || Nposition == Nrotation, ...
            'length of position and rotation must be same length');
      end
      shape = repelem(shape, 1, max([Nposition, Nrotation]));

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

    function shape = ott.shape.TriangularMesh(shape, varargin)
      % Cast to TriangularMesh via PatchMesh
      shape = ott.shape.PatchMesh(shape, varargin{:});
      shape = ott.shape.TriangularMesh(shape);
    end

    function shape = times(a, b)
      % Scale the object
      %
      % Usage
      %   shape = shape * scalar     or     shape = scalar * shape

      if isnumeric(a)
        shape = b.scale(a);
      elseif isnumeric(b)
        shape = a.scale(b);
      else
        error('Either a or b must be numeric');
      end
    end

    function shape = mtimes(a, b)
      % Defers to times
      shape = times(a, b);
    end

    function shape = rdivide(a, b)
      % Scale the shape by a factor
      %
      % Usage
      %   shape = shape ./ scalar

      if isnumeric(b)
        shape = a.scale(1./b);
      else
        error('b must be numeric');
      end
    end

    function shape = mrdivide(a, b)
      % Defers to rdivide
      shape = rdivide(a, b);
    end

    function shape = not(shape)
      % Take the inverse of the shape
      %
      % Usage
      %   shape = ~shape

      shape = ott.shape.Inverse(shape);
    end

    function shape = or(a, b)
      % Create a union of two shapes
      %
      % Usage
      %   shape = shape1 | shape2

      shape = ott.shape.Union([a, b]);
    end

    function shape = and(a, b)
      % Create a intersection of two shapes
      %
      % Usage
      %   shape = shape1 & shape2

      shape = ott.shape.Intersection([a, b]);
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
      bbshape = ott.shape.RectangularPrism(bbsize);

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
      %
      % Parameters
      %   - filename (char | string) -- Filename for file to write to.

      shape = ott.shape.TriangularMesh(shape);
      shape.writeWavefrontObj(filename);
    end

    function R = starRadii(shape, theta, phi)
      % Calculate radii for star shaped particles.
      %
      % Checks that the particle is star shaped, then uses intersect
      % to determine the surface locations.  For intersect, starts each
      % ray from the shape center (i.e., ``shape.position``).
      %
      % Usage
      %   R = shape.starRadii(theta, phi)
      %
      % Parameters
      %   - theta, phi (numeric) -- Angular coordinates

      % Check radii
      if ~shape.starShaped
        warning('ott:shapes:Shape:not_star_shaped', ...
            'Shape may not be star shaped, radii may not be unique');
      end

      % Construct vectors for intersection calculation
      directions = ott.utils.rtp2xyz(ones(size(theta)), theta, phi);

      % Find intersections
      fudge = 1.8;   % Add ~sqrt(3) for bounding box intersection
      locs = shape.intersect(shape.position .* ones(size(directions)), ...
          shape.position + directions.*shape.maxRadius*fudge);

      % Calculate radii
      R = vecnorm(locs);
      
      % Ensure output size matches input
      R = reshape(R, size(theta));
    end

    function [locs, norms] = intersectBoundingBox(shape, x0, x1, varargin)
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
      [locs, norms] = bbshape.intersectAll(x0, x1);
    end

    function [locs, norms, dist] = intersectAll(shape, x0, x1, varargin)
      % Calculate intersection with all faces in a direction
      %
      % Any rays/faces that don't intersect result in nans being returned.
      %
      % Usage
      %   [locs, norms, dist] = shape.intersectAll(shape, x0, x1, ...)
      %
      % Parameters
      %   - x0,x1 (3xM numeric) -- Describes the intersection ray.
      %     x0 is the starting point, x1 is in the forward direction.
      %
      % Returns
      %   - locs (3xNxM numeric) -- intersections with N locations
      %   - norms (3xNxM numeric) -- surface normals at N intersections
      %   - dist (1xNxM numeric) -- distance from vector origin
      %
      % Optional named arguments
      %   - origin (enum) -- Coordinate origin.  'local' or 'global'.
      %
      %   - removeNan (logical) -- If true, removes nan values (i.e.,
      %     faces which are not intersected).  Default: ``false``.
      %     Only works when size(x0) is 3x1.
      %
      % Additional arguments passed to intersectAllInternal.

      p = inputParser;
      p.addParameter('origin', 'global');
      p.addParameter('removeNan', false);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Translate to local coordinates
      switch p.Results.origin
        case 'global'
          x0 = shape.global2local(x0);
          x1 = shape.global2local(x1);
        case 'local'
          % Nothing to do
        otherwise
          error('Unknown origin specified');
      end

      % Defer to shape-specific implementation
      [locs, norms, dist] = shape.intersectAllInternal(x0, x1, unmatched{:});

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
      
      if p.Results.removeNan && size(x0, 2) == 1 && size(x1, 2) == 1
        rmVals = any(isnan(locs), 1);
        locs(:, rmVals) = [];
        dist(:, rmVals) = [];
        norms(:, rmVals) = [];
      end
    end

    function [locs, norms] = intersect(shape, x0, x1, varargin)
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
      
      assert(isnumeric(x0) && size(x0, 1) == 3 && ismatrix(x0), ...
          'x0 must be a 3xN matrix');
      assert(isnumeric(x1) && size(x1, 1) == 3 && ismatrix(x1), ...
          'x1 must be a 3xN matrix');
      assert(size(x0, 2) == size(x1, 2), ...
          'x0 and x1 must be the same size');

      % Translate to local coordinates
      switch p.Results.origin
        case 'global'
          x0 = shape.global2local(x0);
          x1 = shape.global2local(x1);
        case 'local'
          % Nothing to do
        otherwise
          error('Unknown origin specified');
      end

      % Defer to shape-specific implementation
      [locs, norms] = shape.intersectInternal(x0, x1, unmatched{:});

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
      %   - axes (handle) -- Axes hanlde to place plot in.
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

  methods (Sealed)
    function varargout = surf(shape, varargin)
      % Generate a visualisation of the shape
      %
      % Usage
      %   shape.surf(...)
      %   Display visualisation of shape(s).
      %
      %   S = shape.surf(...)
      %   Returns a array of handles to the generated patches or when
      %   no visualisation is enabled, creates a cell array of structures
      %   that can be passed to patch.
      %
      % Optional named parameters
      %   - axes (handle) -- axis to place surface in (default: gca)
      %
      %   - surfOptions (cell) -- options to be passed to patch (default: {})
      %
      %   - showNormals (logical) -- Show surface normals (default: false)
      %
      %   - origin (enum) -- Coordinate origin for drawing.
      %     Can be 'global' or 'local'  Default: 'global'.
      %
      %   - visualise (logical) -- Show the visualisation.
      %     Default: ``true``
      %
      %   - normalScale (numeric) -- Scale for normal vectors.
      %     Default: ``0.1``.
      %
      % Additional parameters passed to :meth:`surfInternal`.

      p = inputParser;
      p.addParameter('surfOptions', {});
      p.addParameter('axes', []);
      p.addParameter('showNormals', false);
      p.addParameter('origin', 'global', ...
          @(x) sum(strcmpi(x, {'local', 'global'}) == 1));
      p.addParameter('visualise', true);
      p.addParameter('normalScale', 0.1);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Get patch from child methods
      pch = cell(1, numel(shape));
      for ii = 1:numel(shape)
        pch{ii} = shape(ii).surfInternal(unmatched{:});

        % Translate patch to desired origin
        if strcmpi(p.Results.origin, 'global')
          pch{ii}.Vertices = shape(ii).local2global(pch{ii}.Vertices.').';
        end
      end

      % Generate visualisation
      if p.Results.visualise

        % Place the surface in the specified axes
        our_axes = p.Results.axes;
        if isempty(our_axes)
          our_axes = gca();
        end

        % Patch doesn't watch for hold, so clear it ourselves
        isholdon = ishold(our_axes);
        if ~isholdon
          cla(our_axes);
        end

        % Convert patch structures to patches
        sp = matlab.graphics.primitive.Patch.empty(1, 0);
        for ii = 1:numel(pch)
          sp(ii) = patch(our_axes, pch{ii}, p.Results.surfOptions{:});
        end
        pch = sp;

        % Clear the orientation/aspect if hold isn't on
        if ~isholdon
          view(our_axes, [60, 30]);
          daspect(our_axes, [1, 1, 1]);
        end

        if p.Results.showNormals

          % Preserve hold-no state
          if ~isholdon
            hold(our_axes, 'on');
          end

          % Add quiver plot for each patch
          for ii = 1:numel(pch)

            % Find mean of each patch face
            mXyz = zeros(3, size(pch(ii).Faces, 1));
            for jj = 1:size(pch(ii).Faces, 2)
              mXyz = mXyz + pch(ii).Vertices(pch(ii).Faces(:, ii), :).';
            end
            mXyz = mXyz ./ size(pch(ii).Faces, 2);

            % Calculate normals
            nxyz = shape(ii).normalsXyz(mXyz);

            % Generate plot of surface normals
            s = p.Results.normalScale;
            quiver3(our_axes, mXyz(1, :), mXyz(2, :), mXyz(3, :), ...
                s.*nxyz(1, :), s.*nxyz(2, :), s.*nxyz(3, :), 0);
          end

          if ~isholdon
            hold(our_axes, 'off');
          end
        end
      end

      if nargout ~= 0
        varargout{1} = pch;
      end
    end

    function shape = scale(shape, sc, varargin)
      % Scale the shape by a factor
      %
      % Usage
      %   shape = shape.scale(sc, ...)
      %
      % Optional named parameters
      %   - origin (enum) -- Either 'local' or 'global'.  If 'global',
      %     scales the position coordinate too (default).

      p = inputParser;
      p.addParameter('origin', 'global');
      p.parse(varargin{:});

      assert(isnumeric(sc) && isscalar(sc), ...
          'sc must be a numeric scalar');

      origin = p.Results.origin;
      assert(any(strcmpi(origin, {'local', 'global'})), ...
          'origin must be local or global');

      for ii = 1:numel(shape)
        shape(ii) = shape(ii).scaleInternal(sc);

        if strcmpi(origin, 'global')
          shape(ii).position = shape(ii).position * sc;
        end
      end
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
      default_object = ott.shape.Empty;
    end
  end
end
