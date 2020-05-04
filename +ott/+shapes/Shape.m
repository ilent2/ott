classdef (Abstract) Shape < ott.utils.RotateHelper
% Shape abstract class for optical tweezers toolbox shapes.
% Inherits from :class:`ott.utils.RotateHelper`.
%
% Properties
%   - maxRadius  -- maximum distance from shape origin
%   - volume     -- volume of shape
%   - position   -- Location of shape ``[x, y, z]``
%   - rotation   -- Orientation of the particle (3x3 rotation matrix)
%
% Methods (abstract)
%   insideRtp       -- Determine if Spherical point is inside shape
%   insideXyz       -- Determine if Cartesian point is inside shape
%   normalsRtp      -- Calculate normals at surface location
%   normalsXyz      -- Calculate normals at surface location
%   get_maxRadius   -- Get max distance from origin
%   get_volume      -- Get shape volume
%
% Methods
%   writeWavefrontObj(shape, ...) write shape to Wavefront OBJ file
%       only implemented if shape supports this action.
%   simple(...) simplified constructor for shape-like objects.
%  - rotate      -- Rotate the particle specifying a rotation matrix
%  - rotate*     -- Rotate the particle around the X,Y,Z axis
%  - intersect   -- Calculate intersection between vectors and surface
%
% See also simple, ott.shapes.Cube, ott.shapes.TriangularMesh.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% TODO: Use rotations where appropriate

  properties
    position       % Location of shape ``[x, y, z]``
    rotation       % Orientation of the shape (3x3 matrix)
  end

  methods (Static)

    function shape = simple(name, parameters)
			%SIMPLE constructs a shape for a bunch of simple particles
      % This method replaces shapesurface.m from OTTv1.
      %
      % Supported shapes [parameters]:
      %   'sphere'          Sphere [ radius ]
      %   'cylinder'        z-axis aligned cylinder [ radius height ]
      %   'ellipsoid'       Ellipsoid [ a b c]
      %   'superellipsoid'  Superellipsoid [ a b c e n ]
      %   'cone-tipped-cylinder'      [ radius height cone_height ]
      %   'cube'            Cube [ width ]
      %   'axisym'          Axis-symetric particle [ rho(:) z(:) ]
      %   'stl'             load STL file [filename]
      %   'obj'             load Wavefront OBJ file [filename]

      switch lower(name)
        case 'sphere'
          shape = ott.shapes.Sphere(parameters(:));
        case 'cylinder'
          shape = ott.shapes.Cylinder(parameters(:));
        case 'ellipsoid'
          shape = ott.shapes.Ellipsoid(parameters(:));
        case 'superellipsoid'
          shape = ott.shapes.Superellipsoid(parameters(:));
        case 'cone-tipped-cylinder'

          radius = parameters(1);
          height = parameters(2);
          cone_height = parameters(3);

          z = [ height/2 + cone_height, height/2, ...
              -height/2, -height/2 - cone_height];
          rho = [ 0.0, radius, radius, 0.0 ];

          shape = ott.shapes.AxisymLerp(rho, z);
        case 'cube'
          shape = ott.shapes.Cube(parameters(:));
        case 'axisym'
          shape = ott.shapes.AxisymLerp(parameters(:));
        case 'obj'
          shape = ott.shapes.WavefrontObj(parameters(:));
        case 'stl'
          shape = ott.shapes.StlLoader(parameters(:));
        otherwise
          error('Unsupported simple particle shape');
      end
    end
  end

  methods (Static, Hidden)
    function xyz = insideXyzParseArgs(origin, varargin)
      % Helper for parsing arguments for insideXyz functions

      % Collect arguments
      p = inputParser;
      p.addRequired('x');
      p.addOptional('y', [], @isnumeric);
      p.addOptional('z', [], @isnumeric);
      p.addParameter('origin', 'world');
      p.parse(varargin{:});

      % Handle xyz or [x, y, z] arguments
      if isempty(p.Results.x)
        error('Must supply either xyz or x, y and z');
      elseif isempty(p.Results.y) && isempty(p.Results.z)
        xyz = p.Results.x;
      else
        x = p.Results.x(:);
        y = p.Results.y(:);
        z = p.Results.z(:);
        [x, y, z] = ott.utils.matchsize(x, y, z);
        xyz = [x y z].';
      end

      % Translate to shape origin
      if strcmpi(p.Results.origin, 'world')
        xyz = xyz - origin(:);
      elseif strcmpi(p.Results.origin, 'shape')
        % Nothing to do
      else
        error('origin must be ''world'' or ''shape''');
      end
    end

    function rtp = insideRtpParseArgs(origin, varargin)
      % Argument parser helper for insideRtp function

      % Collect arguments
      p = inputParser;
      p.addRequired('r');
      p.addOptional('t', [], @isnumeric);
      p.addOptional('p', [], @isnumeric);
      p.addParameter('origin', 'world');
      p.parse(varargin{:});

      % Handle rtp or [r, t, p] arguments
      if isempty(p.Results.r)
        error('Must supply either rtp or r, t and p');
      elseif isempty(p.Results.t) && isempty(p.Results.p)
        rtp = p.Results.r;
      else
        r = p.Results.r(:);
        theta = p.Results.t(:);
        phi = p.Results.p(:);
        [r, theta, phi] = ott.utils.matchsize(r, theta, phi);
        rtp = [r theta phi].';
      end

      % Translate to shape origin
      if strcmpi(p.Results.origin, 'world')

        % Only do work if we need to
        if vecnorm(origin) ~= 0
          xyz = ott.utils.rtp2xyz(rtp);
          xyz = xyz - origin(:);
          rtp = ott.utils.xyz2rtp(xyz);
        end

      elseif strcmpi(p.Results.origin, 'shape')
        % Nothing to do
      else
        error('origin must be ''world'' or ''shape''');
      end

      assert(all(rtp(1, :) >= 0), 'Radii must be positive');
    end
  end

  properties (Dependent)
    maxRadius       % Maximum particle radius (useful for Nmax calculation)
    volume          % Volume of the particle
  end

  methods (Abstract)
    insideRtp(shape, varargin)   % Determine if point is inside shape
    insideXyz(shape, varargin)   % Determine if point is inside shape

    get_maxRadius(shape, varargin)    % Get max distance from origin
    get_volume(shape, varargin);      % Get shape volume
  end

  methods (Access=protected)

    function writeWavefrontObj_helper(shape, filename, verts, faces)
      % Helper to write faces and vertices to a file
      %
      % helper(filename, verts, faces)
      %   verts 3xN array of floats
      %   faces mxN array of integers (starting at 1)

      fp = fopen(filename, 'w');

      % Write documentation
      fprintf(fp, '# Shape description generated by OTT\n');

      % Write verts
      for ii = 1:size(verts, 2)
        fprintf(fp, 'v %f %f %f\n', verts(:, ii));
      end

      % Write faces
      for ii = 1:size(faces, 2)
        fprintf(fp, 'f ');
        fprintf(fp, '%d ', faces(:, ii));
        fprintf(fp, '\n');
      end

      % Close file
      fclose(fp);

    end
  end

  methods
    function shape = Shape(varargin)
      % Construct a new shape instance.
      %
      % This class cannot be instanced directly, use one of the other
      % shape descriptions to create a new shape.
      %
      % Usage
      %   shape = Shape(...)
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

      % Store initial position and rotation
      shape.position = p.Results.position;
      shape.rotation = p.Results.rotation;
    end

    function shape = rotate(shape, R)
      % Apply the rotation matrix to the shapes internal rotation

      % Check at least one output
      shape.nargoutCheck(nargout);

      % Apply rotation
      shape = R * shape.rotation;
    end

    function r = get.maxRadius(shape)
      % Get the particle max radius
      r = shape.get_maxRadius();
    end

    function r = get.volume(shape)
      % Get the particle volume
      r = shape.get_volume();
    end

    function writeWavefrontObj(shape, filename)
      % Write representation of shape to Wavefront OBJ file
      %
      % This is a placeholder, not supported by all shape types
      error('Shape does not support writing to WavefrontOBJ file');
    end

    function surf(shape, varargin)
      % SURF generate a visualisation of the shape
      %
      % SURF() displays a visualisation of the shape in the current figure.
      %
      % SURF(..., 'surfoptions', {varargin}) specifies the options to
      % pass to the surf function.
      error('Shape does not support surf visualisation');
    end

    function [locs, norms] = intersect(shape, vecs)
      % Calculate the intersection between the shape and a vector set.
      %
      % If the vector is inside the shape, finds the intersect between
      % the shapes surface.
      %
      % Usage
      %   [locs, norms] = shape.intersect(vec)
      %
      % Parameters
      %   - vec (utils.Vector) -- A vector or type that can be cast
      %     to a Vector.
      %
      % The default implementation uses insideXyz to test when the
      % rays are inside the shape.  For some shapes there may be
      % faster methods.  Normals are calculated with the normalsXyz method.

      if ~isa(vecs, 'ott.utils.Vector')
        vecs = ott.utils.Vector(vecs);
      end

      dx = shape.maxRadius ./ 100;
      if ~isfinite(dx)
        dx = 1.0e-3;
      end

      % Start by calculating intersects with bounding box
      ints = shape.intersectBoundingBox(vecs);
      ints(1:3, isnan(ints(1, :))) = vecs.origin(:, isnan(ints(1, :)));

      % Calculate distance we want to search for each vector
      search_distance = vecnorm(ints(4:6, :) - ints(1:3, :));

      % Distance from bounding box to intersection
      locs = zeros(1, numel(vecs));
      orgs = ints(1:3, :);
      found = false(1, numel(vecs));
      dirs = vecs.direction ./ vecnorm(vecs.direction);
      
      % Determine which points were already inside
      was_inside = shape.insideXyz(orgs);

      % Find a point between the intersects in the shape
      remaining = locs < search_distance & ~found;
      while any(remaining)

        % March the ray
        locs(remaining) = locs(remaining) + dx;

        % Determine if new point is inside
        inside = shape.insideXyz(orgs(:, remaining) ...
          + locs(remaining).*dirs(:, remaining));
        found(remaining) = inside ~= was_inside(remaining);

        % Determine which points remain
        remaining = locs < search_distance & ~found;
      end

      % Remove logs that weren't found
      locs(~found) = nan;

      % Convert length to location vector
      locs = ints(1:3, :) + locs .* dirs;

      if nargout > 1
        norms = shape.normalXyz(locs);
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
      %   - 'visualise'     Show the visualisation (default: nargout == 0)
      %   - origin (enum) -- Coordinate system origin.  Either 'world'
      %     or 'shape' for world coordinates or shape coordinates.

      p = inputParser;
      p.addOptional('spacing', shape.maxRadius/10);
      p.addParameter('plotoptions', []);
      p.addParameter('visualise', nargout == 0);
      p.addParameter('scale', 1.0);
      p.addParameter('axes', []);
      p.addParameter('origin', 'shape');
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

      % Determine which points are inside
      mask = shape.insideXyz(xx / p.Results.scale, ...
          yy / p.Results.scale, zz / p.Results.scale, ...
          'origin', 'shape');
      xyz = [xx(mask).'; yy(mask).'; zz(mask).'];

      % Translate to world origin
      if strcmpi(p.Results.origin, 'world')
        xyz = xyz + shape.position;
      elseif strcmpi(p.Results.origin, 'shape')
        % Nothing to do
      else
        error('origin must be ''world'' or ''shape''');
      end

      % Visualise the result
      if p.Results.visualise
        
        % Get the axes to use
        our_axes = p.Results.axes;
        if isempty(our_axes)
          our_axes = axes();
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
  end

  methods (Hidden)
    function varargout = surfCommon(shape, p, sz, xx, yy, zz)
      % Helper for doing the final parts of surf
      %
      % Usage
      %   [varargout{1:nargout}] = shape.surfCommon(p, sz, X, Y, Z)
      %
      % Parameters
      %   - p -- inputParser
      %   - sz -- Desired size of shape values
      %   - X, Y, Z -- shape values

      % Apply a rotation to the surface
      if ~isempty(p.Results.rotation)
        XYZ = [xx(:), yy(:), zz(:)].';
        XYZ = p.Results.rotation * XYZ;
        xx = reshape(XYZ(1, :), sz);
        yy = reshape(XYZ(2, :), sz);
        zz = reshape(XYZ(3, :), sz);
      end

      % Apply a translation to the surface
      if ~isempty(p.Results.position)
        xx = xx + p.Results.position(1);
        yy = yy + p.Results.position(2);
        zz = zz + p.Results.position(3);
      end

      % Generate the surface
      if nargout == 0 || ~isempty(p.Results.axes)

        % Place the surface in the specified axes
        our_axes = p.Results.axes;
        if isempty(our_axes)
          our_axes = gca();
        end

        surf(our_axes, xx, yy, zz, p.Results.surfoptions{:});
        axis(our_axes, 'equal');

        % Compute and add surface normals
        if p.Results.show_normals
          nxyz = shape.normalsXyz([xx(:), yy(:), zz(:)].');

          % Preserve hold-no state
          isholdon = ishold();
          if ~isholdon
            hold('on');
          end

          % Generate plot of surface normals
          s = 0.1*shape.maxRadius;
          quiver3(xx(:).', yy(:).', zz(:).', ...
              s.*nxyz(1, :), s.*nxyz(2, :), s.*nxyz(3, :), 0);

          if ~isholdon
            hold('off');
          end
        end
      end

      % Set outputs if requested
      if nargout ~= 0
        varargout = { xx, yy, zz };
      end
    end

    function ints = intersectBoundingBox(shape, vecs)
      % Calculate the bounding box intersections for the vector
      %
      % Usage
      %   ints = intersectBoundingBox(vecs)
      %   vecs is a ott.utils.Vector
      %   ints is a 6xN matrix of intersection locations or nan

      dirs = vecs.direction ./ vecnorm(vecs.direction);
      orgs = vecs.origin;

      R = shape.maxRadius;
      if ~isfinite(R)
        ints = nan(6, size(dirs, 2));
        return;
      end

      % Calculate intersection with planes
      lz = orgs - (orgs(3, :) - R) .* dirs ./ dirs(3, :);
      uz = orgs - (orgs(3, :) + R) .* dirs ./ dirs(3, :);
      ly = orgs - (orgs(2, :) - R) .* dirs ./ dirs(2, :);
      uy = orgs - (orgs(2, :) + R) .* dirs ./ dirs(2, :);
      lx = orgs - (orgs(1, :) - R) .* dirs ./ dirs(1, :);
      ux = orgs - (orgs(1, :) + R) .* dirs ./ dirs(1, :);

      % Find which planes the ray intersects
      lxb = all(lx(2:3, :) < R & lx(2:3, :) > -R, 1);
      uxb = all(ux(2:3, :) < R & ux(2:3, :) > -R, 1);
      lyb = all(ly([1,3], :) < R & ly([1,3], :) > -R, 1);
      uyb = all(uy([1,3], :) < R & uy([1,3], :) > -R, 1);
      lzb = all(lz(1:2, :) < R & lz(1:2, :) > -R, 1);
      uzb = all(uz(1:2, :) < R & uz(1:2, :) > -R, 1);
      brr = [lxb; uxb; lyb; uyb; lzb; uzb];

      ints = zeros(6, size(dirs, 2));

      % Assign intersects
      for ii = 1:size(ints, 2)

        this_ints = zeros(6, 1);
        jj = 1;

        % Get which vector matched
        % This feels like KLUDGE
        if lxb(ii)
          this_ints(jj:jj+2) = lx(:, ii);
          jj = jj + 3;
        end
        if uxb(ii)
          this_ints(jj:jj+2) = ux(:, ii);
          jj = jj + 3;
        end
        if lyb(ii)
          this_ints(jj:jj+2) = ly(:, ii);
          jj = jj + 3;
        end
        if uyb(ii)
          this_ints(jj:jj+2) = uy(:, ii);
          jj = jj + 3;
        end
        if lzb(ii)
          this_ints(jj:jj+2) = lz(:, ii);
          jj = jj + 3;
        end
        if uzb(ii)
          this_ints(jj:jj+2) = uz(:, ii);
          jj = jj + 3;
        end

        % Filter out non-intersects
        if jj == 1
          ints(:, ii) = nan;

        elseif jj == 4

          % This shouldn't happen since we haven't checked direction yet
          warning('Unexpected intersection event');
          ints(1:3, ii) = nan;
          ints(4:6, ii) = this_ints(1:3);

        else

          % Sort the intercepts
          l1 = dot(this_ints(1:3) - orgs(:, ii), dirs(:, ii));
          l2 = dot(this_ints(4:6) - orgs(:, ii), dirs(:, ii));
          if l2 < l1
            this_ints = [this_ints(4:6); this_ints(1:3)];
          end

          % Ray originated from inside, discard other ray
          if l1 < 0 || l2 < 0
            this_ints(1:3) = nan;
          end

          % Store result
          ints(:, ii) = this_ints;
        end

      end
    end
  end

  methods % Getters/Setters
    function shape = set.position(shape, value)
      % Check position values
      assert(numel(value) == 3, 'Position must be 3 element vector');
      assert(isnumeric(value), 'Position must be numeric');
      shape.position = value(:);
    end
  end
end
