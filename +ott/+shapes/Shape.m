classdef Shape
% Shape abstract class for optical tweezers toolbox shapes
%
% Properties
%   - maxRadius  -- maximum distance from shape origin
%   - volume     -- volume of shape
%   - position   -- Location of shape ``[x, y, z]``
%
% Methods (abstract)
%   insideRtp(shape, ...) determine if Spherical point is inside shape
%   insideXyz(shape, ...) determine if Cartesian point is inside shape
%
% Methods
%   writeWavefrontObj(shape, ...) write shape to Wavefront OBJ file
%       only implemented if shape supports this action.
%   simple(...) simplified constructor for shape-like objects.
%
% See also simple, ott.shapes.Cube, ott.shapes.TriangularMesh.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    position = [0;0;0];       % Location of shape ``[x, y, z]``
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

    function shape = set.position(shape, value)
      % Check position values
      assert(numel(value) == 3, 'Position must be 3 element vector');
      assert(isnumeric(value), 'Position must be numeric');
      shape.position = value(:);
    end
  end
end
