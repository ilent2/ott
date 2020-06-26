classdef SphericalLens
% Finite aperture spherical lens approximation.
% Inherits from :class:`ott.shapes.Plane`.
%
% The spherical lens focusses rays which hit the sphere to the centre
% of the sphere.  The radius of the sphere is the lens focal length.
%
% Properties
%   - focal_length      -- Focal length of the lens
%   - front_surface     -- The front surface of the lens (Plane)
%   - lens_surface      -- The lens surface (Sphere)
%   - normal            -- Normal to the lens front surface
%   - position          -- Position of the lens front surface
%
% Methods
%   - surf              -- Draw the lens
%   - intersection      -- Calculates the intersection locations
%   - scatter           -- Calculate the transmitted beam

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  properties (SetAccess=protected)
    focal_length        % Focal length of the lens
  end

  properties (Dependent)
    lens_surface        % Sphere object for the lens surface
    front_surface       % A Plane representing the lens front surface
    focal_position      % Lens focal position
  end

  methods
    function lens = SphericalLens(focal_length, varargin)
      % Construct a new spherical lens
      %
      % Usage
      %   lens = SphericalLens(focal_length, ...)
      %
      % Parameters
      %   - focal_length (numeric) -- Focal length of the lens.
      %     Must be positive.
      %
      % Optional named arguments
      %   - normal (3x1 numeric) -- Normal to the lens front surface.
      %     Default: ``[0;0;1]``.
      %
      %   - position (3x1 numeric) -- Position of the lens front surface.
      %     Default: ``[0;0;0]``.  Ignored if focal_position is set.
      %
      %   - focal_position (3x1 numeric) -- Location for the lens focal
      %     point.  Default: ``[]``.  Directly related to position.

      p = inputParser;
      p.addParameter('position', [0;0;0]);
      p.addParameter('normal', [0;0;1]);
      p.addParameter('focal_position', []);
      p.parse(varargin{:});

      lens = lens@ott.shapes.Plane(p.Results.normal, 0.0);
      lens.position = p.Results.position;
      lens.focal_length = focal_length;
      if ~isempty(p.Results.focal_position)
        lens.focal_position = p.Results.focal_position;
      end
    end

    function varargout = surf(lens, varargin)
      % Generate a visualisation of the lens
      %
      % Usage
      %   lens.surf(...)
      %
      % Optional named arguments
      %   - surf_type (enum) -- Type of surface to draw.  Can either be
      %     'sphere', 'hemi-sphere', or 'plane' (for the front plane).
      %     Default: ``'sphere'``.
      %
      % For other optional arguments, see base class.

      % TODO: This method should change to getGeometry

      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('surf_type', 'sphere');
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      switch p.Results.surf_type
        case 'sphere'
          [varargout{1:nargout}] = lens.lens_surface.surf(unmatched{:});
        case 'hemi-sphere'
          error('Not yet implemented');
        case 'plane'
          [varargout{1:nargout}] = surf@ott.shapes.Plane(lens, unmatched{:});
        otherwise
          error('Unknown surf_type parameter value');
      end
    end

    function varargout = intersect(lens, vecs)
      % Calculate the intersection point on the sphere surface
      %
      % Usage
      %   [locs, norms] = shape.intersect(vec)
      %   Returns a 3xN matrix of intersection locations or nan.

      [varargout{1:nargout}] = lens.lens_surface.intersect(vecs);
    end

    function tbeam = scatter(lens, ibeam)
      % Calculate the transmitted beam
      %
      % Usage
      %   tbeam = lens.scatter(ibeam)
      %
      % Calculates the beam transmitted through the thin lens.

      % Calculate intersection locations
      locs = lens.intersect(ibeam);

      % Get focal position
      fp = lens.focal_position;

      % Calculate new directions
      directions = sign(lens.focal_length)*(fp - locs);
      directions = directions ./ vecnorm(directions);

      % Create transmitted beam
      tbeam = ott.beam.ScatteredRay(ibeam, 'origin', locs, ...
          'direction', directions, 'like', ibeam);
    end
  end

  methods % Getters/setters
    function pln = get.front_surface(lens)

      % Return a copy of ourselves but with the Plane type
      pln = ott.shapes.Plane(lens.normal, 0.0);
      pln.position = lens.position;
      pln.rotation = lens.rotation;
    end

    function sph = get.lens_surface(lens)

      % Return a new sphere instance
      sph = ott.shapes.Sphere(lens.focal_length);
      sph.position = lens.focal_position;

    end

    function lens = set.focal_length(lens, val)
      assert(isnumeric(val) && isscalar(val) && val > 0, ...
        'focal_length must be positive numeric scalar');
      lens.focal_length = val;
    end

    function val = get.focal_position(lens)
      val = lens.position + lens.focal_length .* lens.normal;
    end
    function lens = set.focal_position(lens, val)
      lens.position = val - lens.focal_length .* lens.normal;
      lens.lens_surface.position = lens.focal_position;
    end
  end
end

