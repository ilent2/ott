classdef Ray < ott.optics.beam.PlaneWave
% Geometric optics ray collection.
% Inherits from :class:`ott.optics.beam.PlaneWave`.
%
% Provides the methods for calculating scattering, forces and torques
% from Ray representations of plane wave beams.
%
% Unlike the PlaneWave beam this object inherits from, the Ray set
% power is finite.
%
% Properties
%   - origin        -- Ray origins, n-dimensional array with 3 rows
%   - direction     -- Ray directions, n-dimensional array with 3 rows
%   - field         -- Field parallel and perpendicular to polarisation
%   - polarisation  -- Polarisation, n-dimensional array with 3 rows
%   - power         -- Finite power of the ray set
%
% Methods
%   - focus         -- Focus rays to a point
%   - scatter       -- Calculate how rays are scattered
%   - rotate        -- Rotate the ray and the polarisation
%   - rotate*       -- Rotate the particle around the X,Y,Z axis
%
% Static methods
%   - FromNearfield -- Construct a Ray instance from a near-field beam slice
%   - FromFarfield  -- Construct a Ray instance from beam's far-field
%   - FromParaxial  -- Construct a Ray instance from beam's paraxial far-field
%
% Loosely based on the go.Ray class from OTGO.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

% TODO: Review other features in OTGO (not just in the Ray class)
% TODO: snellslaw

  methods
    function ray = Ray(origin, direction, power, polarisation)
      % Construct a new geometric optics ray instance
      %
      % Usage
      %   ray = Ray(origin, direction, power, polarisation)
      %
      %   ray = Ray(vector, power, polarisation)
      %
      % Parameters
      %   - origin (numeric) -- ray origins, 3 row n-dimensional
      %   - direction (numeric) -- ray directions, 3 row, n-dimensional
      %   - power (numeric) -- ray power, 1 row, n-dimensional
      %   - polarisation (numeric) -- polarisation, 3 row, n-dimensional
      %
      %   - vector (numeric) -- ray directions or ray vectors.
      %     Must be either a :class:`ott.utils.Vector` or a 3 row or
      %     6 row n-dimensional array.

      % Handle arguments and call base constructor
      if nargin == 3
        polarisation = power;
        power = direction;
        vecArgs = {origin};
      else
        vecArgs = {origin, direction};
      end

      % Call the base constructor (must be outside if statement)
      ray = ray@ott.utils.Vector(vecArgs{:});

      % Store polarisation and power
      ray.power = power;
      ray.polarisation = polarisation;
    end

    function ray = focus(ray, location)
      % Focus rays to a point
      %
      % Rotates each ray to face the specified location.
      % For this method to do anything sensible, the origins of
      % the origin rays must be set.
      %
      % Usage
      %   fray = ray.focus(location)
      %
      % Parameters
      %   - location (3x1 numeric) -- Location to focus rays towards.

      rdir = cross(ray, location - ray.origin);
      rmat = ott.utils.rotation_matrix(rdir.direction);
      ray = ray.rotate(rmat);
    end

    function [rray, tray, perp] = scatter(ray, normals, n1, n2, varargin)
      % Calculate rays scattered by a surface using Snell's law
      %
      % Usage
      %   [rray, tray] = iray.scatter(normals, n1, n2)
      %   Calculates the transmitted and reflected rays
      %   for scattering by set of planes described by their normal vectors.
      %
      %   [rray, tray, perp] - iray.scatter(...) as above, but outputs
      %   the vectors perpendicular to the scatter plane.
      %
      % The scattered rays have one additional dimension compared to
      % the incident rays, the size of the dimension corresponds to the
      % number of surface normals.  For example, if the incident ray
      % object is 3x5, the scattered ray objects will be 3x5xN.
      %
      % Parameters
      %   - normals (3xN numeric) -- surface normals
      %   - n1, n2 (numeric) -- Refractive indices in the first
      %     and second medium.  Can be array.
      %
      % Optional named arguments
      %   - individual (logical) -- If true, each ray is scattered by
      %     a separate normal.  Size of rays must match size of normals.
      %     Default: ``false``.
      %
      %   - position (enum|3xN numeric) -- Positions to use for scattered
      %     ray origins.  If the normals correspond to intersections with
      %     a surface, positions is the corresponding position.
      %     If position is [], uses the original ray positions, if
      %     position is 'zero', uses the coordinate origin.
      %     Default: ``[]``.  Use original ray positions.
      %
      %   - internal (logical) -- If the refractive indices should be flipped.
      %     Can be scalar or array with same size as normals.
      %     Default: ``false``.

      p = inputParser;
      p.addParameter('individual', false, @(x) islogical(x) & isscalar(x));
      p.addParameter('position', []);
      p.addParameter('internal', false, @islogical);
      p.parse(varargin{:});

      % Check required inputs
      assert(isnumeric(n1), ...
          'n1 must be a numeric scalar');
      assert(isnumeric(n2), ...
          'n2 must be a numeric scalar');
      assert(isnumeric(normals) && size(normals, 1) == 3, ...
          'normals must be a 3xN numeric array');

      % Normalise normals and ray directions
      normals = normals ./ vecnorm(normals, 2, 1);
      dirs = ray.direction ./ vecnorm(ray.direction, 2, 1);

      % Ensure normals and dirs have correct size/shape (duplicate as needed)
      if p.Results.individual
        szdirs = size(dirs);
        sznorms = size(normals);
        normals = reshape(normals, [3, szdirs(2:end), sznorms(2:end)]);
        normals = repmat(normals, ...
            [1, szdirs(2:end), ones(1, numel(sznorms)-1)]);
        dirs = repmat(dirs, ...
            [1, ones(1, numel(szdirs)-1), sznorms(2:end)]);
      else
        szdirs = size(dirs);
        sznorms = size(normals);
        assert(all(szdirs == sznorms), ...
          'Array sizes must match for individual normals');
      end

      % Calculate incident angle and s-direction
      sdirection = cross(dirs, normals);
      iangle = real(asin(vecnorm(sdirection, 2, 1)));   % [0, pi]

      % Ensure incident angle is between 0 and pi/2
      iangle(iangle > pi/2) = pi - iangle(iangle > pi/2); % [0, pi/2]

      % Handle swapping of internal/external indices
      if p.Results.internal
        if isscalar(p.Results.internal)
          [n2, n1] = deal(n1, n2);
        else
          old_n1 = n1;
          old_n2 = n2;

          if isscalar(n1)
            old_n1 = repmat(n1, size(p.Results.internal));
          end
          if isscalar(n2)
            old_n2 = repmat(n2, size(p.Results.internal));
          end

          n1 = old_n1;
          n2 = old_n2;
          n1(p.Results.internal) = old_n2(p.Results.internal);
          n1(p.Results.internal) = old_n1(p.Results.internal);
        end
      end

      % Calculate Fresnel coefficients in s and p directions
      [Rs,Ts,rs,ts] = Ray.rtcoeffs(iangle,n1,n2);
      [Rp,Tp,rp,tp] = Ray.rtcoeffp(iangle,n1,n2);

      % Calculate reflected directions
      ndirs = dot(dirs, normals).*normals;
      rdirs = 2.*ndirs - dirs;

      % Calculate transmitted direction
      bdirs = dirs - ndirs;
      tdirs = dirs + (n1./n2.*sin(iangle) - bdirs).*bdirs./vecnorm(bdirs);

      % TODO: How do we calculate polarisation?
      % TODO: Calculate power
      % TODO: TIR

      % Calculate transmission angle
      % TODO: Do we need this?
      tangle = asin(n1./n2.*sin(iangle));

      % Construct output rays
      rray = ott.optics.geometric.Ray(ray.position, rdirs, ...
          power, polarisation);
      tray = ott.optics.geometric.Ray(ray.position, tdirs, ...
          power, polarisation);

      % Override position
      if ~isempty(p.Results.position)
        if isnumeric(p.Results.position)
          rray.position = p.Results.position;
          tray.position = p.Results.position;
        elseif strcmpi(p.Results.position, 'none')
          rray.position = zeros(size(rray.position));
          tray.position = zeros(size(tray.position));
        else
          error('Unsupported position value');
        end
      end

      % Generate vector structure for perp (if requested)
      if nargout == 3
        perp = ott.utils.Vector(tray.position, sdirection);
      end
    end
  end

  methods (Hidden)
    function p = getBeamPower(beam)
      % Returns finite power of the ray set
      p = sum(abs(beam.fields(:)));
    end

    function E = efieldInternal(beam, xyz, varargin)
      % Calculate the E-field
      %
      % Usage
      %   E = beam.efieldInternal(xyz, ...)
      %
      % Optional named arguments
      %   - method (enum) -- Method to use when drawing rays
      %     Can be one of
      %     - inf -- Infinite extend plane waves
      %     - ray_invr -- Power scaled with 1/distance from centre
      %     - invr -- Power scaled with 1/distance from origin
      %     Default: ``'ray_invr'``.

      % Change default parameters
      p = inputParser;
      p.addParameter('method', 'ray_invr');
      p.parse(varargin{:});

      E = efieldInternal@ott.optics.beam.PlaneWave(beam, xyz, ...
          'method', p.results.method);
    end
  end
end

