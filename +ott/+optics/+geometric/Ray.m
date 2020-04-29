classdef Ray < ott.utils.Vector
% Geometric optics ray
% Inherits from :class:`ott.utils.Vector`.
%
% Properties
%   - origin        -- Ray origins, n-dimensional array with 3 rows
%   - direction     -- Ray directions, n-dimensional array with 3 rows
%   - power         -- Power, n-dimensional array with 1 row
%   - polarisation  -- Polarisation, n-dimensional array with 3 rows
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

% TODO: It would be nice to keep track of Ray phases

  properties
    power          % Power, n-dimensional array with 1 row
    polarisation   % Polarisation, n-dimensional array with 3 rows
  end

  methods (Static)
    function ray = FromFieldVectors(E, H)
      % Construct a new Ray from E and H field vectors
      %
      % Usage
      %   ray = FromFieldVectors(E, H)

      assert(isa(E, 'ott.utils.FieldVector'), 'E must be a FieldVector');
      assert(isa(H, 'ott.utils.FieldVector'), 'H must be a FieldVector');

      % Get coordinates
      if strcmpi(E.type, 'spherical')
        xyz = ott.utils.rtp2xyz(E.locations);
      elseif strcmpi(E.type, 'cartesian')
        % Nothing to do
      else
        error('Unsupported field vector type');
      end

      Exyz = E.vxyz;
      Hxyz = H.vxyz;

      % Calculate Poynting vector
      S = 0.5.*real(cross(Exyz, conj(Hxyz)));
      Snorm = vecnorm(S, 2, 1);

      % Filter points with no intensity
      mask = Snorm == 0;
      xyz(:, mask) = [];
      S(:, mask) = [];
      Snorm(:, mask) = [];

      % Normalise S afterwards (to avoid nan)
      S = S ./ Snorm;

      % Calculate orthogonal vectors (if present/needed)
      Ex = real(Exyz(:, mask));
      Ey = imag(Eyxz(:, mask));
      Exnorm = vecnorm(Ex, 2, 1) ./ vecnorm(abs(Exyz(:, mask)), 2, 1);
      Eynorm = vecnorm(Ey, 2, 1) ./ vecnorm(abs(Eyxz(:, mask)), 2, 1);
      Snorm = [Snorm, Snorm] .* [Exnorm, Eynorm];

      % Filter properties (and duplicate as required)
      xyz = [xyz(:, Exnorm ~= 0), xyz(:, Eynorm ~= 0)];
      S = [S(:, Exnorm ~= 0), S(:, Eynorm ~= 0)];
      Snorm = Snorm(:, [Exnorm ~= 0, Eynorm ~= 0]);
      pol = [Ex(:, Exnorm ~= 0), Ey(:, Eynorm ~= 0)];

      % Construct ray object
      ray = ott.optics.geometric.Ray(xyz, S, normS, pol);
    end

    function ray = FromFarfield(beam)
      % Constructs a new Ray instance from the beam far-field
      %
      % Uses the ``ehfarfield()`` methods from a :class:`ott.optics.beam.Beam`
      % to sample a slice of the near-field.  The ray direction and intensity
      % is set from the cross product of the E and H fields.
      %
      % Each sample location produces zeros, one or two rays, corresponding
      % to the real and imaginary parts of E.  If either the real or
      % imaginary part is zero, only one ray is produced.
      %
      % Usage
      %   ray = FromFarfield(beam, ...)

      % Calculate fields on a grid
      [xx, yy, zz] = sphere(40);
      xyz = [xx(:), yy(:), zz(:)].';
      rtp = ott.utils.xyz2rtp(xyz);
      [E, H] = beam.ehfarfield(rtp);

      % Generate from field vectors
      ray = ott.optics.geometric.Ray.FromFieldVectors(E, H);
    end

    function ray = FromNearfield(beam)
      % Constructs a new Ray instance from a near-field beam slice
      %
      % Uses the ``ehfield()`` methods from a :class:`ott.optics.beam.Beam`
      % to sample a slice of the near-field.  The ray direction and intensity
      % is set from the cross product of the E and H fields.
      %
      % Each sample location produces zero, one or two rays, corresponding
      % to the real and imaginary parts of E.  If either the real or
      % imaginary part is zero, only one ray is produced.
      %
      % Usage
      %   ray = FromNearfield(beam, ...)

      % Calculate fields on a grid
      [xx, yy, zz] = meshgrid(linspace(-1, 1), linspace(-1, 1), 0);
      xyz = [xx(:), yy(:), zz(:)].';
      [E, H] = beam.ehfield(xyz);

      % Generate from field vectors
      ray = ott.optics.geometric.Ray.FromFieldVectors(E, H);
    end

    function ray = FromParaxial(beam)
      % Constructs a new Ray instance from the beam paraxial far-field
      %
      % Uses the ``ehparaxial()`` methods from a :class:`ott.optics.beam.Beam`
      % to sample a slice of the near-field.  The ray direction and intensity
      % is set from the cross product of the E and H fields.
      %
      % Each sample location produces zero, one or two rays, corresponding
      % to the real and imaginary parts of E.  If either the real or
      % imaginary part is zero, only one ray is produced.
      %
      % Usage
      %   ray = FromNearfield(beam, ...)

      % Calculate fields on a grid
      [xx, yy, zz] = meshgrid(linspace(-1, 1), linspace(-1, 1), 0);
      xyz = [xx(:), yy(:), zz(:)].';
      [E, H] = beam.ehparaxial(xyz);

      % Generate from field vectors
      ray = ott.optics.geometric.Ray.FromFieldVectors(E, H);
    end
  end

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

    % TODO: Rotations, translations, angle, snellslaw
    % TODO: beam2rays beam2focussed
    % TODO: disp method

    function varargout = plot(vec, varargin)
      % Plots the vector set in 3-D.
      %
      % Uses the quiver function to generate a visualisation of the
      % vector set.
      %
      % Usage
      %   h = plot(vec, ...)
      %
      % Optional named arguments
      %   - Scale (numeric) -- rescales the coordinates and components
      %     of the vector before plotting.  Can either be a scalar
      %     or vector ``[S1, S2]`` specifying separate scaling for the
      %     coordinates and components.  Default: ``[1, 1]``.
      %
      % Any unmatched named arguments are applied to the plot handles
      % returned by the quiver function calls.

      % Parse inputs
      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('Scale', [1, 1]);
      p.parse(varargin{:});

      S1 = p.Results.Scale(1);
      S2 = p.Results.Scale(2);

      isholdon = ishold();

      % Generate plot of directions
      h = quiver3(S1*vec.origin(1, :), S1*vec.origin(2, :), ...
          S1*vec.origin(3, :), S2*vec.direction(1, :), ...
          S2*vec.direction(2, :), S2*vec.direction(3, :), 0);

      if ~isholdon
        hold('on');
      end

      % Generate plot of polarisations
      h(2) = quiver3(S1*vec.origin(1, :), S1*vec.origin(2, :), ...
          S1*vec.origin(3, :), S2*vec.polarisation(1, :), ...
          S2*vec.polarisation(2, :), S2*vec.polarisation(3, :), 0);

      if ~isholdon
        hold('off');
      end

      % Apply unmatched arguments to plot handle
      unmatched = [fieldnames(p.Unmatched).'; struct2cell(p.Unmatched).'];
      if ~isempty(unmatched)
        set(h, unmatched{:});
      end

      % Assign outputs
      if nargout > 0
        varargout{1} = h;
      end
    end
  end

  methods
    function ray = set.power(ray, val)

      assert(size(val, 1) == 1, 'power must have 1 rows');

      szpol = size(val);
      szpol = szpol(2:end);
      if numel(szpol) == 1
        szpol = [szpol, 1];
      end
      assert(all(szpol == size(ray)), ...
          'Power should have same number of elements and size as rays');

      assert(isreal(val) && all(val(:) >= 0.0), ...
          'Power must be positive and real');
      ray.power = val;
    end

    function ray = set.polarisation(ray, val)
      
      % Only interested in the direction
      if isa(val, 'ott.utils.Vector')
        val = val.direction;
      end

      assert(size(val, 1) == 3, 'polarisation must have 3 rows');
      assert(all(size(ray.origin) == size(val)), ...
          'origin and polarisation must have same size');
      assert(isnumeric(val) && isreal(val), ...
          'polarisation must be real numeric matrix');

      ray.polarisation = val;
    end
  end
end

