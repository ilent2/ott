classdef (Abstract) ShapeForce
% Base class for the shape-induced force field approximation
%
% The approximation involves calculating the electric field along the
% surface of the particle and estimating the force according to::
%
%     \langle f \rangle = \frac{1}{4} |E|^2 \Delta\epsilon \hat{n}
%
% where :math:`\hat{n}` is the normal to the surface,
% :math:`|E|^2` is the field amplitude and :math:`\Delta\epsilon` is
% the difference in permittivity of the particle and medium.
%
% To be accurate, the above should use the total field (incident+scatted),
% however, a zero-th order approximation is to use just the incident
% field, this approximation was shown to work rather well in
%
%   D. B. Phillips, et al., Shape-induced force fields in optical trapping
%   Nature Photonics volume 8, pages400–405(2014)
%   https://doi.org/10.1038/nphoton.2014.74
%
% Properties
%   - epsilon     -- Difference in permittivity of particle and medium
%   - beam        -- An object with an efield method.
%
% Abstract methods
%   - locations   -- Returns a vector of Cartesian surface locations

  properties
    epsilon     % Difference in permittivity of particle and medium
    beam        % An object with an efield method.
  end

  methods (Abstract)
    location    % Returns a vector of Cartesian surface locations
  end

  methods

    function obj = ShapeForce(epsilon)
      % Base class constructor for shape induced force approximation
      %
      % Usage
      %   obj = obj@ott.optics.shapeforce.ShapeForce(epsilon)
      %
      % Parameters
      %   - epsilon (numeric) -- difference in permittivity

      obj.permittivity = epsilon;
    end

    function f = force(obj, varargin)
      % Calculate the force using the shape-surface approximation
      %
      % Implements::
      %
      %     \langle f \rangle = \frac{1}{4} |E|^2 \Delta\epsilon \hat{n}
      %
      % Usage
      %   f = obj.force(...)
      %
      % Optional named arguments
      %   - position (3xN numeric) -- position of the particle
      %   - rotation (3x3N numeric) -- orientation of the particle

      p = inputParser;
      p.addParameter('position', []);
      p.addParameter('rotation', []);
      p.parse(varargin{:});

      % Check if we have much work to do
      if size(p.Results.position, 2) > 1 || size(p.Results.rotation, 2) > 3
        f = ott.utils.prxfun(@obj.force, ...
            'position', p.Results.position, ...
            'rotation', p.Results.rotation);
        return;
      end

      assert(isempty(p.Results.position) ...
          || all(size(p.Results.position) == [3, 1]), ...
          'Position must be empty or 3x1 vector');
      assert(isempty(p.Results.rotation) ...
          || all(size(p.Results.rotation) == [3, 3]), ...
          'Rotation must be empty or 3x3 matrix');

      % Get locations, normals and area elements
      [xyz, normals, dA] = obj.locations();

      % If present, apply rotation
      if ~isempty(p.Results.rotation)
        xyz = p.Results.rotation * xyz;
        normals = p.Results.rotation * normals;
      end

      % If present, apply translation
      if ~isempty(p.Rseults.position)
        xyz = xyz + p.Results.position;
      end

      % Calculate field intensity at locations
      E2 = abs(obj.beam.efield(xyz)).^2;

      % Calculate force
      f = 0.25 .* E2 .* obj.epsilon .* normals .* dA;

    end
  end

  methods
    function obj = set.epsilon(obj, val)
      assert(isnumeric(val) && isscalar(val), ...
        'Epsilon must be numeric scalar');
      obj.epsilon = val;
    end
  end

end

