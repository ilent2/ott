classdef TwoAxis < ott.optics.ForceTorqueMethod
% Harmonic model for a particle with two alignment directions
% Inherits from :class:`ott.optics.ForceTorqueMethod`.
%
% This class extends the simple harmonic model to particle which have two
% alignment directions, this can be used to model cubes and ellipsoids.
% As with the simple Spherical particle harmonic mode, the force is given by::
%
%   F = -k x
%
% where :math:`x` is the position of the particle and :math:`k` is
% the trap stiffness.  The trap stiffness can either be a scalar,
% 3 element vector or 3x3 tensor.
%
% The torque is given by::
%
%   \tau = - (g_1, g_2, 0) (R_w) g_3 - (f_1, f_2, 0) (R_w) f_3
%
% where :math:``g_1, g_2, f_1, f_2`` are vectors perpendicular to the
% alignment directions :math:``g_3, f_3`` and :math:``R_w`` is the particle
% orientation. The magnitude of the vectors determines the trap stiffness.
%
% Properties
%   posStiffness - Translational stiffness (scalar, 3-vector or 3x3 matrix)
%   rotStiffness1 - Rotational stiffness in direction 1 (3x3 matrix)
%   rotStiffness2 - Rotational stiffness in direction 2 (3x3 matrix)
%
% Methods
%   force - Calculate the optical force on the particle
%   torque - Calculate the optical force on the particle
%   forceTorque - Calculate the optical force and torque on the particle

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    posStiffness   % Translational stiffness (scalar, 3-vector or 3x3 matrix)
    rotStiffness1  % Rotational stiffness in direction 1 (3x3 matrix)
    rotStiffness2  % Rotational stiffness in direction 2 (3x3 matrix)
  end

  methods (Hidden)
    function f = calculateForce(obj, position, rotation, time)
      % Calculate the harmonic optical force

      if isscalar(obj.posStiffness) || isvector(obj.posStiffness)
        f = - obj.posStiffness .* position;
      else
        f = - obj.posStiffness * position;
      end
    end

    function t = calculateTorque(obj, position, rotation, time)
      % Calculate alignment torque to a single axes

      % Change rotation to 3x3xN
      rotation = reshape(rotation, 3, 3, []);

      % First axis
      t = [obj.rotStiffness1(:, 1:2), [0;0;0]] * rotation ...
          * obj.rotStiffness1(:, 3);

      % Second axis
      t = t + [obj.rotStiffness2(:, 1:2), [0;0;0]] * rotation ...
          * obj.rotStiffness2(:, 3);

      % Ensure output is 3xN
      t = reshape(t, 3, []);
    end
  end

  methods
    function obj = TwoAxis(varargin)
      % Construct a harmonic model with two alignment axes
      %
      % Usage
      %   obj = TwoAxis(...)
      %
      % Optional named parameters
      %   - posStiffness (numeric 1x1|3x1|3x3) -- Optical trap
      %     translational stiffness.  Default: ``0``.
      %
      %   - rotStiffness1 (numeric 1x1|2x1|3x3) -- Optical trap
      %     rotational stiffness.  Default: ``zeros(3,3)``.
      %
      %   - rotStiffness2 (numeric 1x1|2x1|3x3) -- Optical trap
      %     rotational stiffness.  Default: ``zeros(3,3)``.

      p = inputParser;
      p.addParameter('posStiffness', 0);
      p.addParameter('rotStiffness1', zeros(3,3));
      p.addParameter('rotStiffness2', zeros(3,3));
      p.parse(varargin{:});

      obj.posStiffness = p.Results.posStiffness;
      obj.rotStiffness1 = p.Results.rotStiffness1;
      obj.rotStiffness2 = p.Results.rotStiffness2;
    end
  end

  methods
    function obj = set.posStiffness(obj, val)
      % Check valid values for stiffness

      assert(isnumeric(val), 'Stiffness must be numeric');
      assert(ndims(val) == 2 ...
          && (isscalar(val) || all(size(val) == [3, 1]) ...
          || all(size(val) == [3, 3])), ...
          'Stiffness must be scalar, 3x1 vector or 3x3 matrix');

      obj.posStiffness = val;
    end

    function obj = set.rotStiffness1(obj, val)
      % Check valid values for stiffness

      assert(isnumeric(val), 'Stiffness must be numeric');
      assert(ndims(val) == 2 && all(size(val) == [3,3]), ...
          'Rot. Stiffness must be 3x3 matrix');

      obj.rotStiffness1 = val;
    end

    function obj = set.rotStiffness2(obj, val)
      % Check valid values for stiffness

      assert(isnumeric(val), 'Stiffness must be numeric');
      assert(ndims(val) == 2 && all(size(val) == [3,3]), ...
          'Rot. Stiffness must be 3x3 matrix');

      obj.rotStiffness2 = val;
    end
  end
end

