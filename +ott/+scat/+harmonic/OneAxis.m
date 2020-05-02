classdef OneAxis < ott.optics.ForceTorqueMethod
% Harmonic model for a particle with one alignment direction.
% Inherits from :class:`ott.optics.ForceTorqueMethod`.
%
% This class extends the simple harmonic model to particle with a single
% alignment direction, which could be used to model elongated particles.
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
%   \tau = - (g_1, g_2, 0) (R_w) g_3
%
% where :math:``g_1, g_2`` are two vectors perpendicular to the alignment
% direction :math:``g_3`` and :math:``R_w`` is the particle orientation.
% The magnitude of the vectors determines the trap stiffness.
%
% Properties
%   posStiffness - Translational stiffness (scalar, 3-vector or 3x3 matrix)
%   rotStiffness - Rotational stiffness (scalar, 2-vector or 3x3 matrix)
%
% Methods
%   force - Calculate the optical force on the particle
%   torque - Calculate the optical force on the particle
%   forceTorque - Calculate the optical force and torque on the particle

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    posStiffness  % Translational stiffness (scalar, 3-vector or 3x3 matrix)
    rotStiffness  % Rotational stiffness (scalar, 2-vector or 3x3 matrix)
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

      if isscalar(obj.rotStiffness)
        t = diag([1, 1, 0]*obj.rotStiffness) * rotation * [0;0;1];
      elseif isvector(obj.rotStiffness)
        t = diag([obj.rotStiffness(:); 0]) * rotation * [0;0;1];
      else
        t = [obj.rotStiffness(:, 1:2), [0;0;0]] * rotation ...
            * obj.rotStiffness(:, 3);
      end

      % Ensure output is 3xN
      t = reshape(t, 3, []);
    end
  end

  methods
    function obj = OneAxis(varargin)
      % Construct a harmonic model with a alignment axis
      %
      % Usage
      %   obj = Sphere(...)
      %
      % Optional named parameters
      %   - posStiffness (numeric 1x1|3x1|3x3) -- Optical trap
      %     translational stiffness.  Default: ``0``.
      %
      %   - rotStiffness (numeric 1x1|2x1|3x3) -- Optical trap
      %     rotational stiffness.  Default: ``0``.

      p = inputParser;
      p.addParameter('posStiffness', 0);
      p.addParameter('rotStiffness', 0);
      p.parse(varargin{:});

      obj.posStiffness = p.Results.posStiffness;
      obj.rotStiffness = p.Results.rotStiffness;
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

    function obj = set.rotStiffness(obj, val)
      % Check valid values for stiffness

      assert(isnumeric(val), 'Stiffness must be numeric');
      assert(ndims(val) == 2 ...
          && (isscalar(val) || all(size(val) == [2, 1]) ...
          || all(size(val) == [3, 3])), ...
          'Stiffness must be scalar, 2x1 vector or 3x3 matrix');

      obj.rotStiffness = val;
    end
  end
end

