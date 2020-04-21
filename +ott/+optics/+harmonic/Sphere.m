classdef Sphere < ott.optics.ForceMethod
% Harmonic model for a spherical particle in harmonic optical trap
% Inherits from :class:`ott.optics.ForceMethod`.
%
% The force on the particle is given by::
%
%   F = -k x
%
% where :math:`x` is the position of the particle and :math:`k` is
% the trap stiffness.  The trap stiffness can either be a scalar,
% 3 element vector or 3x3 tensor.
%
% Properties
%   stiffness - Trap stiffness (scalar | 3-vector | 3x3 matrix)
%
% Methods
%   force - Calculate the optical force on the particle

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    stiffness  % Trap stiffness (scalar | 3-vector | 3x3 matrix)
  end

  methods (Hidden)
    function f = calculateForce(obj, position, rotation, time)
      % Calculate the harmonic optical force

      if isscalar(obj.stiffness) || isvector(obj.stiffness)
        f = - obj.stiffness .* position;
      else
        f = - obj.stiffness * position;
      end
    end
  end

  methods
    function obj = Sphere(stiffness)
      % Construct a simple harmonic model
      %
      % Usage
      %   obj = Sphere(stiffness)
      %
      % Parameters
      %   - stiffness (numeric scalar|3-vector|3x3-matrix) -- Optical
      %     trap stiffness.

      obj.stiffness = stiffness;
    end
  end

  methods
    function obj = set.stiffness(obj, val)
      % Check valid values for stiffness

      assert(isnumeric(val), 'Stiffness must be numeric');
      assert(ndims(val) == 2 ...
          && (isscalar(val) || all(size(val) == [3, 1]) ...
          || all(size(val) == [3, 3])), ...
          'Stiffness must be scalar, 3x1 vector or 3x3 matrix');

      obj.stiffness = val;
    end
  end
end

