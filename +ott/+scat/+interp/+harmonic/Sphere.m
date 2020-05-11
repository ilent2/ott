classdef Sphere < ott.scat.utils.ZeroScattered ...
    & ott.shapes.Sphere
% Harmonic model for a spherical particle in harmonic optical trap
% Inherits from :class:`ott.scat.utils.ZeroScattered` and
% :class:`ott.shapes.Spehre`.
%
% The force on the particle is given by::
%
%   F = -k x
%
% where :math:`x` is the position of the particle and :math:`k` is
% the trap stiffness.  The trap stiffness can either be a scalar,
% 3 element vector or 3x3 tensor.
%
% Torque is always zero for this method.  For a harmonic approximation
% with torque, see :class:`OneAxis` or :class:`TwoAxis`.
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

  methods
    function obj = Sphere(stiffness, varargin)
      % Construct a simple harmonic model
      %
      % Usage
      %   obj = Sphere(stiffness, ...)
      %
      % Parameters
      %   - stiffness (numeric scalar|3-vector|3x3-matrix) -- Optical
      %     trap stiffness.
      %
      % Optional named arguments
      %   - position (3x1 numeric) -- Offset for centre of trap.
      %     Default: ``[0;0;0]``.
      %
      %   - radius (numeric) -- Radius of sphere (for visualisation
      %     purposes only).  Default: ``1.0``.
      %
      % Unmatched arguments are passed to Shape.

      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('radius', 1.0);
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      obj = obj@ott.shapes.Sphere(p.Results.radius, unmatched{:});
      obj.stiffness = stiffness;
    end
  end

  methods (Hidden)
    function f = forceInternal(obj, beam, varargin)
      % Calculate the harmonic optical force
      % Ignores particle rotation but should use beam rotation.

      % TODO: Use beam rotation

      if isscalar(obj.stiffness) || isvector(obj.stiffness)
        f = - obj.stiffness .* (beam.position - obj.position);
      else
        f = - obj.stiffness * (beam.position - obj.position);
      end
    end

    function t = torqueInternal(obj, beam, varargin)
      t = [0;0;0];
    end

    function [f, t] = forcetorqueInternal(particle, beam, varargin)
      % Calculate force and torque

      f = particle.forceInternal(beam, varargin{:});
      t = particle.torqueInternal(beam, varargin{:});
    end
  end

  methods % Setters/getters
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

