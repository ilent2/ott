classdef (Abstract) ForceTorqueMethod < ...
    ott.optics.ForceMethod & ott.optics.TorqueMethod
% Abstract class for systems with force and torque calculation methods
%
% Abstract methods
%   - calculateForceTorque -- Called by the calculation method
%
% Methods
%   - forceTorque -- Public function for force and torque calculation

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  methods (Hidden)
    function f = calculateForceTorque(obj, position, rotation, time)
      % Calculate force and torque acting on a particle
      % This method is called by :meth:`forceTorque`.
      % We use a hidden implementation so we can use common documentation.
      %
      % The default implementation calls the force and torque methods.
      % However, in some instances it may be more efficient to calculate
      % force and torque simultaneously, in that case, this method should
      % be overwritten.
      %
      % Usage
      %   f = calcualteForceTorque(position, rotation, time)
      %   See :meth:`forceTorque` for details about arguments.

      % Calculate force and torque
      force = obj.force(position, rotation, time);
      torque = obj.torque(position, rotation, time);

      % Replicate elements if required
      if size(force, 2) ~= size(torque, 2)
        if size(force, 2) < size(torque, 2)
          force = repelem(force, [1, size(torque, 2)]);
        else
          torque = repelem(torque, [1, size(force, 2)]);
        end
      end

      % Combine for output
      f = [force; torque];
    end
  end

  methods

    function varargout = forceTorque(obj, varargin)
      % Calculate the force and torque acting on a particle
      %
      % Usage
      %   f = obj.forceTorque(...) calculates the force and torque acting
      %   on the particle.  The returned torque is a 6xN matrix,
      %   where N depends on the optional arguments.
      %
      %   [f, t] = obj.forceTorque(...) splits the force and torque
      %   into two 3xN matricies ``f`` and ``t`` respectively.
      %
      % Optional named arguments
      %   - position (numeric 3xN) - particle position.
      %     Default: ``[0; 0; 0]``.
      %
      %   - rotation (numeric 3x3N) - particle orientation.
      %     Default: ``eye(3)``.
      %
      %   - time (numeric 1xN) - simulation time.
      %     Default: ``0.0``.

      % Parse inputs
      [position, rotation, time] = ...
          ott.optics.Common.parsePositionRotationTime(varargin{:});

      % Call the force calculation method to do the work
      f = obj.calculateForceTorque(position, rotation, time);

      if nargout == 1
        varargout{1} = f;
      elseif nargout == 2
        varargout{1} = f(1:3, :);
        varargout{2} = f(4:6, :);
      else
        warning('Too many output arguments requested');
      end
    end
  end
end

