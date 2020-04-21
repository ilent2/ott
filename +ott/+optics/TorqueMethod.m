classdef (Abstract) TorqueMethod < ott.optics.Common
% Abstract class for optical torque calculation methods
%
% Abstract methods
%   - calculateTorque -- Called by the torque calculation method
%
% Methods
%   - torque -- Public function for torque calculation

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  methods (Abstract, Hidden)
    % Calculate torque acting on a particle
    % This method is called by :meth:`torque`.
    % We use a hidden implementation so we can use common documentation.
    %
    % Usage
    %   f = calcualteTorque(position, rotation, time)
    %   See :meth:`torque` for details about arguments.
    calculateTorque;
  end

  methods

    function f = torque(obj, varargin)
      % Calculate the torque acting on a particle
      %
      % Usage
      %   f = obj.torque(...) calculates the torque acting on the particle.
      %   The returned torque is a 3xN matrix, where N depends on the
      %   optional arguments.
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
      f = obj.calculateTorque(position, rotation, time);
    end
  end
end

