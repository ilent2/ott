classdef (Abstract) ForceMethod < ott.optics.Common
% Abstract class for optical force calculation methods
%
% Abstract methods
%   - calculateForce -- Called by the force calculation method
%
% Methods
%   - force -- Public function for force calculation

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  methods (Abstract, Hidden)
    % Calculate force acting on a particle
    % This method is called by :meth:`force`.
    % We use a hidden implementation so we can use common documentation.
    %
    % Usage
    %   f = calcualteForce(position, rotation, time)
    %   See :meth:`force` for details about arguments.
    calculateForce;
  end

  methods

    function f = force(obj, varargin)
      % Calculate the force acting on a particle
      %
      % Usage
      %   f = obj.force(...) calculates the force acting on the particle.
      %   The returned force is a 3xN matrix, where N depends on the
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
      f = obj.calculateForce(position, rotation, time);
    end
  end
end

