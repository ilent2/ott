classdef (Abstract) Common
% Common base class for all optical calculation methods

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  methods (Hidden, Static)
    function [position, rotation, time] = parsePositionRotationTime(varargin)
      % Parses the position, rotation and time optional arguments
      % This function is used by :meth:`ForceMethod.force` and similar.
      %
      % Checks that inputs have correct sizes.  Inputs must be scalar
      % or have matching number of elements.
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
      p = inputParser;
      p.addParameter('position', [0;0;0]);
      p.addParameter('rotation', eye(3));
      p.addParameter('time', 0.0);
      p.parse(varargin{:});

      % Check position
      position = p.Results.position;
      assert(ismatrix(position) && size(position, 1) == 3, ...
          'Position must be 3xN matrix');
      Nposition = size(position, 2);

      % Check rotation
      rotation = p.Results.rotation;
      assert(ismatrix(rotation) && size(rotation, 1) == 3 ...
          && mod(size(rotation, 2), 3) == 0, ...
          'Rotation must be 3x3N matrix');
      Nrotation = size(rotation, 2)/3;

      % Check time
      time = p.Results.time;
      assert(ismatrix(time) && size(time, 1) == 1, ...
          'Time must be 1xN vector');
      Ntime = size(time, 2);

      % Check length of optionals
      N = unique([Nposition, Nrotation, Ntime]);
      assert(numel(N) == 1 || (numel(N) == 2 && min(N) == 1), ...
          'Optional arguments should have length 1 or same length');
    end
  end
end
