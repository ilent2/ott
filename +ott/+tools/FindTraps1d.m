classdef FindTraps1d
% Descriptions and method for finding traps in 1-dimension.
%
% This class contains static methods for finding traps in one dimension.
% Instances of the class simply describe properties of traps, such as
% trap location, depth, stiffness.
%
% Units of properties depend on the method that created the FindTraps1d
% instance (see documentation for static methods).
%
% Properties
%   - position        -- Trap position
%   - stiffness       -- Trap stiffness at equilibrium (force/position units)
%   - depth           -- Optical trap depth (force units)
%   - range           -- Range of optical trap (position units)
%   - minforce        -- Minimum force in trap (force units)
%   - maxforce        -- Maximum force in trap (force units)
%   - minposition     -- Position of minimum force (position units)
%   - maxposition     -- Position of maximum force (position units)
%   - globalStiffness -- Stiffness between min/max force locations
%
% Methods
%   - groupStable     -- Groups stable traps based on trap depth
%   - plot            -- Generate a plot of the traps in the array
%
% Static methods
%   - FromArray       -- Find traps from an array of force/position data
%   - Bisection       -- Uses bisection method to find trap position
%   - Newton          -- uses Newton's method to find trap position

% Copyright 2020 Isaac Lenton (aka ilent2)
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    position(1,1) double        % Trap position
    stiffness(1,1) double       % Trap stiffness at equilibrium
    depth(1,1) double           % Optical trap depth
    range(2,1) double           % Range of optical trap
    minforce(1,1) double        % Minimum force in trap
    maxforce(1,1) double        % Maximum force in trap
    minposition(1,1) double     % Position of minimum force
    maxposition(1,1) double     % Position of maximum force
  end

  properties (Dependent)
    globalStiffness             % Stiffness between min/max force locations
  end

  methods (Static)
    function traps = FromArray(position, force, varargin)
      % Find traps in an array of position/force values
      %
      % Usage
      %   traps = ott.tools.FindTraps1d.FromArray(position, force, ...)
      %
      % Parameters
      %   - position (N numeric) -- Location of each force sample.
      %
      %   - force (N numeric) -- Optical force at each position.
      %
      % Optional named parameters
      %   - keep_unstable (logical) -- If true, keeps both stable and
      %     unstable equilibriums.  Default: ``false``.

      p = inputParser;
      p.addParameter('keep_unstable', false, @islogical);
      p.parse(varargin{:});

      % Check inputs
      assert(isvector(position) && isnumeric(position), ...
        'position must be a numeric vector');
      assert(isvector(force) && isnumeric(force), ...
        'force must be a numeric vector');
      assert(numel(position) == numel(force), ...
        'number of positions must match number of forces');

      % Make sure shapes match
      position = position(:);
      force = force(:);

      % TODO
    end

    function traps = Bisection(position, force, varargin)
      % Find traps using bisection method and a force calculation function
      %
      % The default stopping criteria assume the force/position are in
      % SI units with typical optical tweezers beams/particles.  You may
      % need to adjust these parameters.
      %
      % Usage
      %   traps = ott.tools.FindTraps1d.Bisection(position, force, ...)
      %
      % Parameters
      %   - position (2 numeric) -- Bounds for trap location.  When the
      %     force function is evaluated, one point should give a positive
      %     force and the other a negative force.
      %
      %   - force (function_handle) -- Function to calculate force.
      %     Signature: ``@(position) scalar``.
      %
      % Optional named parameters
      %   - ForceTol (numeric) -- Absolute tolerance for stopping.
      %     Default: ``1e-13`` (should work for 1mW beams in SI units).
      %
      %   - StepTol (numeric) -- Minimum step size tolerance for stopping.
      %     Default: ``1e-10`` (should work for 1um wavelength
      %     beams with nm particles in SI units).

      % TODO
    end

    function traps = Newton(position, force, varargin)
      % Find trap starting with a initial guess at the location.
      %
      % Uses Newtons method and moves the particle in the direction
      % given by the optical force at each step.  The trap should be
      % a stable equilibrium.
      %
      % Usage
      %   traps = ott.tools.FindTraps1d.Newton(position, force, ...)
      %
      % Parameters
      %   - position (numeric) -- Initial guess at position.
      %
      %   - force (function_handle) -- Function for calculating force.
      %     Signature: ``@(position) scalar``.
      %
      % Optional named parameters
      %   - maxIterations (numeric) -- Maximum number of iterations
      %     before stopping.  If solution not found, raises a warning.
      %     Default: ``100``.
      %
      %   - smallStep (numeric) -- Small step used for derivative
      %     calculation.  Default: ``1e-10``.
      %
      %   - ForceTol (numeric) -- Absolute tolerance for stopping.
      %     Default: ``1e-13`` (should work for 1mW beams in SI units).

      % TODO
    end
  end

  methods
    function varargout = plot(traps, varargin)
      % Generate plot of the traps in the array
      %
      % Usage
      %   [...] = traps.plot(...)
      %
      % Additional parameters/outputs passed to/from ``plot``.

      % Plot trap locations
      plot([traps.position], zeros(size(traps)), 'x');

      % Setup hold
      oldHold = ishold();
      hold on;

      % Plot trap ranges
      lineType = [traps.globalStiffness] < 0;
      minForces = [traps.minforce];
      maxForces = [traps.maxforce];
      plot([minForces(lineType), maxForces(lineType)], ...
          zeros(sum(lineType), 2), 'kx-');
      plot([minForces(~lineType), maxForces(~lineType)], ...
          zeros(sum(~lineType), 2), 'ro-');

      % Return hold to initial state
      if ~oldHold
        hold off;
      end
    end

    function [grouped, stats] = groupStable(traps, varargin)
      % Groups stable traps separated by smaller unstable traps.
      %
      % This is useful for calculating the total trap depth of an
      % optical trap which contains local unstable equilibria.
      %
      % Usage
      %   [grouped, stats] = traps.groupStable(...)
      %
      % Returns
      %   - grouped (cell) -- Cell array of grouped traps.
      %   - stats (FindTraps1d) -- Meta traps (with no position or
      %     stiffness) describing each group.

      % Identify trap groups
      %
      % We do two passes to identify the maximum and minimum trapping
      % forces, then we combine the indices to identify unique traps

      idxer = 0;
      group_idx_max = zeros(1, length(traps));
      maxforce = 0;
      for ii = 1:length(traps)
        if traps(ii).minmax_force(1) > maxforce
          maxforce = traps(ii).minmax_force(1);
          idxer = idxer + 1;
        end
        group_idx_max(ii) = idxer;
      end

      idxer = 0;
      minforce = 0;
      group_idx_min = zeros(1, length(traps));
      for ii = length(traps):-1:1
        if traps(ii).minmax_force(2) < minforce
          minforce = traps(ii).minmax_force(2);
          idxer = idxer + 1;
        end
        group_idx_min(ii) = idxer;
      end

      group_idx = group_idx_max + (max(group_idx_min) - group_idx_min);

      % Create list of trap groups
      unique_group_idx = unique(group_idx);
      grouped = cell(1, numel(unique_group_idx));
      stats(numel(grouped)) = ott.tools.FieldTraps1d;
      for ii = 1:numel(unique_group_idx)
        grouped{ii} = traps(group_idx == unique_group_idx(ii));

        stats(ii).position = nan;
        stats(ii).stiffness = nan;
        stats(ii).range(1) = min([grouped{ii}.range]);
        stats(ii).range(2) = max([grouped{ii}.range]);
        stats(ii).minposition = grouped{ii}(1).minposition;
        stats(ii).maxposition = grouped{ii}(end).maxposition;
        stats(ii).minforce = grouped{ii}(1).minforce;
        stats(ii).maxforce = grouped{ii}(end).maxforce;
        stats(ii).depth = min(abs(...
            [stats(ii).minforce, stats(ii).maxforce]));
      end
    end

    function traps = filterShallow(traps, varargin)
      % Filter traps that are shallow
      %
      % Usage
      %   traps = traps.filterShallow(...)

      % TODO
    end

    function traps = filterUnstable(traps, varargin)
      % Filter traps that are unstable
      %
      % Usage
      %   traps = traps.fitlerUnstable(...)

      % TODO
    end
  end

  methods % Getters/setters
    function val = get.globalStiffness(traps)
      val = (traps.maxforce - traps.minforce) ...
          ./ (traps.maxposition - traps.minposition);
    end
  end
end

