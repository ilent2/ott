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
    range(1,2) double           % Range of optical trap
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
      % Find traps in an array of position/force values.
      %
      % It may be useful to pre-filter the arrays, such as applying a
      % low pass filter to the forces to remove high frequency noise or
      % rounding small values to zero.
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

      % Find rough equilibrium positions

      % Find first equilibrium position
      if force(1) > 0
        last = find(force < 0, 1);
      elseif force(1) == 0
        last = find(force ~= 0, 1);
      else
        last = find(force > 0, 1);
      end

      % Find reaming equilibriums
      eqs = [];
      while ~isempty(last)
        eqs(end+1) = last; %#ok<AGROW>

        if force(last) > 0
          last = find(force(last:end) < 0, 1) + last - 1;
        else
          last = find(force(last:end) > 0, 1) + last - 1;
        end
      end

      % Find more precise equilibrium positions and trap stiffness
      % This is based on ott.find_equilibrium

      peqs = zeros(1, length(eqs));
      pstiff = zeros(1, length(eqs));

      for ii = 1:length(eqs)

        % For extended zeros we might have several points equal to zero.
        % For the polyfit, we want to use a point at each end of the equilibrium
        % and all points in between.
        first = find(force(1:eqs(ii)-1)*sign(force(eqs(ii))) < 0, 1, 'last');
        if isempty(first)
          first = 1;
        end

        % Equilibrium is between eqs(ii)-1 and eqs(ii)
        eqguess = (position(eqs(ii)) + position(first))/2;

        % Fit polynomial to points aroung equilibrium
        eqrange = max([first-1,1]):min([eqs(ii)+1,length(position)]);
        z = position(eqrange);

        % Scale position and force before polyfit
        zmin = min(z);
        zmax = max(z);
        z = 2 * (z - zmin) / (zmax - zmin) - 1;
        zzero = 2 * (eqguess - zmin) / (zmax - zmin) - 1;

        % Find equilibrium: fit local points to 3rd order polynomail
        % Requires small distance between positions
        if length(z) < 4
          pz=polyfit(z, force(eqrange), 2);
          dpz=[2*pz(1),1*pz(2)];
        else
          pz=polyfit(z, force(eqrange), 3);
          dpz=[3*pz(1),2*pz(2),1*pz(3)];
        end
        root_z=roots(pz);

        % Ignore non-real roots
        real_z=root_z(imag(root_z)==0);
        if numel(real_z) == 0
          error('No real roots');
        end

        % Keep only one root closest to position(eqs(ii))
        zeqs_idx = abs(real_z-zzero) == min(abs(real_z-zzero));
        real_z = real_z(zeqs_idx);
        real_z = real_z(1);

        % Get the equilibrium
        peqs(ii) = real_z;

        % Calculate stiffness (using derivative of 3rd order polynomial)
        pstiff(ii) = polyval(dpz, peqs(ii));

        % Inverse scaling of position and force after fitting
        peqs(ii) = (peqs(ii) + 1)/2*(zmax - zmin) + zmin;
        pstiff(ii) = pstiff(ii)*2/(zmax - zmin);

      end

      % Calculate other properties needed for traps

      traps = ott.tools.FindTraps1d.empty(1, 0);

      for ii = 1:length(eqs)

        % Calculate trap depth (part 1)
        frange = [1, length(force)];
        if ii ~= 1
          frange(1) = eqs(ii-1);
        end
        if ii ~= length(eqs)
          frange(2) = eqs(ii+1)-1;
        end
        [fn, fnidx] = min(force(frange(1):frange(2)));
        [fx, fxidx] = max(force(frange(1):frange(2)));

        idx = length(traps) + 1;

        % Store stiffness and equilibrium
        traps(idx).position = peqs(ii);
        traps(idx).stiffness = pstiff(ii);

        % Calculate trap range
        traps(idx).range = [-Inf, Inf];
        if ii ~= 1
          traps(idx).range(1) = peqs(ii-1);
        end
        if ii ~= length(eqs)
          traps(idx).range(2) = peqs(ii+1);
        end

        % Calculate trap depth (part 2)
        if fnidx < fxidx
          [fx, fn] = deal(fn, fx);
          [fxidx, fnidx] = deal(fnidx, fxidx);
        end
        traps(idx).minforce = fx;
        traps(idx).maxforce = fn;
        traps(idx).minposition = position(frange(1)+fxidx-1);
        traps(idx).maxposition = position(frange(1)+fnidx-1);
        traps(idx).depth = min(abs(...
            [traps(idx).maxforce, traps(idx).minforce]));
      end

      % Remove unstable traps
      if ~p.Results.keep_unstable
        traps = traps.filterUnstable();
      end
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
      %     beams with nm particles in SI units).  This is also used
      %     for calculating the ``stiffness`` property.

      p = inputParser;
      p.addParameter('ForceTol', 1e-10);
      p.addParameter('StepTol', 1e-13);
      p.parse(varargin{:});

      stepTol = p.Results.StepTol;
      forceTol = p.Results.ForceTol;

      % Check inputs
      assert(isa(force, 'function_handle'), ...
          'Force must be a function handle');
      assert(isnumeric(position) && numel(position) == 2, ...
          'Position must be 2 element numeric vector');
      assert(isnumeric(stepTol) && isscalar(stepTol) && stepTol > 0, ...
          'StepTol must be positive numeric scalar');
      assert(isnumeric(forceTol) && isscalar(forceTol) && forceTol > 0, ...
          'ForceTol must be positive numeric scalar');

      % Evaluate initial forces
      f1 = force(position(1));
      f2 = force(position(2));

      % Check forces
      assert(isnumeric(f1) && isscalar(f1), ...
          'Value returned by force function should be numeric scalar');

      % Check sign on initial points
      if sign(f1) == sign(f2)
        warning('ott:tools:FindTraps1d:Bisection:bounds_sign', ...
            'Initial forces have same sign, may be no stable equilibrium');
      end

      while diff(position) > stepTol && ~any(abs([f1, f2]) < forceTol)

        % Calculate new point
        midp = mean(position);
        midf = force(midp);

        if sign(midf) == sign(f1)
          f1 = midf;
          position(1) = midp;
        else
          f2 = midf;
          position(2) = midp;
        end
      end

      % Find best trap
      [f1, I] = min(abs([f1, f2]));
      traps = ott.tools.FindTraps1d();
      traps.position = position(I);
      traps.stiffness = (f1 - force(position(I) - stepTol))./stepTol;
      traps.depth = nan;
      traps.range = [nan, nan];
      traps.minforce = nan;
      traps.maxforce = nan;
      traps.minposition = nan;
      traps.maxposition = nan;
    end

    function traps = Newton(position, force, varargin)
      % Find trap starting with a initial guess at the location.
      %
      % Uses Newtons method and moves the particle in the direction
      % given by the optical force at each step.  Searches for a stable
      % trap, if you want to find an unstable trap, modify the force
      % function to give negative forces.
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

      p = inputParser;
      p.addParameter('maxIterations', 100);
      p.addParameter('smallStep', 1e-10);
      p.addParameter('ForceTol', 1e-13);
      p.parse(varargin{:});

      maxIter = p.Results.maxIterations;
      smallStep = p.Results.smallStep;
      forceTol = p.Results.ForceTol;

      % Check inputs
      assert(isnumeric(position) && isscalar(position), ...
          'position must be numeric scalar');
      assert(isa(force, 'function_handle'), ...
          'force must be a function handle');
      assert(isnumeric(maxIter) && isscalar(maxIter) && maxIter > 0, ...
          'maxIter must be positive numeric scalar');
      assert(isnumeric(smallStep) && isscalar(smallStep), ...
          'smallStep must be numeric scalar');
      assert(isnumeric(forceTol) && isscalar(forceTol) && forceTol > 0, ...
          'forceTol must be positive numeric scalar');

      % Get initial force
      fv = force(position);
      dvf = (force(position + smallStep) - fv) ./ smallStep;

      % Check force
      assert(isnumeric(fv) && isscalar(fv), ...
          'force must return a numeric scalar');

      % Newton's iteration
      for ii = 1:maxIter
        
        % Calculate step
        position = position - fv ./ dvf;
        
        % Calculate new values
        fv = force(position);
        dvf = (force(position + smallStep) - fv) ./ smallStep;

        % Check for stopping condition
        if abs(fv) <= forceTol
          break;
        end
      end

      % Final check of force value
      if abs(fv) > forceTol
        warning('ott:tools:FindTraps1d:Newton:final_force_check', ...
            'Final force not within force tolerance');
      end

      % Setup output
      traps = ott.tools.FindTraps1d();
      traps.position = position;
      traps.stiffness = dvf;
      traps.depth = nan;
      traps.range = [nan, nan];
      traps.minforce = nan;
      traps.maxforce = nan;
      traps.minposition = nan;
      traps.maxposition = nan;
    end
  end

  methods
    function h = plot(traps, varargin)
      % Generate plot of the traps in the array.
      %
      % Squares mark equilibriums.  Solid lines (with crosses) mark
      % stable traps.  Dashed lines (with circles) mark unstable traps.
      %
      % Usage
      %   h = traps.plot(...)
      %
      % Returns
      %   - h -- An array of plot handles.  The first handle is the
      %     eqilibrium markers.  Remaining handles are trap lines.

      p = inputParser;
      p.parse(varargin{:});

      % Plot trap locations
      h = plot([traps.position], zeros(size(traps)), 's');

      % Setup hold
      oldHold = ishold();
      hold on;

      % Plot trap ranges
      lineType = [traps.globalStiffness] < 0;
      minPositions = [traps.minposition];
      maxPositions = [traps.maxposition];
      h = [h; plot([minPositions(lineType); maxPositions(lineType)], ...
          zeros(2, sum(lineType)), 'kx-')];
      h = [h; plot([minPositions(~lineType); maxPositions(~lineType)], ...
          zeros(2, sum(~lineType)), 'ro--')];

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

      p = inputParser;
      p.parse(varargin{:});

      % Identify trap groups
      %
      % We do two passes to identify the maximum and minimum trapping
      % forces, then we combine the indices to identify unique traps

      idxer = 0;
      group_idx_max = zeros(1, length(traps));
      omaxforce = 0;
      for ii = 1:length(traps)
        if traps(ii).minforce > omaxforce
          omaxforce = traps(ii).minforce;
          idxer = idxer + 1;
        end
        group_idx_max(ii) = idxer;
      end

      idxer = 0;
      ominforce = 0;
      group_idx_min = zeros(1, length(traps));
      for ii = length(traps):-1:1
        if traps(ii).maxforce < ominforce
          ominforce = traps(ii).maxforce;
          idxer = idxer + 1;
        end
        group_idx_min(ii) = idxer;
      end

      group_idx = group_idx_max + (max(group_idx_min) - group_idx_min);

      % Create list of trap groups
      unique_group_idx = unique(group_idx);
      grouped = cell(1, numel(unique_group_idx));
      stats(numel(grouped)) = ott.tools.FindTraps1d;
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
      %
      % Optional named parameters
      %   - AbsTol (numeric | []) -- Absolute tolerance for removing
      %     shallow traps.  Default: ``[]``.
      %
      %   - RelTol (numeic | []) -- Relative tolerance for removing shallow
      %     traps.  Relative to maximum trap depth.  Default: ``1e-2``.

      p = inputParser;
      p.addParameter('AbsTol', []);
      p.addParameter('RelTol', []);
      p.parse(varargin{:});

      notShallow = true(size(traps));

      if ~isempty(p.Results.AbsTol)
        notShallow = notShallow & [traps.depth] >= p.Results.AbsTol;
      end

      if ~isempty(p.Results.RelTol)
        maxDepth = max([traps.depth]);
        notShallow = notShallow & [traps.depth] >= p.Results.RelTol*maxDepth;
      end

      traps = traps(notShallow);
    end

    function traps = filterUnstable(traps, varargin)
      % Filter traps that are unstable
      %
      % Usage
      %   traps = traps.fitlerUnstable(...)

      p = inputParser;
      p.parse(varargin{:});

      stable = [traps.stiffness] < 0 | [traps.globalStiffness] < 0;
      traps = traps(stable);
    end
  end

  methods % Getters/setters
    function val = get.globalStiffness(traps)
      val = (traps.maxforce - traps.minforce) ...
          ./ (traps.maxposition - traps.minposition);
    end
  end
end

