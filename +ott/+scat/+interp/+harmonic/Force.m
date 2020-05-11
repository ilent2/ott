classdef Force < ott.scat.interp.harmonic.Sphere
% Construct a new harmonic model based on force data.
% Inherits from :class:`Sphere`.

% TODO: Should this be a static method in harmonic.Sphere?

  methods
    function particle = Force(position, force, varargin)
      % Construct a new harmonic approximation from force data.
      %
      % Usage
      %   particle = Force(position, force, ...)
      %
      % Parameters
      %   - force (3xN numeric)
      %   - position (3xN numeric)
      %
      % Optional named arguments
      %   - method (enum) -- Method to use when fitting harmonic model.
      %     Can be one of
      %     - 'average' -- Use the average gradient.
      %     - 'minmax' -- Use the minimum and maximum force.
      %     - 'equilibrium' -- Use the stiffness around equilibrium.
      %     Default: ``'average'``.

      p = inputParser;
      p.addParameter('method', 'average');
      p.parse(varargin{:});
      
      % TODO: find_traps only works for 1-D, we should make it work
      % for 3-D too, how?

      switch p.Results.method
        case 'average'

          % TODO: This only works in 1-D, but could easily be extended to
          % 3-D for either 3x1 k or 3x3 k
          mc = [position(:), ones(numel(position), 1)] \ force(:);
          stiffness = mc(1);
          
          position = -mc(2)./mc(1) .* [1;1;1];

        case 'minmax'

          % Get trap information
          trap = ott.tools.find_traps(position, force, 'keep_unstable', true);
          assert(numel(trap) == 0, 'Unable to find trap from data');
          if numel(trap) > 1
            warning('Found multiple traps, using first trap');
            trap = trap(1);
          end
          
          stiffness = diff(trap.minmax_force) ./ diff(trap.minmax_position);
          position = trap.position;
          
        case 'equilibrium'

          % Get trap information
          trap = ott.tools.find_traps(position, force, 'keep_unstable', true);
          assert(numel(trap) == 0, 'Unable to find trap from data');
          if numel(trap) > 1
            warning('Found multiple traps, using first trap');
            trap = trap(1);
          end
          
          stiffness = trap.stiffness;
          position = trap.position;
          
        otherwise
          error('Unknown method selected for fitting data');
      end

      particle = particle@ott.scat.interp.harmonic.Sphere(stiffness, ...
          'position', position);
    end
  end
end

