classdef Dda
% Minimal implementation of the discrete dipole approximation.
%
% This class calculates the dipole polarizations from the interaction
% matrix `A` and the incident electric field `E`::
%
%   p = A \ E
%
% The class implements methods for calculating the interaction matrix
% with optimisation for mirror, rotational symmetry and low memory.
%
% Properties
%   - locations       -- Location of each dipole
%   - polarizability  -- Dipole polarizabilities
%   - xySymmetry      -- True if using z-mirror symmetry
%   - zRotSymmetry    -- Order of z-rotational symmetry (0 - infinite)
%   - ndipoles        -- Number of dipoles in the array
%
% Methods
%   - solve           -- Solve for dipole polarizations
%   - interaction_matrix -- Calculate the interaction matrix
%
% Static methods
%   - FromShape               -- Construct instance from geometric shape

  properties (SetAccess=protected)
    locations        % Location of each dipole
    polarizability   % Dipole polarizabilities
  end

  properties
    xySymmetry       % True if using z-mirror symmetry
    zRotSymmetry     % Order of z-rotational symmetry (0 - infinite)
  end

  properties (Dependent)
    ndipoles         % Number of dipoles in the array
  end

  methods (Static)
    function dda = FromShape(shape, varargin)
      % Construct a DDA instance from a geometric shape.
      %
      % Usage
      %   dda = ott.tmatrix.dda.Dda.FromShape(shape, ...)
      %
      % Optional named arguments
      %   - spacing -- (numeric) -- Dipole spacing in wavelength units.
      %     Default: ``1/20``.
      %
      %   - polarizability -- (function_handle | 1x1 | 3x3 numeric)
      %     Particle polarizability.  Must be a function handle for
      %     inhomogeneous materials or either a scalar or 3x3 matrix for
      %     homogeneous materials.
      %     Default: ``@(xyz, spacing, ri) polarizability.LDR(spacing, ri)``
      %
      %   - index_relative -- (function_handle | numeric) Method to calculate
      %     relative refractive index or homogeneous value.  Ignored if
      %     polarizability is not a function handle.
      %
      % Unmatched parameters are passed to :meth:`FromPolarizability`.

      p = inputParser;
      p.addParameter('spacing', 1/20, @isnumeric);
      p.addParameter('polarizability', ...
          @(~, s, r) ott.tmatrix.dda.polarizability.LDR(s, r));
      p.addParameter('index_relative', 1.0);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      spacing = p.Results.spacing;
      assert(isscalar(spacing) && spacing > 0 && isreal(spacing), ...
          'spacing must be positive numeric real value');

      % Should we avoid 0 for voxels?
      use_even = mod(shape.zRotSymmetry, 2) == 0 || shape.xySymmetry;

      % Calculate voxel locations
      voxels = shape.voxels(spacing, 'even', use_even);

      % Calculate polarizability
      polarizability = p.Results.polarizability;
      if isnumeric(polarizability)
        assert(ismatrix(polarizability) && (isscalar(polarizability) ...
            || all(size(polarizability) == [3, 3])), ...
            'polarizability must be function handle or 3x3 numeric matrix');
      else

        ri = p.Results.index_relative;
        if isnumeric(ri)
          assert(isnumeric(ri) && isscalar(ri), ...
              'relative index must be function handle or numeric scalar');
        else
          ri = ri(voxels);
        end

        % Calculate polarizability
        polarizability = polarizability(voxels, spacing, ri);
      end

      % Reshape polarizability
      if numel(polarizability) == size(voxels, 2)
        polarizability = eye(3) .* reshape(polarizability, 1, 1, []);
        polarizability = reshape(polarizability, 3, []);
      elseif isscalar(polarizability)
        polarizability = repmat(eye(3)*polarizability, 1, size(voxels, 2));
      elseif size(polarizability, 2) == 3
        polarizability = repmat(polarizability, 1, size(voxels, 2));
      else
        error('polarizability must be scalar, 3x3 or 3x3N matrix');
      end

      % Hand over to FromPolarzability
      dda = ott.tmatrix.dda.Dda(voxels, polarizability, ...
          'xySymmetry', shape.xySymmetry, ...
          'zRotSymmetry', shape.zRotSymmetry, ...
          unmatched{:});
    end
  end

  methods
    function dda = Dda(varargin)
      % Construct new DDA solver instance.
      %
      % Usage
      %   dda = Dda(locations, polarizabilities, ...)
      %
      % Parameters
      %   - voxels -- (3xN numeric) Voxel locations in Cartesian coordinates.
      %
      %   - polarizabilities -- (3x3N numeric) Array of 3x3 dipole
      %     polarizability tensors.
      %
      % Optional named arguments
      %   - xySymmetry (logical) -- If the particle has
      %     z-mirror symmetry.  Default: ``false``.
      %     If true, locations with z<0 are removed.
      %
      %   - zRotSymmetry (numeric) -- Order of the particle
      %     z-rotational symmetry.  Default: ``1``.
      %     If ``0``, uses fourth order rotational symmetry (might change
      %     in a future release).
      %     If zRotSymmetry is 2, removes x<0 points.  If zRotSymmetry is
      %     mod 4, removes x<0 | y<0 points.  No other options supported
      %     for now.

      p = inputParser;
      p.addOptional('locations', [], @isnumeric);
      p.addOptional('polarizability', [], @isnumeric);
      p.addParameter('xySymmetry', false);
      p.addParameter('zRotSymmetry', 1);
      p.parse(varargin{:});

      xyz = p.Results.locations;
      alpha = p.Results.polarizability;

      assert(ismatrix(xyz) && size(xyz, 1) == 3, ...
          'voxels must be 3xN numeric matrix');
      assert(ismatrix(alpha) && size(alpha, 1) == 3 ...
          && mod(size(alpha, 2), 3) == 0, ...
          'polarizability must be 3x3N numeric matrix');
      assert(size(xyz, 2)*3 == size(alpha, 2), ...
          'number of polarizabilities must match number of voxels');

      % Filter xySymmetry points
      if p.Results.xySymmetry
        rmidx = xyz(3, :) < 0;
        alpha(:, repelem(rmidx, 1, 3)) = [];
        xyz(:, rmidx) = [];
      end

      % Filter zRotSymmetry
      if p.Results.zRotSymmetry ~= 1
        if mod(p.Results.zRotSymmetry, 4) == 0
          rmidx = xyz(1, :) < 0 | xyz(2, :) < 0;
          alpha(:, repelem(rmidx, 1, 3)) = [];
          xyz(:, rmidx) = [];
        elseif p.Results.zRotSymmetry == 2
          rmidx = xyz(1, :) < 0;
          alpha(:, repelem(rmidx, 1, 3)) = [];
          xyz(:, rmidx) = [];
        else
          error('zRotSymemtry must be either mod 4 or 2 for now');
        end
      end

      dda.locations = xyz;
      dda.polarizability = alpha;
      dda.xySymmetry = p.Results.xySymmetry;
      dda.zRotSymmetry = p.Results.zRotSymmetry;
    end

    function dipoles = solve(dda, Eincident, varargin)
      % Solve for dipole polarizations
      %
      % Usage
      %   dipoles = dda.solve(Eincident, ...)
      %
      % Parameters
      %   - Eincident -- (3NxM numeric) Incident electric field at each
      %     dipole.  Format: [x1;y1;z1;x2;y2;z2;...].
      %
      % Optional named arguments
      %   - solver -- (function_handle) Solver to use.  Good solvers to try
      %     include ``gmres``, ``bicgstab`` and ``\``.
      %     Default: ``@(A, E) A \ E``.
      %
      %   - multiBeam -- (logical) If true, passes the whole beam into the
      %     solver, otherwise iterates over each beam.  Default: ``true``.
      %
      %   - parity (enum) -- Parity of incident beam (even or odd).
      %     Only used when using z_mirror.  Default: ``'even'``.
      %
      %   - rorder (numeric) -- Rotational order of incident beam.
      %     Only used when using z_rotation.  Default: ``0``.

      p = inputParser;
      p.addParameter('solver', @(A, E) A \ E, @(x) isa(x, 'function_handle'));
      p.addParameter('multiBeam', true, @islogical);
      p.addParameter('parity', 'even');
      p.addParameter('rorder', 1);
      p.parse(varargin{:});

      assert(isnumeric(Eincident) && ismatrix(Eincident) ...
          && size(Eincident, 1) == 3*dda.ndipoles, ...
          'Eindicent must be 3NxM numeric matrix');

      % Get or calculate interaction matrix
      A = dda.interaction_matrix('rorder', p.Results.rorder, ...
          'parity', p.Results.parity);

      % Calculate polarizations
      if p.Results.multiBeam
        P = p.Results.solver(A, Eincident);
      else
        P = zeros(3*dda.ndipoles, size(Eincident, 2));
        for ii = 1:size(Eincident, 2)
          P(:, ii) = p.Results.solver(A, Eincident(:, ii));
        end
      end

      % Package output
      dipoles = ott.tmatrix.dda.Dipole(dda.locations, P, ...
          'xySymmetry', dda.xySymmetry, 'zRotSymmetry', dda.zRotSymmetry, ...
          'parity', p.Results.parity, 'rorder', p.Results.rorder);
    end

    function A = interaction_matrix(dda, varargin)
      % Calculate the interaction matrix assuming memory is limited
      %
      % This method is rather time consuming and involves a lot of
      % redundant calculations if called repeatedly with different rorder
      % or parity parameters.  For a faster but more memory intensive
      % method see :class:`DdaHighMemory`.
      %
      % Usage
      %   A = dda.interaction_matrix(rorder, parity)
      %
      % Optional parameters
      %   - parity (enum) -- Parity of incident beam (even or odd).
      %     Only used when using xySymmetry.  Default: ``'even'``.
      %
      %   - rorder (numeric) -- Rotational order of incident beam.
      %     Only used when using zRotSymmetry.  Default: ``0``.

      p = inputParser;
      p.addOptional('parity', 'even');
      p.addOptional('rorder', 0);
      p.parse(varargin{:});

      A = ott.tmatrix.dda.Dipole.build_field_matrix(...
          dda.locations, dda.locations, ...
          @dda.interaction_efield_column, ...
          'low_memory', true, ...
          'xySymmetry', dda.xySymmetry, ...
          'zRotSymmetry', dda.zRotSymmetry, ...
          'rorder', p.Results.rorder, ...
          'parity', p.Results.parity);
    end
  end

  methods (Hidden)
    function D = interaction_efield_column(dda, ii, jj, kk, vxyz, nvec, M)
      % Build column of interaction matrix
      
      % Calculate interaction term
      D = ott.tmatrix.dda.Dipole.enearfield_matrix_column(vxyz, nvec, M);

      % Replace self interaction term with inverse polarizability
      if jj == 1 && kk == 1
        D((1:3) + (ii-1)*3, :) = ...
            inv(dda.polarizability(:, (1:3) + (ii-1)*3));
      end
    end
  end

  methods % Getters/Setters
    function dda = set.locations(dda, val)
      assert(isnumeric(val) && ismatrix(val) && size(val, 1) == 3, ...
          'location must be 3xN numeric matrix');
      dda.locations = val;
    end

    function dda = set.polarizability(dda, val)
      assert(isnumeric(val) && ismatrix(val) ...
          && size(val, 1) == 3 && mod(size(val, 2), 3) == 0, ...
          'polarizability must be a 3x3N matrix');
      dda.polarizability = val;
    end

    function n = get.ndipoles(dda)
      n = size(dda.locations, 2);
    end

    function dda = set.xySymmetry(dda, val)
      assert(islogical(val) && isscalar(val), ...
        'xySymmetry must be logical scalar');
      dda.xySymmetry = val;
    end

    function dda = set.zRotSymmetry(dda, val)
      assert(isnumeric(val) && isscalar(val) && round(val) == val ...
          && val >= 0, 'zRotSymmetry must be non-negative integer');
      dda.zRotSymmetry = val;
    end
  end
end

