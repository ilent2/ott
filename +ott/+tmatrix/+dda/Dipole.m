classdef Dipole
% Describes an array of radiating dipoles.
%
% This class stores the dipole locations and polarizations.  The scatted
% fields can be calculated using::
%
%   Es = F p
%
% where `F` describes the locations where fields should be calculated and
% `p` describes the polarization of each dipole.  The class provides
% methods for calculating `F` and `Es`.
%
% Methods
%   - setDipoles      -- Set the dipole data (location/polarization)
%   - efield          -- Calculate E near-fields
%   - hfield          -- Calculate H near-fields
%   - efarfield       -- Calculate E far-fields
%   - hfarfield       -- Calculate H far-fields
%   - efarfield_matrix    -- Calculate far-field matrix for field calculation
%   - hfarfield_matrix    -- Calculate far-field matrix for field calculation
%   - enearfield_matrix   -- Calculate near-field matrix for field calculation
%   - hnearfield_matrix   -- Calculate near-field matrix for field calculation
%   - mtimes              -- Apply field matrix and calculate fields
%
% Properties
%   - location        -- Dipole locations
%   - polarization    -- Dipole polarization
%   - xySymmetry      -- True if using z-mirror symmetry
%   - zRotSymmetry    -- Order of z-rotational symmetry (0 - infinite)
%   - parity          -- Parity of incident beam
%   - rorder          -- Rotational order of incident beam
%   - ndipoles        -- Number of dipoles in the array
%   - nbeams          -- Number of beams in the array

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (SetAccess=protected)
    locations         % (3xN numeric) Dipole locations
    polarization     % (3NxM numeric) Dipole polarizations
  end

  properties
    xySymmetry       % True if using z-mirror symmetry
    zRotSymmetry     % Order of z-rotational symmetry (0 - infinite)
    parity           % Parity of incident beam
    rorder           % Rotational order of incident beam
  end

  properties (Dependent)
    ndipoles         % Number of dipoles in the array
    nbeams           % Number of beams in the array
  end

  methods
    function beam = Dipole(varargin)
      % Construct a new dipole array
      %
      % Usage
      %   beam = Dipole(locations, polarization, ...)
      %   Parameters can also be passed as named arguments.
      %
      % Parameters
      %   - locations (3xN numeric) -- Locations of dipoles
      %   - polarization (3NxM) -- Dipole polarizations sorted
      %     packaged in [x1;y1;z1; x2;y2;z2; ...] order.
      %
      % Optional named parameters
      %   - parity (enum) -- Parity of incident beam (even or odd).
      %     Only used when using z_mirror.  Default: ``'even'``.
      %
      %   - rorder (numeric) -- Rotational order of incident beam.
      %     Only used when using z_rotation.  Default: ``0``.
      %
      %   - xySymmetry (logical) -- If the particle has
      %     z-mirror symmetry.  Default: ``false``.
      %
      %   - zRotSymmetry (numeric) -- Order of the particle is
      %     z-rotational symmetric.  Default: ``1``.
      %     If ``0``, uses fourth order rotational symmetry (might change
      %     in a future release).

      p = inputParser;
      p.addOptional('locations', [], @isnumeric);
      p.addOptional('polarization', [], @isnumeric);
      p.addParameter('parity', 'even');
      p.addParameter('rorder', 0);
      p.addParameter('xySymmetry', false);
      p.addParameter('zRotSymmetry', 1);
      p.parse(varargin{:});

      beam = beam.setDipoles(p.Results.locations, p.Results.polarization);
      beam.parity = p.Results.parity;
      beam.xySymmetry = p.Results.xySymmetry;
      beam.rorder = p.Results.rorder;
      beam.zRotSymmetry = p.Results.zRotSymmetry;
    end

    function beam = setDipoles(beam, locations, polarization)
      % Set the dipole position and polarization data
      %
      % Usage
      %   beam = beam.setDipoles(location, polarization)
      %
      % Parameters
      %   - location (3xN numeric) -- Locations of dipoles
      %   - polarization (3NxM) -- Dipole polarizations sorted
      %     packaged in [x1;y1;z1; x2;y2;z2; ...] order.

      ott.utils.nargoutCheck(beam, nargout);

      assert(isnumeric(locations) && ismatrix(locations) ...
          && size(locations, 1) == 3, ...
          'location must be 3xN numeric matrix');
      assert(isnumeric(polarization) && ismatrix(polarization) ...
          && size(polarization, 1) == numel(locations), ...
          'polarization must be 3NxM numeric matrix');

      beam.locations = locations;
      beam.polarization = polarization;
    end

    function F = efarfield_matrix(beam, rtp, varargin)
      % Evaluates the electric far-field matrix for a set of points
      %
      % Can be applied to the beam to evaluate the scattered fields.
      %
      %   Es = F * beam
      %
      % Usage
      %   F = beam.farfield_matrix(rtp, ...)
      %
      % Parameters
      %   - rtp (2xN|3xN numeric) -- Far-field coordinates.  Can either
      %     be [theta; phi] or [r; theta; phi] in which case `r` is
      %     ignored.
      %
      % Optional named arguments
      %   - low_memory (logical) -- If true, evaluates the low-memory
      %     version of F.  Default: ``false``.

      % Convert to Cartesian coordinates
      [~, rtp] = ott.utils.rtpFarfield(rtp);
      xyz = ott.utils.rtp2xyz(rtp);

      % Construct matrix
      F = beam.field_matrix_internal(xyz, ...
          @beam.efarfield_matrix_column, varargin{:});
    end

    function F = hfarfield_matrix(beam, rtp, varargin)
      % Evaluates the magnetic far-field matrix for a set of points
      %
      % See :meth:`efarfield_matrix` for usage and parameters.

      % Convert to Cartesian coordinates
      [~, rtp] = ott.utils.rtpFarfield(rtp);
      xyz = ott.utils.rtp2xyz(rtp);

      % Construct matrix
      F = beam.field_matrix_internal(xyz, ...
          @beam.hfarfield_matrix_column, varargin{:});
    end

    function F = enearfield_matrix(beam, xyz, varargin)
      % Evaluates the electric near-field matrix for a set of points
      %
      % Can be applied to the beam to evaluate the scattered fields.
      %
      %   Es = F * beam
      %
      % Usage
      %   F = beam.nearfield_matrix(xyz, ...)
      %
      % Parameters
      %   - xyz (3xN numeric) -- Near-field locations.
      %
      % Optional named arguments
      %   - low_memory (logical) -- If true, evaluates the low-memory
      %     version of F.  Default: ``false``.

      assert(isnumeric(xyz) && ismatrix(xyz) && size(xyz, 1) == 3, ...
          'xyz must be 3xN numeric matrix');

      F = beam.field_matrix_internal(xyz, ...
          @beam.enearfield_matrix_column, varargin{:});
    end

    function F = hnearfield_matrix(beam, xyz, varargin)
      % Evaluates the magnetic near-field matrix for a set of points
      %
      % See :meth:`enearfield_matrix` for parameters and usage.

      assert(isnumeric(xyz) && ismatrix(xyz) && size(xyz, 1) == 3, ...
          'xyz must be 3xN numeric matrix');

      F = beam.field_matrix_internal(xyz, ...
          @beam.hnearfield_matrix_column, varargin{:});
    end

    function Es = mtimes(F, beam)
      % Matrix multiplication for calculating scattered fields
      %
      % Evaluates::
      %
      %   Es = F * p
      %
      % Usage
      %   Es = F * beam
      %   Does not reshape the output.
      %
      % Parameters
      %   - F (numeric) -- Field matrix calculated using
      %     :meth:`enearfield_matrix`, :meth:`efarfield_matrix`,
      %     :meth:`hnearfield_matrix`, or :meth:`hfarfield_matrix`.

      assert(isa(beam, 'ott.tmatrix.dda.Dipole'), ...
          'Second argument to mtimes must be Dipole beam');

      assert(size(F, 2) == size(beam.polarization, 1), ...
          'Matrix dimensions must agree');

      % Convert from mirror/rotsym F to full F
      if ~ismatrix(F)
        F = beam.combine_rotsym_matrix(F, beam.rorder, beam.parity);
      end

      Es = F * beam.polarization;
    end

    function E = efield(beam, xyz, varargin)
      % Calculate the E-field
      %
      % Evaluates::
      %
      %     Es = F * p
      %
      % Usage
      %   E = beam.efield(xyz)
      %
      % Parameters
      %   - xyz -- (3xN numeric) Cartesian coordinates.
      %
      % Unmatched parameters are passed to :meth:`enearfield_matrix`.

      % Calculate the near-field matrix
      F = beam.enearfield_matrix(xyz, varargin{:});

      E = F * beam;

      % Package output
      E = ott.utils.FieldVector(xyz, reshape(E, 3, []), 'cartesian');
    end

    function H = hfield(beam, xyz, varargin)
      % Calculate the H-field
      %
      % Usage
      %   E = beam.efield(xyz)
      %
      % Parameters
      %   - xyz -- (3xN numeric) Cartesian coordinates.
      %
      % Unmatched parameters are passed to :meth:`hnearfield_matrix`.

      % Calculate the near-field matrix
      F = beam.hnearfield_matrix(xyz, varargin{:});

      H = F * beam;

      % Package output
      H = ott.utils.FieldVector(xyz, reshape(H, 3, []), 'cartesian');
    end

    function E = efarfield(beam, rtp, varargin)
      % Calculate the E-field
      %
      % Usage
      %   E = beam.efarfield(rtp)
      %
      % Parameters
      %   - rtp -- (3xN | 2xN numeric) Spherical coordinates.
      %     Either [radius; theta; phi] or [theta; phi].  Radius is ignored.
      %
      % Unmatched parameters are passed to :meth:`efarfield_matrix`.

      % Ensure rtp is 3xN
      [~, rtp] = ott.utils.rtpFarfield(rtp);

      % Calculate the near-field matrix
      F = beam.efarfield_matrix(rtp, varargin{:});

      E = F * beam;

      % Package output
      xyz = ott.utils.rtp2xyz(rtp);
      E = ott.utils.FieldVector(xyz, reshape(E, 3, []), 'cartesian');
    end

    function H = hfarfield(beam, rtp, varargin)
      % Calculate the H-field
      %
      % Usage
      %   H = beam.hfarfield(rtp)
      %
      % Parameters
      %   - rtp -- (3xN | 2xN numeric) Spherical coordinates.
      %     Either [radius; theta; phi] or [theta; phi].  Radius is ignored.
      %
      % Unmatched parameters are passed to :meth:`hfarfield_matrix`.

      % Ensure rtp is 3xN
      [~, rtp] = ott.utils.rtpFarfield(rtp);

      % Calculate the near-field matrix
      F = beam.hfarfield_matrix(rtp, varargin{:});

      H = F * beam;

      % Package output
      xyz = ott.utils.rtp2xyz(rtp);
      H = ott.utils.FieldVector(xyz, reshape(H, 3, []), 'cartesian');
    end
  end

  methods (Hidden, Static)
    function D = enearfield_matrix_column(dipole_xyz, target_xyz, M)
      % Calculate a column of the near-field matrix

      assert(size(dipole_xyz, 1) == 3, 'dipole_xyz must be 3xN array');
      assert(size(target_xyz, 1) == 3, 'target_xyz must be 3xN array');
      assert(all(size(M) == [3,3]), 'M must be 3x3 matrix');

      n_targets = size(target_xyz, 2);
      k = 2*pi;

      % Calculate distance and unit vector from dipole to targets
      r_vec = target_xyz - dipole_xyz;
      r_jk = vecnorm(r_vec);
      r_hat = r_vec./r_jk;

      % Compute outer product of r_hat * r_hat
      rr = reshape(r_hat, 3, 1, n_targets) .* reshape(r_hat, 1, 3, n_targets);

      % Reshape arrays for multiplication
      rr = reshape(rr, [3, 3, n_targets]);
      r_jk = reshape(r_jk, [1, 1, n_targets]);

      % Calculate Greens function between dipole and target location
      % This matches Eq 17 from https://doi.org/10.1364/JOSAA.18.001944
      F = -exp(1i*k*r_jk)./r_jk .* ...
        (k^2*(rr - eye(3)) + (1i*k*r_jk - 1)./r_jk.^2.*(3*rr - eye(3)));

      % Apply coordinate transformations for dipoles/targets (D = F * M)
      D = sum(reshape(F, [3, 3, 1, n_targets]) ...
          .* reshape(M, [1, size(M)]), 2);
      D = reshape(D, 3, 3, n_targets);

      % Change back to 3Nx3 result
      D = permute(D, [1, 3, 2]);
      D = reshape(D, [3*n_targets, 3]);
    end

    function D = hnearfield_matrix_column(dipole_xyz, target_xyz, M)
      % Calculate a column of the magnetic near-field matrix

      assert(size(dipole_xyz, 1) == 3, 'dipole_xyz must be 3xN array');
      assert(size(target_xyz, 1) == 3, 'target_xyz must be 3xN array');
      assert(all(size(M) == [3,3]), 'M must be 3x3 matrix');

      n_targets = size(target_xyz, 2);
      k = 2*pi;

      % Calculate distance and unit vector from dipole to targets
      r_vec = target_xyz - dipole_xyz;
      r_jk = vecnorm(r_vec);

      % Compute the cross-product matrix
      rcross = zeros(3, 3, n_targets);
      rcross(1, 2, :) = -r_vec(3, :);
      rcross(1, 3, :) = r_vec(2, :);
      rcross(2, 3, :) = -r_vec(1, :);
      rcross(2, 1, :) = r_vec(3, :);
      rcross(3, 1, :) = -r_vec(2, :);
      rcross(3, 2, :) = r_vec(1, :);

      % Reshape arrays for multiplication
      r_jk = reshape(r_jk, [1, 1, n_targets]);

      % Calculate Greens function between dipole and target location
      % Eq 9 from https://doi.org/10.1364/JOSAA.18.001944
      F = k^2 .* exp(1i*k*r_jk)./r_jk.^2 .* rcross .* (1 - 1./(1i.*k.*r_jk));

      % Apply coordinate transformations for dipoles/targets (D = F * M)
      D = sum(reshape(F, [3, 3, 1, n_targets]) ...
          .* reshape(M, [1, size(M)]), 2);
      D = reshape(D, 3, 3, n_targets);

      % Change back to 3Nx3 result
      D = permute(D, [1, 3, 2]);
      D = reshape(D, [3*n_targets, 3]);
    end

    function D = hfarfield_matrix_column(dipole_xyz, n_vec, M)
      % Calculate a column of the far-field matrix

      assert(size(dipole_xyz, 1) == 3, 'dipole_xyz must be 3xN array');
      assert(size(n_vec, 1) == 3, 'n_vec must be 3xN array');

      k = 2*pi;
      k_vec = k .* n_vec;
      k_rd = sum(k_vec .* dipole_xyz, 1);

      n_targets = size(n_vec, 2);

      % Compute the cross-product matrix
      rcross = zeros(3, 3, n_targets);
      rcross(1, 2, :) = -n_vec(3, :);
      rcross(1, 3, :) = n_vec(2, :);
      rcross(2, 3, :) = -n_vec(1, :);
      rcross(2, 1, :) = n_vec(3, :);
      rcross(3, 1, :) = -n_vec(2, :);
      rcross(3, 2, :) = n_vec(1, :);

      % Reshape arrays for multiplication
      k_rd = reshape(k_rd, [1, 1, n_targets]);

      % Far-field limit of dipole approximation (only 1/r term)
      F = k^2 .* exp(-1i .* k_rd) .* rcross;

      % Apply coordinate transformations for dipoles/targets (D = F * M)
      D = sum(reshape(F, [3, 3, 1, n_targets]) ...
          .* reshape(M, [1, size(M)]), 2);
      D = reshape(D, 3, 3, n_targets);

      % Change back to 3Nx3 result
      D = permute(D, [1, 3, 2]);
      D = reshape(D, [3*n_targets, 3]);
    end

    function D = efarfield_matrix_column(dipole_xyz, n_vec, M)
      % Calculate a column of the magnetic far-field matrix

      assert(size(dipole_xyz, 1) == 3, 'dipole_xyz must be 3xN array');
      assert(size(n_vec, 1) == 3, 'n_vec must be 3xN array');

      k = 2*pi;
      k_vec = k .* n_vec;
      k_rd = sum(k_vec .* dipole_xyz, 1);

      n_targets = size(n_vec, 2);

      nn = reshape(n_vec, 3, 1, n_targets) .* reshape(n_vec, 1, 3, n_targets);

      % Reshape arrays for multiplication
      nn = reshape(nn, [3, 3, n_targets]);
      k_rd = reshape(k_rd, [1, 1, n_targets]);

      % Far-field limit of dipole approximation (only 1/r term)
      %
      % Matches equation 5 from
      % https://www.sciencedirect.com/science/article/pii/S0022407315002307
      F = -k^2 .* exp(-1i .* k_rd) .* (nn - eye(3));

      % Apply coordinate transformations for dipoles/targets (D = F * M)
      D = sum(reshape(F, [3, 3, 1, n_targets]) ...
          .* reshape(M, [1, size(M)]), 2);
      D = reshape(D, 3, 3, n_targets);

      % Change back to 3Nx3 result
      D = permute(D, [1, 3, 2]);
      D = reshape(D, [3*n_targets, 3]);
    end

    function F = combine_rotsym_matrix(F, zorder, parity)
      % Warning: Function usage/definition may change without notice

      % Combine z rotational symmetry slices
      if size(F, 3) ~= 1
        midx = 1:size(F, 3);
        phase_factor = exp(1i*zorder*2*pi*(midx-1)./size(F, 3));
        F = sum(F .* reshape(phase_factor, 1, 1, [], 1), 3);
      end

      % Calculate even and odd parity reflections
      if size(F, 4) ~= 1
        if strcmpi(parity, 'odd')
          F = F(:, :, :, 1) - F(:, :, :, 2);
        else
          F = F(:, :, :, 1) + F(:, :, :, 2);
        end
      end
    end

    function F = build_field_matrix(vxyz, txyz, func, varargin)
      % Build field/interaction matrix
      %
      % This function includes a check that the generated matrix will
      % fit into available physical memory.  This check could raise
      % two warnings:
      %
      %   - 'ott:utils:arrayMaxSize:cant_get_memory' is raised if the
      %     memory command fails.
      %   - 'ott:tmatrix:dda:Dipole:memory_may_exceed_physical' is raised
      %     if the memory exceeds available physical memory.
      %
      % You may want to set the error flag for these warnings with::
      %
      %   warning('error', MSGID)
      %
      % Usage
      %   F = ott.tmatrix.dda.Dipole.build_field_matrix(vxyz, txyz, func, ...)
      %
      % Parameters
      %   - vxyz -- (3xN numeric) Dipole locations
      %
      %   - txyz -- (3xN numeric) Target locations
      %
      %   - func -- (function handle) Field function with the signature
      %     ``func(col_idx, rot_idx, mirror_idx, dxyz, txyz, M)``
      %     where ``_idx`` are the column, rotation and mirror indices,
      %     ``dxyz`` and ``txyz`` are the dipole locations, and ``M`` is
      %     a rotation matrix.
      %
      % Optional named parameters
      %   - low_memory -- (logical) If matrix should be low memory variant.
      %
      %   - xySymmetry -- (logical) If matrix has xy symmetry.
      %
      %   - parity -- (enum) Either 'even' or 'odd' for beam parity.
      %
      %   - zRotSymmetry -- (numeric) Degree of z rotational symmetry.
      %
      %   - rorder -- (numeric) Rotational order of beam.

      p = inputParser;
      p.addParameter('low_memory', false);
      p.addParameter('xySymmetry', false);
      p.addParameter('zRotSymmetry', 1);
      p.addParameter('parity', 'even');
      p.addParameter('rorder', 1);
      p.parse(varargin{:});

      % Calculate amount of work to do
      rotsym_wk = p.Results.zRotSymmetry;
      mirror_wk = p.Results.xySymmetry + 1;
      if rotsym_wk == 0
        rotsym_wk = 4;
      end

      % Calculate matrix sizes
      if p.Results.low_memory
        rotsym_sz = 1;
        mirror_sz = 1;
      else
        rotsym_sz = rotsym_wk;
        mirror_sz = mirror_wk;
      end
      ncols = size(vxyz, 2);
      nrows3 = 3*size(txyz, 2);

      % Do a memory check before trying the following
      if prod([nrows3, 3, mirror_sz, rotsym_sz, ncols]) ...
          > ott.utils.arrayMaxSize()
        warning('ott:tmatrix:dda:Dipole:memory_may_exceed_physical', ...
            'Required memory may exceed available physical memory');
      end

      % Allocate memory for F
      F = zeros(nrows3, 3, mirror_sz, rotsym_sz, ncols);

      % Calculate phase factors for mirror/rotational symmetry
      if p.Results.low_memory
        midx = 1:rotsym_wk;
        rphase_factor = exp(1i*p.Results.rorder*2*pi*(midx-1)/rotsym_wk);
        if strcmpi(p.Results.parity, 'even')
          mphase_factor = [1, 1];
        else
          mphase_factor = [1, -1];
        end
      end

      % Get Spherical Coordinates of voxels
      vrtp = ott.utils.xyz2rtp(vxyz);

      for ii = 1:size(vxyz, 2)        % Loop over each dipole

        % Calculate cartesian to spherical conversion for each dipole
        M_cart2sph = ott.utils.cart2sph_mat(vrtp(2, ii), vrtp(3, ii));

        for jj = 1:rotsym_wk     % duplicates for z mirror symmetry

          % Calculate spherical to cartesian for each mirror version
          M_sph2cart = ott.utils.sph2cart_mat(...
              vrtp(2, ii), vrtp(3, ii) + 2*pi*(jj-1)/rotsym_wk);

          dipole_xyz = ott.utils.rtp2xyz(vrtp(1, ii), ...
            vrtp(2, ii), vrtp(3, ii) + 2*pi*(jj-1)/rotsym_wk);

          for kk = 1:mirror_wk     % duplicates for z rotational symmetry

            % Apply mirror symmetry rotation
            if kk == 1
              M_dipole = M_sph2cart * M_cart2sph;
            else
              M_dipole = diag([1, 1, -1]) * M_sph2cart * M_cart2sph;
              dipole_xyz(3) = -dipole_xyz(3);
            end

            % Calculate columns of F
            Fc = func(ii, jj, kk, dipole_xyz, txyz, M_dipole);

            if p.Results.low_memory
              F(:, :, 1, 1, ii) = F(:, :, 1, 1, ii) ...
                  + Fc .* rphase_factor(jj) .* mphase_factor(kk);
            else
              F(:, :, kk, jj, ii) = Fc;
            end
          end
        end
      end

      % Convert from 3xN*3*M*L*O to 3xN*3xM*L*O
      F = permute(F, [1, 2, 5, 4, 3]);
      F = reshape(F, [nrows3, 3*ncols, rotsym_sz, mirror_sz]);
    end
  end

  methods (Hidden)
    function F = field_matrix_internal(beam, locs, field_func, varargin)
      % Calculate the field-matrices

      p = inputParser;
      p.addParameter('low_memory', false);
      p.parse(varargin{:});

      % Don't need row/col indices for our field funcs, so discard them
      func = @(~, ~, ~, a, b, c) field_func(a, b, c);

      F = beam.build_field_matrix(beam.locations, locs, func, ...
          'low_memory', p.Results.low_memory, ...
          'xySymmetry', beam.xySymmetry, ...
          'zRotSymmetry', beam.zRotSymmetry, ...
          'rorder', beam.rorder, ...
          'parity', beam.parity);
    end
  end

  methods % Getters/setters
    function beam = set.polarization(beam, val)
      assert(isnumeric(val) && ismatrix(val) && mod(size(val, 1), 3) == 0, ...
          'polarization must be a 3NxM numeric matrix');
      beam.polarization = val;
    end

    function beam = set.locations(beam, val)
      assert(isnumeric(val) && ismatrix(val) && size(val, 1) == 3, ...
          'location must be 3xN numeric matrix');
      beam.locations = val;
    end

    function n = get.ndipoles(beam)
      n = size(beam.locations, 2);
    end

    function n = get.nbeams(beam)
      n = size(beam.polarization, 2);
    end

    function beam = set.xySymmetry(beam, val)
      assert(islogical(val) && isscalar(val), ...
        'z_mirror must be logical scalar');
      beam.xySymmetry = val;
    end

    function beam = set.zRotSymmetry(beam, val)
      assert(isnumeric(val) && isscalar(val) && round(val) == val ...
          && val >= 0, 'zRotSymmetry must be non-negative integer');
      beam.zRotSymmetry = val;
    end

    function beam = set.parity(beam, val)
      assert(any(strcmpi(val, {'even', 'odd'})), ...
          'Parity must be ''even'' or ''odd''');
      beam.parity = val;
    end

    function beam = set.rorder(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
          'rorder must be numeric scalar');
      beam.rorder = val;
    end
  end
end

