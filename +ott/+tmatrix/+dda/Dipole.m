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
    location         % (3xN numeric) Dipole locations
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
      %   beam = Dipole(location, polarization, ...)
      %   Parameters can also be passed as named arguments.
      %
      % Parameters
      %   - location (3xN numeric) -- Locations of dipoles
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
      p.addOptional('location', [], @isnumeric);
      p.addOptional('polarization', [], @isnumeric);
      p.addParameter('parity', 'even');
      p.addParameter('rorder', 0);
      p.addParameter('xySymmetry', false);
      p.addParameter('zRotSymmetry', 1);
      p.parse(varargin{:});

      beam = beam.setDipoles(p.Results.location, p.Results.polarization);
      beam.parity = p.Results.parity;
      beam.xySymmetry = p.Results.xySymmetry;
      beam.rorder = p.Results.rorder;
      beam.zRotSymmetry = p.Results.zRotSymmetry;
    end

    function beam = setDipoles(beam, location, polarization)
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

      assert(isnumeric(location) && ismatrix(location) ...
          && size(location, 1) == 3, ...
          'location must be 3xN numeric matrix');
      assert(isnumeric(polarization) && ismatrix(polarization) ...
          && size(polarization, 1) == numel(location), ...
          'polarization must be 3NxM numeric matrix');

      beam.location = location;
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
      r_hat = r_vec./r_jk;

      % Compute the cross-product matrix
      rcross = zeros(3, 3, n_targets);
      rcross(1, 2, :) = -r_hat(3, :);
      rcross(1, 3, :) = r_hat(2, :);
      rcross(2, 3, :) = -r_hat(1, :);
      rcross(2, 1, :) = r_hat(3, :);
      rcross(3, 1, :) = -r_hat(2, :);
      rcross(3, 2, :) = r_hat(1, :);

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

    function D = efarfield_matrix_column(dipole_xyz, n_vec, M)
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

    function D = hfarfield_matrix_column(dipole_xyz, n_vec, M)
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
  end

  methods (Hidden)
    function F = field_matrix_internal(beam, locs, field_func, varargin)
      % Calculate the field-matrices

      p = inputParser;
      p.addParameter('low_memory', false);
      p.parse(varargin{:});

      % Calculate amount of work we have to do
      full_rotsym_sz = beam.zRotSymmetry;
      if full_rotsym_sz == 0
        full_rotsym_sz = 4;
      end
      if beam.xySymmetry
        full_mirror_sz = 2;
      else
        full_mirror_sz = 1;
      end

      % Calculate rotsym and mirror sizes for F
      if p.Results.low_memory
        rotsym_sz = 1;
        mirror_sz = 1;
      else
        rotsym_sz = full_rotsym_sz;
        mirror_sz = full_mirror_sz;
      end

      % Get sizes of xyz and dipoles
      npts = size(locs, 2);
      nrows = 3 * npts;

      % Allocate memory for F
      F = zeros(nrows, 3, mirror_sz, rotsym_sz, beam.ndipoles);

      % Calculate phase factor for rotational symmetry
      if p.Results.low_memory
        midx = 1:full_rotsym_sz;
        rphase_factor = exp(1i*beam.rorder*2*pi*(midx-1)/full_rotsym_sz);
        if strcmpi(beam.parity, 'even')
          mphase_factor = [1, 1];
        else
          mphase_factor = [1, -1];
        end
      end

      % Get Spherical Coordinates
      rtp = ott.utils.xyz2rtp(beam.location);

      for ii = 1:beam.ndipoles        % Loop over each dipole

        % Calculate cartesian to spherical conversion for each dipole
        M_cart2sph = ott.utils.cart2sph_mat(rtp(2, ii), rtp(3, ii));

        for jj = 1:full_rotsym_sz     % duplicates for z mirror symmetry

          % Calculate spherical to cartesian for each mirror version
          M_sph2cart = ott.utils.sph2cart_mat(...
              rtp(2, ii), rtp(3, ii) + 2*pi*(jj-1)/full_rotsym_sz);

          dipole_xyz = ott.utils.rtp2xyz(rtp(1, ii), ...
            rtp(2, ii), rtp(3, ii) + 2*pi*(jj-1)/full_rotsym_sz);

          for kk = 1:full_mirror_sz     % duplicates for z rotational symmetry

            % Apply mirror symmetry rotation
            if kk == 1
              M_dipole = M_sph2cart * M_cart2sph;
            else
              M_dipole = diag([1, 1, -1]) * M_sph2cart * M_cart2sph;
              dipole_xyz(3) = -dipole_xyz(3);
            end

            % Calculate columns of F
            Fc = field_func(dipole_xyz, locs, M_dipole);

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
      F = reshape(F, [nrows, 3*beam.ndipoles, rotsym_sz, mirror_sz]);
    end
  end

  methods % Getters/setters
    function beam = set.polarization(beam, val)
      assert(isnumeric(val) && ismatrix(val) && mod(size(val, 1), 3) == 0, ...
          'polarization must be a 3NxM numeric matrix');
      beam.polarization = val;
    end

    function beam = set.location(beam, val)
      assert(isnumeric(val) && ismatrix(val) && size(val, 1) == 3, ...
          'location must be 3xN numeric matrix');
      beam.location = val;
    end

    function n = get.ndipoles(beam)
      n = size(beam.location, 2);
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

