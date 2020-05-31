classdef Dipole < ott.beam.properties.DipoleArray ...
    & ott.beam.properties.AnyArrayType ...
    & ott.beam.Beam
% Describes the field produced by a polarisable dipole.
% Inherits from :class:`ott.beam.Beam`, :class:`ott.beam.properties.Dipole`
% and :class:`ott.beam.utils.AnyArrayType`.
%
% Stores the dipole locations and polarisabiities.  The scattered fields
% are evaluated using::
%
%   Es = F p
%
% This class can be used to describe coherent arrays of dipoles and
% arrays of beams from coherent dipoles collections.
%
% This class supports rotational symmetry about the z-axis and
% and mirror symmetry about the xy-plane.
%
% Methods
%   - efarfield_matrix   -- Generate electric far-field F-matrix
%   - hfarfield_matrix   -- Generate magnetic far-field F-matrix
%   - enearfield_matrix  -- Generate electric near-field F-matrix
%   - hnearfield_matrix  -- Generate magnetic near-field F-matrix
%   - mtimes             -- Calculate scattered fields using F-matrix
%   - setDipoles         -- Set the dipole locations and polarizations
%   - size               -- Size of the beam array
%
% Properties
%   - location        -- Dipole locations
%   - polarization    -- Dipole polarization
%   - z_mirror        -- True if using z-mirror symmetry
%   - z_rotation      -- Order of z-rotational symmetry (0 - infinite)
%   - parity          -- Parity of incident beam
%   - rorder          -- Rotational order of incident beam
%   - ndipoles        -- Number of dipoles in the array `size(beam, 1)`
%   - nbeams          -- Number of beams in the array `size(beam, 2)`

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    z_mirror          % True if using z-mirror symmetry
    z_rotation        % Order of z-rotational symmetry (0 - infinite)
    parity            % Parity of incident beam
    rorder            % Rotational order of incident beam
  end

  properties (Dependent)
    power             % Power of scattered field
  end

  methods (Static)
    function args = likeProperties(other, args)
      % Construct an array of like-properties
      if isa(other, 'ott.beam.Dipole')
        args = ott.utils.addDefaultParameter('z_mirror', other.z_mirror, args);
        args = ott.utils.addDefaultParameter('parity', other.parity, args);
        args = ott.utils.addDefaultParameter('rorder', other.rorder, args);
        args = ott.utils.addDefaultParameter(...
            'z_rotation', other.z_rotation, args);
      end
      args = ott.beam.properties.DipoleArray.likeProperties(other, args);
    end

    function beam = like(other, varargin)
      % Create a beam like another beam
      %
      % Usage
      %   beam = Dipole.like(other, ...)
      %
      % See constructor for arguments.

      args = ott.beam.Dipole.likeProperties(other, varargin);
      beam = ott.beam.Dipole(args{:});
    end
  end

  methods
    function beam = Dipole(varargin)
      % Construct a beam representing a dipole's scattered field
      %
      % Usage
      %   beam = Dipole(location, polarization, ...)
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
      %   - z_mirror (logical) -- If the particle/beam has z-mirror symmetry.
      %     Default: ``false``.
      %
      %   - z_rotation (numeric) -- Order of the particle/beam has
      %     xy-rotational symmetry.  Default: ``1``.
      %     If ``0``, uses fourth order rotational symmetry (might change
      %     in a future release).
      %
      %   - array_type (enum) -- Type of beam array.  Default: ``array``.

      p = inputParser;
      p.addOptional('location', [], @isnumeric);
      p.addOptional('polarization', [], @isnumeric);
      p.addParameter('parity', 'even');
      p.addParameter('rorder', 1);
      p.addParameter('z_mirror', false);
      p.addParameter('z_rotation', 1);
      p.addParameter('array_type', 'array');
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Construct base
      beam = beam@ott.beam.properties.AnyArrayType(p.Results.array_type);
      beam = beam@ott.beam.properties.DipoleArray(unmatched{:}, ...
          'polarization', p.Results.polarization, ...
          'location', p.Results.location);

      % Store optional parameters
      beam.parity = p.Results.parity;
      beam.z_mirror = p.Results.z_mirror;
      beam.rorder = p.Results.rorder;
      beam.z_rotation = p.Results.z_rotation;
    end

    function beam = defaultArrayType(beam, array_type, elements)
      % Construct a new array for this type

      if beam.contains('array')
        assert(strcmpi(array_type, 'array'), ...
            'type must be array for beams containing generic arrays');
      elseif beam.contains('incoherent')
        assert(strcmpi(array_type, {'array', 'incoherent'}), ...
            'type must be array/incoherent for incoherent content');
      end

      if nargin == 2
        if strcmpi(array_type, 'coherent')
          beam = @(arg) ott.beam.Dipole(arg);
        else
          beam = @(arg) ott.beam.Array(array_type, arg);
        end
      else
        if strcmpi(array_type, 'coherent')
          beam = ott.beam.Dipole(elements);
        else
          beam = ott.beam.Array(array_type, elements);
        end
      end
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

      assert(isnumeric(rtp) && ismatrix(rtp) ...
          && size(rtp, 1) == 2 || size(rtp, 1) == 3, ...
          'rtp must be 2xN or 3xN numeric matrix');

      % Convert to unit vectors
      if size(rtp, 1) == 3
        xyz = ott.utils.rtp2xyz([ones(1, size(rtp, 2)), rtp(2:3, :)]);
      else
        xyz = ott.utils.rtp2xyz([ones(1, size(rtp, 2)), rtp]);
      end

      % Construct matrix
      F = beam.field_matrix_internal(xyz, ...
          @beam.efarfield_matrix_column, varargin{:});
    end

    function F = hfarfield_matrix(beam, rtp, varargin)
      % Evaluates the magnetic far-field matrix for a set of points
      %
      % See :meth:`efarfield_matrix` for usage and parameters.

      assert(isnumeric(rtp) && ismatrix(rtp) ...
          && size(rtp, 1) == 2 || size(rtp, 1) == 3, ...
          'rtp must be 2xN or 3xN numeric matrix');

      % Convert to unit vectors
      if size(rtp, 1) == 3
        xyz = ott.utils.rtp2xyz([ones(1, size(rtp, 2)); rtp(2:3, :)]);
      else
        xyz = ott.utils.rtp2xyz([ones(1, size(rtp, 2)); rtp]);
      end

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
      % Parametesr
      %   - F (numeric) -- Field matrix calculated using
      %     :meth:`nearfield_matrix` or :meth:`farfield_matrix`.

      assert(isa(beam, 'ott.beam.Dipole'), ...
          'Second argument to mtimes must be Dipole beam');

      assert(size(F, 2) == size(beam.polarization, 1), ...
          'Matrix dimensions must agree');

      % Convert from mirror/rotsym F to full F
      if ~ismatrix(F)
        % Combine z rotational symmetry slices
        if size(F, 3) ~= 1
          midx = 1:beam.z_rotation;
          phase_factor = exp(1i*beam.rorder*2*pi*(midx-1)/beam.z_rotation);
          F = sum(F .* reshape(phase_factor, [1, 1, beam.z_rotation, 1]), 3);
        end

        % Calculate even and odd parity reflections
        if size(F, 4) ~= 1
          if strcmpi(beam.parity, 'odd')
            F = F(:, :, :, 1) - F(:, :, :, 2);
          else
            F = F(:, :, :, 1) + F(:, :, :, 2);
          end
        end
      end

      % This could be moved to the setter instead
      polarization = beam.polarization;
      if size(polarization, 1) == 3
        polarization = repmat(polarization, beam.ndipoles, 1);
      end

      % Calculate fields
      Es = F * polarization;
    end
  end

  methods (Hidden)
    function D = enearfield_matrix_column(beam, dipole_xyz, target_xyz, M)
      % Calculate a column of the near-field matrix

      assert(size(dipole_xyz, 1) == 3, 'dipole_xyz must be 3xN array');
      assert(size(target_xyz, 1) == 3, 'target_xyz must be 3xN array');
      assert(all(size(M) == [3,3]), 'M must be 3x3 matrix');

      n_targets = size(target_xyz, 2);
      k = beam.wavenumber;

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

    function D = hnearfield_matrix_column(beam, dipole_xyz, target_xyz, M)
      % Calculate a column of the magnetic near-field matrix

      assert(size(dipole_xyz, 1) == 3, 'dipole_xyz must be 3xN array');
      assert(size(target_xyz, 1) == 3, 'target_xyz must be 3xN array');
      assert(all(size(M) == [3,3]), 'M must be 3x3 matrix');

      n_targets = size(target_xyz, 2);
      k = beam.wavenumber;

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

    function D = efarfield_matrix_column(beam, dipole_xyz, n_vec, M)
      % Calculate a column of the far-field matrix

      assert(size(dipole_xyz, 1) == 3, 'dipole_xyz must be 3xN array');
      assert(size(n_vec, 1) == 3, 'n_vec must be 3xN array');

      k = beam.wavenumber;
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

    function D = hfarfield_matrix_column(beam, dipole_xyz, n_vec, M)
      % Calculate a column of the magnetic far-field matrix

      assert(size(dipole_xyz, 1) == 3, 'dipole_xyz must be 3xN array');
      assert(size(n_vec, 1) == 3, 'n_vec must be 3xN array');

      k = beam.wavenumber;
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

    function F = field_matrix_internal(beam, locs, field_func, varargin)
      % Calculate the field-matrices

      p = inputParser;
      p.addParameter('low_memory', false);
      p.parse(varargin{:});

      % Calculate amount of work we have to do
      full_rotsym_sz = beam.z_rotation;
      if full_rotsym_sz == 0
        full_rotsym_sz = 4;
      end
      if beam.z_mirror
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
      ndipoles = beam.ndipoles;

      % Allocate memory for F
      F = zeros(nrows, 3, mirror_sz, rotsym_sz, ndipoles);

      % Calculate phase factor for rotational symmetry
      if p.Results.low_memory
        midx = 1:beam.z_rotation;
        rphase_factor = exp(1i*beam.rorder*2*pi*(midx-1)/beam.z_rotation);
      end
      mphase_factor = [1, -1];

      % Get Spherical Coordinates
      rtp = ott.utils.xyz2rtp(beam.location);

      for ii = 1:ndipoles        % Loop over each dipole

        % Calculate cartesian to spherical conversion for each dipole
        M_cart2sph = ott.utils.cart2sph_mat(rtp(2, ii), rtp(3, ii));

        for jj = 1:full_rotsym_sz     % duplicates for z mirror symmetry

          % Calculate spherical to cartesian for each mirror version
          M_sph2cart = ott.utils.sph2cart_mat(...
              rtp(2, ii), rtp(3, ii) + 2*pi*(jj-1)/beam.z_rotation);

          dipole_xyz = ott.utils.rtp2xyz(rtp(1, ii), ...
            rtp(2, ii), rtp(3, ii) + 2*pi*(jj-1)/beam.z_rotation);

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
      F = reshape(F, [nrows, 3*ndipoles, rotsym_sz, mirror_sz]);
    end

    function E = efieldInternal(beam, xyz, varargin)
      % Calculate the E-field
      %
      % Evaluates::
      %
      %     Es = F * p

      % Calculate the near-field matrix
      F = beam.enearfield_matrix(xyz);

      E = F * beam;

      % Package output
      E = beam.packageArrayOutput(E, xyz, 'cartesian');
    end

    function H = hfieldInternal(beam, xyz, varargin)
      % Calculate the H-field

      % Calculate the near-field matrix
      F = beam.hnearfield_matrix(xyz);

      H = F * beam;

      % Package output
      H = beam.packageArrayOutput(H, xyz, 'cartesian');
    end

    function E = efarfieldInternal(beam, rtp, varargin)
      % Calculate the E-field

      % Calculate the near-field matrix
      F = beam.efarfield_matrix(rtp);

      E = F * beam;

      % Package output
      rtp(1, :) = 1.0;
      xyz = ott.utils.rtp2xyz(rtp);
      E = beam.packageArrayOutput(E, xyz, 'cartesian');
    end

    function H = hfarfieldInternal(beam, rtp, varargin)
      % Calculate the H-field

      % Calculate the near-field matrix
      F = beam.hfarfield_matrix(rtp);

      H = F * beam;

      % Package output
      rtp(1, :) = 1.0;
      xyz = ott.utils.rtp2xyz(rtp);
      H = beam.packageArrayOutput(H, xyz, 'cartesian');
    end

    function E = packageArrayOutput(beam, E, cords, cords_type)
      % Package output for array types
      E = reshape(E, 3, beam.nbeams, []);
      if strcmpi(beam.array_type, 'coherent') || beam.nbeams == 1
        E = squeeze(sum(E, 2));
        E = ott.utils.FieldVector(cords, E, cords_type);
      else
        E = mat2cell(E, 3, ones(1, size(E, 2)), size(E, 3));
        E = cellfun(@squeeze, E, 'UniformOutput', false);
        E = cellfun(@(x) ott.utils.FieldVector(cords, x, cords_type), ...
            E, 'UniformOutput', false);
      end
    end
  end

  methods % Getters/setters
    % z_mirror          % True if using z-mirror symmetry
    % z_rotation        % Order of z-rotational symmetry (0 - infinite)
    % parity            % Parity of incident beam
    % rorder            % Rotational order of incident beam

    function beam = set.z_mirror(beam, val)
      assert(islogical(val) && isscalar(val), ...
        'z_mirror must be logical scalar');
      beam.z_mirror = val;
    end

    function beam = set.z_rotation(beam, val)
      assert(isnumeric(val) && isscalar(val) && round(val) == val ...
          && val >= 0, 'z_rotation must be non-negative integer');
      beam.z_rotation = val;
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

    function p = get.power(~)
      % TODO: Implement this, numerical integration?
      %   Are there cases with simpler solutions?
      error('Not yet implemented');
    end
    function beam = set.power(beam, val)
      % TODO: Is this something we want to support?  Can it be done by
      % scaling the dipole magnitudes?
      error('Not yet implemented');
    end
  end
end
