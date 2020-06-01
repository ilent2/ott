classdef NearfieldPm < ott.beam.vswf.Pointmatch
% Point-matching for VSWF in the near-field.
% Inherits from :class:`Pointmatch`.
%
% Performs point-matching in the near-field, that is, solves the system::
%
%   M \vec{ab} = \vec{E}
%
% where :math:`\vec{ab}` are the VSWF coefficients, :math:`\vec{E}` is
% the near-field and :math:`M` describes the near-field corresponding
% to each VSWF coefficient.  This works best if the system of equations
% is over-determined (i.e., more near-field point than VSWF coefficients).
%
% Static methods
%   - coefficient_matrix      -- Generate near-field coefficient matrix.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Static)
    function cm = coefficient_matrix(nn, mm, krtp)
      % Generate the coefficient matrix for the specified VSWF and point.
      %
      % Usage
      %   cm = NearfieldPm.coefficient_matrix(nn, mm, krtp)
      %   If inputs are cell arrays, output is a cell array of matrices.
      %
      % Parameters
      %   - nn (M-numeric|cell) -- Vector of VSWF n-indices or cell array
      %     with vectors for each coefficient matrix to generate.
      %
      %   - mm (M-numeric|cell) -- Vector of VSWF m-indices or cell array.
      %
      %   - krtp (3xN numeric|cell) -- Locations to evacuate fields.
      %     Spherical coordinates.  Radius is scaled by medium wave-number.

      % Much work to be done
      if iscell(nn) || iscell(mm) || iscell(krtp)

        % Convert any non-cells to cell
        if ~iscell(nn), nn = {nn}; end;
        if ~iscell(mm), mm = {mm}; end;
        if ~iscell(krtp), krtp = {krtp}; end;

        % Check lengths
        Ncells = max([numel(nn), numel(mm), numel(krtp)]);
        assert(numel(nn) == 1 || numel(nn) == Ncells, ...
            'nn must have same length as other cells or length 1');
        assert(numel(mm) == 1 || numel(mm) == Ncells, ...
            'nn must have same length as other cells or length 1');
        assert(numel(krtp) == 1 || numel(krtp) == Ncells, ...
            'nn must have same length as other cells or length 1');

        % Grow any short cells
        if numel(nn) == 1, nn = repmat(nn, 1, Ncells); end;
        if numel(mm) == 1, mm = repmat(mm, 1, Ncells); end;
        if numel(krtp) == 1, krtp = repmat(krtp, 1, Ncells); end;

        % Call ourselves for each cell
        cm = cellfun(@ott.beam.vswf.NearfieldPm.coefficient_matrix, ...
            nn, mm, krtp, 'UniformOutput', false);
      end

      % Check size of inputs
      assert(isnumeric(krtp) && ismatrix(krtp) && size(krtp, 1) == 3, ...
          'krtp must be 3xN numeric matrix');
      assert(numel(nn) == numel(mm), 'size of nn and mm must match');
      assert(isnumeric(nn) && isvector(nn), 'nn must be numeric vector');
      assert(isnumeric(mm) && isvector(mm), 'nn must be numeric vector');

      cm = zeros(3*size(krtp, 2), numel(nn));
      for ii = 1:numel(nn)

        % Calculate field for specified mode (regular basis)
        [M, N] = ott.utils.vswfcart(nn(ii), mm(ii), krtp(1, :), ...
            krtp(2, :), krtp(3, :), 'regular');

        % Store the fields
        if rem(nn(n),2) == 0
          % Even n
          cm(:, ii) = M(:);
        else
          % Odd n
          cm(:, ii) = N(:);
        end

      end
    end
  end

  methods
    function bsc = NearfieldPm(varargin)
      % Construct a BSC instance using far-field point matching.
      %
      % Usage
      %   bsc = NearfieldPm(nn, mm, rtp, E_farfield)
      %
      %   bsc = NearfieldPm(nn, mm, coefficient_matrix, E_farfield)
      %
      %   bsc = NearfieldPm(coefficient_matrix, E_farfield)
      %   Assumes coefficient matrix is a full matrix (all nn/mm values).
      %
      %   bsc = NearfieldPm(...) Constructs an empty Bsc object.
      %
      % Parameters
      %   - nn, mm (N-numeric|H-cell) -- VSWF mode indices.  Either vectors
      %     or cell arrays of vectors for each BSC to match.
      %
      %   - E_farfield (numeric|cell) -- Far-field values to match.
      %     Must be in `[Etheta; Ephi]` format.  Can either be
      %     a MxL matrix, H cells with MxL matrices.  For each coefficient
      %     matrix, solves the point matching problem for L BSC vectors.
      %
      %   - rtp (2xM numeric|H-cell) -- Coordinates of E-field samples.
      %     Spherical coordinates (radius/theta/phi).
      %
      %   - coefficient_matrix (MxN numeric|H-cell) -- A coefficient matrix
      %     generated using `NearfieldPm.coefficient_matrix`.
      %
      % Unmatched arguments are passed to base class.

      bsc = bsc@ott.beam.vswf.Pointmatch(varargin{:});
    end
  end

  methods (Hidden)
    function cm = build_coefficients(beam, nn, mm, rtp)
      % Build the coefficient matrix, scales radial direction

      % Scale radial direction by wavenumber
      if iscell(rtp)
        for ii = 1:numel(rtp)
          rtp{ii}(1, :) = rtp{ii}(1, :) .* beam.wavenumber;
        end
      else
        rtp(1, :) = rtp(1, :) .* beam.wavenumber;
      end

      cm = beam.coefficient_matrix(nn, mm, rtp);
    end

    function [a, b] = unpack_coefficients(beam, fab, nn, mm)
      % Unpack near-field coefficients

      % Get maximum n
      maxn = max(cellfun(@(x) max(x), nn));
      numrows = ott.utils.combined_index(maxn, maxn);
      numcols = sum(cellfun(@(x) size(x, 2), fab));

      % Unpack coefficients
      a = sparse(numrows, numcols);
      b = sparse(numrows, numcols);
      offset = 0;
      for ii = 1:numel(fab)

        colidx = (1:size(fab{ii}, 2)) + offset;
        offset = max(colidx);

        for jj = 1:numel(nn{ii})

          rowidx = ott.utils.combined_index(nn{ii}(jj), mm{ii}(jj));

          % TODO: Should this be sign(mm) or something else (mm == 0?)
          if rem(nn{ii}(jj), 2) == 0
            a(rowidx, colidx) = fab{ii}(jj, :);
            b(rowidx, colidx) = fab{ii}(jj, :) * sign(mm{ii}(jj));
          else
            a(rowidx, colidx) = fab{ii}(jj, :) * sign(mm{ii}(jj));
            b(rowidx, colidx) = fab{ii}(jj, :);
          end
        end
      end
    end
  end
end
