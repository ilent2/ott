classdef FarfieldPm < ott.beam.vswf.Pointmatch
% Point-matching for VSWF in the far-field.
% Inherits from :class:`Pointmatch`.
%
% Performs point-matching in the far-field, that is, solves the system::
%
%   M \vec{ab} = \vec{E}
%
% where :math:`\vec{ab}` are the VSWF coefficients, :math:`\vec{E}` is
% the far-field and :math:`M` describes the far-field corresponding
% to each VSWF coefficient.  This works best if the system of equations
% is over-determined (i.e., more far-field points than VSWF coefficients).
%
% Static methods
%   - coefficient_matrix      -- Generate far-field coefficient matrix.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Static)
    function cm = coefficient_matrix(nn, mm, tp)
      % Generate the coefficient matrix for the specified VSWF and point.
      %
      % Usage
      %   cm = NearfieldPm.coefficient_matrix(nn, mm, tp)
      %   If inputs are cell arrays, output is a cell array of matrices.
      %
      % Parameters
      %   - nn (M-numeric|cell) -- Vector of VSWF n-indices or cell array
      %     with vectors for each coefficient matrix to generate.
      %
      %   - mm (M-numeric|cell) -- Vector of VSWF m-indices or cell array.
      %
      %   - tp (2xN numeric|cell) -- Locations to evacuate fields.
      %     Spherical coordinates for far-field point (theta/phi).

      % Much work to be done
      if iscell(nn) || iscell(mm) || iscell(tp)

        % Convert any non-cells to cell
        if ~iscell(nn), nn = {nn}; end;
        if ~iscell(mm), mm = {mm}; end;
        if ~iscell(tp), tp = {tp}; end;

        % Check lengths
        Ncells = max([numel(nn), numel(mm), numel(tp)]);
        assert(numel(nn) == 1 || numel(nn) == Ncells, ...
            'nn must have same length as other cells or length 1');
        assert(numel(mm) == 1 || numel(mm) == Ncells, ...
            'nn must have same length as other cells or length 1');
        assert(numel(tp) == 1 || numel(tp) == Ncells, ...
            'nn must have same length as other cells or length 1');

        % Grow any short cells
        if numel(nn) == 1, nn = repmat(nn, 1, Ncells); end;
        if numel(mm) == 1, mm = repmat(mm, 1, Ncells); end;
        if numel(tp) == 1, tp = repmat(tp, 1, Ncells); end;

        % Call ourselves for each cell
        cm = cellfun(@ott.beam.vswf.NearfieldPm.coefficient_matrix, ...
            nn, mm, tp, 'UniformOutput', false);
      end

      % Check size of inputs
      assert(isnumeric(tp) && ismatrix(tp) && size(tp, 1) == 2, ...
          'tp must be 2xN numeric matrix');
      assert(numel(nn) == numel(mm), 'size of nn and mm must match');
      assert(isnumeric(nn) && isvector(nn), 'nn must be numeric vector');
      assert(isnumeric(mm) && isvector(mm), 'nn must be numeric vector');

      cm = zeros(2*size(tp, 2), 2*numel(nn));

      unn = unique(nn);

      for ii = 1:numel(unn)

        % Find which nn and mm correspond to this unique-nn
        n = unn(ii);
        ci = find(nn == n);

        % Evaluate VSWF in far-field
        [~,dtY,dpY]= ott.utils.spharm(n,mm(ci),tp(1, :), tp(2, :));

        % Store coefficients
        cm(:, ci) = [dpY;-dtY] * 1i^(n+1)/sqrt(n*(n+1));
        cm(:, ci+numel(nn)) = [dtY;dpY]*1i^(n)/sqrt(n*(n+1));
      end
    end

    function bsc = empty(varargin)
      % Create an empty FarfieldPm beam collection
      %
      % Usage
      %   beam = FarfieldPm.empty()
      %
      % Additional arguments are passed to constructor.

      bsc = ott.beam.vswf.FarfieldPm([], [], [], [], varargin{:});
    end
  end

  methods
    function bsc = FarfieldPm(varargin)
      % Construct a BSC instance using far-field point matching.
      %
      % Usage
      %   bsc = FarfieldPm(nn, mm, tp, E_farfield)
      %
      %   bsc = FarfieldPm(nn, mm, coefficient_matrix, E_farfield)
      %
      %   bsc = FarfieldPm(coefficient_matrix, E_farfield)
      %   Assumes coefficient matrix is a full matrix (all nn/mm values).
      %
      %   bsc = FarfieldPm(...) Constructs an empty Bsc object.
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
      %   - tp (2xM numeric|H-cell) -- Coordinates of E-field samples.
      %     Spherical coordinates (theta/phi).
      %
      %   - coefficient_matrix (MxN numeric|H-cell) -- A coefficient matrix
      %     generated using `FarfieldPm.coefficient_matrix`.
      %
      % Unmatched arguments are passed to base class.

      bsc = bsc@ott.beam.vswf.Pointmatch(varargin{:});
    end
  end

  methods (Hidden)
    function cm = build_coefficients(beam, nn, mm, tp)
      % Build the coefficient matrix: no scaling, just call static method.
      cm = beam.coefficient_matrix(nn, mm, tp);
    end
  end
end
