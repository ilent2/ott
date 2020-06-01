classdef Pointmatch < ott.beam.vswf.Bsc
% Generates Bsc using point matching.
% Inherits from :class:`Bsc`.
%
% Hidden methods
%   - unpack_coefficients   -- Called by constructor to build a/b vectors
%   - build_coefficients    -- Called by constructor for coefficient matrix

% Based on bsc_pointmatch_focalplane and bsc_pointmatch_farfield
% from version 1 of the optical tweezers toolbox.

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  methods
    function beam = Pointmatch(varargin)
      % Construct a new Bsc using point matching.
      %
      % Usage
      %   bsc = Pointmatch(nn, mm, locations, Efield)
      %   Builds coefficient matrix for modes at specified locations
      %   and applies point matching with Efield.  Only supported
      %   when sub-classed.
      %
      %   bsc = Pointmatch(nn, mm, coefficient_matrix, Efield)
      %
      %   bsc = Pointmatch(coefficient_matrix, Efield)
      %   Assumes coefficient matrix is a full matrix (all nn/mm values)
      %
      %   bsc = Pointmatch(...)
      %
      % Parameters
      %   - nn, mm (N-numeric|H-cell) -- VSWF mode indices.  Either vectors
      %     or cell arrays of vectors for each BSC to match.
      %
      %   - Efield (numeric|cell) -- Field values to match.
      %     Format depends on coefficient matrix format.  Number of rows
      %     must match number of rows in coefficient matrix.
      %     Can either be a MxL matrix, H cells with MxL matrices.
      %     For each coefficient matrix, solves the point matching
      %     problem for L BSC vectors.
      %
      %   - coefficient_matrix (MxN numeric|H-cell) -- A coefficient matrix.
      %     The number of rows must match the number of rows in Efield,
      %     otherwise :methd:`beam.build_coefficients` is called (only
      %     implemented for sub-classes).
      %
      % Unmatched arguments are passed to base class.

      p = inputParser;
      p.addOptional('arg1', [], @isnumeric);
      p.addOptional('arg2', [], @isnumeric);
      p.addOptional('arg3', [], @isnumeric);
      p.addOptional('arg4', [], @isnumeric);
      p.KeepUnmatched = true;
      p.parse(varargin{:});

      % Construct base
      unmatched = ott.utils.unmatchedArgs(p);
      beam = beam@ott.beam.vswf.Bsc(unmatched{:});

      arg1 = p.Results.arg1;
      arg2 = p.Results.arg2;
      arg3 = p.Results.arg3;
      arg4 = p.Results.arg4;

      % Check number of arguments
      num_args = ~isempty(arg1) + ~isempty(arg2) + ~isempty(arg3) + ~isempty(arg4);
      assert(num_args == 0 || num_args == 2 || num_args == 4, ...
          'Must provide either 0, 2 or 4 positional arguments');

      % Get coefficient matrix and Efield
      if num_args == 4
        nn = arg1;
        mm = arg2;
        Efield = arg4;
        if ~iscell(nn), nn = {nn}; end;
        if ~iscell(mm), mm = {mm}; end;
        if iscell(arg3) && size(arg3{1}, 1) == size(Efield, 1) ...
          || ~iscell(arg3) && size(arg3, 1) == size(Efield, 1)
          coefficient_matrix = arg3;
        else
          coefficient_matrix = beam.build_coefficients(arg1, arg2, arg3);
        end
      elseif num_args == 2
        coefficient_matrix = arg1;
        if iscell(coefficient_matrix)
          lens = cellfun(@(x) size(x, 1)/2, coefficient_matrix);
          [nn, mm] = cellfun(@(x) ott.utils.combined_index(1:x), ...
              num2cell(lens), 'UniformOutput', false);
        else
          lens = size(x, 1)/2;
          [nn, mm] = ott.utils.combined_index(1:lens);
          nn = {nn};
          mm = {mm};
        end
        Efield = arg2;
      else
        coefficient_matrix = [];
        Efield = [];
        nn = {[]};
        mm = {[]};
      end

      % Solve system of equations
      if iscell(coefficient_matrix) && iscell(Efield)
        assert(numel(coefficient_matrix) == numel(Efield), ...
          'Efield and coefficient_matrix cell arrays must have same length');

        fab = cellfun(@(cm, e) cm \ e, coefficient_matrix, Efield, ...
            'UniformOutput', false);

      elseif iscell(coefficient_matrix)
        assert(isnumeric(Efield), 'Efield must be numeric');

        fab = cellfun(@(cm) cm \ Efield, coefficient_matrix, ...
            'UniformOutput', false);

      else
        assert(isnumeric(Efield) && isnumeric(coefficient_matrix), ...
            'Efield and coefficient_matrix must be numeric');

        fab = { coefficient_matrix \ Efield };

      end

      % Unpack result
      [a, b] = beam.unpack_coefficients(fab, nn, mm);
      beam = beam.setCoefficients(a, b);
    end
  end

  methods (Hidden)
    function cm = build_coefficients(beam, nn, mm, locs)
      % Build the coefficient matrix.
      %
      % This should be overloaded by your sub-class.  This implementation
      % simply raises and error when called.

      error(['Building coefficient matrix not supported', newline, ...
          'Use an overloaded class instead']);
    end

    function [a, b] = unpack_coefficients(beam, fab, nn, mm)
      % Unpack beam shape coefficients.
      %
      % This default implementation assumes there is one row for each
      % a/b coefficient in the form `[a; b]`.  Overload this method
      % if your `build_coefficients` function uses a different packing.
      
      % Check if any work to do
      allempty = all(cellfun(@(x) isempty(x), nn));
      if allempty
        a = [];
        b = [];
        return;
      end

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

        rowidx = ott.utils.combined_index(nn{ii}, mm{ii});

        a(rowidx, colidx) = fab{ii}(1:end/2, :);
        b(rowidx, colidx) = fab{ii}(1+end/2:end, :);

      end
    end
  end
end
