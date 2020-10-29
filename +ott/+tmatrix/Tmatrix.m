classdef Tmatrix < matlab.mixin.Heterogeneous
% Class representing T-matrix of a scattering particle or lens.
% This class can either be instantiated directly or used as a base
% class for defining custom T-matrix types.
%
% This class is the base class for all other T-matrix objects, you
% should inherit from this class when defining your own T-matrix
% creation methods. This class doesn't inherit from ``double`` or ``single``,
% instead the internal array type can be set at creation allowing the
% use of different data types such as ``sparse`` or ``gpuArray``.
%
% Properties
%   - data        -- The T-matrix this class encapsulates
%   - type        -- Type of T-matrix (total, scattered or internal)
%   - Nmax        -- Size of the T-matrix data (number of multipoles)
%   - total       -- Total-field instance of the T-matrix
%   - scattered   -- Scattered-field instance of the T-matrix
%
% Methods
%   - Tmatrix     -- Class constructor
%   - issparse    -- Returns true if the internal data is sparse
%   - full        -- Convert internal data to full
%   - sparse      -- Convert internal data to sparse
%   - makeSparse  -- Make the data sparse (with additional options)
%   - setNmax     -- Resize data to desired Nmax
%   - shrinkNmax  -- Reduce Nmax while preserving column/row power
%   - gpuArray    -- Make T-matrix a gpuArray
%   - gather      -- Apply gather to T-matrix data
%   - setType     -- Set the T-matrix type property (doesn't change data)
%   - columnCheck -- Calculate and check T-matrix column power
%   - mergeCols   -- Merge columns of tmatrices
%
% Mathematical and matrix operations
%   - times       -- Scalar multiplication of data
%   - mtimes      -- Scalar and matrix multiplication of data
%   - rdivide     -- Scalar division of data
%   - mrdivide    -- Scalar division of data
%   - uminus      -- Negation of data
%   - minus       -- Subtraction of data
%   - plus        -- Addition of data
%   - real        -- Extract real part of T-matrix
%   - imag        -- Extract imaginary part of T-matrix
%   - diag        -- Extract the diagonal of the T-matrix
%
% Static methods
%   - FromShape   -- Take a guess at a suitable T-matrix method
%   - SmartCylinder -- Smart method selection for cylindrical particles
%
% Casts
%   - ott.bsc.Bsc -- Convert each T-matrix column to beam vector
%   - ott.tmatrix.Tmatrix -- Downcast T-matrix superclass to base class
%
% See also :meth:`Tmatrix` and :class:`ott.tmatrix.Mie`.

% Copyright 2020 Isaac Lenton (aka ilent2)
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    type          % Type of T-matrix (total, scattered or internal)
    data          % The matrix this class encapsulates
  end

  properties (Dependent)
    Nmax          % Current size of T-matrix
    total         % Total version of the T-matrix
    scattered     % Scattered version of the T-matrix
  end

  methods (Static)
    function tmatrix = FromShape(shape, relative_index)
      % Take a guess at a suitable T-matrix method.
      %
      % The T-matrix can only be calculated easily and accurately for
      % a few very specific cases; as such, this function defaults to
      % DDA for most particles which may result in very slow calculations
      % (or failures due to memory limitations).  The resulting T-matrices
      % may not be accurate and it is recommended to inspect the fields
      % and compare results with another method.
      %
      % This method does the following
      %   - spheres -- Uses :class:`Mie`.
      %   - spheroids -- Uses :class:`Smarties`
      %   - cylinders -- Uses :meth:`SmartCylinder`.  This method may also
      %     work well for other semi-elongated rotationally symmetry shapes.
      %   - rotationally symmetric -- Uses :class:`Ebcm`.
      %   - star shaped -- Uses :class:`Pointmatch`.
      %   - otherwise -- Uses :class:`Dda`.
      %
      % For many types of particles it would be better to use another
      % method (such as geometric optics for very large particles).
      % This method may change in future releases when other methods
      % are added or when limits of existing methods are further explored.
      %
      % Usage
      %   tmatrix = ott.tmatrix.Tmatrix.FromShape(shape, relative_index)
      %
      % Parameters
      %   - shape (ott.shape.Shape) -- Shape to generate T-matrix for.
      %     Shape dimensions should be in units of wavelength.
      %
      %   - relative_index (numeric) -- Relative refractive index.

      if isa(shape, 'ott.shape.Sphere')
        tmatrix = ott.tmatrix.Mie.FromShape(shape, relative_index);
        return;
      end

      if shape.zRotSymmetry == 0 && isa(shape, 'ott.shape.Ellipsoid')
        tmatrix = ott.tmatrix.Smarties.FromShape(shape, relative_index);
        return;
      end

      if isa(shape, 'ott.shape.Superellipsoid') ...
          && shape.isEllipsoid && shape.zRotSymmetry == 0
        tmatrix = ott.tmatrix.Smarties.FromShape(shape, relative_index);
        return;
      end

      if isa(shape, 'ott.shape.Cylinder')
        tmatrix = ott.tmatrix.Tmatrix.SmartCylinder(shape, relative_index);
        return;
      end

      if shape.zRotSymmetry == 0
        tmatrix = ott.tmatrix.Ebcm.FromShape(shape, relative_index);
        return;
      end

      if shape.starShaped
        tmatrix = ott.tmatrix.Pointmatch.FromShape(shape, relative_index);
        return;
      end

      tmatrix = ott.tmatrix.Dda.FromShape(shape, relative_index);
    end

    function tmatrix = SmartCylinder(shape, relative_index, varargin)
      % Constructs a T-matrix for a cylinder using smart method selection
      %
      % Either uses DDA, EBCM or Pointmatch depending on the cylinder's
      % aspect ratio and material.  Pointmatching is the preferred method,
      % otherwise defaults to EBCM or DDA if error tolerance is not met.
      %
      % Uses the results from
      %
      %   Qi et al.,
      %   Optics Letters Vol. 39, Issue 16, pp. 4827-4830 (2014)
      %   https://doi.org/10.1364/OL.39.004827
      %
      % Usage
      %   tmatrix = SmartCylinder(shape, relative_index, ...)
      %
      % Parameters
      %   - shape (ott.shape.Shape) -- Shape to generate T-matrix for.
      %     Shape dimensions should be in units of wavelength.
      %     If the shape is not a cylinder, attempts to cast to Cylinder.
      %
      %   - relative_index (numeric) -- Relative refractive index.
      %
      % Optional named arguments
      %   - tolerance (enum) -- Error tolerance, can either be
      %     'one' or 'ten' for approximately 1% and 10% contours
      %     from the paper. Default: ``'ten'``.

      p = inputParser;
      p.addParameter('tolerance', 'ten');
      p.parse(varargin{:});

      assert(isscalar(shape), 'shape must be a single shape');
      if ~isa(shape, 'ott.shape.Cylinder')
        shape = ott.shape.Cylinder(shape);
      end

      % EBCM 1% data
      ebcm1 = {};
      ebcm1.x = [73, 169, 198, 228, 261, 391, 586, 718, 718, ...
          657, 523, 457, 262, 73];
      ebcm1.y = [409, 406, 418, 423, 397, 412, 400, 375, 223, ...
          193, 195, 165, 204, 390];
      ebcm1.x = (ebcm1.x - 73) * 2.0 / (718 - 73);
      ebcm1.y = -(ebcm1.y - 438) * 6.0 / (438 - 9);

      % PM 1% data
      pm1 = {};
      pm1.x = [297, 355, 394, 718, 718, 591, 525, 391, 361, 297];
      pm1.y = [943, 933, 946, 894, 868, 846, 874, 864, 913, 913];
      pm1.x = (pm1.x - 73) * 2.0 / (718 - 73);
      pm1.y = -(pm1.y - 985) * 6.0 / (985 - 555);

      % EBCM 10% data
      ebcm10 = {};
      ebcm10.x = [73, 193, 718, 718, 525, 328, 229, 160, 73];
      ebcm10.y = [430, 426, 381, 37, 94, 177, 214, 274, 375];
      ebcm10.x = (ebcm10.x - 73) * 2.0 / (718 - 73);
      ebcm10.y = -(ebcm10.y - 438) * 6.0 / (438 - 9);

      % PM 10% data
      pm10 = {};
      pm10.x = [130, 160, 328, 397, 462, 522, 589, 718, 718, ...
          654, 589, 522, 328, 265, 130];
      pm10.y = [961, 970, 967, 951, 946, 946, 925, 912, 753, ...
          784, 798, 798, 865, 874, 948];
      pm10.x = (pm10.x - 73) * 2.0 / (718 - 73);
      pm10.y = -(pm10.y - 985) * 6.0 / (985 - 555);

      % Paper parameters were in physical units
      lambda = 1.064;
      diameter = lambda * shape.radius * 2;
      len = shape.height;

      switch p.Results.tolerance
        case 'ten'
          if inpolygon(diameter, len, pm10.x, pm10.y)
            method = 'pm';
          elseif inpolygon(diameter, len, ebcm10.x, ebcm10.y)
              method = 'ebcm';
          else
              method = 'dda';
          end
        case 'one'
          if inpolygon(diameter, len, pm1.x, pm1.y)
            method = 'pm';
          elseif inpolygon(diameter, len, ebcm1.x, ebcm1.y)
            method = 'ebcm';
          else
            method = 'dda';
          end
        otherwise
          error('tolerance must be ''one'' or ''ten''');
      end

      switch method
        case 'dda'
          tmatrix = ott.tmatrix.Dda.FromShape(shape, relative_index);
        case 'ebcm'
          tmatrix = ott.tmatrix.Ebcm.FromShape(shape, relative_index);
        case 'pm'
          tmatrix = ott.tmatrix.Pointmatch.FromShape(shape, relative_index);
        otherwise
          error('Internal error');
      end
    end
  end

  methods
    function tmatrix = Tmatrix(varargin)
      % Construct a new T-matrix object.
      %
      % Usage
      %   tmatrix = Tmatrix(...)
      %   New empty T-matrix.  Leaves the data uninitialised.
      %
      %   tmatrix = Tmatrix(data, ...)
      %   Initializes the data with the matrix `data`.
      %
      % Parameters
      %   - data (NxM numeric | cell) -- The T-matrix data.
      %     Typically a sparse or full matrix.  Data must be empty or
      %     valid T-matrix size.  If cell, describes T-matrix array
      %     and elements must be a cell array of NxM matrices.
      %
      % Optional named arguments
      %   - type (enum) -- Type of T-matrix.  Must be 'internal',
      %     'scattered' or 'total'.  Default: ``'scattered'``.
      %
      % Example
      %   The following example creates an identity T-matrix which
      %   represents a particle which doesn't scatter light::
      %
      %     data = eye(16);
      %     tmatrix = ott.scat.vswf.Tmatrix(data, 'type', 'total');

      p = inputParser;
      p.addOptional('data', []);
      p.addParameter('type', 'scattered');
      p.parse(varargin{:});

      if iscell(p.Results.data)
        % T-matrix array
        tmatrix = repmat(tmatrix, size(p.Results.data));
        for ii = 1:numel(tmatrix)
          tmatrix(ii).data = p.Results.data{ii};
          tmatrix(ii) = tmatrix(ii).setType(p.Results.type);
        end
      else
        % Single T-matrix
        tmatrix.data = p.Results.data;
        tmatrix = tmatrix.setType(p.Results.type);
      end
    end

    function beam = ott.bsc.Bsc(tmatrix)
      % Create beams for each column in T-matrix
      %
      % Usage
      %   bsc = ott.bsc.Bsc(tmatrix)
      %
      % Parameters
      %   - tmatrix (instance of ott.tmatrix.Tmatrix) -- The T-matrix.

      beam = ott.bsc.Bsc.empty();
      for ii = 1:numel(tmatrix)
        oa = tmatrix(ii).data(1:end/2, :);
        ob = tmatrix(ii).data(end/2+1:end, :);
        beam = [beam, ott.bsc.Bsc(oa, ob)]; %#ok<AGROW>
      end
    end

    function tmatrix = ott.tmatrix.Tmatrix(tmatrix)
      % Downcast T-matrix to base class
      %
      % Usage
      %   tmatrix = ott.tmatrix.Tmatrix(tmatrix)

      for ii = 1:numel(tmatrix)
        tmatrix(ii) = ott.tmatrix.Tmatrix(tmatrix(ii).data, ...
            'type', tmatrix(ii).type);
      end
    end

    %
    % Sparsity functions
    %

    function b = issparse(tmatrix)
      % Returns true if the data is sparse
      %
      % Usage
      %   b = issparse(tmatrix)

      b = issparse(tmatrix.data);
    end

    function tmatrix = full(tmatrix)
      % Convert the data to a full matrix
      %
      % Usage
      %   tmatrix = full(tmatrix)

      ott.utils.nargoutCheck(tmatrix, nargout);

      tmatrix.data = full(tmatrix.data);
    end

    function tmatrix = sparse(tmatrix)
      % Convert the data to a sparse matrix
      %
      % This function doesn't change the data.  For a method that removes
      % near-zeros elements, see :meth:`makeSparse`.
      %
      % Usage
      %   tmatrix = sparse(tmatrix)

      ott.utils.nargoutCheck(tmatrix, nargout);

      tmatrix.data = sparse(tmatrix.data);
    end

    function tmatrix = makeSparse(tmatrix, varargin)
      % Make the T-matrix data sparse by removing near-zero power elements
      %
      % Treats each column as a beam shape vector and applies
      % :meth:`ott.bsc.Bsc.makeSparse` to each column.
      %
      % Usage
      %   tmatrix = tmatrix.makeSparse(...)
      %
      % Optional named arguments
      %   - AbsTol (numeric) -- Absolute tolerance for removing elements.
      %     Default: ``[]``.
      %
      %   - RelTol (numeric) -- Relative tolerance for removing elements.
      %     Power is relative to power in each column.
      %     Default: ``1.0e-15``.
      %
      % If both AbsTol and RelTol are specified, only elements satisfying
      % both conditions are kept.

      ott.utils.nargoutCheck(tmatrix, nargout);

      p = inputParser;
      p.addParameter('AbsTol', [], @isnumeric);
      p.addParameter('RelTol', 1.0e-15, @isnumeric);
      p.parse(varargin{:});

      for ii = 1:numel(tmatrix)

        % Cast to beam vector
        beam = ott.bsc.Bsc(tmatrix(ii));

        % Make beam vector sparse
        beam = beam.makeSparse('AbsTol', p.Results.AbsTol, ...
            'RelTol', p.Results.RelTol);

        % Re-form T-matrix
        tmatrix(ii) = ott.tmatrix.Tmatrix(beam);

      end

    end

    %
    % gpuArray functions
    %

    function tmatrix = gpuArray(tmatrix)
      % Copies the tmatrix data to the GPU
      %
      % Usage
      %   tmatrix = gpuArray(tmatrix)

      ott.utils.nargoutCheck(tmatrix, nargout);

      tmatrix.data = gpuArray(tmatrix.data);
    end

    function tmatrix = gather(tmatrix)
      % Apply `gather` to data.
      %
      % If the data is a gpuArray, returns a copy
      % of the data in the local workspace with data transferred from the GPU.
      %
      % Usage
      %   tmatrix = gather(tmatrix)

      ott.utils.nargoutCheck(tmatrix, nargout);

      tmatrix.data = gather(tmatrix.data);
    end

    %
    % Nmax functions
    %

    function tmatrix = setNmax(tmatrix, nmax, varargin)
      % Resize the T-matrix, with additional options
      %
      % Usage
      %   tmarix = tmatrix.setNmax(nmax, ...)   or    tmatrix.Nmax = nmax
      %   Set the Nmax, a optional warning is issued if truncation occurs.
      %
      % Parameters
      %   - Nmax (1 | 2 numeric) -- Nmax for both dimensions or vector
      %     with Nmax for `[rows, cols]`.
      %
      % Optional named arguments
      %   - AbsTol (numeric) -- Absolute tolerance for removing elements.
      %     Default: ``[]``.
      %
      %   - RelTol (numeric) -- Relative tolerance for removing rows.
      %     Power is relative to power in each column.
      %     Default: ``1.0e-15``.
      %
      %   - ColTol (numeric) -- Absolute tolerance for removing columns.
      %     Default: ``1.0e-15``.
      %
      %   - powerloss (enum) -- Action to take when column power is lost.
      %     Can be one of 'ignore', 'warn' or 'error'.
      %     Default: ``'warn'``.

      ott.utils.nargoutCheck(tmatrix, nargout);

      p = inputParser;
      p.addParameter('AbsTol', [], @isnumeric);
      p.addParameter('ColTol', 1.0e-15, @isnumeric);
      p.addParameter('RelTol', 1.0e-15, @isnumeric);
      p.addParameter('powerloss', 'warn');
      p.parse(varargin{:});

      assert(numel(nmax) == 1 || numel(nmax) == 2, ...
          'nmax must be 1 or 2 element vector');
      assert(isnumeric(nmax) && all(nmax == floor(nmax)) && all(nmax >= 0), ...
          'nmax must be positive numeric integer or zero');

      % Convert the input to row/column sizes
      if length(nmax) == 2
        nmax1 = nmax(1);
        nmax2 = nmax(2);
      else
        nmax1 = nmax(1);
        nmax2 = nmax(1);
      end

      for ii = 1:numel(tmatrix)

        % Check if we have work to do
        if all([nmax1, nmax2] == tmatrix(ii).Nmax)
          continue;
        end

        total_orders1 = ott.utils.combined_index(nmax1, nmax1);
        total_orders2 = ott.utils.combined_index(nmax2, nmax2);

        midpoint1 = size(tmatrix(ii).data, 1)/2;
        midpoint2 = size(tmatrix(ii).data, 2)/2;

        % Convert the T-matrix to scattered if we are growing the size
        % The conversion back to total ensures the new elements have the
        % correct values.
        %
        % No change necessary for scattered or internal T-matrices.
        old_type = tmatrix(ii).type;
        if (total_orders1 > midpoint1 || total_orders2 > midpoint2) ...
            && strcmpi(old_type, 'total')
          tmatrix(ii) = tmatrix(ii).scattered;
        end

        % Split T-matrix into quadrants
        A11 = tmatrix(ii).data(1:midpoint1, 1:midpoint2);
        A12 = tmatrix(ii).data(1:midpoint1, (midpoint2+1):end);
        A21 = tmatrix(ii).data((midpoint1+1):end, 1:midpoint2);
        A22 = tmatrix(ii).data((midpoint1+1):end, (midpoint2+1):end);

        % Resize rows
        if total_orders1 > midpoint1
          A11(total_orders1, :) = 0;
          A12(total_orders1, :) = 0;
          A21(total_orders1, :) = 0;
          A22(total_orders1, :) = 0;
        else

          % Calculate power in modes to be removed
          Arm = [A11(total_orders1+1:end, :), A12(total_orders1+1:end, :)];
          Brm = [A21(total_orders1+1:end, :), A22(total_orders1+1:end, :)];
          pw = abs(Arm).^2 + abs(Brm).^2;
          pw0 = sum(abs(tmatrix(ii).data).^2, 1);

          % Check RelTol
          if ~isempty(p.Results.RelTol) ...
              && any(sum(pw, 1)./pw0 > p.Results.RelTol)
            tmatrix.setNmaxWarning(p.Results.powerloss, 'reltol', ...
                'Relative tolerance for column power not satisfied');
          end

          % Check AbsTol
          if ~isempty(p.Results.AbsTol) ...
              && any(pw(:) > p.Results.AbsTol)
            tmatrix.setNmaxWarning(p.Results.powerloss, 'abstol', ...
                'Absolute tolerance for truncation not satisfied');
          end

          A11 = A11(1:total_orders1, :);
          A12 = A12(1:total_orders1, :);
          A21 = A21(1:total_orders1, :);
          A22 = A22(1:total_orders1, :);
        end

        % Resize columns
        if total_orders2 > midpoint2
          A11(:, total_orders2) = 0;
          A12(:, total_orders2) = 0;
          A21(:, total_orders2) = 0;
          A22(:, total_orders2) = 0;
        else

          % Calculate power in modes to be removed
          Arm = [A11(:, total_orders2+1:end), A12(:, total_orders2+1:end)];
          Brm = [A21(:, total_orders2+1:end), A22(:, total_orders2+1:end)];
          pw = abs(Arm).^2 + abs(Brm).^2;

          % Check ColTol
          if ~isempty(p.Results.ColTol) ...
              && any(sum(pw, 1) > p.Results.ColTol)
            tmatrix.setNmaxWarning(p.Results.powerloss, 'coltol', ...
                'Column power exceeds ColTol');
          end

          % Check AbsTol
          if ~isempty(p.Results.AbsTol) ...
              && any(pw(:) > p.Results.AbsTol)
            tmatrix.setNmaxWarning(p.Results.powerloss, 'abstol', ...
                'Absolute tolerance for truncation not satisfied');
          end

          A11 = A11(:, 1:total_orders2);
          A12 = A12(:, 1:total_orders2);
          A21 = A21(:, 1:total_orders2);
          A22 = A22(:, 1:total_orders2);
        end

        % Re-form T-matrix from quadrants
        tmatrix(ii).data = [ A11 A12; A21 A22 ];

        % If we were originally total field, convert back
        if strcmpi(old_type, 'total')
          tmatrix(ii) = tmatrix(ii).total;
        end
      end
    end

    function tmatrix = shrinkNmax(tmatrix, varargin)
      % Shrink the size of the T-matrix while preserving power
      %
      % Converts to a scattered or internal T-matrix and then removes
      % columns with no significant power and rows by passing each
      % column to :meth:`ott.bsc.Bsc.shrinkNmax`.
      %
      % Usage
      %   tmatrix = tmatrix.shrinkNmax(...)
      %
      % Optional named arguments
      %   - AbsTol (numeric) -- Absolute tolerance for removing elements.
      %     Default: ``[]``.
      %
      %   - RelTol (numeric) -- Relative tolerance for removing rows.
      %     Power is relative to power in each column.
      %     Default: ``1.0e-15``.
      %
      %   - ColTol (numeric) -- Absolute tolerance for removing columns.
      %     Default: ``1.0e-15``.

      ott.utils.nargoutCheck(tmatrix, nargout);

      p = inputParser;
      p.addParameter('AbsTol', [], @isnumeric);
      p.addParameter('ColTol', 1.0e-15, @isnumeric);
      p.addParameter('RelTol', 1.0e-15, @isnumeric);
      p.parse(varargin{:});

      for ii = 1:numel(tmatrix)

        % Calculate power of all elements
        pw = abs(tmatrix(ii).data).^2;

        colNmax = 0;

        % Find column Nmax that satisfies AbsTol
        if ~isempty(p.Results.AbsTol)
          pwa = [pw(:, 1:end/2); pw(:, end/2+1:end)];
          last_idx = find(any(pwa > p.Results.AbsTol, 1), 1, 'last');
          if ~isempty(last_idx)
            colNmax = max(ott.utils.combined_index(last_idx) + 1, colNmax);
          end
        end

        % Find column Nmax that satisfies ColTol
        if ~isempty(p.Results.ColTol)
          pwb = sum(pw, 1);
          pwb = [pwb(:, 1:end/2); pwb(:, end/2+1:end)];
          last_idx = find(any(pwb > p.Results.ColTol, 1), 1, 'last');
          if ~isempty(last_idx)
            colNmax = max(ott.utils.combined_index(last_idx) + 1, colNmax);
          end
        end

        % Change column Nmax
        tmatrix(ii).Nmax(2) = colNmax;

        % Convert to Bsc and apply shrinkNmax to remaining columns
        bsc = ott.bsc.Bsc(tmatrix(ii).data);
        bsc = bsc.shrinkNmax('RelTol', p.Results.RelTol, ...
            'AbsTol', p.Results.AbsTol);
        tmatrix(ii) = ott.tmatrix.Tmatrix(bsc);
      end
    end

    %
    % Mathematical operations
    %

    function tmatrix = plus(a, b)
      % Apply plus operation on T-matrix data
      %
      % Usage
      %   tmatrix = tmatrix1 + tmatrix2;
      %
      %   tmatrix = other + tmatrix;
      %   tmatrix = tmatrix + other;

      if isa(a, 'ott.tmatrix.Tmatrix') && isa(b, 'ott.tmatrix.Tmatrix')

        % Ensure both T-matrices have same Nmax
        Nmax1 = a.Nmax;
        Nmax2 = b.Nmax;
        oNmax = max(Nmax1, Nmax2);
        a.Nmax = oNmax;
        b.Nmax = oNmax;

        % Add data
        tmatrix = ott.tmatrix.Tmatrix(a.data + b.data);

      elseif isa(a, 'ott.tmatrix.Tmatrix')
        tmatrix = ott.tmatrix.Tmatrix(a.data + b);
      else
        tmatrix = ott.tmatrix.Tmatrix(a + b.data);
      end
    end

    function tmatrix = uminus(tmatrix)
      % Unary minus of data
      %
      % Usage
      %   tmatrix = -tmatrix

      ott.utils.nargoutCheck(tmatrix, nargout);

      tmatrix.data = -tmatrix.data;
    end

    function tmatrix = minus(a, b)
      % Minus operation on tmatrix
      %
      % Usage
      %   tmatrix = tmatrix1 - tmatrix2
      %
      %   tmatrix = other - tmatrix
      %   tmatrix = tmatrix - other
      %
      % Note: Uses the uminus operation of the second argument.

      tmatrix = a + (-b);
    end

    function tmatrix = rdivide(tmatrix, o)
      % Scalar division of T-matrix data
      %
      % Usage
      %   tmatrix = tmatrix ./ scalar

      tmatrix = ott.tmatrix.Tmatrix(tmatrix.data ./ o);

    end

    function tmatrix = mrdivide(tmatrix, o)
      % Scalar division of T-matrix data
      %
      % Usage
      %   tmatrix = tmatrix / scalar

      tmatrix = ott.tmatrix.Tmatrix(tmatrix.data / o);

    end

    function S = mtimes(a, b)
      % Matrix and scalar multiplication of T-matrix data
      %
      % Usage
      %   beam = tmatrix * beam
      %   Calculate how a beam is scattered by the T-matrix, increase
      %   beam or T-matrix Nmax if required.
      %
      %   tmatrix = tmatrix1 * tmatrix2
      %   Multiply T-matrix data, increasing Nmax if required.
      %
      %   vector = vector * tmatrix
      %   vector = tmatrix * vector
      %
      %   tmatrix = scalar * tmatrix
      %   tmatrix = tmatrix * scalar

      if isa(a, 'ott.tmatrix.Tmatrix') && isa(b, 'ott.tmatrix.Tmatrix')

        % Ensure both T-matrices have compatible Nmax
        Nmax1 = a.Nmax(2);
        Nmax2 = b.Nmax(1);
        oNmax = max(Nmax1, Nmax2);
        a.Nmax(2) = oNmax;
        b.Nmax(1) = oNmax;

        % Multiply data
        S = ott.tmatrix.Tmatrix(a.data * b.data);

      elseif isa(a, 'ott.tmatrix.Tmatrix') && isa(b, 'ott.bsc.Bsc')

        % Ensure T-matrix and beam have compatible Nmax
        Nmax1 = a.Nmax(2);
        Nmax2 = b.Nmax;
        oNmax = max(Nmax1, Nmax2);
        a.Nmax(2) = oNmax;
        b.Nmax = oNmax;

        S = a.data * b;

      elseif isa(b, 'ott.tmatrix.Tmatrix') && isa(a, 'ott.bsc.Bsc')
        error('ott:tmatrix:Tmatrix:mtimes:beam_order', ...
            'Cannot multiply beam by T-matrix, check multiplication order');
      elseif isa(a, 'ott.tmatrix.Tmatrix') && isscalar(b)
        S = ott.tmatrix.Tmatrix(a.data * b);
      elseif isa(b, 'ott.tmatrix.Tmatrix') && isscalar(a)
        S = ott.tmatrix.Tmatrix(a * b.data);
      elseif isa(a, 'ott.tmatrix.Tmatrix')
        S = a.data * b;
      elseif isa(b, 'ott.tmatrix.Tmatrix')
        S = a * b.data;
      else
        error('ott:tmatrix:Tmatrix:mtimes:unknown_args', ...
            'Unable to multiply arguments');
      end
    end

    function tmatrix = times(a, b)
      % Element-wise multiplication of T-matrix data
      %
      % Usage
      %   tmatrix = tmatrix1 .* tmatrix2;
      %
      %   tmatrix = other .* tmatrix;
      %   tmatrix = tmatrix .* other;

      if isa(a, 'ott.tmatrix.Tmatrix') && isa(b, 'ott.tmatrix.Tmatrix')

        % Ensure both T-matrices have same Nmax
        Nmax1 = a.Nmax;
        Nmax2 = b.Nmax;
        oNmax = max(Nmax1, Nmax2);
        a.Nmax = oNmax;
        b.Nmax = oNmax;

        % Multiply data
        tmatrix = ott.tmatrix.Tmatrix(a.data .* b.data);

      elseif isa(a, 'ott.tmatrix.Tmatrix')
        tmatrix = ott.tmatrix.Tmatrix(a.data .* b);
      else
        tmatrix = ott.tmatrix.Tmatrix(a .* b.data);
      end
    end

    function tmatrix = real(tmatrix)
      % Extract real part of T-matrix
      %
      % Usage
      %   tmatrix = real(tmatrix);

      tmatrix.data = real(tmatrix.data);
    end

    function tmatrix = imag(tmatrix)
      % Extract imaginary part of T-matrix
      %
      % Usage
      %   tmatrix = imag(tmatrix);

      tmatrix.data = imag(tmatrix.data);
    end

    function tmatrix = abs(tmatrix)
      % Calculate absolute value of T-matrix data
      %
      % Usage
      %   tmatrix = abs(tmatrix);

      tmatrix.data = abs(tmatrix.data);
    end

    function varargout = diag(tmatrix, k)
      % Extract the T-matrix diagonal
      %
      % Usage
      %   d = diag(tmatrix)
      %   Returns the diagonal in vector format
      %
      %   d = diag(tmatrix, K)
      %   Extracts the K-th diagonal of the T-matrix.
      %
      %   [dA, dB] = diag(...)
      %   Returns the diagonal of the upper and lower parts separately.

      if nargin == 1
        k = 0;
      end

      dA = diag(tmatrix.data(1:end/2, 1:end/2), k);
      dB = diag(tmatrix.data(end/2+1:end, end/2+1:end), k);

      if nargout == 1
        varargout{1} = [dA; dB];
      else
        varargout{1} = dA;
        varargout{2} = dB;
      end
    end

    %
    % Other T-matrix methods
    %

    function tmatrix = setType(tmatrix, val)
      % Set the T-matrix type paramter (without raising a warning)
      %
      % Usage
      %   tmatrix = tmatrix.setType(val);

      % Check output arguments
      ott.utils.nargoutCheck(tmatrix, nargout);

      S = warning('off', 'ott:scat:vswf:Tmatrix:type_change');
      tmatrix.type = val;
      warning(S);
    end

    function varargout = columnCheck(tmatrix, varargin)
      % Check the power in each column
      %
      % For a non-absorbing total-field T-matrix, the power in each column
      % should add up to unity (power should be conserved).
      %
      % Usage
      %   tmatrix.columnCheck(...)
      %   Raises a warning if the power drops bellow a threshold.
      %
      %   column_power = tmatrix.columnCheck(...)
      %   Returns the power in each column.
      %
      % Optional named arguments
      %   - threshold (numeric) -- Threshold for power loss warnings/errors.
      %     Default: ``1.0e-3``.
      %
      %   - action (enum) -- Action to take if power lost/gained.
      %     Can be 'warn', 'error' or 'none'.  If no outputs present,
      %     the default behaviour is 'warn'.  Otherwise 'none'.

      default_action = 'none';
      if nargout == 0
        default_action = 'warn';
      end

      p = inputParser;
      p.addParameter('action', default_action);
      p.addParameter('threshold', 1.0e-3);
      p.parse(varargin{:});

      % Calculate total T-matrix column power
      column_power = sum(abs(tmatrix.total.data).^2, 1);

      % Raise warnings or errors
      if any(abs(1.0 - column_power) > p.Results.threshold)
        switch p.Results.action
          case 'warn'
            warning('ott:tmatrix:Tmatrix:columnCheck:action_warn', ...
                'T-matrix power may not be conserved');
          case 'error'
            error('ott:tmatrix:Tmatrix:columnCheck:action_error', ...
                'T-matrix power may not be conserved');
          case 'none'
            % Nothing to do
          otherwise
            error('ott:tmatrix:Tmatrix:columnCheck:unknown_action', ...
                'Unknown action parameter value');
        end
      end

      % Assign outputs
      if nargout ~= 0
        varargout{1} = column_power;
      end
    end

    function tmatrix = mergeCols(tmatrix, tmatrix2, ci)
      % Merge columns of two T-matrices
      %
      % Usage
      %   tmatrix = tmatrix.mergeCols(tmatrix2, ci)
      %   Keeps all cols of the first T-matrix except thouse replaced
      %   by tmatrix2 (specified by ci).
      %
      % Parameters
      %   - tmatrix1, tmatrix2 -- T-matrices to merge
      %
      %   - ci -- (N numeric) Combined indices of T-matrix columns
      %     to replace.  For each ``ci`` entry, keeps 2 columns from
      %     ``tmatrix2`` (i.e., both TE and TM modes).

      ott.utils.nargoutCheck(tmatrix, nargout);

      % Check sizes match
      oNmax = max(tmatrix.Nmax, tmatrix2.Nmax);
      tmatrix.Nmax = oNmax;
      tmatrix2.Nmax = oNmax;

      % Replace columns
      midpoint = size(tmatrix.data, 2)/2;
      tmatrix.data(:, [ci, ci+midpoint]) = ...
          tmatrix2.data(:, [ci, ci+midpoint]);
    end
  end

  methods (Static, Hidden)
    function setNmaxWarning(warn_action, warn_field, msg)
      % Helper function for setNmax warnings

      if strcmpi(warn_action, 'warn')
        warning(['ott:tmatrix:Tmatrix:setNmaxWarning:' warn_field], msg);
      elseif strcmpi(warn_action, 'error')
        error(['ott:tmatrix:Tmatrix:setNmaxWarning:' warn_field], msg);
      elseif strcmpi(warn_action, 'ignore')
        % Nothing to do
      else
        error('ott:tmatrix:Tmatrix:setNmaxWarning:unknown_action', ...
          'powerloss should be one of ignore, warn or error');
      end
    end
  end

  methods % Getters/setters
    function tmatrix = set.data(tmatrix, val)

      % Check data type and size
      assert(isnumeric(val) && ismatrix(val), ...
        'tmatrix.data must be numeric matrix');
      assert(isempty(val) || ...
          (all(mod(size(val), 2) == 0) && all(size(val)./2 >= 3) && ...
          all(sqrt(size(val)./2+1) == floor(sqrt(size(val)./2+1)))), ...
          'tmatrix.data dimensions must be empty, 6, 16, 30, 48, ...');

      tmatrix.data = val;
    end

    function tmatrix = set.type(tmatrix, val)

      % Check type
      assert(any(strcmpi(val, {'internal', 'total', 'scattered'})), ...
          'type must be ''internal'' ''total'' or ''scattered''');

      % Warn user they may be doing the wrong thing
      warning('ott:scat:vswf:Tmatrix:type_change', ...
        ['Changing the type property doesnt change the type', newline, ...
        'Use tmatrix.total or tmatrix.scattered instead']);

      tmatrix.type = val;
    end

    function tmatrix = get.total(tmatrix)
      % Convert to a total T-matrix (if possible)

      switch tmatrix.type
        case 'internal'
          error('Cannot convert from internal to total T-matrix');
        case 'scattered'
          tmatrix.data = 2.0*tmatrix.data + speye(size(tmatrix.data));
          tmatrix = tmatrix.setType('total');
        case 'total'
          % Nothing to do
        otherwise
          error('Internal error: T-matrix has invalid type');
      end
    end

    function tmatrix = get.scattered(tmatrix)
      % Convert to a scattered T-matrix (if possible)

      switch tmatrix.type
        case 'internal'
          error('Cannot convert from internal to total T-matrix');
        case 'scattered'
          % Nothing to do
        case 'total'
          tmatrix.data = 0.5*(tmatrix.data - speye(size(tmatrix.data)));
          tmatrix = tmatrix.setType('scattered');
        otherwise
          error('Internal error: T-matrix has invalid type');
      end
    end

    function nmax = get.Nmax(tmatrix)
      %get.Nmax calculate Nmax from the current T-matrix data
      nmax1 = ott.utils.combined_index(size(tmatrix.data, 1)/2);
      nmax2 = ott.utils.combined_index(size(tmatrix.data, 2)/2);

      % Support non-square T-matrices
      nmax = [nmax1 nmax2];
    end
    function tmatrix = set.Nmax(tmatrix, nmax)
      %set.Nmax resizes the T-matrix
      tmatrix = tmatrix.setNmax(nmax);
    end
  end
end

