classdef Tmatrix < ott.scat.utils.Particle ...
    & ott.scat.utils.BeamForce & matlab.mixin.Heterogeneous ...
    & ott.utils.RotationPositionProp
% Class representing the T-matrix of a scattering particle or lens.
% This class can either be instantiated directly or used as a base
% class for defining custom T-matrix types.
%
% This class is the base class for all other T-matrix object, you
% should inherit from this class when defining your own T-matrix
% creation methods. This class doesn't inherit from ``double`` or ``single``,
% instead the internal array type can be set at creation allowing the
% use of different data types such as ``sparse`` or ``gpuArray``.
%
% Properties
%   - data        -- The T-matrix this class encapsulates
%   - type        -- Type of T-matrix (total, scattered or internal)
%   - position    -- Position of the particle
%   - rotation    -- Rotation of the particle
%   - Nmax        -- Size of the T-matrix data (number of multipoles)
%   - total       -- Total-field instance of the T-matrix
%   - scattered   -- Scattered-field instance of the T-matrix
%
% Methods
%   - setType     -- Set the T-matrix type property (doesn't change data)
%   - columnCheck -- Calculate and check T-matrix column power
%   - real        -- Extract real part of T-matrix
%   - imag        -- Extract imaginary part of T-matrix
%   - issparse    -- Returns true if the internal data is sparse
%   - full        -- Convert internal data to full
%   - sparse      -- Convert internal data to sparse
%   - rotate*     -- Methods for rotating particle
%
% Static methods
%   - simple      -- Construct a simple particle T-matrix
%
% See also Tmatrix, simple, :class:`+ott.scat.vswf.Mie`.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  % TODO: Use position/rotation (should this be in the base class?)
  % TODO: Review need for simple method.  Should we have shape casts?
  %   Or perhaps we should have a T-matrix/Shape class which stores
  %   the shape data?  Having a shape would make visualisation easier.
  % TODO: Review `defaultMethod`
  % TODO: Remove parse_wavenumber, parser_k_medium, parser_k_particle
  % TODO: Update defaults to use SMARTIES and DDA and others?
  % TODO: Rotation/translatiosn between T-matrices (in mtimes)

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
    function method = defaultMethod(shape, varargin)
      % Determine the appropriate method for a particular shape
      % Returns one of 'mie', smarties', 'dda', 'ebcm', or 'pm'.
      %
      % DEFAULTMETHOD(shape) determine the default method to use for
      % a ott.shapes.Shape obejct.
      %
      % DEFAULTMETHOD(name, parameters) determine the default method
      % for a shape described by its name and parameters.
      %
      % Supported shape names [parameters]:
      %   'sphere'          Spherical (or layered sphere) [ radius ]
      %   'cylinder'        z-axis aligned cylinder [ radius height ]
      %   'ellipsoid'       Ellipsoid [ a b c]
      %   'superellipsoid'  Superellipsoid [ a b c e n ]
      %   'cone-tipped-cylinder'      [ radius height cone_height ]
      %   'cube'            Cube [ width ]
      %   'axisym'          Axis-symetric particle [ rho(:) z(:) ]

      p = inputParser;
      p.addOptional('parameters', []);
      p.addParameter('method_tol', []);

      % Things required for k_medium
      p.addParameter('k_medium', []);
      p.addParameter('wavelength_medium', []);
      p.addParameter('index_medium', []);
      p.addParameter('wavelength0', []);

      p.parse(varargin{:});

      % Parse k_medium
      k_medium = ott.scat.vswf.Tmatrix.parser_k_medium(p, 2.0*pi);

      % Get a shape object from the inputs
      if ischar(shape) && ~isempty(p.Results.parameters)
        shape = ott.shapes.Shape.simple(shape, p.Results.parameters);
      elseif ~isa(shape, 'ott.shapes.Shape') || ~isempty(p.Results.parameters)
        error('Must input either Shape object or string and parameters');
      end

      if isa(shape, 'ott.shapes.Sphere') ...
          || (isa(shape, 'ott.shapes.Ellipsoid') && shape.isSphere) ...
          || (isa(shape, 'ott.shapes.Superellipsoid') && shape.isSphere)
        method = 'mie';
      elseif isa(shape, 'ott.shapes.Ellipsoid') ...
          || (isa(shape, 'ott.shapes.Superellipsoid') && shape.isEllipsoid)
        % TODO: Where does SMARTIES fail?
        method = 'smarties';
      elseif isa(shape, 'ott.shapes.Superellipsoid')
        % TODO: Where does PM fail?
        method = 'pm';
      elseif isa(shape, 'ott.shapes.Cube') ...
          || isa(shape, 'ott.shapes.RectangularPrism')
        % TODO: Where does PM fail?
        method = 'pm';
      elseif isa(shape, 'ott.shapes.Cylinder') ...
          || isa(shape, 'ott.shapes.AxisymLerp')

        if isa(shape, 'ott.shapes.Cylinder')
          parameters = [ shape.radius, shape.height ];
        elseif isa(shape, 'ott.shapes.AxisymLerp')
          parameters = [ max(shape.rho), max(shape.z) - min(shape.z) ];
        end

        method = ott.scat.vswf.Tmatrix.cylinder_preferred_method(...
            parameters, k_medium, p.Results.method_tol);

        if strcmp(method, 'other')
          method = 'dda';
        end

      else
        error('ott:Tmatrix:simple:no_shape', 'Unsupported particle shape');
      end

    end

    function tmatrix = simple(shape, varargin)
      % Constructs a T-matrix for different simple particle shapes.
      % This method creates an instance of one of the other T-matrix
      % classes and is here only as a helper method.
      %
      % Usage
      %   SIMPLE(shape) constructs a new simple T-matrix for the given
      %   :class:`+ott.+shapes.Shape` object.
      %
      %   SIMPLE(name, parameters) constructs a new T-matrix for the
      %   shape described by the name and parameters.
      %
      % Supported shape names [parameters]
      %   - 'sphere'       -- Spherical (or layered sphere) [ radius ]
      %   - 'cylinder'     -- z-axis aligned cylinder [ radius height ]
      %   - 'ellipsoid'    -- Ellipsoid [ a b c]
      %   - 'superellipsoid' -- Superellipsoid [ a b c e n ]
      %   - 'cone-tipped-cylinder'  --  [ radius height cone_height ]
      %   - 'cube'         -- Cube [ width ]
      %   - 'axisym'       -- Axis-symetric particle [ rho(:) z(:) ]
      %
      % Optional named arguments
      %   - method (enum) -- Allows you to choose the preferred method
      %     to use for T-matrix calculation.  Supported methods are
      %     'mie', smarties', 'dda', 'ebcm', and 'pm'.
      %     Default: `''`.
      %
      %   - method_tol (numeric) -- Specifies the error tolerances,
      %     a number between (0, 1] to use for method selection.
      %     Smaller values correspond to more accurate methods.
      %     Default: `[]`.
      %
      % Example
      %   The following example creates a T-matrix for a cube with side
      %   length of 1 micron using a guess at the best available method.
      %   Illumination wavelength is 1064 nm, relative index 1.5/1.33::
      %
      %     tmatrix = ott.scat.vswf.Tmatrix.simple('cube', 1.0e-6, ...
      %       'wavelength0', 1064e-9, ...
      %       'index_medium', 1.33, 'index_particle', 1.5);

      % Parse inputs
      p = inputParser;
      p.KeepUnmatched = true;
      p.addOptional('parameters', []);
      p.addParameter('method', '');
      p.addParameter('method_tol', []);

      % Things required for k_medium
      p.addParameter('k_medium', []);
      p.addParameter('wavelength_medium', []);
      p.addParameter('index_medium', []);
      p.addParameter('wavelength0', []);

      p.parse(varargin{:});

      % Parse k_medium
      k_medium = ott.scat.vswf.Tmatrix.parser_k_medium(p, 2.0*pi);

      % Get a shape object from the inputs
      if ischar(shape) && ~isempty(p.Results.parameters)
        shape = ott.shapes.Shape.simple(shape, p.Results.parameters);
        varargin = varargin(2:end);
      elseif ~isa(shape, 'ott.shapes.Shape') || ~isempty(p.Results.parameters)
        error('Must input either Shape object or string and parameters');
      end

      % Call the appropriate class to do the work
      method = p.Results.method;
      if isempty(method)
        method = ott.scat.vswf.Tmatrix.defaultMethod(shape, ...
          'method_tol', p.Results.method_tol, ...
          'k_medium', k_medium);
      end
      switch method
        case 'mie'
          tmatrix = ott.scat.vswf.Mie.simple(shape, varargin{:});
        case 'smarties'
          tmatrix = ott.scat.vswf.Smarties.simple(shape, varargin{:});
        case 'pm'
          tmatrix = ott.scat.vswf.Pm.simple(shape, varargin{:});
        case 'dda'
          tmatrix = ott.scat.vswf.Dda.simple(shape, varargin{:});
        case 'ebcm'
          tmatrix = ott.scat.vswf.Ebcm.simple(shape, varargin{:});
        otherwise
          error('Unsupported method specified');
      end
    end

    function method = cylinder_preferred_method(parameters, k, tolerance)
      %CYLINDER_PREFERRED_METHOD selects a method for cylinder calculation
      % Uses the results of Qi et al., 2014 to choose either PM or
      % EBCM to find the T-matrix for a cylinder shaped particle.
      %
      % CYLINDER_PREFERRED_METHOD(parameters, k, tolerance) calculates
      %   the prefered method, either 'pm', 'ebcm' or 'other'.
      %
      %   parameters is a vector with [ particle_radius, height ].
      %
      %   k is the wavenumber in the medium.
      %
      %   tolerance specifies the error tolerance, Qi et al. report
      %   results for tolerances of 0.1 and 0.01.  If tolerance is [],
      %   defaults to 0.01 as the tolerance.

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

      % Conversion factor, paper uses 1064nm illumination
      k = k * 1.064 / 2.0 / pi;
      diameter = 2.0 * parameters(1) * k;
      len = parameters(2) * k;

      if isempty(tolerance)
        tolerance = 0.01;
      end

      if tolerance >= 0.1
        if inpolygon(diameter, len, pm10.x, pm10.y)
            method = 'pm';
        elseif inpolygon(diameter, len, ebcm10.x, ebcm10.y)
            method = 'ebcm';
        else
            method = 'other';
        end
      elseif tolerance >= 0.01
        if inpolygon(diameter, len, pm1.x, pm1.y)
          method = 'pm';
        elseif inpolygon(diameter, len, ebcm1.x, ebcm1.y)
          method = 'ebcm';
        else
          method = 'other';
        end
      else
        method = 'other';
      end
    end

    function [km, kp] = parser_wavenumber(p, default)
      % Parses both k_medium and k_particle, provides n_relative support
      %
      % default is the default value for k_medium;

      % Run the original parsers with default arguments
      km = ott.optics.vswf.tmatrix.Tmatrix.parser_k_medium(p, []);
      kp = ott.optics.vswf.tmatrix.Tmatrix.parser_k_particle(p, []);

      % If we don't yet have any information, set km from default
      if isempty(km) && isempty(kp)
        km = default;
      end

      % Support for index_relative
      if ~isempty(p.Results.index_relative)
        if isempty(km) && ~isempty(kp)
          km = kp ./ p.Results.index_relative;
        elseif ~isempty(km) && isempty(kp)
          kp = km .* p.Results.index_relative;
        else
          error('index_relative specified but both indices already known');
        end
      elseif isempty(kp)
        error('Unable to determine particle wavenumber from inputs');
      end
    end

    function k_medium = parser_k_medium(p, default)
      %PARSER_K_MEDIUM helper to get k_medium from a parser object

      if ~isempty(p.Results.k_medium)
        k_medium = p.Results.k_medium;
      elseif ~isempty(p.Results.wavelength_medium)
        k_medium = 2.0*pi/p.Results.wavelength_medium;
      elseif ~isempty(p.Results.index_medium)
        if isempty(p.Results.wavelength0)
          error('wavelength0 must be specified to use index_medium');
        end
        k_medium = p.Results.index_medium*2.0*pi/p.Results.wavelength0;
      elseif nargin == 2
        k_medium = default;
      else
        error('Unable to determine k_medium from inputs');
      end
    end

    function k_particle = parser_k_particle(p, default)
      %PARSER_K_PARTICLE helper to get k_particle from a parser object

      if ~isempty(p.Results.k_particle)
        k_particle = p.Results.k_particle;
      elseif ~isempty(p.Results.wavelength_particle)
        k_particle = 2.0*pi/p.Results.wavelength_particle;
      elseif ~isempty(p.Results.index_particle)
        if isempty(p.Results.wavelength0)
          error('wavelength0 must be specified to use index_particle');
        end
        k_particle = p.Results.index_particle ...
            * 2.0*pi/p.Results.wavelength0;
      elseif nargin == 2
        k_particle = default;
      else
        error('Unable to determine k_particle from inputs');
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
      %   - data (numeric) -- The T-matrix data.  Typically a sparse or
      %     full matrix.  Data must be empty or valid T-matrix size.
      %
      % Optional named arguments
      %   - type (enum) -- Type of T-matrix.  Must be 'internal',
      %     'scattered' or 'total'.  Default: ``'scattered'``.
      %
      %   - position (3x1 numeric) -- Position of particle.
      %     Default: ``[0;0;0]``.
      %
      %   - rotation (3x3 numeric) -- Rotation of particle.
      %     Default: ``eye(3)``.
      %
      % Example
      %   The following example creates an identity T-matrix which
      %   represents a particle which doesn't scatter light::
      %
      %     data = eye(16);
      %     tmatrix = ott.scat.vswf.Tmatrix(data, 'type', 'total');

      p = inputParser;
      p.addOptional('data', [], @isnumeric);
      p.addParameter('type', 'scattered');
      p.addParameter('position', [0;0;0]);
      p.addParameter('rotation', eye(3));
      p.parse(varargin{:});

      % Store properties
      tmatrix.data = p.Results.data;
      tmatrix = tmatrix.setType(p.Results.type);
      tmatrix.position = p.Results.position;
      tmatrix.rotation = p.Results.rotation;
    end

    function tmatrix = setNmax(tmatrix, nmax, varargin)
      % Resize the T-matrix, with additional options.
      %
      % Usage
      %   tmatrix = tmatrix.setNmax(nmax, ...)  or  tmatrix.Nmax = nmax
      %   Change the size of the T-matrix.
      %
      % Parameters
      %   - Nmax (numeric) -- Nmax should be a scalar or 2-element vector
      %     with row/column Nmax values.
      %
      % Optional named arguments
      %   - tolerance (numeric) -- Tolerance to use for power warnings.
      %     Default: ``1.0e-6``.
      %
      % 	- powerloss (enum) -- behaviour for when power loss is detected.
      %     Can be 'ignore', 'warn' or 'error'.  Default: ``'warn'``.

      % TODO: What about internal T-matrices?

      p = inputParser;
      p.addParameter('tolerance', 1.0e-6);
      p.addParameter('powerloss', 'warn');
      p.parse(varargin{:});

      % Convert the input to row/column sizes
      if length(nmax) == 2
        nmax1 = nmax(1);
        nmax2 = nmax(2);
      else
        nmax1 = nmax(1);
        nmax2 = nmax(1);
      end

      % Check if we need to do anything
      if all([nmax1, nmax2] == tmatrix.Nmax)
        return;
      end

      total_orders1 = ott.utils.combined_index(nmax1, nmax1);
      total_orders2 = ott.utils.combined_index(nmax2, nmax2);

      midpoint1 = size(tmatrix.data, 1)/2;
      midpoint2 = size(tmatrix.data, 2)/2;

      % Convert the T-matrix to scattered if we are growing the size
      % The conversion back to total ensures the new elements have the
      % correct values.
      old_type = tmatrix.type;
      if total_orders1 > midpoint1 || total_orders2 > midpoint2 ...
          && strcmpi(old_type, 'total')
        tmatrix = tmatrix.scattered;
      end

      % Split T-matrix into quadrants
      A11 = tmatrix.data(1:midpoint1, 1:midpoint2);
      A12 = tmatrix.data(1:midpoint1, (midpoint2+1):end);
      A21 = tmatrix.data((midpoint1+1):end, 1:midpoint2);
      A22 = tmatrix.data((midpoint1+1):end, (midpoint2+1):end);

      % Resize rows
      if total_orders1 > midpoint1

        if issparse(tmatrix)
          [row_index,col_index,a] = find(A11);
          A11 = sparse(row_index,col_index,a,total_orders1,midpoint2);

          [row_index,col_index,a] = find(A12);
          A12 = sparse(row_index,col_index,a,total_orders1,midpoint2);

          [row_index,col_index,a] = find(A21);
          A21 = sparse(row_index,col_index,a,total_orders1,midpoint2);

          [row_index,col_index,a] = find(A22);
          A22 = sparse(row_index,col_index,a,total_orders1,midpoint2);
        else
          A11 = padarray(A11, [total_orders1 - midpoint1, 0], 0, 'post');
          A12 = padarray(A12, [total_orders1 - midpoint1, 0], 0, 'post');
          A21 = padarray(A21, [total_orders1 - midpoint1, 0], 0, 'post');
          A22 = padarray(A22, [total_orders1 - midpoint1, 0], 0, 'post');
        end

      elseif total_orders1 < midpoint1

        A11 = A11(1:total_orders1, :);
        A12 = A12(1:total_orders1, :);
        A21 = A21(1:total_orders1, :);
        A22 = A22(1:total_orders1, :);

      end

      % Resize cols
      if total_orders2 > midpoint2

        if issparse(tmatrix)
          [row_index,col_index,a] = find(A11);
          A11 = sparse(row_index,col_index,a,total_orders1,total_orders2);

          [row_index,col_index,a] = find(A12);
          A12 = sparse(row_index,col_index,a,total_orders1,total_orders2);

          [row_index,col_index,a] = find(A21);
          A21 = sparse(row_index,col_index,a,total_orders1,total_orders2);

          [row_index,col_index,a] = find(A22);
          A22 = sparse(row_index,col_index,a,total_orders1,total_orders2);
        else
          A11 = padarray(A11, [0, total_orders2 - midpoint2], 0, 'post');
          A12 = padarray(A12, [0, total_orders2 - midpoint2], 0, 'post');
          A21 = padarray(A21, [0, total_orders2 - midpoint2], 0, 'post');
          A22 = padarray(A22, [0, total_orders2 - midpoint2], 0, 'post');
        end

      elseif total_orders2 < midpoint2

        A11 = A11(:, 1:total_orders2);
        A12 = A12(:, 1:total_orders2);
        A21 = A21(:, 1:total_orders2);
        A22 = A22(:, 1:total_orders2);

      end

      % Get existing T-matrix magnitude
      % TODO: This should be applied to tmatrix.scattered or 'internal'
      if total_orders1 < midpoint1 || total_orders2 < midpoint2
        magA = full(sum(sum(abs(tmatrix.scattered.data).^2)));
      end

      % Recombined T-matrix from quadrants
      tmatrix.data = [ A11 A12; A21 A22 ];

      % Check for change in T-matrix magnitude
      if total_orders1 < midpoint1 || total_orders2 < midpoint2

        % TODO: This should be applied to tmatrix.scattered or 'internal'
        magB = full(sum(sum(abs(tmatrix.scattered.data).^2)));
        apparent_error = abs( magA - magB )/magA;

        if apparent_error > p.Results.tolerance
          if strcmpi(p.Results.powerloss, 'warn')
            warning('ott:Tmatrix:setNmax:truncation', ...
                ['Apparent error of ' num2str(apparent_error)]);
          elseif strcmpi(p.Results.powerloss, 'error')
            error('ott:Tmatrix:setNmax:truncation', ...
                ['Apparent error of ' num2str(apparent_error)]);
          elseif strcmpi(p.Results.powerloss, 'ignore')
            % Nothing to do
          else
            error('Unrecognized option for powerloss');
          end
        end
      end

      % If we were originally total field, convert back
      tmatrix = tmatrix.(old_type);
    end

    function S = mtimes(A, B)
      % Provide T-matrix multiplication overload
      %
      % Usage
      %   sbeam = tmatrix * ibeam
      %   Scatters an incident beam by a T-matrix.  The scattered beam
      %   is a instance of :class:`ott.beam.vswf.Scattered`.
      %
      %   S = tmatrix1 * tmatrix2
      %   Multiplies T-matrix data, increasing Nmax if required.
      %   `S` has the same type/properties as `tmatrix1`.
      %
      %   S = tmatrix * M    or    S = M * tmatrix
      %   Normal matrix multiplication.  `S` is numeric if
      %   `M` is non-scalar, otherwise `S` has the same type as `tmatrix`.

      if isa(A, 'ott.scat.vswf.Tmatrix') && isa(B, 'ott.scat.vswf.Tmatrix')

        % Ensure T-matrices have same size
        ab_Nmax = max(A.Nmax(2), B.Nmax(1));
        A.Nmax = [A.Nmax(1), ab_Nmax];
        B.Nmax = [ab_Nmax, B.Nmax(2)];

        S = A;
        S.data = S.data * B.data;

      elseif isa(B, 'ott.beam.abstract.Beam')

        % Calculate scattered beam
        S = A.scatter(B);

      elseif isa(A, 'ott.beam.abstract.Beam')
        error('Cannot multiply beam by T-matrix, check multiplication order');

      elseif isa(A, 'ott.scat.vswf.Tmatrix')
        if isscalar(B)
          S = A;
          S.data = S.data * B;
        else
          S = A.data * B;
        end

      elseif isa(B, 'ott.scat.vswf.Tmatrix')
        if isscalar(A)
          S = B;
          S.data = A * S.data;
        else
          S = A * B.data;
        end

      else
        error('Atleast one input must be T-matrix');
      end
    end

    function tmatrix = uminus(tmatrix)
      % Unary minus (negates elements of T-matrix)
      %
      % Usage
      %   tmatrix = -tmatrix;

      tmatrix.data = -tmatrix.data;
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
      %   - threshold (numeric) -- Power loss (or gain) threshold.
      %     Default: ``1.0e-3``.
      %
      %   - action (enum) -- Action to take if power lost/gained.
      %     Can be 'warn', 'error' or 'none'.  If no outputs present,
      %     the default behaviour is 'warn'.  Otherwise 'none'.

      p = inputParser;
      p.addParameter('threshold', 1.0e-3);
      if nargout == 0
        p.addParameter('action', 'warn');
      else
        p.addParameter('action', 'none');
      end
      p.parse(varargin{:});

      % Calculate total T-matrix column power
      column_power = sum(abs(tmatrix.total.data).^2, 1);

      % Raise warnings or errors
      if any(abs(1.0 - column_power) > p.Results.threshold)
        switch p.Results.action
          case 'warn'
            warning('T-matrix power may not be conserved');
          case 'error'
            error('T-matrix power may not be conserved');
          case 'none'
            % Nothing to do
          otherwise
            error('Unknown action parameter value');
        end
      end

      % Assign outputs
      if nargout ~= 0
        varargout{1} = column_power;
      end
    end

    function b = issparse(tmatrix)
      % Returns true if the T-matrix data is sparse
      %
      % Usage
      %   b = issparse(tmatrix)

      b = issparse(tmatrix.data);
    end

    function tmatrix = full(tmatrix)
      % Convert the T-matrix data to a full matrix
      %
      % Usage
      %   tmatrix = full(tmatrix)

      tmatrix.data = full(tmatrix.data);
    end

    function tmatrix = sparse(tmatrix)
      % Convert the T-matrix data to a sparse matrix
      %
      % Usage
      %   tmatrix = sparse(tmatrix)

      tmatrix.data = sparse(tmatrix.data);
    end

    function tmatrix = setType(tmatrix, val)
      % Set the T-matrix type paramter (without raising a warning)
      %
      % Usage
      %   tmatrix = tmatrix.setType(val);

      % Check output arguments
      tmatrix.nargoutCheck(nargout);

      S = warning('off', 'ott:scat:vswf:Tmatrix:type_change');
      tmatrix.type = val;
      warning(S);
    end

    function varargout = scatter(tmatrix, beam, varargin)
      % Calculate the beam scattered by a T-matrix.
      %
      % Usage
      %   scattered_beam = tmatrix.scatter(incident_beam, ...)
      %
      % Parameters and outputs
      %   - incident_beam (ott.beam.abstract.Beam) -- A beam object
      %     which is convertible to a VSWF beam.
      %
      %   - scattered_beam (ott.beam.vswf.Scattered) -- The scattered beam.
      %
      % Optional named arguments
      %   - position (3xN numeric) -- Translation applied to the beam
      %     before the beam is scattered by the particle.  Default: ``[]``.
      %     Only used with T-matrix input.
      %
      %   - rotation (3x3N numeric) -- Rotation applied to the beam,
      %     calculates the scattered beam and applies the inverse rotation,
      %     effectively rotating the particle.  Default: ``[]``.
      %     Only used with T-matrix input.
      %
      %   - array_type (enum) -- Array type for scattered beam with respect
      %     to rotation and position arguments.  Default: ``'array'``.
      %
      %   - store_tmatrix (logical) -- If false, the T-matrix is not
      %     stored in the scattered beam. Default: ``true``.
      %
      %   - store_beam (logical) -- If false, the incident beam is not
      %     in the scattered beam. Default: ``true``.
      %
      % If both position and rotation are present, the translation is
      % applied first, followed by the rotation.
      % The length of position and rotation must match or be scalar.

      % Only add additional documentation
      [varargout{1:nargout}] = scatter@ott.scat.utils.Particle(...
          tmatrix, beam, varargin{:});
    end
  end

  methods (Hidden)
    function sbeam = scatterInternal(tmatrix, beam, varargin)
      % Calculate the beam scattered by a T-matrix.
      %
      % Operates on a single T-matrix
      % TODO: Handling for multiple particles in scatter?

      % Setup input parser
      p = inputParser;
      p.addParameter('store_beam', true);
      p.addParameter('store_tmatrix', true);
      p.parse(varargin{:});

      % Ensure the beam is a Bsc
      % TODO: For some beams this could (should?) be done earlier
      if ~isa(beam, 'ott.beam.vswf.Bsc')

        % Calculate suggested Nmax from T-matrix and translation
        %
        % TODO: There are two optimal cases we could implement
        %   * For beams with an almost exact representation we should
        %     use the minimum Nmax.
        %   * For PlaneWave beams, we should pre-calculate the
        %     rotations and translatiosn and only convert at the end.
        %     So Nmax matches the T-matrix Nmax.
        maxPosition = max(vecnorm(beam.position));
        particleRadius = ott.utils.nmax2ka(tmatrix.Nmax(2));
        Nmax = ott.utils.ka2nmax(maxPosition*beam.wavenumber + particleRadius);

        beam = ott.beam.vswf.Bsc(beam, 'suggestedNmax', Nmax);
      end

      % Pre-combine coherent beams
      % TODO: For some beams this could be done earlier
      if strcmpi(beam.array_type, 'coherent')
        beam = sum(beam);
        beam.array_type = p.Results.array_type;
      end

      % Determine the maximum tmatrix.Nmax(2) and check type
      maxNmax1 = tmatrix.Nmax(1);
      maxNmax2 = tmatrix.Nmax(2);

      % If the T is scattered, we can save time by throwing away columns
      if strcmpi(tmatrix(1).type, 'scattered')
        maxNmax2 = min(maxNmax2, beam.Nmax);
        tmatrix.Nmax = [maxNmax1, maxNmax2];
      end

      % Apply translation to the beam
      if beam.position ~= [0;0;0]

        % Requires scattered beam, convert if needed
        if ~strcmpi(tmatrix.type, 'scattered')
          tmatrix = tmatrix.scattered;
          maxNmax2 = min(maxNmax2, beam.Nmax);
          tmatrix.Nmax = [maxNmax1, maxNmax2];
        end

        % Apply translation
        % We need Nmax+1 terms for the force calculation
        beam = beam.translateXyz(beam.position, 'Nmax', maxNmax2+1);
      end

      % Apply rotation to the beam
      rbeam = beam;
      if beam.rotation ~= eye(3)
        rbeam = rbeam.rotate(rbeam.rotation, 'Nmax', maxNmax1);
      end

      % Ensure the Nmax for the inner dimension matches
      if strcmpi(tmatrix.type, 'scattered')
        % T-matrix is already done
        rbeam = rbeam.setNmax(maxNmax2, 'powerloss', 'ignore');
      else
        tmatrix = tmatrix.setNmax([maxNmax1, rbeam.Nmax], ...
            'powerloss', 'ignore');

        if ~strcmpi(tmatrix.type, 'internal')
          warning('ott:Bsc:scatter', ...
              'It may be more optimal to use a scattered T-matrix');
        end
      end

      % Calculate the resulting beams
      sbeam = ott.beam.vswf.Scattered(tmatrix.data * rbeam);

      % Store the incident beam and T-matrices if requested
      if p.Results.store_beam
        sbeam.incident_beam = beam;
      end
      if p.Results.store_tmatrix
        sbeam.tmatrix = tmatrix;
      end

      % TODO: Set the sbeam position/rotation

      % Apply the inverse rotation
      % TODO: Should we just set the beam coordinate rotation?
      if beam.rotation ~= eye(3)
        % Faster to recompute rotation than to use wigner-D
        sbeam = sbeam.rotate(beam.rotation.');
      end

      % Assign a type to the resulting beam
      switch tmatrix.type
        case 'total'
          sbeam = sbeam.setType('total');
          sbeam.basis = 'regular';
        case 'scattered'
          sbeam = sbeam.setType('scattered');
          sbeam.basis = 'outgoing';
        case 'internal'
          sbeam = sbeam.setType('internal');
          sbeam.basis = 'regular';

          % Wavelength has changed, update it
          sbeam = sbeam.setWavenumber(...
              tmatrix.wavenumber_particle, 'medium');

        otherwise
          error('Unrecognized T-matrix type');
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
