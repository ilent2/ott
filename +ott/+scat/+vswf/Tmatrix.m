classdef Tmatrix < ott.scat.Particle ...
    & ott.scat.utils.BeamForce ...
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
% Distance units are all relative to the medium wavelength, usually this
% will be the wavelength of the medium surrounding the particle.
% Most methods have a `wavelength` parameter which can be used to
% specify the scale of the distance parameters in the input.
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
  % TODO: Rotation/translatiosn between T-matrices (in mtimes)
  % TODO: Add a shrink method and makeSparse (similar to Bsc)

  properties
    type          % Type of T-matrix (total, scattered or internal)
    data          % The matrix this class encapsulates
  end

  properties (Dependent)
    Nmax          % Current size of T-matrix
    total         % Total version of the T-matrix
    scattered     % Scattered version of the T-matrix
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
      if strcmpi(old_type, 'total')
        tmatrix = tmatrix.total;
      end
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

      % Ensure the beam is a Bsc
      if ~isa(beam, 'ott.beam.vswf.Bsc')

        % Calculate suggested Nmax from T-matrix and translation
        %
        % TODO: There are two optimal cases we could implement
        %   * For beams with an almost exact representation we should
        %     use the minimum Nmax.
        %   * For PlaneWave beams, we should pre-calculate the
        %     rotations and translatiosn and only convert at the end.
        %     So Nmax matches the T-matrix Nmax.
        %
        % TODO: Are there any cases when this should be converted later?
        maxPosition = max(vecnorm(beam.position));
        particleRadius = ott.utils.nmax2ka(tmatrix.Nmax(2));
        Nmax = ott.utils.ka2nmax(maxPosition*beam.wavenumber + particleRadius);

        beam = ott.beam.vswf.Bsc(beam, 'suggestedNmax', Nmax);
      end

      % Pre-combine coherent beams
      % TODO: Are there any cases when this should be done later?
      if strcmpi(beam.array_type, 'coherent')
        beam = sum(beam);
        beam.array_type = p.Results.array_type;
      end
      
      % Only add additional documentation
      [varargout{1:nargout}] = scatter@ott.scat.utils.Particle(...
          tmatrix, beam, varargin{:});
    end

    function shape = NmaxSphere(tmatrix, varargin)
      % Get a sphere representing the particle's Nmax
      %
      % Usage
      %   shape = tmatrix.NmaxSphere(...)
      %
      % Optional named parameters
      %   - wavelength (numeric) -- Wavelength of medium (default: 1.0)
      %
      %   - nmaxType (enum) -- Which Nmax to use: either 'rows',
      %     'cols' or 'both'.  Default: ``'cols'``.

      p = inputParser;
      p.addParameter('wavelength', 1.0);
      p.addParameter('nmaxType', 'cols');
      p.parse(varargin{:});

      switch p.Results.nmaxType
        case 'both'
          oNmax = tmatrix.Nmax;
        case 'cols'
          oNmax = tmatrix.Nmax(2);
        case 'rows'
          oNmax = tmatrix.Nmax(1);
        otherwise
          error('nmaxType must be both, cols or rows');
      end

      radius = ott.utils.nmax2ka(oNmax)./(2*pi).*p.Results.wavelength;
      shape = ott.shapes.Sphere(radius, ...
          'position', tmatrix.position*p.Results.wavelength, ...
          'rotation', tmatrix.rotation);
    end
  end

  methods (Hidden)
    function shape = getGeometry(tmatrix, wavelength, varargin)
      % Get a shape representing the T-matrix
      %
      % The default method returns a Nmax sphere

      shape = tmatrix.NmaxSphere('wavelength', wavelength, varargin{:});
    end

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


      % Determine the maximum tmatrix.Nmax(2) and check type
      maxNmax1 = tmatrix.Nmax(1);
      maxNmax2 = tmatrix.Nmax(2);

      % If the T is scattered, we can save time by throwing away columns
      if strcmpi(tmatrix(1).type, 'scattered')
        maxNmax2 = min(maxNmax2, beam.Nmax);
        tmatrix.Nmax = [maxNmax1, maxNmax2];
      end

      % Apply translation to the beam
      if any(beam.position ~= 0)

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
      if any(beam.rotation ~= eye(3))
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
      if any(beam.rotation ~= eye(3))
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
