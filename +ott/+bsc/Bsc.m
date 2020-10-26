classdef Bsc < matlab.mixin.Heterogeneous ...
    & ott.utils.RotateHelper & ott.utils.TranslateHelper
% Class representing vector spherical wave function beam shape coefficients.
% Inherits from :class:`ott.utils.RotateHelper`
% and :class:`ott.utils.TranslateHelper`.
%
% Unlike the previous version of the toolbox, the `Bsc` class does not track
% the position/rotation of the beam.  This is to avoid suggesting that the
% beam actually has a well defined position/rotation.
% Other notable differences include the use of Hetrogeneous arrays,
% support for arbitrary beam data data-types, and moving user friendly
% features such as nice units to :class:`ott.beam.Beam`.
%
% `Bsc` rotations have units of radians and displacements have units of
% medium wavelength.
%
% Properties
%   - a         -- Beam shape coefficients `a` vector
%   - b         -- Beam shape coefficients `b` vector
%   - Nmax      -- (Dependent) Truncation number for VSWF coefficients
%   - power     -- (Dependent) Power of the beam shape coefficients
%
% Static methods
%   - FromDenseBeamVectors -- Construct beam from dense beam vectors.
%   - VisualisationData    -- Calculate visualisation data.
%
% Methods
%   - Bsc        -- Class constructor
%   - issparse   -- Check if the beam data is sparse
%   - full       -- Make the beam data full
%   - sparse     -- Make the beam data sparse
%   - makeSparse -- Make the beam data sparse (with additional options)
%   - setNmax    -- Resize beam data to desired Nmax
%   - shrinkNmax -- Reduce Nmax while preserving beam power
%   - gpuArray   -- Make beam a gpuArray
%   - gather     -- Apply gather to beam data
%   - rotate*    -- (Inherited) Functions for rotating the beam
%   - translate* -- (Inherited) Functions for translating the beam
%   - translateZ -- Translation along the theta=0 (z) axis
%   - getCoefficients -- Get a/b vectors with additional options
%   - setCoefficients -- Set a/b vectors with additional options
%
% Mathematical operations
%   - sum       -- Combine array of beams using summation
%   - times     -- Scalar multiplication of beam vectors
%   - mtimes    -- Scalar and matrix multiplication of beam vectors
%   - rdivide   -- Scalar division of beam vectors
%   - mrdivide  -- Scalar division of beam vectors
%   - uminus    -- Negation of beam vectors
%   - minus     -- Subtraction of beam vectors
%   - plus      -- Addition of beam vectors
%   - real      -- Extract real part of BSC data
%   - imag      -- Extract imag part of BSC data
%   - abs       -- Calculate absolute value of BSC data
%
% Field calculation methods
%   - efield     -- Calculate electric field around the origin
%   - hfield     -- Calculate magnetic field around the origin
%   - ehfield    -- Calculate electric and magnetic fields around the origin
%   - efarfield  -- Calculate electric fields in the far-field
%   - hfarfield  -- Calculate magnetic fields in the far-field
%   - ehfarfield -- Calculate electric and magnetic fields in the far-field
%   - eparaxial  -- Calculate electric fields in the paraxial far-field
%   - hparaxial  -- Calculate magnetic fields in the paraxial far-field
%   - ehparaxial -- Calculate electric and magnetic paraxial far-fields
%
% Force and torque related methods
%   - intensityMoment -- Calculate moment of beam intensity in the far-field
%   - force           -- Calculate the change in momentum between two beams
%   - torque          -- Calculate change in angular momentum between beams
%   - spin            -- Calculate change in spin momentum between beams
%   - forcetorque     -- Calculate the force and the torque between beams
%
% Field visualisation methods
%   - visNearfield      -- Generate a visualisation around the origin
%   - visFarfield       -- Generate a visualisation at the far-field
%   - visFarfieldSlice  -- Visualise the field on a angular slice
%   - visFarfieldSphere -- Visualise the filed on a sphere
%
% Casts
%   - ott.bsc.Bsc     -- Downcast BSC superclass to base class
%   - ott.tmatrix.Tmatrix -- Create T-matrix from beam array

% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    a          % Beam shape coefficients `a` vector
    b          % Beam shape coefficients `b` vector
  end

  properties (Dependent)
    Nmax       % Truncation number for VSWF coefficients
    power      % Power of the beam shape coefficients
  end

  methods (Static)
    function bsc = FromDenseBeamVectors(a, b, n, m)
      % Construct a new Bsc instance from dense beam vectors.
      %
      % Usage
      %   bsc = Bsc.FromDenseBeamVectors(a, b, n, m)
      %
      % Parameters
      %   - a,b (numeric) -- Dense beam vectors.  Can be Nx1 or NxM for
      %     M beams.  Each row corresponds to n/m indices.
      %
      %   - n,m (numeric) -- Mode indices for beam shape coefficients.
      %     Must have same number of rows as a/b.

      % Check inputs
      assert(all(size(a) == size(b)), 'size of a must match size of b');
      assert(numel(n) == numel(m), 'length of m and n must match');
      assert(numel(n) == size(a, 1), 'length of n must match rows of a');

      % Generate combine indices
      ci = ott.utils.combined_index(n, m);

      % Calculate total_order and number of beams
      nbeams = size(a, 2);
      Nmax = max(n);
      total_orders = ott.utils.combined_index(Nmax, Nmax);

      % Replicate ci for 2-D beam vectors (multiple beams)
      [ci, cinbeams] = ndgrid(ci, 1:nbeams);

      if ~isempty(ci)
        fa = sparse(ci, cinbeams, a, total_orders, nbeams);
        fb = sparse(ci, cinbeams, b, total_orders, nbeams);
      else
        fa = sparse(0, 0);
        fb = sparse(0, 0);
      end

      bsc = ott.bsc.Bsc(fa, fb);
    end

    function data = VisualisationData(field_type, field)
      % Helper to generate the visualisation data output.
      %
      % Usage
      %   VisualisationData(field_type, fieldvector)
      %   Takes a field_type string and a :class:`ott.utils.FieldVector`
      %   with the coordinates and field vectors.
      %
      % Parameters
      %   - fieldvector (3xN FieldVector) -- Field vector values
      %     and coordinate locations.
      %
      %   - field_type -- (enum) Type of field to calculate.
      %     Supported types include:
      %
      %     - irradiance  -- :math:`\sqrt{Ex^2 + Ey^2 + Ez^2}`
      %     - E2          -- :math:`Ex^2 + Ey^2 + Ez^2`
      %     - Sum(Abs(E)) -- :math:`|Ex| + |Ey| + |Ez|`
      %
      %     - Re(Er), Re(Et), Re(Ep), Re(Ex), Re(Ey), Re(Ez)
      %     - Abs(Er), Abs(Et), Abs(Ep), Abs(Ex), Abs(Ey), Abs(Ez)
      %     - Arg(Er), Arg(Et), Arg(Ep), Arg(Ex), Arg(Ey), Arg(Ez)

      assert(size(field, 1) == 3, 'fieldvector must be 3xN matrix');
      assert(isa(field, 'ott.utils.FieldVector'), ...
          'fieldvector must be a instance of FieldVector');

      % Generate the requested field
      switch field_type
        case 'irradiance'
          data = sqrt(sum(abs(field.vxyz).^2, 1));
        case 'E2'
          data = sum(abs(field.vxyz).^2, 1);
        case 'sum(Abs(E))'
          data = sum(abs(field.vxyz), 1);

        case 'Re(Er)'
          data = real(field.vrtp(1, :));
        case 'Re(Et)'
          data = real(field.vrtp(2, :));
        case 'Re(Ep)'
          data = real(field.vrtp(3, :));

        case 'Re(Ex)'
          data = real(field.vxyz(1, :));
        case 'Re(Ey)'
          data = real(field.vxyz(2, :));
        case 'Re(Ez)'
          data = real(field.vxyz(3, :));

        case 'Abs(Er)'
          data = abs(field.vrtp(1, :));
        case 'Abs(Et)'
          data = abs(field.vrtp(2, :));
        case 'Abs(Ep)'
          data = abs(field.vrtp(3, :));

        case 'Abs(Ex)'
          data = abs(field.vxyz(1, :));
        case 'Abs(Ey)'
          data = abs(field.vxyz(2, :));
        case 'Abs(Ez)'
          data = abs(field.vxyz(3, :));

        case 'Arg(Er)'
          data = angle(field.vrtp(1, :));
        case 'Arg(Et)'
          data = angle(field.vrtp(2, :));
        case 'Arg(Ep)'
          data = angle(field.vrtp(3, :));

        case 'Arg(Ex)'
          data = angle(field.vxyz(1, :));
        case 'Arg(Ey)'
          data = angle(field.vxyz(2, :));
        case 'Arg(Ez)'
          data = angle(field.vxyz(3, :));

        case 'Er'
          data = field.vrtp(1, :);
        case 'Et'
          data = field.vrtp(2, :);
        case 'Ep'
          data = field.vrtp(3, :);

        case 'Ex'
          data = field.vxyz(1, :);
        case 'Ey'
          data = field.vxyz(2, :);
        case 'Ez'
          data = field.vxyz(3, :);

        otherwise
          error('ott:bsc:Bsc:GetVisualisationData:unknown_field_type', ...
              'Unknown value for field_type.');
      end
    end
  end

  methods
    function beam = Bsc(varargin)
      % Construct a new beam object
      %
      % Usage
      %   beam = Bsc() Constructs an empty Bsc beam.
      %
      %   beam = Bsc(a, b, ...) constructs beam from a/b coefficients.
      %
      % Parameters
      %   - a,b (numeric) -- Vectors of VSWF coefficients

      p = inputParser;
      p.addOptional('a', zeros(0, 1), @isnumeric);
      p.addOptional('b', zeros(0, 1), @isnumeric);
      p.parse(varargin{:});

      if ~ismatrix(p.Results.a) || size(p.Results.a, 2) > 1 ...
          || ~ismatrix(p.Results.b) || size(p.Results.b, 2) > 1
        % Construct array of beams
        sza = size(p.Results.a);
        szb = size(p.Results.b);
        assert(numel(sza) == numel(szb) && all(sza == szb), ...
          'dimensions of ''a'' and ''b'' must match');
        beam = repmat(beam, 1, sza(2:end));
        for ii = 1:numel(beam)
          beam(ii).a = p.Results.a(:, ii);
          beam(ii).b = p.Results.b(:, ii);
        end
      else
        % Construct single beam
        beam.a = p.Results.a;
        beam.b = p.Results.b;
      end
    end

    function beam = ott.bsc.Bsc(other)
      % Cast to Bsc class
      %
      % Usage
      %   bsc = ott.bsc.Bsc(other)
      %
      % Parameters
      %   - other (instance of ott.bsc.Bsc) -- The class to cast to Bsc.

      beam = ott.bsc.Bsc(other.a, other.b);
    end

    function tmatrix = ott.tmatrix.Tmatrix(beam)
      % Create T-matrix from beam array
      %
      % Usage
      %   tmatrix = ott.tmatrix.Tmatrix(beam)
      %
      % Each beam in the beam array becomes a column of the T-matrix.
      % Uses :meth:`getCoefficients` to get the T-matrix data.

      data = beam.getCoefficients();
      tmatrix = ott.tmatrix.Tmatrix(data);
    end

    %
    % Field calculation functions
    %

    function [E, data] = efarfield(beam, rtp, varargin)
      % Calculate E far-field
      %
      % Usage
      %   [E, data] = beam.efarfield(rtp, ...)
      %
      % Parameters
      %   - rtp (3xN | 2xN numeric) -- Spherical coordinates for field
      %     calculation. Packaged [r; theta; phi] or [theta; phi].
      %
      % Optional named parameters
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.
      %
      %   - basis (enum) -- Basis for field calculation.  Must be either
      %     'incoming' or 'outgoing'.  Default: ``'incoming'``.

      p = inputParser;
      p.addParameter('data', ott.utils.VswfData(), ...
          @(x) isa(x, 'ott.utils.VswfData'));
      p.addParameter('basis', 'incoming');
      p.parse(varargin{:});

      % Ensure rtp size is 3xN
      [~, rtp] = ott.utils.rtpFarfield(rtp);

      ci = beam.getCombinedIndex();
      [n, ~] = ott.utils.combined_index(ci);
      Nn = 1./sqrt(n.*(n+1));
      [oa, ob] = beam.getCoefficients(ci);

      if strcmpi(p.Results.basis, 'incoming')
        ai = oa .* (1i).^(n+1) .* Nn;
        bi = ob .* (1i).^n .* Nn;
      elseif strcmpi(p.Results.basis, 'outgoing')
        ai = oa .* (-1i).^(n+1) .* Nn;
        bi = ob .* (-1i).^n .* Nn;
      else
        error('ott:bsc:Bsc:efarfield:unknown_basis', ...
            'Unknown value for basis parameter');
      end

      % Get or calculate spherical harmonics
      [~, Ytheta, Yphi, data] = p.Results.data.evaluateYtp(...
          ci, rtp(2, :), rtp(3, :));

      % Re-arrange a/b for multiplication
      ai = permute(full(ai), [1, 3, 2]);
      bi = permute(full(bi), [1, 3, 2]);

      % Calculate field components
      Etheta = sum(ai .* Yphi + bi .* Ytheta, 1);
      Ephi = sum(-ai .* Ytheta + bi .* Yphi, 1);

      % Package output
      Ertp = [zeros(size(Etheta)); Etheta; Ephi];
      E = ott.utils.FieldVector(rtp, Ertp, 'spherical');
    end

    function varargout = hfarfield(beam, rtp, varargin)
      % Calculate H far-field (uses :meth:`efarfield`)
      %
      % Usage
      %   [H, data] = beam.hfarfield(rtp, ...)
      %
      % See :meth:`efarfield` for further details.
      
      % Ensure rtp size is 3xN
      [~, rtp] = ott.utils.rtpFarfield(rtp);

      % Calculate E as normal
      [varargout{1:nargout}] = beam.efarfield(rtp, varargin{:});

      % Swap theta and phi components
      H = varargout{1};
      H = ott.utils.FieldVector(rtp, ...
          -1i .* H.vrtp([1, 3, 2], :), 'spherical');
      varargout{1} = H;
    end

    function [E, H, data] = ehfarfield(beam, rtp, varargin)
      % Calculate E and H far-fields
      %
      % Usage
      %   [E, H, data] = beam.ehfarfield(rtp, ...)
      %   E and H are of type :class:`ott.utils.FieldVector`.
      %
      % Parameters
      %   - rtp (3xN | 2xN numeric) -- Spherical coordinates for field
      %     calculation. Packaged [r; theta; phi] or [theta; phi].
      %
      % Optional named parameters
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.
      
      % Ensure rtp size is 3xN
      [~, rtp] = ott.utils.rtpFarfield(rtp);

      % Calculate E as normal
      [E, data] = beam.efarfield(rtp, varargin{:});

      % Swap theta and phi components
      H = ott.utils.FieldVector(rtp, ...
          -1i .* E.vrtp([1, 3, 2], :), 'spherical');
    end

    function varargout = eparaxial(beam, xy, varargin)
      % Calculate E paraxial far-field
      %
      % Usage
      %   [E, data] = beam.eparaxial(xy, ...)
      %   E is of type :class:`ott.utils.FieldVector`.
      %
      % Parameters
      %   - xy (2xN numeric) -- Paraxial coordinates for field
      %     calculation. Packaged [x; y].
      %
      % Optional named parameters
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.
      %
      %   - mapping (enum) -- Mapping for paraxial projection.
      %     See :func:`ott.utils.paraxial2rtp` for details.
      %     Default: ``'sin'``.
      %
      %   - direction (enum | 2 numeric) -- Mapping direction.
      %     See :func:`ott.utils.paraxial2rtp` for details.
      %     Default: ``'pos'``.
      %
      %   - basis (enum) -- Basis for field calculation.  Must be either
      %     'incoming' or 'outgoing'.  Default: ``'incoming'``.

      p = inputParser;
      p.addParameter('data', ott.utils.VswfData(), ...
          @(x) isa(x, 'ott.utils.VswfData'));
      p.addParameter('mapping', 'sin');
      p.addParameter('basis', 'incoming');
      p.addParameter('direction', 'pos');
      p.parse(varargin{:});

      rtp = ott.utils.paraxial2rtp(xy, ...
          p.Results.mapping, p.Results.direction);
      [varargout{1:nargout}] = beam.efarfield(rtp, ...
        'data', p.Results.data, 'basis', p.Results.basis);
    end

    function varargout = hparaxial(beam, xy, varargin)
      % Calculate H paraxial far-field
      %
      % Usage
      %   [H, data] = beam.hparaxial(xy, ...)
      %   H is of type :class:`ott.utils.FieldVector`.
      %
      % Parameters
      %   - xy (2xN numeric) -- Paraxial coordinates for field
      %     calculation. Packaged [x; y].
      %
      % Optional named parameters
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.
      %
      %   - mapping (enum) -- Mapping for paraxial projection.
      %     See :func:`ott.utils.paraxial2rtp` for details.
      %     Default: ``'sin'``.
      %
      %   - direction (enum | 2 numeric) -- Mapping direction.
      %     See :func:`ott.utils.paraxial2rtp` for details.
      %     Default: ``'pos'``.
      %
      %   - basis (enum) -- Basis for field calculation.  Must be either
      %     'incoming' or 'outgoing'.  Default: ``'incoming'``.

      p = inputParser;
      p.addParameter('data', ott.utils.VswfData(), ...
          @(x) isa(x, 'ott.utils.VswfData'));
      p.addParameter('mapping', 'sin');
      p.addParameter('direction', 'pos');
      p.addParameter('basis', 'incoming');
      p.parse(varargin{:});

      rtp = ott.utils.paraxial2rtp(xy, ...
          p.Results.mapping, p.Results.direction);
      [varargout{1:nargout}] = beam.hfarfield(rtp, ...
        'data', p.Results.data, 'basis', p.Results.basis);
    end

    function varargout = ehparaxial(beam, xy, varargin)
      % Calculate E and H paraxial far-fields
      %
      % Usage
      %   [E, H, data] = beam.ehparaxial(xy, ...)
      %   E and H are of type :class:`ott.utils.FieldVector`.
      %
      % Parameters
      %   - xy (2xN numeric) -- Paraxial coordinates for field
      %     calculation. Packaged [x; y].
      %
      % Optional named parameters
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.
      %
      %   - mapping (enum) -- Mapping for paraxial projection.
      %     See :func:`ott.utils.paraxial2rtp` for details.
      %     Default: ``'sin'``.
      %
      %   - direction (enum | 2 numeric) -- Mapping direction.
      %     See :func:`ott.utils.paraxial2rtp` for details.
      %     Default: ``'pos'``.

      p = inputParser;
      p.addParameter('data', ott.utils.VswfData(), ...
          @(x) isa(x, 'ott.utils.VswfData'));
      p.addParameter('mapping', 'sin');
      p.addParameter('direction', 'pos');
      p.parse(varargin{:});

      rtp = ott.utils.paraxial2rtp(xy, ...
          p.Results.mapping, p.Results.direction);
      [varargout{1:nargout}] = beam.ehfarfield(rtp, 'data', p.Results.data);
    end

    function varargout = efield(beam, xyz, varargin)
      % Calculate E field from Cartesian coordinates.
      %
      % Usage
      %   [E, data] = beam.efield(xyz, ...)
      %   E is of type :class:`ott.utils.FieldVector`.
      %
      % Parameters
      %   - xyz (3xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of wavelength.  Packaged [x; y; z].
      %
      % Optional named parameters
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.

      rtp = ott.utils.xyz2rtp(xyz);
      [varargout{1:nargout}] = beam.efieldRtp(rtp, varargin{:});
    end

    function varargout = hfield(beam, xyz, varargin)
      % Calculate H field from Cartesian coordinates.
      %
      % Usage
      %   [H, data] = beam.hfield(xyz, ...)
      %   H is of type :class:`ott.utils.FieldVector`.
      %
      % Parameters
      %   - xyz (3xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of wavelength.  Packaged [x; y; z].
      %
      % Optional named parameters
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.

      rtp = ott.utils.xyz2rtp(xyz);
      [varargout{1:nargout}] = beam.hfieldRtp(rtp, varargin{:});
    end

    function varargout = ehfield(beam, xyz, varargin)
      % Calculate E and H fields from Cartesian coordinates.
      %
      % Usage
      %   [E, H, data] = beam.ehfield(xyz, ...)
      %   E and H are of type :class:`ott.utils.FieldVector`.
      %
      % Parameters
      %   - xyz (3xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of wavelength.  Packaged [x; y; z].
      %
      % Optional named parameters
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.

      rtp = ott.utils.xyz2rtp(xyz);
      [varargout{1:nargout}] = beam.ehfieldRtp(rtp, varargin{:});
    end

    function [E, data] = efieldRtp(beam, rtp, varargin)
      % Calculate E-field around beam focus
      %
      % Usage
      %   [E, data] = beam.efieldRtp(rtp, ...)
      %
      % Parameters
      %   - rtp (3xN numeric) -- Spherical coordinates for field calculation.
      %     Units of wavelength.  Packaged [r; theta; phi].
      %
      % Optional named parameters
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.
      %
      %   - basis (enum) -- Vector spherical wave function basis to
      %     visualise.  Must be one of 'incoming', 'outgoing' or 'regular'.
      %     Default: ``'regular'``.

      p = inputParser;
      p.addParameter('data', ott.utils.VswfData(), ...
          @(x) isa(x, 'ott.utils.VswfData'));
      p.addParameter('basis', 'regular');
      p.parse(varargin{:});

      assert(ismatrix(rtp) && isnumeric(rtp) && size(rtp, 1) == 3, ...
          'rtp must be 3xN numeric matrix');

      % Multiply r by 2*pi
      rtp(1, :) = rtp(1, :) * 2*pi;

      ci = beam.getCombinedIndex();
      [n, ~] = ott.utils.combined_index(ci);
      Nn = 1./sqrt(n.*(n+1));

      % Get or calculate harmonics/Bessel data
      data = p.Results.data;
      [Y, Ytheta, Yphi, data] = data.evaluateYtp(ci, rtp(2, :), rtp(3, :));
      [hn, dhn, data] = data.evaluateBessel(n, rtp(1, :), p.Results.basis);

      % Get beam coefficients
      [ai, bi] = beam.getCoefficients(ci);
      ai = permute(full(ai), [1, 3, 2]);
      bi = permute(full(bi), [1, 3, 2]);

      % Calculate field components
      Er = sum(Nn.*n.*(n+1)./rtp(1, :).*hn.*Y.*bi, 1);
      Etheta = sum(Nn .* (ai .* Yphi .* hn + bi .* Ytheta .* dhn), 1);
      Ephi = sum(Nn .* (-ai .* Ytheta .* hn + bi .* Yphi .* dhn), 1);

      % Package output
      Ertp = [Er; Etheta; Ephi];
      E = ott.utils.FieldVector(rtp, Ertp, 'spherical');
    end

    function [H, data] = hfieldRtp(beam, rtp, varargin)
      % Calculate H-field around beam focus
      %
      % Usage
      %   [E, data] = beam.hfieldRtp(rtp, ...)
      %
      % Parameters
      %   - rtp (3xN numeric) -- Spherical coordinates for field calculation.
      %     Units of wavelength.  Packaged [r; theta; phi].
      %
      % Optional named parameters
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.
      %
      %   - basis (enum) -- Vector spherical wave function basis to
      %     visualise.  Must be one of 'incoming', 'outgoing' or 'regular'.
      %     Default: ``'regular'``.

      p = inputParser;
      p.addParameter('data', ott.utils.VswfData(), ...
          @(x) isa(x, 'ott.utils.VswfData'));
      p.addParameter('basis', 'regular');
      p.parse(varargin{:});

      assert(ismatrix(rtp) && isnumeric(rtp) && size(rtp, 1) == 3, ...
          'rtp must be 3xN numeric matrix');

      % Multiply r by 2*pi
      rtp(1, :) = rtp(1, :) * 2*pi;

      ci = beam.getCombinedIndex();
      [n, ~] = ott.utils.combined_index(ci);
      Nn = 1./sqrt(n.*(n+1));

      % Get or calculate harmonics/Bessel data
      data = p.Results.data;
      [Y, Ytheta, Yphi, data] = data.evaluateYtp(ci, rtp(2, :), rtp(3, :));
      [hn, dhn, data] = data.evaluateBessel(n, rtp(1, :), p.Results.basis);

      % Get beam coefficients
      [ai, bi] = beam.getCoefficients(ci);
      ai = permute(ai, [1, 3, 2]);
      bi = permute(bi, [1, 3, 2]);

      % Calculate field components
      Hr = sum(Nn.*n.*(n+1)./rtp(1, :).*hn.*Y.*ai, 1);
      Htheta = sum(Nn .* (bi .* Yphi .* hn + ai .* Ytheta .* dhn), 1);
      Hphi = sum(Nn .* (-bi .* Ytheta .* hn + ai .* Yphi .* dhn), 1);

      % Package output
      Hrtp = -1i*[Hr; Htheta; Hphi];
      H = ott.utils.FieldVector(rtp, Hrtp, 'spherical');
    end

    function [E, H, data] = ehfieldRtp(beam, rtp, varargin)
      % Calculate E and H fields in spherical coordinates.
      %
      % Usage
      %   [E, H, data] = beam.ehfieldRtp(rtp, ...)
      %   E and H are of type :class:`ott.utils.FieldVector`.
      %
      % Parameters
      %   - rtp (3xN numeric) -- Spherical coordinates for field calculation.
      %     Units of wavelength.  Packaged [r; theta; phi].
      %
      % Optional named parameters
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.

      p = inputParser;
      p.addParameter('data', ott.utils.VswfData(), ...
          @(x) isa(x, 'ott.utils.VswfData'));
      p.parse(varargin{:});

      data = p.Results.data;
      [E, data] = beam.efieldRtp(rtp, 'data', data);
      [H, data] = beam.hfieldRtp(rtp, 'data', data);
    end

    %
    % Visualisation functions
    %

    function varargout = visNearfield(beam, varargin)
      % Create a visualisation of the beam
      %
      % Usage
      %   beam.visualise(...) displays an image of the beam in the current
      %   figure window.
      %
      %   [im, data] = beam.visualise(...) returns a image of the beam.
      %   If the beam object is an array, returns an image for each beam.
      %   Also returns a :class:`ott.utils.VswfData` structure for fast
      %   repeated calculations.
      %
      % Optional named arguments
      %   - size (2 numeric) -- Number of rows and columns in image.
      %     Default: ``[80, 80]``.
      %
      %   - field (enum) -- Type of visualisation to generate, see
      %     :meth:`VisualisationData` for valid options.
      %     Default: ``'irradiance'``.
      %
      %   - axis (enum|cell) -- Describes the slice to visualise.
      %     Can either be 'x', 'y' or 'z' for a plane perpendicular to
      %     the x, y and z axes respectively; or a cell array with 2
      %     or 3 unit vectors (3x1) for x, y, [z] directions.
      %     Default: ``'z'``.
      %
      %   - offset (numeric) -- Plane offset along axis (default: 0.0)
      %
      %   - range (numeric|cell) -- Range of points to visualise.
      %     Can either be a cell array { x, y }, two scalars for
      %     range [-x, x], [-y, y] or 4 scalars [ x0, x1, y0, y1 ].
      %     Default: ``[1, 1].*beam.wavelength``.
      %
      %   - mask (function_handle) Describes regions to remove from the
      %     generated image.  Function should take one argument for the
      %     3xN field xyz field locations and return a logical array mask.
      %     Default: ``[]``.
      %
      %   - plot_axes (axes handle) -- Axes to place the visualisation in.
      %     If empty, no visualisation is generated.
      %     Default: ``gca()`` if ``nargout == 0`` otherwise ``[]``.
      %
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.
      %
      %   - basis (enum) -- Vector spherical wave function basis to
      %     visualise.  Must be one of 'incoming', 'outgoing' or 'regular'.
      %     Default: ``'regular'``.

      p = inputParser;
      p.addParameter('field', 'irradiance');
      p.addParameter('basis', 'regular');
      p.addParameter('size', []);
      p.addParameter('axis', 'z');
      p.addParameter('offset', 0.0);
      p.addParameter('range', [1, 1].*ott.utils.nmax2ka(beam.Nmax)./(2*pi));
      p.addParameter('mask', []);
      p.addParameter('plot_axes', []);
      p.addParameter('data', ott.utils.VswfData(), ...
          @(x) isa(x, 'ott.utils.VswfData'));
      p.parse(varargin{:});

      % Get range and size of data
      default_sz = [80, 80];
      [xrange, yrange, sz] = beam.visualiseGetRange(p.Results.range, ...
          p.Results.size, default_sz);

      % Generate grid of positions
      assert(isscalar(p.Results.offset) && isnumeric(p.Results.offset), ...
          'offset must be numeric scalar');
      [xx, yy, zz] = meshgrid(xrange, yrange, p.Results.offset);

      % Change the grid to the user requested orientation
      [xyz, labels] = beam.visualiseGetXyz(xx, yy, zz, p.Results.axis);

      % Calculate mask
      if isempty(p.Results.mask)
        mask = true(sz);
      else
        mask = p.Results.mask(xyz);
      end

      % Calculate which points are needed
      our_xyz = xyz(:, mask(:));

      % Calculate fields
      [E, data] = beam.efield(our_xyz, 'data', p.Results.data, ...
          'basis', p.Results.basis);

      % Unpack masked image (use nans for other values)
      Ef = nan(size(xyz));
      Ef(:, mask(:)) = E.vxyz;
      Ef = ott.utils.FieldVector(xyz, Ef, 'cartesian');

      % Generate visualisation data
      imout = beam.VisualisationData(p.Results.field, Ef);

      % Display the visualisation
      beam.visualiseShowPlot(...
          nargout, p.Results.plot_axes, imout, ...
          {xrange, yrange}, labels);

      % Assign output
      if nargout >= 1
        varargout{1} = imout;
        if nargout >= 2
          varargout{2} = data;
        end
      end
    end

    function varargout = visFarfield(beam, varargin)
      % Create a visualisation of the beam by projecting the far-field
      % onto a plane.
      %
      % Usage
      %   beam.visualiseFarfield(...) displays an image of the beam
      %   in the current axes.
      %
      %   im = beam.visualise(...) returns a image of the beam.
      %   If the beam object is an array, returns an image for each beam.
      %
      % Optional named arguments
      %   - size (2 numeric) -- Number of rows and columns in image.
      %     Default: ``[80, 80]``.
      %
      %   - field (enum) -- Type of visualisation to generate, see
      %     :meth:`VisualisationData` for valid options.
      %     Default: ``'irradiance'``.
      %
      %   - range (numeric|cell) -- Range of points to visualise.
      %     Can either be a cell array { x, y }, two scalars for
      %     range [-x, x], [-y, y] or 4 scalars [ x0, x1, y0, y1 ].
      %     Default: ``[1, 1]``.
      %
      %   - direction (enum|2 numeric|3x3 numeric) -- Direction to visualise.
      %     Either 'pos', 'neg' for the positive and negative z hemispheres,
      %     [roty, rotz] for the rotation angles (in degrees), or
      %     a 3x3 rotation matrix to apply to the coordinates.
      %     Default: ``'pos'``.
      %
      %   - mapping (enum) -- Mapping from theta-phi to far-field.
      %     Must be one of 'sin', 'tan' or 'theta'.
      %     Default: ``'sin'``.
      %
      %   - plot_axes (axes handle) -- Axes to place the visualisation in.
      %     If empty, no visualisation is generated.
      %     Default: ``gca()`` if ``nargout == 0`` otherwise ``[]``.
      %
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.
      %
      %   - basis (enum) -- Vector spherical wave function basis to
      %     visualise.  Must be one of 'incoming' or 'outgoing'.
      %     Default: ``'incoming'``.

      % Parse arguments
      p = inputParser;
      p.addParameter('field', 'irradiance');
      p.addParameter('basis', 'incoming');
      p.addParameter('size', []);
      p.addParameter('direction', 'pos');
      p.addParameter('range', [1, 1]);
      p.addParameter('mapping', 'sin');
      p.addParameter('plot_axes', []);
      p.addParameter('data', ott.utils.VswfData(), ...
          @(x) isa(x, 'ott.utils.VswfData'));
      p.parse(varargin{:});

      % Generate grid of coordinates
      default_sz = [80, 80];
      [xrange, yrange, ~] = beam.visualiseGetRange(p.Results.range, ...
          p.Results.size, default_sz);
      [xx, yy] = meshgrid(xrange, yrange);
      nxy = [xx(:), yy(:)].';

      % Calculate fields
      [E, data] = beam.eparaxial(nxy, 'data', p.Results.data, ...
          'basis', p.Results.basis, 'mapping', p.Results.mapping);

      % Generate visualisation data
      imout = beam.VisualisationData(p.Results.field, E);

      % Display the visualisation
      beam.visualiseShowPlot(...
          nargout, p.Results.plot_axes, imout, ...
          {xrange, yrange}, {'Direction 1', 'Direction 2'});

      % Assign output
      if nargout >= 1
        varargout{1} = imout;
        if nargout >= 2
          varargout{2} = data;
        end
      end
    end

    function varargout = visFarfieldSphere(beam, varargin)
      % Generate a spherical surface visualisation of the far-field
      %
      % Usage
      %   beam.visualiseFarfieldSphere(...)
      %   Generate a visualisation of the far-field in the current axes.
      %
      %   [I, XYZ, data] = beam.visualiseFarfieldSphere(...)
      %   Outputs the field value and three coordinate matrices that
      %   can be passed to ``surf(XYZ{1}, XYZ{2}, XYZ{3})``
      %
      % Optional named arguments
      %   - field (enum) -- Type of visualisation to generate, see
      %     :meth:`VisualisationData` for valid options.
      %     Default: ``'irradiance'``.
      %
      %   - normalise (logical) -- If the field value should be normalized.
      %     Default: ``false``.
      %
      %   - type (enum) -- Type of surface visualisation.
      %     Can be either 'sphere' or '3dpolar'.
      %     Default: ``'sphere'``.
      %
      %   - npts (numeric) -- Number of points for sphere surface.
      %     Passed to ``sphere(npts)``.  Default: ``100``.
      %
      %   - plot_axes (axes handle) -- Axes to place the visualisation in.
      %     If empty, no visualisation is generated.
      %     Default: ``gca()`` if ``nargout == 0`` otherwise ``[]``.
      %
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.
      %
      %   - basis (enum) -- Vector spherical wave function basis to
      %     visualise.  Must be one of 'incoming' or 'outgoing'.
      %     Default: ``'incoming'``.

      p = inputParser;
      p.addParameter('field', 'irradiance');
      p.addParameter('basis', 'incoming');
      p.addParameter('normalise', false);
      p.addParameter('type', 'sphere');
      p.addParameter('npts', 100);
      p.addParameter('plot_axes', []);
      p.addParameter('data', ott.utils.VswfData(), ...
          @(x) isa(x, 'ott.utils.VswfData'));
      p.parse(varargin{:});

      % Generate grid
      [X,Y,Z] = sphere(p.Results.npts);

      % Calculate fields
      rtp = ott.utils.xyz2rtp(X, Y, Z);
      [E, data] = beam.efarfield(rtp, 'data', p.Results.data, ...
          'basis', p.Results.basis);

      % Generate visualisation data
      imout = beam.VisualisationData(p.Results.field, E);

      % Generate visualisation
      if nargout == 0 || ~isempty(p.Results.plot_axes)

        % Get the default axes
        our_axes = p.Results.plot_axes;
        if isempty(our_axes)
          our_axes = gca();
        end

        % Check that we only have one visualisation to show
        if ~ismatrix(imout)
          warning('ott:beam:Beam:non_matrix_plot', ...
              'Multiple beams generated, only showing first');
        end

        I = imout(:, :, 1);
        if p.Results.normalise
          I = I ./ max(abs(I(:)));
        end

        switch p.Results.type
          case 'sphere'
            surf(our_axes, X, Y, Z, I, ...
                'facecolor','interp','edgecolor','none');
          case '3dpolar'
            surf(our_axes, abs(I).*X,abs(I).*Y,abs(I).*Z,I,...
              'facecolor','interp','edgecolor','none');
          otherwise
            error('Unknown visualisation type');
        end

        zlabel(our_axes, 'Z');
        xlabel(our_axes, 'X');
        ylabel(our_axes, 'Y');
        view(our_axes, 50, 20);
        axis(our_axes, 'equal');
      end

      % Assign output data
      if nargout >= 1
        varargout{1} = imout;
        if nargout >= 2
          varargout{2} = {X, Y, Z};
          if nargout >= 3
            varargout{3} = data;
          end
        end
      end
    end

    function varargout = visFarfieldSlice(beam, varargin)
      % Generate a 2-D slice through the far-field
      %
      % Usage
      %   beam.visualiseFarfieldSlice(phi, ...)
      %   Generates a 2-D slice at angle phi around the z-axis.
      %   If `phi` is not prsent, uses default of ``[0, pi/2]``.
      %   Plots into the current axes.
      %
      %   [im, theta, data] = beam.visualiseFarfieldSlice(...)
      %   Outputs the calculated values and corresponding angles.
      %
      % Optional named arguments
      %   - field (enum) -- Type of visualisation to generate, see
      %     :meth:`VisualisationData` for valid options.
      %     Default: ``'irradiance'``.
      %
      %   - normalise (logical) -- If the field value should be normalized.
      %     Default: ``false``.
      %
      %   - npts (numeric) -- Number of points for sphere surface.
      %     Passed to ``sphere(npts)``.  Default: ``100``.
      %
      %   - plot_axes (axes handle) -- Axes to place the visualisation in.
      %     If empty, no visualisation is generated.
      %     Default: ``gca()`` if ``nargout == 0`` otherwise ``[]``.
      %
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.
      %
      %   - basis (enum) -- Vector spherical wave function basis to
      %     visualise.  Must be one of 'incoming' or 'outgoing'.
      %     Default: ``'incoming'``.

      p = inputParser;
      p.addOptional('phi', [0, pi/2], @isnumeric);
      p.addParameter('field', 'irradiance');
      p.addParameter('basis', 'incoming');
      p.addParameter('normalise', false);
      p.addParameter('npts', 100);
      p.addParameter('plot_axes', []);
      p.addParameter('data', ott.utils.VswfData(), ...
          @(x) isa(x, 'ott.utils.VswfData'));
      p.parse(varargin{:});

      % Generate grid
      ptheta = linspace(0, 2*pi, p.Results.npts);
      [r, theta, phi] = ott.utils.matchsize(0, ptheta(:), p.Results.phi);
      rtp = [r(:), theta(:), phi(:)].';

      % Calculate fields
      [E, data] = beam.efarfield(rtp, 'data', p.Results.data, ...
          'basis', p.Results.basis);

      % Generate visualisation data
      imout = beam.VisualisationData(p.Results.field, E);

      % Generate visualisation
      if nargout == 0 || ~isempty(p.Results.plot_axes)

        % Get the default axes
        our_axes = p.Results.plot_axes;
        if isempty(our_axes)
          our_axes = gca();
        end

        % Check that we only have one visualisation to show
        if ~ismatrix(imout)
          warning('ott:beam:Beam:non_matrix_plot', ...
              'Multiple beams generated, only showing first');
        end

        I = imout(:, :, 1);
        if p.Results.normalise
          I = I ./ max(abs(I(:)));
        end

        % Generate plot
        polaraxes(our_axes);
        polarplot(ptheta, I);
      end

      % Assign outputs
      if nargout >= 1
        varargout{1} = imout;
        if nargout >= 2
          varargout{2} = theta;
          if nargout >= 3
            varargout{3} = data;
          end
        end
      end
    end

    %
    % Sparsity functions
    %

    function b = issparse(beam)
      % Returns true if the VSWF data is sparse
      %
      % Usage
      %   b = issparse(beam)

      b = false(size(beam));
      for ii = 1:numel(beam)
        b(ii) = issparse(beam(ii).a) & issparse(beam(ii).b);
      end
    end

    function beam = full(beam)
      % Convert the VSWF data to a full matrix
      %
      % Usage
      %   beam = full(beam)

      ott.utils.nargoutCheck(beam, nargout);

      for ii = 1:numel(beam)
        beam(ii).a = full(beam(ii).a);
        beam(ii).b = full(beam(ii).b);
      end
    end

    function beam = sparse(beam)
      % Convert the VSWF data to a sparse matrix
      %
      % This function doesn't change the data.  For a method that removes
      % near-zeros elements, see :meth:`makeSparse`.
      %
      % Usage
      %   beam = sparse(beam)

      ott.utils.nargoutCheck(beam, nargout);

      for ii = 1:numel(beam)
        beam(ii).a = sparse(beam(ii).a);
        beam(ii).b = sparse(beam(ii).b);
      end
    end

    function beam = makeSparse(beam, varargin)
      % Make the beam sparse by removing near-zero power elements
      %
      % Usage
      %   beam = beam.makeSparse(...)
      %
      % Optional named arguments
      %   - AbsTol (numeric) -- Absolute tolerance for removing elements.
      %     Default: ``[]``.
      %
      %   - RelTol (numeric) -- Relative tolerance for removing elements.
      %     Power is relative to power in each beam.
      %     Default: ``1.0e-15``.
      %
      % If both AbsTol and RelTol are specified, only elements satisfying
      % both conditions are kept.

      ott.utils.nargoutCheck(beam, nargout);
      
      p = inputParser;
      p.addParameter('AbsTol', [], @isnumeric);
      p.addParameter('RelTol', 1.0e-15, @isnumeric);
      p.parse(varargin{:});
      
      for ii = 1:numel(beam)

        [oa, ob] = beam(ii).getCoefficients();
        pw = abs(oa).^2 + abs(ob).^2;

        non_zero = true(size(pw));
        if ~isempty(p.Results.AbsTol)
          non_zero = non_zero & pw > p.Results.AbsTol;
        end
        if ~isempty(p.Results.RelTol)
          non_zero = non_zero & pw > p.Results.RelTol*sum(pw, 1);
        end

        oa(~non_zero) = 0;
        ob(~non_zero) = 0;

        beam(ii).a = sparse(oa);
        beam(ii).b = sparse(ob);
      end
    end

    %
    % gpuArray functions
    %

    function beam = gpuArray(beam)
      % Copies the beam shape coefficient data to the GPU
      %
      % Usage
      %   beam = gpuArray(beam)

      ott.utils.nargoutCheck(beam, nargout);

      beam.a = gpuArray(beam.a);
      beam.b = gpuArray(beam.b);

    end

    function beam = gather(beam)
      % Apply `gather` to beam shape coefficient data.
      %
      % If the beam shape coefficient data is a gpuArray, returns a copy
      % of the beam in the local workspace with data transferred from the GPU.
      %
      % Usage
      %   beam = gather(beam)

      ott.utils.nargoutCheck(beam, nargout);

      beam.a = gather(beam.a);
      beam.b = gather(beam.b);

    end

    %
    % Nmax functions
    %

    function beam = setNmax(beam, nmax, varargin)
      % Resize the beam, with additional options
      %
      % Usage
      %   beam = beam.setNmax(nmax, ...)   or    beam.Nmax = nmax
      %   Set the Nmax, a optional warning is issued if truncation occurs.
      %
      % Optional named arguments
      %   - AbsTol (numeric) -- Absolute tolerance for removing elements.
      %     Default: ``[]``.
      %
      %   - RelTol (numeric) -- Relative tolerance for removing elements.
      %     Power is relative to power in each beam.
      %     Default: ``1.0e-15``.
      %
      %   - powerloss (enum) -- Action to take when beam power is lost.
      %     Can be one of 'ignore', 'warn' or 'error'.
      %     Default: ``'warn'``.

      ott.utils.nargoutCheck(beam, nargout);

      p = inputParser;
      p.addParameter('AbsTol', [], @isnumeric);
      p.addParameter('RelTol', 1.0e-15, @isnumeric);
      p.addParameter('powerloss', 'warn');
      p.parse(varargin{:});

      total_orders = ott.utils.combined_index(nmax, nmax);
      if size(beam.a, 1) > total_orders

        % Check AbsTol
        if ~isempty(p.Results.AbsTol) ...
            && ~strcmpi(p.Results.powerloss, 'ignore')

          pw = abs(beam.a).^2 + abs(beam.b).^2;
          last_idx = find(pw > p.Results.AbsTol, 1, 'last');
          if last_idx < total_orders
            if strcmpi(p.Results.powerloss, 'warn')
              warning('ott:beam:vswf:Bsc:setNmax:truncation', ...
                  ['Apparent errors of ' num2str(aapparent_error) ...
                      ', ' num2str(bapparent_error) ]);
            elseif strcmpi(p.Results.powerloss, 'error')
              error('ott:beam:vswf:Bsc:setNmax:truncation', ...
                  ['Apparent errors of ' num2str(aapparent_error) ...
                      ', ' num2str(bapparent_error) ]);
            else
              error('ott:beam:vswf:Bsc:setNmax:truncation', ...
                'powerloss should be one of ignore, warn or error');
            end
          end
        end

        mag0 = beam.power;

        beam.a = beam.a(1:total_orders, :);
        beam.b = beam.b(1:total_orders, :);

        % Check RelTol
        if ~strcmpi(p.Results.powerloss, 'ignore')

          mag1 = beam.power;
          apparent_error = abs( mag1 - mag0 )/mag0;

          if apparent_error > p.Results.RelTol
            if strcmpi(p.Results.powerloss, 'warn')
              warning('ott:beam:vswf:Bsc:setNmax:truncation', ...
                  ['Apparent errors of ' num2str(aapparent_error) ...
                      ', ' num2str(bapparent_error) ]);
            elseif strcmpi(p.Results.powerloss, 'error')
              error('ott:beam:vswf:Bsc:setNmax:truncation', ...
                  ['Apparent errors of ' num2str(aapparent_error) ...
                      ', ' num2str(bapparent_error) ]);
            else
              error('ott:beam:vswf:Bsc:setNmax:truncation', ...
                'powerloss should be one of ignore, warn or error');
            end
          end
        end
      elseif size(beam.a, 1) < total_orders
        beam.a(total_orders, :) = 0;
        beam.b(total_orders, :) = 0;
      end
    end

    function nbeam = shrinkNmax(beam, varargin)
      % Reduces the size of the beam while preserving power
      %
      % Usage
      %   beam = beam.shrinkNmax(...)
      %
      % Optional named arguments
      %   - AbsTol (numeric) -- Absolute tolerance for removing elements.
      %     Default: ``[]``.
      %
      %   - RelTol (numeric) -- Relative tolerance for removing elements.
      %     Power is relative to power in each beam.
      %     Default: ``1.0e-15``.
      %
      % If both AbsTol and RelTol are specified, only elements satisfying

      ott.utils.nargoutCheck(beam, nargout);

      p = inputParser;
      p.addParameter('AbsTol', [], @isnumeric);
      p.addParameter('RelTol', 1.0e-15, @isnumeric);
      p.parse(varargin{:});

      [oa, ob] = beam.getCoefficients();

      pw = abs(oa).^2 + abs(ob).^2;

      % Find last element that satisfies AbsTol
      if ~isempty(p.Results.AbsTol)
        last_idx = find(pw > p.Results.AbsTol, 1, 'last');
        minNmax = ott.utils.combined_index(last_idx) + 1;
      else
        minNmax = 0;
      end

      pw0 = [beam.power];

      for ii = minNmax:beam.Nmax

        total_orders = ott.utils.combined_index(ii, ii);
        na = oa(1:total_orders, :);
        nb = ob(1:total_orders, :);
        nbeam = ott.bsc.Bsc(na, nb);

        pw1 = nbeam.power;
        apparent_error = abs( pw1 - pw0 )/pw0;
        if apparent_error < p.Results.RelTol
          break;
        end
      end
    end

    %
    % Coefficient access
    %

    function beam = setCoefficients(beam, a, b)
      % Set the `a` and `b` coefficients.
      %
      % The beam coefficients can also be set directly by accessing the
      % `a` and `b` properties of the beam.
      %
      % Usage
      %   beam = beam.setCoefficients(a, b)
      %
      %   beam = beam.setCoefficients(ab)
      %   Set coefficients from a [a; b] vector.
      %
      %   beam = beam.setCoefficients(other_beam)
      %   Get the beam coefficients from another beam

      ott.utils.nargoutCheck(beam, nargout);

      if nargin == 2
        if isa(a, 'ott.beam.vswf.Bsc')
          b = a.b;
          a = a.a;
        else
          assert(isnumeric(a) && ismatrix(a), 'ab must be numeric matrix');
          assert(mod(size(a, 1), 2) == 0, 'ab must be 2NxM in size');
          b = a(end/2+1:end, :);
          a = a(1:end/2, :);
        end
      end

      assert(all(size(a) == size(b)), 'size of a and b must match');
      assert(nargout == 1, 'Expected one output from function');

      assert(isempty(a) || ...
          (size(a, 1) >= 3 && ...
          sqrt(size(a, 1)+1) == floor(sqrt(size(a, 1)+1))), ...
        'number of multipole terms must be empty, 3, 8, 15, 24, ...');

      beam.a = a;
      beam.b = b;
    end

    function varargout = getCoefficients(beam, ci)
      % Gets the beam `a` and `b` coefficients
      %
      % The `a` and `b` coefficients can also be retrieved by accessing the
      % `a` and `b` properties of the beam directly.
      % This function is useful for getting specific beam shape coefficients
      % or packaging the coefficients in specific ways.
      %
      % Usage
      %   ab = beam.getCoefficients() gets the beam coefficients packed
      %   into a single vector, suitable for multiplying by a T-matrix.
      %
      %   [a, b] = beam.getCoefficients() get the coefficients in two
      %   beam vectors.
      %
      %   [...] = beam.getCoefficients(ci) behaves as above but only returns
      %   the requested beam cofficients a(ci) and b(ci) in a dense format.

      % Get size of each beam vector
      na = zeros(1, numel(beam));
      nb = zeros(1, numel(beam));
      for ii = 1:numel(beam)
        na(ii) = numel(beam(ii).a);
        nb(ii) = numel(beam(ii).b);
      end

      % If ci omitted, return all a and b
      if nargin == 1
        ci = 1:max([na(:); nb(:)]);
      end

      % Get data from all beams
      oa = zeros(numel(ci), numel(beam), 'like', beam(1).a);
      ob = zeros(numel(ci), numel(beam), 'like', beam(1).b);
      for ii = 1:numel(beam)
        oa(ci <= na(ii), ii) = beam(ii).a(ci <= na(ii));
        ob(ci <= nb(ii), ii) = beam(ii).b(ci <= nb(ii));
      end

      % Package output
      if nargout == 1
        varargout{1} = [oa; ob];
      else
        [varargout{1:2}] = deal(oa, ob);
      end
    end

    function ci = getCombinedIndex(beam, varargin)
      % Get combined index for non-zero rows in all beams
      %
      % Usage
      %   ci = beam.getCombinedIndex(...)
      %
      % Optional named arguments
      %   - full (logical) -- If true, gets the combined index for every
      %     row (including empty rows).  Default: ``false``.

      p = inputParser;
      p.addParameter('full', false);
      p.parse(varargin{:});

      if p.Results.full

        % Get size of each beam vector
        na = zeros(1, numel(beam));
        nb = zeros(1, numel(beam));
        for ii = 1:numel(beam)
          na(ii) = numel(beam(ii).a);
          nb(ii) = numel(beam(ii).b);
        end

        % Construct full list of cis
        ci = (1:max([na; nb])).';

      else

        % Find rows with non-zero elements (works best with sparse)
        [oa, ob] = beam.getCoefficients();
        [rowIdx, ~] = find(oa ~= 0 | ob ~= 0);
        ci = unique(rowIdx);
      end
    end

    function [n, m] = getModeIndices(beam, varargin)
      % Get mode indices (useful for creating a dense beam vector)
      %
      % Usage
      %   [n, m] = beam.getModeIndices()
      %
      % Optional named arguments
      %   - full (logical) -- If true, gets the indices for every
      %     row (including empty rows).  Default: ``false``.

      p = inputParser;
      p.addParameter('full', false);
      p.parse(varargin{:});

      ci = beam.getCombinedIndex('full', p.Results.full);
      [n, m] = ott.utils.combined_index(ci);
    end

    %
    % Mathematical operations
    %

    function beam = plus(beam1, beam2)
      % Add beam shape coefficients of two beams
      %
      % Usage
      %   beam = beam1 + beam2;

      assert(isa(beam1, 'ott.bsc.Bsc'), 'beam1 must be a ott.bsc.Bsc');
      assert(isa(beam2, 'ott.bsc.Bsc'), 'beam2 must be a ott.bsc.Bsc');

      ott.utils.nargoutCheck(beam1, nargout);

      ci1 = beam1.getCombinedIndex('full', true);
      ci2 = beam2.getCombinedIndex('full', true);
      if numel(ci1) > numel(ci2)
        ci = ci1;
      else
        ci = ci2;
      end

      [a1, b1] = beam1.getCoefficients(ci);
      [a2, b2] = beam2.getCoefficients(ci);
      beam = ott.bsc.Bsc(a1+a2, b1+b2);
    end

    function beam = uminus(beam)
      % Unary minus of beam shape coefficient data
      %
      % Usage
      %   beam = -beam

      ott.utils.nargoutCheck(beam, nargout);

      beam.a = -beam.a;
      beam.b = -beam.b;
    end

    function beam = minus(beam1, beam2)
      % Minus operation on two beams
      %
      % Usage
      %   beam = beam1 - beam2

      assert(isa(beam1, 'ott.bsc.Bsc'), 'beam1 must be a ott.bsc.Bsc');
      assert(isa(beam2, 'ott.bsc.Bsc'), 'beam2 must be a ott.bsc.Bsc');

      ott.utils.nargoutCheck(beam1, nargout);

      beam = beam1 + (-beam2);
    end

    function beam = rdivide(beam, other)
      % Right array divide
      %
      % Usage
      %   beam = beam ./ other
      %
      % `other` can be a scalar or matrix with a compatible size.
      % Applies the operation to the whole coefficients vector, i.e.::
      %
      %   [a;b] = [a;b] ./ other

      ott.utils.nargoutCheck(beam, nargout);
      beam = beam.setCoefficients(beam.getCoefficients ./ other);
    end

    function beam = mrdivide(beam, other)
      % Right matrix divide
      %
      % Usage
      %   beam = beam / other
      %
      % Applies the operation to the whole coefficients vector, i.e.::
      %
      %   [a;b] = [a;b] / other

      ott.utils.nargoutCheck(beam, nargout);
      beam = beam.setCoefficients(beam.getCoefficients / other);
    end

    function beam = times(beam, rv)
      % Provides element-wise multiplication of BSC coefficients.
      %
      % Usage
      %   beam = beam .* row_vec
      %   Multiplies each beam by the elements of `row_vec`.

      assert(isa(beam, 'ott.bsc.Bsc'), 'first argument must be a Bsc');
      assert(isnumeric(rv) && ismatrix(rv) && size(rv, 1) == 1 ...
          && size(rv, 2) == numel(beam), ...
          'second argument must be row vector matching N-beams');

      beam.a = beam.a .* rv;
      beam.b = beam.b .* rv;
    end

    function beam = mtimes(a,b)
      % Provides beam matrix and scalar multiplication
      %
      % Usage
      %   beam = scalar * beam   or   beam = beam * scalar
      %   Scalar multiplication of beam shape coefficients.
      %
      %   beam = matrix * beam
      %   Matrix multiplication.  Matrix can either have the same
      %   number of columns as the `a` and `b` beam shape coefficients
      %   or half as many rows.  In other words, the resulting beam
      %   is either `[M*a; M*b]` or `M*[a;b]`.

      if isa(a, 'ott.bsc.Bsc')
        [oa, ob] = a.getCoefficients();
        beam = ott.bsc.Bsc(oa * b, ob * b);
      else
        [oa, ob] = b.getCoefficients();
        if size(a, 2) == 2*size(oa, 1)
          ab = a * [oa; ob];
          beam = ott.bsc.Bsc(ab(1:size(ab, 1)/2, :), ...
              ab(1+size(ab, 1)/2:end, :));
        else
          beam = ott.bsc.Bsc(a * oa, a * ob);
        end
      end
    end

    function beam = sum(beamin, dim)
      % Sum beam coefficients (similar to the builtin sum function).
      %
      % Usage
      %   beam = sum(beam, dim)
      %
      % Parameters
      %   - dim -- (Optional) dimension to sum over.  Default is the
      %     first non-singleton dimension.

      % Handle default value for dimension
      if nargin < 2
        dim = find(size(beamin) > 1, 1);
      end

      % Select the first row in our dimension
      subs = [repmat({':'},1,dim-1), 1, ...
        repmat({':'},1,ndims(beamin)-dim)];
      S = struct('type', '()', 'subs', {subs});
      beam = subsref(beamin, S);

      % Add each beam
      for ii = 2:size(beamin, dim)
        subs = [repmat({':'},1,dim-1), ii, ...
          repmat({':'},1,ndims(beamin)-dim)];
        S = struct('type', '()', 'subs', {subs});
        beam = beam + subsref(beamin, S);
      end
    end

    function beam = real(beam)
      % Extract real part of beam shape coefficients
      %
      % Usage
      %   beam = real(beam)

      [oa, ob] = beam.getCoefficients();
      beam = ott.bsc.Bsc(real(oa), real(ob));
    end

    function beam = imag(beam)
      % Extract imaginary part of beam shape coefficients
      %
      % Usage
      %   beam = imag(beam)

      [oa, ob] = beam.getCoefficients();
      beam = ott.bsc.Bsc(imag(oa), imag(ob));
    end

    function beam = abs(beam)
      % Calculate absolute value of BSC data
      %
      % Usage
      %   beam = abs(beam)

      [oa, ob] = beam.getCoefficients();
      beam = ott.bsc.Bsc(abs(oa), abs(ob));
    end

    %
    % Force calculation methods
    %

    function [moment, ints, data] = intensityMoment(beam, varargin)
      % Calculate moment of the beam intensity in the far-field.
      %
      % By comparing the moment of the incident beam to the scattered
      % beam, this method should produce a similar estimate for the force
      % as :meth:`force`.
      %
      % Usage
      %   [moment, int, data] = beam.intensityMoment(...)
      %
      % Optional named arguments
      %   - theta_range (2 numeric) -- Range of angles to integrate over.
      %     Default: ``[0, pi]``.
      %
      %   - ntheta (numeric) -- Number of theta points.  (Default: 100)
      %   - nphi (numeric) -- Number of phi points.  (Default: 100)
      %
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.
      %
      %   - basis (enum) -- Vector spherical wave function basis to
      %     visualise.  Must be one of 'incoming' or 'outgoing'.
      %     Default: ``'incoming'``.

      p = inputParser;
      p.addParameter('basis', 'incoming');
      p.addParameter('theta_range', [0, pi]);
      p.addParameter('ntheta', 100);
      p.addParameter('nphi', 100);
      p.addParameter('data', ott.utils.VswfData(), ...
          @(x) isa(x, 'ott.utils.VswfData'));
      p.parse(varargin{:});

      % Setup grid
      [theta, phi] = ott.utils.angulargrid(p.Results.ntheta, p.Results.nphi);
      dtheta = diff(theta(1:2));
      dphi = diff(phi(1:2));

      % Truncate the theta range
      keep = theta > p.Results.theta_range(1) & theta < p.Results.theta_range(2);
      theta = theta(keep);
      phi = phi(keep);

      rtp = [ones(numel(theta), 1), theta(:), phi(:)].';

      % Calculate Cartesian coordinates
      % negate z, So integrals match sign convention used in :meth:`force`.
      uxyz = ott.utils.rtp2xyz(rtp);
      uxyz(3, :) = -uxyz(3, :);

      % Calculate field and E2
      [E, data] = beam.efarfield(rtp, 'data', p.Results.data, ...
          'basis', p.Results.basis);
      Eirr = beam.VisualisationData('E2', E);

      % Calculate intensity
      ints = sum(Eirr .* sin(theta.') .* dtheta .* dphi, 2);

      % Calculate moment in Cartesian coordinates
      Eirr_xyz = uxyz .* Eirr;
      moment = sum(Eirr_xyz .* sin(theta.') .* dtheta .* dphi, 2);
    end

    function [force, torque, spin] = forcetorque(ibeam, sbeam)
      % Calculate change in momentum between beams
      %
      % Usage
      %   [f, t, s] = ibeam.forcetorque(sbeam) calculates the force,
      %   torque and spin between the incident beam ``ibeam`` and
      %   scattered beam ``sbeam``.
      %   Outputs 3x[N...] matrix depending on the number and shape of beams.
      %
      % To convert the force to SI units, divide by the speed in the medium
      % (assuming the beam power is in SI units).
      %
      % To convert torque/spin to SI units, divide by the angular frequency
      % (assuming the beam power is in SI units).
      %
      % The scattered beam must be a total field beam (incoming-outgoing).
      %
      % This uses mathematical result of Farsund et al., 1996, in the form of
      % Chricton and Marsden, 2000, and our standard T-matrix notation S.T.
      % E_{inc}=sum_{nm}(aM+bN);

      % Dispatch to other methods to calculate quantities
      force = ibeam.force(sbeam);
      if nargout > 1
        torque = ibeam.torque(sbeam);
        if nargout > 2
          spin = ibeam.spin(sbeam);
        end
      end
    end

    function varargout = force(ibeam, sbeam)
      % Calculate change in linear momentum between beams.
      % For details on usage/arguments see :meth:`forcetorque`.
      %
      % Usage
      %   force = ibeam.force(sbeam)
      %   [fx, fy, fz] = ibeam.force(sbeam)

      % Get the abpq terms for the calculation
      [ai, bi, p, q, n, m, ...
        anp1, bnp1, pnp1, qnp1, ...
        amp1, bmp1, pmp1, qmp1, ...
        anp1mp1, bnp1mp1, pnp1mp1, qnp1mp1, ...
        anp1mm1, bnp1mm1, pnp1mm1, qnp1mm1] = ...
      ibeam.find_abpq_force_terms(sbeam);

      % Calculate the Z force
      Az=m./n./(n+1).*imag(-(ai).*conj(bi)+conj(q).*(p));
      Bz=1./(n+1).*sqrt(n.*(n-m+1).*(n+m+1).*(n+2)./(2*n+3)./(2*n+1)) ... %.*n
          .*imag(anp1.*conj(ai)+bnp1.*conj(bi)-(pnp1).*conj(p) ...
          -(qnp1).*conj(q));
      fz=2*sum(Az+Bz);

      % Calculate the XY force
      Axy=1i./n./(n+1).*sqrt((n-m).*(n+m+1)) ...
          .*(conj(pmp1).*q - conj(amp1).*bi - conj(qmp1).*p + conj(bmp1).*ai);
      Bxy=1i./(n+1).*sqrt(n.*(n+2))./sqrt((2*n+1).*(2*n+3)).* ... %sqrt(n.*)
          ( sqrt((n+m+1).*(n+m+2)) .* ( p.*conj(pnp1mp1) + q.* ...
          conj(qnp1mp1) -ai.*conj(anp1mp1) -bi.*conj(bnp1mp1)) + ...
          sqrt((n-m+1).*(n-m+2)) .* (pnp1mm1.*conj(p) + qnp1mm1.* ...
          conj(q) - anp1mm1.*conj(ai) - bnp1mm1.*conj(bi)) );

      fxy=sum(Axy+Bxy);
      fx=real(fxy);
      fy=imag(fxy);

      % Ensure things are full
      fx = full(fx);
      fy = full(fy);
      fz = full(fz);

      % Package output
      if nargout == 3
        varargout{1:3} = {fx, fy, fz};
      else
        varargout{1} = [fx(:) fy(:) fz(:)].';
      end
    end

    function varargout = torque(ibeam, sbeam)
      % Calculate change in angular momentum between beams
      % For details on usage/arguments see :meth:`forcetorque`.
      %
      % Usage
      %   torque = ibeam.torque(sbeam)
      %   [tx, ty, tz] = ibeam.torque(sbeam)

      % Get the abpq terms for the calculation
      [ai, bi, p, q, n, m, ~, ~, ~, ~, ...
        amp1, bmp1, pmp1, qmp1] = ...
      ibeam.find_abpq_force_terms(sbeam);

      tz=sum(m.*(ai.*conj(ai)+bi.*conj(bi)-p.*conj(p)-q.*conj(q)));

      txy=sum(sqrt((n-m).*(n+m+1)).*(ai.*conj(amp1)+...
        bi.*conj(bmp1)-p.*conj(pmp1)-q.*conj(qmp1)));
      tx=real(txy);
      ty=imag(txy);

      % Ensure things are full
      tx = full(tx);
      ty = full(ty);
      tz = full(tz);

      % Package output
      if nargout == 3
        varargout{1:3} = {tx, ty, tz};
      else
        varargout{1} = [tx(:) ty(:) tz(:)].';
      end
    end

    function varargout = spin(ibeam, sbeam)
      % Calculate change in spin between beams
      % For details on usage/arguments see :meth:`forcetorque`.
      %
      % Usage
      %   torque = ibeam.torque(sbeam)
      %   [tx, ty, tz] = ibeam.torque(sbeam)

      % Get the abpq terms for the calculation
      [ai, bi, p, q, n, m, ...
        anp1, bnp1, pnp1, qnp1, ...
        amp1, bmp1, pmp1, qmp1, ...
        anp1mp1, bnp1mp1, pnp1mp1, qnp1mp1, ...
        anp1mm1, bnp1mm1, pnp1mm1, qnp1mm1] = ...
      ibeam.find_abpq_force_terms(sbeam);

      Cz=m./n./(n+1).*(-(ai).*conj(ai)+conj(q).*(q)-(bi).*conj(bi)+conj(p).*(p));
      Dz=-2./(n+1).*sqrt(n.*(n-m+1).*(n+m+1).*(n+2)./(2*n+3)./(2*n+1)) ...
            .*real(anp1.*conj(bi)-bnp1.*conj(ai)-(pnp1).*conj(q) ...
            +(qnp1).*conj(p));

      sz = sum(Cz+Dz);

      Cxy=1i./n./(n+1).*sqrt((n-m).*(n+m+1)).* ...
          (conj(pmp1).*p - conj(amp1).*ai + conj(qmp1).*q - conj(bmp1).*bi);
      Dxy=1i./(n+1).*sqrt(n.*(n+2))./sqrt((2*n+1).*(2*n+3)).* ...
            ( (sqrt((n+m+1).*(n+m+2)) .* ...
            ( p.*conj(qnp1mp1) - q.* conj(pnp1mp1) - ...
            ai.*conj(bnp1mp1) +bi.*conj(anp1mp1))) + ...
            (sqrt((n-m+1).*(n-m+2)) .* ...
            (pnp1mm1.*conj(q) - qnp1mm1.*conj(p) ...
            - anp1mm1.*conj(bi) + bnp1mm1.*conj(ai))) );

      sxy=sum(Cxy+Dxy);
      sy=real(sxy);
      sx=imag(sxy);

      % Ensure things are full
      sx = full(sx);
      sy = full(sy);
      sz = full(sz);

      % Package output
      if nargout == 3
        varargout{1:3} = {sx, sy, sz};
      else
        varargout{1} = [sx(:) sy(:) sz(:)].';
      end
    end

    %
    % Translation functions
    %

    function [beam, A, B] = translateZ(beam, z, varargin)
      % Apply translation along z axis
      %
      % Usage
      %   beam = beam.translateZ(z, ...)
      %
      %   [beam, A, B] = beam.translateZ(z, ...)
      %   Returns A, B for repeated translation calculations.
      %
      % Optional named arguments
      %   - Nmax (numeric) -- Requested minimum Nmax for translated beam.
      %     Default: ``beam.Nmax``.
      %
      %   - basis (enum) -- Type of translation, must be either
      %     'incoming', 'outgoing' or 'regular'.  Default is 'regular'
      %     which should be used for most types of translations.
      %     'outgoing' should be applied when the beam is a scattered beam.

      p = inputParser;
      p.addParameter('basis', 'regular');
      p.addParameter('Nmax', beam.Nmax);
      p.parse(varargin{:});

      ott.utils.nargoutCheck(beam, nargout);

      % Determine beam type
      switch p.Results.basis
        case 'incoming'
          translation_type = 'sbesselh2';
        case 'outgoing'
          translation_type = 'sbesselh1';
        case 'regular'
          translation_type = 'sbesselj';
      end

      % Calculate tranlsation matrices
      [A, B] = ott.utils.translate_z(p.Results.Nmax, ...
          z, 'type', translation_type);

      % Apply translation to beam
      beam = [A, B; B, A] * beam;

    end

  end

  methods (Hidden)
    function [beam, D] = rotateInternal(beam, R, varargin)
      % Apply a rotation to the beam shape coefficients.
      %
      % This function is called by the rotate* functions and is
      % typically not called directly.
      %
      % Usage
      %   [beam, D] = beam.rotateInternal(R)
      %   Returns the rotated beam and the wigner rotation matrix.
      %
      % Optional named arguments
      %   - Nmax (numeric) -- Requested minimum Nmax for rotated beam.
      %     Default: ``beam.Nmax``.

      p = inputParser;
      p.addParameter('Nmax', beam.Nmax);
      p.parse(varargin{:});

      ott.utils.nargoutCheck(beam, nargout);

      assert(isnumeric(R) && ismatrix(R) && size(R, 1) == 3 ...
          && mod(size(R, 2), 3) == 0, ...
          'R must be 3x3N matrix');

      Nrots = size(R, 2)/3;
      Nbeams = numel(beam);

      assert(Nrots == 1 || Nbeams == 1 || Nrots == Nbeams, ...
          'Number of rotations must match number of beams or be scalar');

      Nwork = max([Nrots, Nbeams]);
      if Nwork > 1
        D = cell(1, Nwork);
        for ii = 1:Nwork
          D = ott.utils.wigner_rotation_matrix(...
              max([beam.Nmax, p.Results.Nmax]), R(:, (1:3) + (ii-1)*3));
        end
      else
        D = {ott.utils.wigner_rotation_matrix(...
            max([beam.Nmax, p.Results.Nmax]), R)};
      end

      % Apply wigner matrices to beam
      for ii = 1:numel(beam)
        nend = max(numel(beam(ii).a), numel(beam(ii).b));
        beam(ii) = D{ii}(:, 1:nend) * beam(ii);
      end

    end

    function [beam, Az, Bz, D] = translateXyzInternal(beam, xyz, varargin)
      % Apply translation to the beam shape coefficients.
      %
      % Applies rotations by first rotating the beam z-axis to align
      % to the translation axis, then translating along z, before
      % rotating back to the original orientation.
      %
      % This function uses units of wavelengths.
      %
      % This function is called by the rotate* functions and is
      % typically not called directly.
      %
      % Usage
      %   beam = beam.translateXyzInternal(xyz)
      %   Applies the rotation to the beam.
      %
      %   [beam, Az, Bz, D] = beam.translateXyzInternal(xyz)
      %   Returns Az, Bz and Wigner rotation matrix for repeated operations.
      %
      % Optional named arguments
      %   - Nmax (numeric) -- Requested minimum Nmax for translated beam.
      %     The Nmax limit is applied during the translation. The first
      %     rotation step uses the full Nmax, the last uses the new Nmax.
      %     Ignored when multiple outputs are requested.
      %     Default: ``beam.Nmax``.
      %
      %   - basis (enum) -- Type of translation, must be either
      %     'incoming', 'outgoing' or 'regular'.  Default is 'regular'
      %     which should be used for most types of translations.
      %     'outgoing' should be applied when the beam is a scattered beam.

      p = inputParser;
      p.addParameter('Nmax', beam.Nmax);
      p.addParameter('basis', 'regular');
      p.parse(varargin{:});

      ott.utils.nargoutCheck(beam, nargout);

      assert(isnumeric(xyz) && ismatrix(xyz) && size(xyz, 1) == 3, ...
          'P must be 3xN numeric matrix');

      Nrots = size(xyz, 2);
      Nbeams = size(beam, 2);

      assert(Nrots == 1 || Nbeams == 1 || Nrots == Nbeams, ...
          'Number of positions must match number of beams or be scalar');

      oNmax = p.Results.Nmax;
      if nargout > 1
        oNmax = beam.Nmax;
      end

      % Check for work to do
      hasWork = false;
      for ii = 1:numel(Nrots)
        hasWork = hasWork | any(xyz(:, ii) ~= [0;0;0]);
      end

      Az = cell(1, Nrots);
      Bz = cell(1, Nrots);
      D = cell(1, Nrots);

      if hasWork

        ibsc = beam;

        for ii = 1:Nrots

          rtp = ott.utils.xyz2rtp(xyz(:, ii));
          R = ott.utils.rotz(rtp(3)*180/pi) * ott.utils.roty(rtp(2)*180/pi);

          if nargout > 1
            [newbeam, D{ii}] = ibsc.rotate(R);
            [newbeam, Az{ii}, Bz{ii}] = newbeam.translateZ(...
                rtp(1), 'Nmax', oNmax, 'basis', p.Results.basis);
            newbeam = D{ii}' * newbeam;
          else
            newbeam = ibsc.rotate(R);
            newbeam = newbeam.translateZ(rtp(1), ...
                'Nmax', oNmax, 'basis', p.Results.basis);
            newbeam = newbeam.rotate(R.');
          end

          beam(ii) = newbeam;
        end
      else
        beam = repmat(beam, 1, Nrots);
        Az = repmat({1}, 1, Nrots);
        Bz = repmat({0}, 1, Nrots);
        D = repmat({1}, 1, Nrots);
      end

      if Nrots == 1 && nargout > 1
        Az = Az{1};
        Bz = Bz{1};
        D = D{1};
      end
    end

    function [a, b, p, q, n, m, ...
      anp1, bnp1, pnp1, qnp1, ...
      amp1, bmp1, pmp1, qmp1, ...
      anp1mp1, bnp1mp1, pnp1mp1, qnp1mp1, ...
      anp1mm1, bnp1mm1, pnp1mm1, qnp1mm1] = ...
    find_abpq_force_terms(ibeam, sbeam)
      % Find the terms required for force/torque calculations
      %
      % From forcetorque function in OTTv1

      % Ensure beams are the same size
      if ibeam.Nmax > sbeam.Nmax
        sbeam.Nmax = ibeam.Nmax;
      elseif ibeam.Nmax < sbeam.Nmax
        ibeam.Nmax = sbeam.Nmax;
      end

      % Get the relevant beam coefficients
      [a, b] = ibeam.getCoefficients();
      [p, q] = sbeam.getCoefficients();
      [n, m] = ott.utils.combined_index((1:size(a, 1)).');

      nmax=ibeam.Nmax;

      b=1i*b;
      q=1i*q;

      addv=zeros(2*nmax+3,1);

      at=[a;repmat(addv, 1, size(a, 2))];
      bt=[b;repmat(addv, 1, size(b, 2))];
      pt=[p;repmat(addv, 1, size(p, 2))];
      qt=[q;repmat(addv, 1, size(q, 2))];

      ci=ott.utils.combined_index(n,m);

      %these preserve order and number of entries!
      np1=2*n+2;
      cinp1=ci+np1;
      cinp1mp1=ci+np1+1;
      cinp1mm1=ci+np1-1;
      cimp1=ci+1;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %this is for m+1... if m+1>n then we'll ignore!
      kimp=(m>n-1);

      anp1=at(cinp1, :);
      bnp1=bt(cinp1, :);
      pnp1=pt(cinp1, :);
      qnp1=qt(cinp1, :);

      anp1mp1=at(cinp1mp1, :);
      bnp1mp1=bt(cinp1mp1, :);
      pnp1mp1=pt(cinp1mp1, :);
      qnp1mp1=qt(cinp1mp1, :);

      anp1mm1=at(cinp1mm1, :);
      bnp1mm1=bt(cinp1mm1, :);
      pnp1mm1=pt(cinp1mm1, :);
      qnp1mm1=qt(cinp1mm1, :);

      amp1=at(cimp1, :);
      bmp1=bt(cimp1, :);
      pmp1=pt(cimp1, :);
      qmp1=qt(cimp1, :);

      amp1(kimp, :)=0;
      bmp1(kimp, :)=0;
      pmp1(kimp, :)=0;
      qmp1(kimp, :)=0;

      a=a(ci, :);
      b=b(ci, :);
      p=p(ci, :);
      q=q(ci, :);

    end
  end

  methods (Static, Hidden)
    function visualiseShowPlot(vis_nargout, our_axes, imout, ...
        ranges, labels)
      % Helper function for generating the visualisation
      % Called by :meth:`visualise`.
      %
      % Usage
      %   beam.visualiseShowPlot(vis_nargout, our_axes, imout, ranges, labels)

      if vis_nargout == 0 || ~isempty(our_axes)

        % Get the default axes
        if isempty(our_axes)
          our_axes = gca();
        end

        % Check that we only have one visualisation to show
        if ~ismatrix(imout)
          warning('ott:beam:Beam:non_matrix_plot', ...
              'Multiple beams generated, only showing first');
        end

        % Generate the plot
        imagesc(our_axes, ranges{1}, ranges{2}, imout(:, :, 1), ...
            'AlphaData', ~isnan(imout(:, :, 1)));
        xlabel(our_axes, labels{1});
        ylabel(our_axes, labels{2});
        axis(our_axes, 'image');
      end
    end

    function [xrange, yrange, sz] = visualiseGetRange(range, sz, default_sz)
      % Get the xrange and yrange dimensions of the plot
      %
      % Usage
      %   [xrange, yrange, sz] = visualiseGetRange(range, sz, default_sz)

      if iscell(range)

        if ~isempty(sz)
          warning('Ignoring size parameter');
        end

        xrange = range{1};
        yrange = range{2};
        sz = [length(yrange), length(xrange)];

      elseif length(range) == 2

        if isempty(sz)
          sz = default_sz;
        end
        assert(isnumeric(sz) && numel(sz) == 2, ...
            'size must be 2 element numeric vector');

        xrange = linspace(-1, 1, sz(2))*range(1);
        yrange = linspace(-1, 1, sz(1))*range(2);

      elseif length(range) == 4

        if isempty(sz)
          sz = default_sz;
        end
        assert(isnumeric(sz) && numel(sz) == 2, ...
            'size must be 2 element numeric vector');

        xrange = linspace(range(1), range(2), sz(2));
        yrange = linspace(range(3), range(4), sz(1));

      else
        error('ott:Bsc:visualise:range_error', ...
            'Incorrect type or size of range arguments');
      end
    end

    function [xyz, labels] = visualiseGetXyz(xx, yy, zz, our_axis)
      % Helper function for :meth:`visualise`
      %
      % Usage
      %   [xyz, labels] = beam.visualiseGetXyz(xx, yy, zz, axis)
      %
      %   [~, labels] = beam.visualiseGetXyz([], [], [], axis)

      % TODO: Should we split this function?

      % Change the grid to the user requested orientation
      if ischar(our_axis)
        switch our_axis
          case 'x'
            xyz = [zz(:), yy(:), xx(:)].';
            labels = {'Z', 'Y'};
          case 'y'
            xyz = [yy(:), zz(:), xx(:)].';
            labels = {'Z', 'X'};
          case 'z'
            xyz = [xx(:), yy(:), zz(:)].';
            labels = {'X', 'Y'};
          otherwise
            error('Unknown axis name specified');
        end
      elseif iscell(our_axis) && numel(our_axis) >= 2
        dir1 = our_axis{1}(:);
        dir2 = our_axis{2}(:);
        if numel(our_axis) == 3
          dir3 = our_axis{3}(:);
        else
          dir3 = cross(dir1(:), dir2(:));
        end

        assert(isnumeric(dir1) && numel(dir1) == 3 ...
            && isnumeric(dir2) && numel(dir2) == 3 ...
            && isnumeric(dir3) && numel(dir3) == 3, ...
            'direction vectors must be 3 element numeric');

        labels = {'Direction 1', 'Direction 2'};

        if ~isempty(xx) && ~isempty(yy) && ~isempty(zz)
          xyz = dir1.*xx(:).' + dir2.*yy(:).' + dir3.*zz(:).';
        else
          xyz = [];
        end
      else
        error('axis must be character or a 2 or 3 element cell array');
      end
    end
  end

  methods % Getters/setters
    function nmax = get.Nmax(beam)
      % Calculates Nmax from the current size of the beam coefficients
      nmax = ott.utils.combined_index(size(beam.a, 1));
    end
    function beam = set.Nmax(beam, nmax)
      % Resizes the beam vectors (a,b)
      beam = beam.setNmax(nmax);
    end

    function p = get.power(beam)
      % Calculate beam power from a/b
      [oa, ob] = beam.getCoefficients();
      p = gather(full(sum(abs(oa).^2 + abs(ob).^2, 1)));
    end
    function beam = set.power(beam, val)
      % Set the beam power
      beam = sqrt(val ./ beam.power) * beam;
    end

    function beam = set.a(beam, val)
      % Validate new a values
      if isempty(val)
        val = zeros(0, 1, 'like', val);
      end
      assert(iscolumn(val), ...
          'a must be a column vector');
      beam.a = val;
    end
    function beam = set.b(beam, val)
      % Validate new b values
      if isempty(val)
        val = zeros(0, 1, 'like', val);
      end
      assert(iscolumn(val), ...
          'b must be a column vector');
      beam.b = val;
    end
  end
end
