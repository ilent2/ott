classdef PlaneBasis < ott.beam.vswf.Bsc ...
    & ott.beam.properties.PlaneWaveArray
% Bsc array using a plane-wave basis set.
% Inherits from :class:`Bsc` and :class:`PlaneWaveArray`.
%
% Inherited properties
%   - field         -- Field parallel and perpendicular to polarisation
%   - directionSet  -- Set of direction vectors describing orientation
%   - origin        -- Position used to calculate beam phase offset
%   - position      -- Beam position
%   - rotation      -- Beam rotation
%   - a             -- Beam shape coefficients a vector
%   - b             -- Beam shape coefficients b vector
%   - basis         -- VSWF beam basis (incoming, outgoing or regular)
%   - absdz         -- Absolute cumulative distance the beam has moved
%   - Nmax          -- Truncation number for VSWF coefficients
%   - power         -- The power of the beam (may be infinite)
%   - omega         -- Beam optical frequency
%   - medium        -- Medium where beam is propagating
%
% Methods
%   - applyTranslation  -- Apply phase shift to plane waves.
%   - applyZTranslation -- Apply phase shift to plane waves.
%   - setData           -- Set data and update VSWF coefficients
%
% Static methods
%   - empty         -- Construct an empty beam array.
%   - likeProperties -- Form argument list of like-properties
%   - like          -- Construct object like another beam
%   - FromDirection -- Construct array from direction vectors
%
% All other methods inherited from base.

% Based on code by Alexander Stilgoe.
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% TODO: Rotation should leave type unchanged

  methods (Static)
    function bsc = empty(varargin)
      % Construct an empty beam array
      %
      % Usage
      %   bsc = PlaneBasis.empty()

      bsc = ott.beam.vswf.PlaneBasis();
    end

    function args = likeProperties(other, args)
      % Construct an array of like-properties
      args = ott.beam.vswf.Bsc.likeProperties(other, args);
      args = ott.beam.properties.PlaneWaveArray.likeProperties(other, args);
    end

    function beam = like(other, varargin)
      % Create a beam like another beam
      %
      % Usage
      %   beam = PlaneBasis.like(other, ...)
      %
      % See constructor for arguments.

      args = ott.beam.vswf.PlaneBasis.likeProperties(other, varargin);
      beam = ott.beam.vswf.PlaneBasis(args{:});
    end

    function beam = FromDirection(varargin)
      % Construct beam from direction/polarisation vectors.
      %
      % Usage
      %   beam = FromDirection(origin, direction, polarisation, field, ...)
      %   Parameters can also be passed as named arguments.
      %
      % Parameters
      %   - origin (3xN numeric) -- Origin (for phase offset) of wave.
      %   - direction (3xN numeric) -- Propagation direction of wave.
      %   - polarisation (3xN numeric) -- Primary polarisation direction.
      %   - field (2xN numeric) -- Field in two polarisation directions.
      %
      % Additional parameters passed to base.  See class constructor
      % for additional information on parameters.

      % TODO: This method is duplicated quite a lot
      %   There are other methods we would like too, perhaps there
      %   is a better way to do this?

      p = inputParser;
      p.addOptional('origin', [], @isnumeric);
      p.addOptional('direction', [], @isnumeric);
      p.addOptional('polarisation', [], @isnumeric);
      p.addOptional('field', [], @isnumeric);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Construct direction set
      directionSet = ott.beam.properties.PlaneWave.DirectionSet(...
          p.Results.direction, p.Results.polarisation);

      % Construct beam
      beam = ott.beam.vswf.PlaneBasis(...
          'origin', p.Results.origin, ...
          'directionSet', directionSet, ...
          'field', p.Results.field, ...
          'Nmax', p.Results.Nmax, ...
          unmatched{:});
    end
  end

  methods
    function beam = PlaneBasis(varargin)
      % Construct a new VSWF plane wave beam basis set.
      %
      % Stores properties and calls :meth:`updateCoefficients` to
      % generate the coefficients.  Setting `Nmax = 0` results in
      % empty beam vectors.
      %
      % Usage
      %   beam = PlaneBasis(origin, directionSet, field, ...)
      %
      % See also :meth:`FromDirection` for construction from
      % direction/polarisation vectors.
      %
      % Parameters
      %   - origin (3xN numeric) -- Plane wave origins.
      %
      %   - directionSet (3x3N numeric) -- Array formed by combining
      %     direction/polarisation vectors into rotation matrices.  The
      %     direction vector should be the last column of the matrix.
      %
      %   - field (2xN numeric) -- Field parallel and perpendicular to
      %     plane wave polarisation direction.
      %
      % Optional named arguments
      %   - array_type (enum) -- Beam array type.  Can be
      %     'coherent', 'incoherent' or 'array'.  Default: ``'coherent'``.
      %
      %   - basis (enum) -- VSWF basis: incoming, outgoing or regular.
      %     Default: ``'regular'``.
      %
      %   - Nmax (numeric) -- Truncation for VSWF coefficients.
      %     For a simpler interface for creating beams without explicit
      %     `Nmax`, see :class:`vswf.PlaneWave`.  Default: ``0``.

      p = inputParser;
      p.addOptional('origin', [], @isnumeric);
      p.addOptional('directionSet', [], @isnumeric);
      p.addOptional('field', [], @isnumeric);
      p.addParameter('Nmax', 0, @isnumeric);
      p.addParameter('array_type', 'coherent');
      p.addParameter('basis', 'regular');
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.vswf.Bsc('basis', p.Results.basis, ...
          'array_type', p.Results.array_type);
      beam = beam@ott.beam.properties.PlaneWaveArray(...
          'origin', p.Results.origin, ...
          'directionSet', p.Results.directionSet, ...
          'field', p.Results.field, ...
          unmatched{:});

      % Generate coefficients
      beam = beam.updateCoefficients(p.Results.Nmax);
    end

    function beam = updateCoefficients(beam, Nmax)
      % Update the VSWF coefficients
      %
      % Usage
      %   beam = beam.updateCoefficients(Nmax)
      %
      % Parameters
      %   - Nmax (numeric) -- Truncation for VSWF coefficients.
      %     For a simpler interface for creating beams without explicit
      %     `Nmax`, see :class:`vswf.PlaneWave`.

      assert(isnumeric(Nmax) && isscalar(Nmax) ...
          && Nmax >= 0 && round(Nmax) == Nmax, ...
          'Nmax should be an single positive integer');

      % Get theta/phi field/coordinates
      [rtpv1, rtp] = ott.utils.xyzv2rtpv(beam.polarisation1, beam.direction);
      [rtpv2, ~] = ott.utils.xyzv2rtpv(beam.polarisation2, beam.direction);
      Ertp = rtpv1 .* beam.field(1, :) + rtpv2 .* beam.field(2, :);

      ablength = ott.utils.combined_index(Nmax, Nmax);

      a = zeros(ablength, size(rtp, 2));
      b = zeros(ablength, size(rtp, 2));

      for n = 1:Nmax
        iter= (n-1)*(n+1)+1:n*(n+2);
        leniter=2*n+1;

        %expand theta and phi components of field to match spherical harmonics
        ET=repmat(Ertp(2, :), [1,leniter]);
        EP=repmat(Ertp(3, :), [1,leniter]);

        %power normalisation.
        Nn = 1/sqrt(n*(n+1));

        %Generate the farfield components of the VSWFs
        [~,dtY,dpY] = ott.utils.spharm(n, -n:n, rtp(2, :), rtp(3, :));

        %equivalent to dot((1i)^(n+1)*C,E);
        a(iter,:) = 4*pi*Nn*(-1i)^(n+1)*(conj(dpY).*ET - conj(dtY).*EP).';
        %equivalent to dot((1i)^(n)*B,E);
        b(iter,:) = 4*pi*Nn*(-1i)^(n)  *(conj(dtY).*ET + conj(dpY).*EP).';
      end

      beam = beam.setCoefficients(a, b);
      beam = beam.makeSparse();
    end

    function varargout = size(varargin)
      % Get the number of beams contained in this object
      %
      % Usage
      %   sz = size(beam)   or    sz = beam.size()
      %   For help on arguments, see builtin ``size``.
      %
      % The leading dimension is always 1.  May change in future.

      [varargout{1:nargout}] = ...
          size@ott.beam.properties.PlaneWaveArray(varargin{:});
    end

    function varargout = applyTranslation(beam, varargin)
      % Apply translation by applying phase shift to the beam.
      %
      % Calculates the component of the translation in the beam
      % direction and applies a phase shift to the beam shape coefficients.
      %
      % Leaves the beam type unchanged.
      %
      % Usage
      %   bsc = bsc.applyTranslation(...)
      %   Applies ``bsc.position`` to the beam shape coefficients.
      %
      %   bsc = bsc.applyTranslation(P, ...)
      %   Applies a specific rotation.  P should be a 3xN matrix.
      %
      %   [bsc, Az, Bz, D] = bsc.applyRotation(...)
      %   Additionally returns the translation matrices.
      %
      % Optional named arguments
      %   - Nmax (numeric) -- Requested minimum Nmax for translated beam.
      %     Only used when AB is passed as an input instead of P.
      %     The Nmax limit is applied during the translation. The first
      %     rotation step uses the full Nmax, the last uses the new Nmax.
      %     Ignored when multiple outputs are requested.
      %     Default: ``bsc.Nmax``.
      %
      %   - AB (2xN|3xN cell) -- Cell array of translation matrices to
      %     apply.  For 2xN cell array, assumes {A; B}, for 3xN cell
      %     array {Az; Bz; D} where D is the rotation to the z axis.

      p = inputParser;
      p.addOptional('P', []);
      p.addParameter('Nmax', beam.Nmax);
      p.addParameter('AB', {});
      p.parse(varargin{:});

      ott.utils.nargoutCheck(beam, nargout);

      if isempty(p.Results.AB)
        P = p.Results.P;
        if isempty(P)
          P = beam.position;
        else
          assert(isnumeric(P) && ismatrix(P) && size(P, 1) == 3, ...
              'P must be 3xN numeric matrix');
        end

        ibeam = beam;
        beam = ott.beam.vswf.Array.empty([1, size(P, 2)]);

        for ii = 1:size(P, 2)
          dz = dot(repmat(P(:, ii), 1, ...
              size(ibeam.direction, 2)), ibeam.direction);
          dz = exp(1i.*dz.*ibeam.wavenumber);

          % Apply translation
          beam(ii) = ibeam .* dz;
          beam(ii).position = [0;0;0];
        end

        if numel(beam) == 1
          beam = beam(1);
        end

        varargout{1} = beam;
      else
        [varargout{1:nargout}] = applyTranslation@ott.beam.vswf.Bsc(...
            beam, varargin{:});
      end
    end

    function beam = applyZTranslation(beam, z, varargin)
      % Apply Z translation by applying phase shift to beam.
      %
      % Calls :meth:`applyTranslation`.
      %
      % Usage
      %   bsc = bsc.applyZTranslation(z, ...)
      %
      % Optional named arguments
      %   - Nmax (numeric) -- Ignored.

      p = inputParser;
      p.addParameter('Nmax', []);
      p.parse(varargin{:});

      xyz = [0*z(:), 0*z(:), z(:)].';
      beam = beam.applyTranslation(xyz);
    end
  end

  methods (Hidden)
    function beam = catInternal(dim, beam, varargin)
      % Concatenate beams
      beam = catInternal@ott.beam.vswf.Bsc(dim, beam, varargin{:});
      beam = catInternal@ott.beam.properties.PlaneWaveArray(dim, ...
          beam, varargin{:});
    end

    function beam = plusInternal(beam1, beam2)
      %PLUS add two beams together
      beam = plusInternal@ott.beam.vswf.Bsc(beam1, beam2);
    end

    function beam = subsrefInternal(beam, subs)
      % Get the subscripted beam
      beam = subsrefInternal@ott.beam.vswf.Bsc(beam, subs);
      beam = subsrefInternal@ott.beam.properties.PlaneWaveArray(beam, subs);
    end

    function beam = subsasgnInternal(beam, subs, rem, other)
      % Assign to the subscripted beam
      beam = subsasgnInternal@ott.beam.vswf.Bsc(beam, subs, rem, other);
      beam = subsasgnInternal@ott.beam.properties.PlaneWaveArray(...
          beam, subs, rem, other);
    end
  end
end
