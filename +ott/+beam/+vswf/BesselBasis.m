classdef BesselBasis < ott.beam.vswf.Bsc ...
    & ott.beam.properties.BesselArray
% Bsc array using Bessel beam basis set.
% Inherits from :class:`BesselArray`.
%
% Properties
%   - angle         -- Far-field angle of Bessel beam (radians)
%   - field         -- Field in theta and phi directions
%   - lmode         -- Azimuthal angular momentum number
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
% Static methods
%   - empty           -- Construct an empty beam array.
%   - likeProperties  -- Form argument list of like-properties
%   - like            -- Construct object like another beam
%   - FromCartesianField      -- Construct from Cartesian polarisation
%   - CartesianFieldWeights   -- Calculate weights for Cartesian field
%
% Additional methods inherited from base.

% Based on code by Alexander Stilgoe.
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  methods (Static)
    function bsc = empty(varargin)
      % Construct an empty beam array
      %
      % Usage
      %   bsc = BesselBasis.empty()

      bsc = ott.beam.vswf.BesselBasis();
    end

    function args = likeProperties(other, args)
      % Construct an array of like-properties
      args = ott.beam.vswf.Bsc.likeProperties(other, args);
      args = ott.beam.properties.BesselArray.likeProperties(other, args);
    end

    function beam = like(other, varargin)
      % Create a beam like another beam
      %
      % Usage
      %   beam = BesselBasis.like(other, ...)
      %
      % See constructor for arguments.

      args = ott.beam.vswf.BesselBasis.likeProperties(other, varargin);
      beam = ott.beam.vswf.BesselBasis(args{:});
    end

    function [bsc, weights] = FromCartesianField(varargin)
      % Construct beam from Cartesian polarisation vector.
      %
      % Usage
      %   bsc = FromParaxialPolarisation(angle, polarisation, ...)
      %   Applies the polarisation weights and returns a :class:`Bsc`.
      %
      %   [bsc, weights] = FromParaxialPolarisation(angle, ...)
      %   Returns a :class:`BesselBasis` and the weights to convert
      %   it to a :class:`Bsc`.  The resulting beam has twice as many
      %   beam shape coefficients, two for each polarisation.
      %
      % Parameters
      %   - angle (N numeric) -- Bessel incoming angle.
      %
      %   - polarisation (2xN numeric) -- Paraxial polarisation [x;y].
      %     If omitted, polarisation weights are not computed or applied.
      %
      % Returns
      %   - bsc -- The beam shape coefficients.
      %   - weights (2xN numeric) -- Weights vector that can be
      %     applied to generate Bessel beams.   Needs 'polarisation'.
      %     See :meth:`CartesianFieldWeights`.
      %
      % Optional named parameters
      %   - lmode (N numeric) -- Orbital angular momentum number.
      %     Default: ``0``.
      %
      % Other inputs are passed to constructor.

      p = inputParser;
      p.addOptional('angle', [], @isnumeric);
      p.addOptional('polarisation', [], @isnumeric);
      p.addParameter('lmode', 0);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      theta = p.Results.angle(:).';
      lmode = p.Results.lmode(:).';

      Ntheta = numel(theta);
      Nlmode = numel(lmode);
      Nwork = max([Ntheta, Nlmode]);

      assert(Ntheta == 1 || Ntheta == Nwork, ...
          'theta must be scalar or inputs must have matching size');
      assert(Nlmode == 1 || Nlmode == Nwork, ...
          'lmode must be scalar or inputs must have matching size');

      % Calculate weights (before resizing theta/lmode)
      if ~isempty(p.Results.polarisation)
        weights = ott.beam.properties.Bessel.CartesianFieldWeights(...
            theta, p.Results.polarisation);
      else
        weights = [];
      end

      % Duplicate as needed
      if Ntheta == 1, theta = repmat(theta, 1, Nwork); end
      if Nlmode == 1, lmode = repmat(lmode, 1, Nwork); end

      % Parameters for basis
      lmode = reshape([lmode-1; lmode+1], 1, []);
      theta = reshape([theta; theta], 1, []);
      Etp = repmat([1,1;-1i,1i], 1, Nwork);

      % Construct beam
      bsc = ott.beam.vswf.BesselBasis('lmode', lmode, ...
          'field', Etp, 'angle', theta, unmatched{:});

      % Apply weights if needed
      if nargout == 1
        assert(~isempty(p.Results.polarisation), ...
            'polarisation needed to calculate weights');
        bsc = bsc.applyWeights(weights);
      end
    end
  end

  methods
    function beam = BesselBasis(varargin)
      % Construct a new Bessel basis set
      %
      % Usage
      %   beam = BesselBasis(angle, field, ...)
      %
      % Parameters
      %   - angle (N numeric) -- Incoming angle of Bessel beams (radians).
      %
      %   - field (2xN numeric) -- Field in spherical coordinates.
      %     Format [Etheta; Ephi].  For an alternative interface using
      %     Cartesian coordinates, see :meth:`FromCartesianField`.
      %
      % Optional named parameters
      %   - lmode (N numeric) -- Orbital angular momentum mode number.
      %     Default: ``0``.
      %
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
      p.addOptional('angle', [], @isnumeric);
      p.addOptional('field', [], @isnumeric);
      p.addParameter('Nmax', 0, @isnumeric);
      p.addParameter('lmode', 0, @isnumeric);
      p.addParameter('array_type', 'coherent');
      p.addParameter('basis', 'regular');
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.vswf.Bsc('basis', p.Results.basis, ...
          'array_type', p.Results.array_type);
      beam = beam@ott.beam.properties.BesselArray(...
          'angle', p.Results.angle, ...
          'lmode', p.Results.lmode, ...
          'field', p.Results.field, ...
          unmatched{:});

      % Generate coefficients
      beam = beam.updateCoefficients(p.Results.Nmax);
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
          size@ott.beam.properties.BesselArray(varargin{:});
    end

    function beam = setData(beam, angle, field, lmode)
      % Set beam data
      %
      % Raises a warning if the angle is near zero or pi.
      %
      % Usage
      %   beam = beam.setData(angle, field, lmode)
      %
      % Parameters
      %   - angle (N numeric) -- Incoming angle (radians).
      %
      %   - field (2xN numeric) -- Field in spherical coordinates.
      %
      %   - lmode (N numeric) -- Azimuthal mode number.

      % Check theta
      if any(abs(mod(angle, pi)) < 1e-6)
        warning('ott:beam:vswf:Bessel:angle_zero', ...
          'angle = 0+/-pi is not a Bessel beam, may have unexpected results');
      end

      beam = setData@ott.beam.properties.Bessel(beam, angle, field, lmode);
    end

    function beam = updateCoefficients(beam, Nmax)
      % Re-calculate VSWF coefficients for specified Nmax
      %
      % Usage
      %   beam = beam.updateCoefficients(Nmax)
      %
      % Parameters
      %   - Nmax (numeric) -- Truncation for VSWF coefficients.
      %     For a simpler interface for creating beams without explicit
      %     `Nmax`, see :class:`vswf.Bessel`.

      assert(isnumeric(Nmax) && isscalar(Nmax) ...
          && Nmax >= 0 && round(Nmax) == Nmax, ...
          'Nmax should be an single positive integer');

      %% calculate the mode indices we are going to find.
      nTheta = length(beam.angle);
      a = zeros((Nmax*(Nmax+2)), nTheta);
      b = zeros((Nmax*(Nmax+2)), nTheta);

      indexes=[1:length(beam.angle)].';
      lmode = beam.lmode.';
      theta = beam.angle.';
      Etheta = beam.field(1, :).';
      Ephi = beam.field(2, :).';

      for n = 1:Nmax

        ci_index_m=find(abs(lmode)<=n);
        indt=n+lmode(ci_index_m)+1;

        %power normalisation.
        Nn = 1/sqrt(n*(n+1));

        %Generate the farfield components of the VSWFs
        [~,dtY,dpY] = ott.utils.spharm(n,theta, 0);

        %slow indexing.
        szA=sub2ind(size(a),(n-1)*(n+1)+indt,indexes(ci_index_m));
        szY=sub2ind(size(dtY),indexes(ci_index_m),indt);

        dtY=dtY(:);
        dpY=dpY(:);

        %equivalent to dot((1i)^(n+1)*C,E);
        a(szA) = 4*pi*Nn*(-1i)^(n+1) ...
            *(conj(dpY(szY)).*Etheta(ci_index_m) ...
            - conj(dtY(szY)).*Ephi(ci_index_m));
        %equivalent to dot((1i)^(n)*B,E);
        b(szA) = 4*pi*Nn*(-1i)^(n)  ...
            *(conj(dtY(szY)).*Etheta(ci_index_m) ...
            + conj(dpY(szY)).*Ephi(ci_index_m));
      end

      % Make the beam vector and store the coefficients
      beam = beam.setCoefficients(a, b);
      beam = beam.makeSparse();
    end

    function beam = applyWeights(beam, weights, varargin)
      % Apply weights to beam and construct a BSC instance
      %
      % Effectively applies
      %
      %   beam = beam * weights(:)
      %
      % Usage
      %   beam = beam.applyWeights(beam, weights)
      %
      % Parameters
      %   - weights (2xM|N numeric) -- Weights to apply.
      %     If the weights vector is 2xM, applies pairs of weights
      %     to consecutive pairs of beams, producing M beams.
      %
      % Optional named arguments
      %   - combine (logical) -- Combine beams.  Default: ``true``.
      %     The combined beam is a Bsc instance not a BesselBasis.

      p = inputParser;
      p.addParameter('combine', true);
      p.parse(varargin{:});

      % Apply weights
      beam = beam .* weights(:).';

      if p.Results.combine
        if size(weights, 1) == 2
          % Combine pairs
          idx = 1:2:numel(beam);
          beam = subsref(beam, substruct('()', {idx})) ...
            + subsref(beam, substruct('()', {idx+1}));
        else
          % Combine all
          beam = sum(beam);
        end
      end
    end

    function varargout = applyZTranslation(beam, varargin)
      % Apply Z translation by applying phase shift to beam.
      %
      % Usage
      %   bsc = bsc.applyZTranslation(...)
      %   Apply the z-component of the translation.
      %
      %   bsc = bsc.applyZTranslation(z, ...)
      %   Apply a specific z translation to the beam.
      %
      % Optional named arguments
      %   - Nmax (numeric) -- Ignored.

      p = inputParser;
      p.addOptional('z', beam.position(3), @isnumeric);
      p.addParameter('Nmax', []);
      p.parse(varargin{:});

      if nargout == 1

        ott.utils.nargoutCheck(beam, nargout);

        z = p.Results.z;
        assert(isnumeric(z) && isvector(z), ...
            'z must be numeric vector');

        ibeam = beam;
        beam = ott.beam.vswf.Array.empty([1, numel(z)]);

        for ii = 1:numel(z)
          dz = z(ii) .* cos(ibeam.angle);
          dz = exp(1i.*dz.*ibeam.wavenumber);

          % Apply translation
          beam(ii) = ibeam .* dz;
          beam(ii).position = beam(ii).position - [0;0;z];
        end

        if numel(beam) == 1
          beam = beam(1);
        end

        varargout{1} = beam;
      else
        [varargout{1:nargout}] = applyZTranslation@ott.beam.vswf.Bsc(...
            beam, varargin{:});
      end
    end
  end

  methods (Hidden)
    function beam = catInternal(dim, beam, varargin)
      % Concatenate beams
      beam = catInternal@ott.beam.vswf.Bsc(dim, beam, varargin{:});
      beam = catInternal@ott.beam.properties.BesselArray(dim, ...
          beam, varargin{:});
    end

    function beam = plusInternal(beam1, beam2)
      %PLUS add two beams together
      beam = plusInternal@ott.beam.vswf.Bsc(beam1, beam2);
    end

    function beam = subsrefInternal(beam, subs)
      % Get the subscripted beam
      beam = subsrefInternal@ott.beam.vswf.Bsc(beam, subs);
      beam = subsrefInternal@ott.beam.properties.BesselArray(beam, subs);
    end

    function beam = subsasgnInternal(beam, subs, rem, other)
      % Set the subscripted beam
      beam = subsasgnInternal@ott.beam.vswf.Bsc(beam, subs, rem, other);
      beam = subsasgnInternal@ott.beam.properties.BesselArray(...
          beam, subs, rem, other);
    end
  end
end
