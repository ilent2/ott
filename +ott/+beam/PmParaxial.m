classdef PmParaxial < ott.beam.BscBeam ...
    & ott.beam.properties.Mapping ...
    & ott.beam.properties.Polarisation ...
    & ott.beam.properties.Profile
% Construct a beam using paraxial far-field point matching.
% Inherits from :class:`BscBeam` and :class:`+properties.Mapping`,
% :class:`+properties.Polarisation` and :class:`+properties.Profile`.
%
% Properties
%   - mapping     -- Paraxial to far-field mapping
%   - polbasis    -- (enum) Polarisation basis ('polar' or 'cartesian')
%   - polfield    -- (2 numeric) Field in theta/phi or x/y directions
%   - truncation_angle -- Maximum angle for beam.  Default: []
%   - data        -- Internal BSC data
%
% Methods
%   - recalculate     -- Update the internal data for new Nmax
%
% Static methods
%   - InterpProfile   -- Generate a beam profile using interpolation
%   - BeamProfile     -- Construct beam profile from another beam
%
% Supported casts
%   - ott.bsc.Bsc     -- Construct bsc instance

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    truncation_angle
  end

  properties (Hidden, SetAccess=protected)
    truncation_angleInternal
  end

  methods (Static)
    function beam = InterpProfile(X, Y, A, varargin)
      % Construct a beam using interpolation
      %
      % Usage
      %   beam = InterpProfile(X, Y, A, ...)
      %
      % Parameters
      %   - X, Y (NxM numeric) -- Paraxial coordinates.  Should be in
      %     ``ndgrid`` format (uses ``griddedInterpolant``).
      %
      %   - A (NxM | 2xNxM numeric) -- Complex field amplitude array.
      %
      % All other parameters passed to constructor.

      if ismatrix(A)
        F = griddedInterpolant(X, Y, A);
        profile = @(x, y) reshape(F(x, y), [1, size(x)]);
      else
        assert(ndims(A) == 3, 'A must be NxM or 2xNxM array');
        F1 = griddedInterpolant(X, Y, squeeze(A(1, :, :)));
        F2 = griddedInterpolant(X, Y, squeeze(A(2, :, :)));
        profile = @(x, y) [reshape(F1(x, y), [1, size(x)]); ...
                     reshape(F2(x, y), [1, size(x)])];
      end
      beam = ott.beam.PmParaxial(profile, varargin{:});
    end

    function beam = BeamProfile(beam, varargin)
      % Construct a beam using another beam as the profile
      %
      % Usage
      %   beam = BeamProfile(other_beam, ...)
      %
      % Optional named arguments
      %   - mapping (enum) -- Mapping for paraxial projection.
      %     See :func:`ott.utils.paraxial2rtp` for details.
      %     Default: ``'sin'``.
      %
      %   - direction (enum | 2 numeric) -- Mapping direction for
      %     other beam.  See :meth:`eparaxial` for details.
      %     Default: ``'pos'``.
      %
      % Additional arguments passed to constructor.

      p = inputParser;
      p.addParameter('mapping', 'sin');
      p.addParameter('direction', 'pos');
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);
      
      profile = @(x, y) reshape(beam.eparaxial([x(:), y(:)].', ...
          'mapping', p.Results.mapping, ...
          'direction', p.Results.direction).vxyz(1:2, :), [2, size(x)]);

      beam = ott.beam.PmParaxial(profile, 'mapping', p.Results.mapping, ...
          unmatched{:});
    end
  end

  methods
    function beam = PmParaxial(varargin)
      % Construct beam using paraxial point matching.
      %
      % The beam always comes from the same direction.  To change the
      % direction, modify the beam rotation property.
      %
      % Usage
      %   beam = PmParaxial(profile, ...)
      %
      % Optional named parameters
      %   - profile (function_handle) -- Profile for
      %     paraxial point matching.  Default: ``@(x,y) ones([1,size(x)])``.
      %     The profile function must return with a 1xNxM or 2xNxM array.
      %
      %   - mapping (enum) -- Mapping method for paraxial far-field.
      %     Can be either 'sin', 'tan' (small angle) or 'theta'.
      %     For a discussion of this parameter, see Documentation
      %     (:ref:`conception-angular-scaling`).  Default: ``'sin'``.
      %
      %   - truncation_angle (numeric) -- Maximum angle for beam.
      %     Must be between 0 and pi/2.  Default: ``pi/2``.
      %
      %   - polfield (2 numeric) -- Field in the [theta, phi] or [x, y]
      %     directions.  Default: ``[1, 1]``.
      %
      %   - polbasis (enum) -- Polarisation basis.  Can be either
      %     'polar' or 'cartesian'.  Default: ``'cartesian'``.
      %
      %   - Nmax (numeric) -- Initial beam Nmax.  Default: ``0``.
      %     This parameter automatically grows when the beam is used.
      %     See also :meth:`recalculate` and :meth:`getData`.

      p = inputParser;
      p.addOptional('profile', @(x,~) zeros([1, size(x)]));
      p.addParameter('truncation_angle', pi/2);
      p.addParameter('polfield', [1, 1]);
      p.addParameter('polbasis', 'cartesian');
      p.addParameter('mapping', 'sin');
      p.addParameter('Nmax', 0);
      p.parse(varargin{:});

      beam.profile = p.Results.profile;
      beam.mapping = p.Results.mapping;
      beam.truncation_angle = p.Results.truncation_angle;
      beam.polfield = p.Results.polfield;
      beam.polbasis = p.Results.polbasis;
      beam = beam.recalculate(p.Results.Nmax);
    end

    function beam = recalculate(beam, Nmax)
      % Re-calculate the beam data.
      %
      % This function can be called to pre-compute the beam data for
      % the current beam parameters.  It is automatically called when
      % the beam is used, however, calling the function explicitly can
      % speed up run-time.
      %
      % Usage
      %   beam = beam.recalcualte(Nmax)
      %
      % Parameters
      %   - Nmax (numeric) -- Truncation number for VSWF expansion.

      nTheta = 2*(Nmax+1);
      nPhi = 2*(Nmax+1);

      % Generate coordinates for far-field mapping
      theta = linspace(0, pi, nTheta);
      phi = linspace(0, 2*pi, nPhi);
      [theta, phi] = meshgrid(theta, phi);

      % Find which coordinates we are interested in
      idx = theta < beam.truncation_angle;

      % Get mapping coordinates in paraxial coordinates
      switch beam.mapping
        case 'sin'
          Xt = sin(theta(idx)).*cos(phi(idx));
          Yt = sin(theta(idx)).*sin(phi(idx));
        case 'tan'
          Xt = tan(theta(idx)).*cos(phi(idx));
          Yt = tan(theta(idx)).*sin(phi(idx));
        case 'theta'
          Xt = theta(idx).*cos(phi(idx));
          Yt = theta(idx).*sin(phi(idx));
        otherwise
          error('Unknown mapping parameter');
      end

      % Calculate field
      E = zeros([2, size(theta)]);
      data = beam.profile(Xt, Yt);
      if size(data, 1) == 1
        data = repmat(data, 2, 1, 1);
      end
      E(:, idx) = data;
      E = E .* beam.wavenumber;

      % Flip theta direction (so beam comes from negative direction)
      % This is to match the convention used by other beams.
      theta = pi - theta;

      % Convert Cartesian fields to polar fields
      convertFields = @(E) deal(...
          sign(cos(theta)).*cos(phi).*squeeze(E(1, :, :)) ...
             + sign(cos(theta)).*sin(phi).*squeeze(E(2, :, :)), ...
          -sin(phi).*squeeze(E(1, :, :)) + cos(phi).*squeeze(E(2, :, :)));
      switch beam.polbasis
        case 'cartesian'
          E = E .* beam.polfield(:);
          [Et, Ep] = convertFields(E);
        case 'polar'
          [Et, Ep] = convertFields(E);
          Et = Et .* beam.polfield(1);
          Ep = Ep .* beam.polfield(2);
        otherwise
          error('Unknown polarisation basis');
      end

      % Do point matching
      ci = 1:ott.utils.combined_index(Nmax, Nmax);
      Efield = [Et(:), Ep(:)].';
      rtp = [theta(:), phi(:)].';
      beam.data = ott.bsc.Bsc.PmFarfield(rtp, Efield, ci);
    end
  end

  methods % Getters/setters
    function beam = set.truncation_angle(beam, val)
      assert(isnumeric(val) && isscalar(val) && val > 0 && val <= pi/2, ...
          'truncation_angle must be numeric scalar in [0, pi/2]');
      beam.truncation_angleInternal = val;
      beam.data = [];
    end
    function val = get.truncation_angle(beam)
      val = beam.truncation_angleInternal;
    end
  end
end

