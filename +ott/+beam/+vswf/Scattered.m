classdef Scattered < ott.beam.vswf.Bsc & ott.beam.Scattered
% A Bsc instance describing a scattered beam.
% Inherits from :class:`Bsc`.
%
% An instance of this class is created by the Bsc object when
% a particle scatters a beam.  By default, this object keeps track of the
% incident beam and the T-matrix, allowing easy total-field and
% scattered-field calculation.
%
% Scattered beams can either have a total-field type or scattered-field
% type, depending on if the beam shape coefficients describe only the
% scattered fields or the total fields.
%
% Properties
%   - incident_beam     -- The incident beam that was scattered (or [])
%   - tmatrix           -- The T-matrix which scattered the beam (or [])
%   - type              -- Type of beam (scattered or total)
%
% Dependent properties
%   - total_beam        -- Instance of the beam with total type
%   - scattered_beam    -- Instance of the beam with scattered type
%
% Methods
%   - totalField        -- Calculate the total field representation
%   - scatteredField    -- Calculated the scattered field representation

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    tmatrix            % The T-matrix which scattered the beam (or [])
  end

  methods (Static)
    function [sbeam, tbeam] = FromTmatrix(tbeam, tmatrix, varargin)
      % Calculate the scattered beam using a T-matrix
      %
      % Usage
      %   beam = Scattered(incident_beam, tmatrix, ...)
      %   Calculate the scattered beam using a T-matrix.
      %
      % Optional named arguments
      %   - store_tmatrix (logical) -- If false, the T-matrix is not stored.
      %     Default: ``true``.
      %
      %   - store_beam (logical) -- If false, the incident beam is not stored.
      %     Default: ``true``.
      %
      %   - position (3xN numeric) -- Translation applied to the beam
      %     before the beam is scattered by the particle.  Default: ``[]``.
      %     Only used with T-matrix input.
      %
      %   - rotation (3x3N numeric) -- Rotation applied to the beam,
      %     calculates the scattered beam and applies the inverse rotation,
      %     effectively rotating the particle.  Default: ``[]``.
      %     Only used with T-matrix input.
      %
      % If both position and rotation are present, the translation is
      % applied first, followed by the rotation.
      % If both position and rotation are arrays, they must have the same
      % number of locations.
      %
      % Additional arguments are passed to class constructor.

      % Setup input parser
      p = inputParser;
      p.addParameter('store_beam', true);
      p.addParameter('store_tmatrix', true);
      p.addParameter('position', []);
      p.addParameter('rotation', []);
      p.parse(varargin{:});
      
      % Ensure the beam is a Bsc
      if ~isa(tbeam, 'ott.beam.vswf.Bsc')
        
        % Calculate suggested Nmax from T-matrix and translation
        %
        % TODO: There are two optimal cases we could implement
        %   * For beams with an almost exact representation we should
        %     use the minimum Nmax.
        %   * For PlaneWave beams, we should pre-calculate the
        %     rotations and translatiosn and only convert at the end.
        %     So Nmax matches the T-matrix Nmax.
        maxPosition = max(vecnorm(p.Results.position));
        particleRadius = ott.utils.nmax2ka(tmatrix.Nmax(2));
        Nmax = ott.utils.ka2nmax(maxPosition*tbeam.wavenumber + particleRadius);
        
        tbeam = ott.beam.vswf.Bsc(tbeam, 'suggestedNmax', Nmax);
      end

      % Pre-combine coherent beams
      if strcmpi(tbeam.array_type, 'coherent')
        tbeam = sum(tbeam);
      end

      % Check if we need to calculate multiple scatters
      if size(p.Results.position, 2) > 1 || size(p.Results.rotation, 2) > 3
        [sbeam, tbeam] = ott.utils.prxfun(...
            @(varargin) tbeam.scatter(tmatrix, varargin{:}), 1, ...
            'position', p.Results.position, ...
            'rotation', p.Results.rotation, ...
            'zeros', @(x) repmat(tbeam, x));
        return;
      end

      % Determine the maximum tmatrix.Nmax(2) and check type
      maxNmax1 = 0;
      maxNmax2 = 0;
      tType = tmatrix(1).type;
      for ii = 1:numel(tmatrix)
        maxNmax1 = max(maxNmax1, tmatrix(ii).Nmax(1));
        maxNmax2 = max(maxNmax2, tmatrix(ii).Nmax(2));
        if ~strcmpi(tmatrix(ii).type, tType)
          error('T-matrices must be same type');
        end
      end

      % If the T is scattered, we can save time by throwing away columns
      if strcmpi(tmatrix(1).type, 'scattered')
        maxNmax2 = min(maxNmax2, tbeam.Nmax);
      end

      % Ensure all T-matrices are the same size
      for ii = 1:numel(tmatrix)
        tmatrix(ii).Nmax = [maxNmax1, maxNmax2];
      end

      % Apply translation to the beam
      if ~isempty(p.Results.position)

        % Requires scattered beam, convert if needed
        if ~strcmpi(tmatrix(1).type, 'scattered')
          maxNmax2 = min(maxNmax2, tbeam.Nmax);
          for ii = 1:numel(tmatrix)
            tmatrix(ii).type = 'scattered';
            tmatrix(ii).Nmax = [maxNmax1, maxNmax2];
          end
        end

        % Apply translation
        % We need Nmax+1 terms for the force calculation
        tbeam = tbeam.translateXyz(p.Results.position, 'Nmax', maxNmax2+1);
      end

      % Apply rotation to the beam
      rbeam = tbeam;
      if ~isempty(p.Results.rotation)
        [rbeam, D] = rbeam.rotate(p.Results.rotation, ...
            'Nmax', maxNmax1);
      end

      % Ensure the Nmax for the inner dimension matches
      if strcmpi(tmatrix(1).type, 'scattered')
        % T-matrix is already done
        rbeam = rbeam.setNmax(maxNmax2, 'powerloss', 'ignore');
      else
        for ii = 1:numel(tmatrix)
          tmatrix(ii) = tmatrix(ii).setNmax([maxNmax1, rbeam.Nmax], ...
              'powerloss', 'ignore');
        end
        if ~strcmpi(tmatrix(1).type, 'internal')
          ott.warning('ott:Bsc:scatter', ...
              'It may be more optimal to use a scattered T-matrix');
        end
      end

      % Calculate the resulting beams
      sbeam = ott.beam.abstract.Empty();
      for ii = 1:numel(tmatrix)
        sbeam = [sbeam, ott.beam.vswf.Scattered(tmatrix(ii).data * rbeam)];
      end
      
      % Store the incident beam and T-matrices if requested
      if p.Results.store_beam
        sbeam.incident_beam = tbeam;
      end
      if p.Results.store_tmatrix
        sbeam.tmatrix = tmatrix;
      end

      % Apply the inverse rotation
      if ~isempty(p.Results.rotation)

        % This seems to take a long time
        %sbeam = sbeam.rotate('wigner', D');

        sbeam = sbeam.rotate(inv(p.Results.rotation));
      end

      % Assign a type to the resulting beam
      switch tmatrix(1).type
        case 'total'
          sbeam.type = 'total';
          sbeam.basis = 'regular';
        case 'scattered'
          sbeam.type = 'scattered';
          sbeam.basis = 'outgoing';
        case 'internal'
          sbeam.type = 'internal';
          sbeam.basis = 'regular';

          % Wavelength has changed, update it
          sbeam = sbeam.setWavenumber(tmatrix(1).wavenumber_particle, 'medium');

        otherwise
          error('Unrecognized T-matrix type');
      end
    end
  end

  methods
    function beam = Scattered(varargin)
      % Calculate the scattered beam from a beam and T-matrix
      %
      % Usage
      %   beam = Scattered(...) construct an empty scattered beam.
      %
      %   beam = Scattered(a, b, ...) construct beam from a/b coefficients.
      %
      %   beam = Scattered(bsc, ...) specify an existing Bsc beam.
      %
      % Optional named arguments
      %   - incident_beam (Bsc) -- Incident beam object.  Default ``[]``.
      %
      %   - tmatrix (Tmatrix) -- T-matrix describing scattering.
      %     Default ``[]``.
      %
      %   - type (enum) -- Type of scattered beam.
      %     Type can be one of 'internal', 'total' or 'scattered'.
      %     Default: ``'total'``.
      %
      %   - like -- Another beam object to use for default parameters.
      %     Used for default 'type' and passed to base class.
      %     Default: ``[]``.
      %
      % Unmatched parameters are passed to the base class.

      % Setup input parser
      p = inputParser;
      p.KeepUnmatched = true;
      p.addOptional('a', [], @(x) isnumeric(x) || isa(x, 'ott.beam.vswf.Bsc'));
      p.addOptional('b', [], @isnumeric);
      p.addParameter('incident_beam', []);
      p.addParameter('tmatrix', []);
      p.addParameter('type', []);
      p.addParameter('like', []);
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Get default value for total
      default_type = 'total';
      if ~isempty(p.Results.like)
        if isa(p.Results.like, 'ott.beam.abstract.Scattered')
          default_type = p.Results.like.type;
        end
      end

      % Get Bsc arguments
      bsc_coeffs = {};
      if ~isempty(p.Results.a)
        bsc_coeffs = [bsc_coeffs, {p.Results.a}];
        if ~isempty(p.Results.b)
          bsc_coeffs = [bsc_coeffs, {p.Results.b}];
        end
      end

      % Get type (for Scattered constructor)
      if isempty(p.Results.type)
        type = default_type;
      else
        type = p.Results.type;
      end

      % Call base class
      beam = beam@ott.beam.Scattered(type);
      beam = beam@ott.beam.vswf.Bsc(bsc_coeffs{:}, ...
          'like', p.Results.like, unmatched{:});

      beam.incident_beam = p.Results.incident_beam;
      beam.tmatrix = p.Results.tmatrix;
    end

    function [sbeam, tbeam] = calculateScatteredBeam(tbeam, ...
        tmatrix, varargin)
      % Calculate the scattered beam
      %
      % Usage
      %   [sbeam, tbeam] = beam.calculateScatteredBeam(tmatrix, ...)
      %   Calcualtes the scattered beam ``sbeam`` and the translated
      %   and possibly rotated beam ``tbeam``.
      %
      % For optional arguments, see Bsc.scatter or Scattered constructor.

      p = inputParser;
      p.addParameter('position', []);
      p.addParameter('rotation', []);
      p.addParameter('combine', []);
      p.parse(varargin{:});

      % Apply combine parameter
      assert(any(strcmpi(p.Results.combine, {'coherent', 'incoherent'})) ...
          || isempty(p.Results.combine), ...
          'combine must be empty, ''cohernet'' or ''incoherent''');
      if strcmpi(p.Results.combine, 'coherent')
        tbeam = sum(tbeam);
      end

      % Check if we need to calculate multiple scatters
      if size(p.Results.position, 2) > 1 || size(p.Results.rotation, 2) > 3
        [sbeam, tbeam] = ott.utils.prxfun(...
            @(varargin) tbeam.scatter(tmatrix, varargin{:}), 1, ...
            'position', p.Results.position, ...
            'rotation', p.Results.rotation, ...
            'zeros', @(x) repmat(tbeam, x));

        % Join outputs (convert from beam arrays to combined beams)
        sbeam = sbeam(1).append(sbeam(2:end));
        tbeam = tbeam(1).append(tbeam(2:end));

        return;
      end

      % Determine the maximum tmatrix.Nmax(2) and check type
      maxNmax1 = 0;
      maxNmax2 = 0;
      tType = tmatrix(1).type;
      for ii = 1:numel(tmatrix)
        maxNmax1 = max(maxNmax1, tmatrix(ii).Nmax(1));
        maxNmax2 = max(maxNmax2, tmatrix(ii).Nmax(2));
        if ~strcmpi(tmatrix(ii).type, tType)
          error('T-matrices must be same type');
        end
      end

      % If the T is scattered, we can save time by throwing away columns
      if strcmpi(tmatrix(1).type, 'scattered')
        maxNmax2 = min(maxNmax2, tbeam.Nmax);
      end

      % Ensure all T-matrices are the same size
      for ii = 1:numel(tmatrix)
        tmatrix(ii).Nmax = [maxNmax1, maxNmax2];
      end

      % Apply translation to the beam
      if ~isempty(p.Results.position)

        % Requires scattered beam, convert if needed
        if ~strcmpi(tmatrix(1).type, 'scattered')
          maxNmax2 = min(maxNmax2, tbeam.Nmax);
          for ii = 1:numel(tmatrix)
            tmatrix(ii).type = 'scattered';
            tmatrix(ii).Nmax = [maxNmax1, maxNmax2];
          end
        end

        % Apply translation
        % We need Nmax+1 terms for the force calculation
        tbeam = tbeam.translateXyz(p.Results.position, 'Nmax', maxNmax2+1);
      end

      % Apply rotation to the beam
      rbeam = tbeam;
      if ~isempty(p.Results.rotation)
        [rbeam, D] = rbeam.rotate(p.Results.rotation, ...
            'Nmax', maxNmax1);
      end

      % Ensure the Nmax for the inner dimension matches
      if strcmpi(tmatrix(1).type, 'scattered')
        % T-matrix is already done
        rbeam = rbeam.set_Nmax(maxNmax2, 'powerloss', 'ignore');
      else
        for ii = 1:numel(tmatrix)
          tmatrix(ii) = tmatrix(ii).set_Nmax([maxNmax1, rbeam.Nmax], ...
              'powerloss', 'ignore');
        end
        if ~strcmpi(tmatrix(1).type, 'internal')
          ott.warning('ott:Bsc:scatter', ...
              'It may be more optimal to use a scattered T-matrix');
        end
      end

      % Calculate the resulting beams
      sbeam = ott.beam.vswf.Bsc();
      for ii = 1:numel(tmatrix)
        sbeam = [sbeam, tmatrix(ii).data * rbeam];
      end

      % Apply the inverse rotation
      if ~isempty(p.Results.rotation)

        % This seems to take a long time
        %sbeam = sbeam.rotate('wigner', D');

        sbeam = sbeam.rotate(inv(p.Results.rotation));
      end

      % Assign a type to the resulting beam
      switch tmatrix(1).type
        case 'total'
          sbeam.typeInternal = 'total';
          sbeam.basis = 'regular';
        case 'scattered'
          sbeam.typeInternal = 'scattered';
          sbeam.basis = 'outgoing';
        case 'internal'
          sbeam.typeInternal = 'internal';
          sbeam.basis = 'regular';

          % Wavelength has changed, update it
          sbeam.wavenumber = tmatrix(1).wavenumber_particle;

        otherwise
          error('Unrecognized T-matrix type');
      end
    end
  end

  methods % Getters/setters
    function beam = set.tmatrix(beam, val)
      assert(isempty(val) || isa(val, 'ott.scat.vswf.Tmatrix'), ...
          'tmatrix must be a valid ott.scat.vswf.Tmatrix');
      beam.tmatrix = val;
    end
  end
end
