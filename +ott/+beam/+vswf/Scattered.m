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
      %   - array_type (enum) -- Array type for scattered beam with respect
      %     to rotation and position arguments.  Default: ``'array'``.
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
      p.addParameter('array_type', 'array');
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
        tbeam.array_type = p.Results.array_type;
      end

      % Check if we need to calculate multiple scatters
      if size(p.Results.position, 2) > 1 || size(p.Results.rotation, 2) > 3
        [sbeam, tbeam] = ott.utils.prxfun(...
            @(varargin) tmatrix.scatter(tbeam, varargin{:}), 1, ...
            'position', p.Results.position, ...
            'rotation', p.Results.rotation, ...
            'zeros', {@(sz) ott.beam.vswf.Scattered.empty(sz, 'like', tbeam), ...
            @(sz) repmat(tbeam, sz)});
        sbeam.array_type = p.Results.array_type;
        tbeam.array_type = p.Results.array_type;
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
      sbeam = ott.beam.vswf.Scattered.empty(numel(tmatrix));
      for ii = 1:numel(tmatrix)
        sbeam(ii) = ott.beam.vswf.Scattered(tmatrix(ii).data * rbeam);
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
    
    function beam = empty(varargin)
      % Construct a new empty scattered beam
      %
      % Usage
      %   beam = Scattered.empty()
      %
      %   beam = Scattered.empty(sz, ...)
      %
      % Parameters
      %   - sz (numeric) -- Number of beams in empty array.
      %     Can either be [1, num] or [num].
      %
      % Additional arguments are passed to constructor.
      
      p = inputParser;
      p.addOptional('sz', 0, @isnumeric);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);
      
      sz = p.Results.sz;
      
      assert(isnumeric(sz) && (sz(1) == 1 || isscalar(sz)), ...
        'sz must be numeric scalar or first element must be 1');
      
      if numel(sz) == 2
        sz = sz(2);
      end
      
      a = zeros(0, sz);
      b = zeros(0, sz);
      
      beam = ott.beam.vswf.Scattered(a, b, unmatched{:});
    end
  end

  methods
    function bsc = Scattered(varargin)
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
      %   - incident_beam (Bsc) -- Incident beam object.
      %     Default ``ott.beam.vswf.Bsc.empty()``.
      %
      %   - tmatrix (Tmatrix) -- T-matrix describing scattering.
      %     Default ``ott.scat.vswf.Tmatrix.empty()``.
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
      p.addParameter('incident_beam', ott.beam.vswf.Bsc.empty());
      p.addParameter('tmatrix', ott.scat.vswf.Tmatrix.empty());
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
      bsc = bsc@ott.beam.Scattered(type);
      bsc = bsc@ott.beam.vswf.Bsc(bsc_coeffs{:}, ...
          'like', p.Results.like, unmatched{:});

      bsc.incident_beam = p.Results.incident_beam;
      bsc.tmatrix = p.Results.tmatrix;
    end
  end
  
  methods (Hidden)
    function beam = catInternal(dim, beam, varargin)
      % Concatenate beams
      
      % Combine a/b coefficients (Bsc)
      beam = catInternal@ott.beam.vswf.Bsc(dim, beam, varargin{:});
      
      % Combine incident_beam and tmatrix properties
      other_tmatrix = {};
      other_ibeam = {};
      for ii = 1:length(varargin)
        other_tmatrix{ii} = varargin{ii}.tmatrix;
        other_ibeam{ii} = varargin{ii}.incident_beam;
      end

      beam.tmatrix = cat(dim, beam.tmatrix, other_tmatrix{:});
      beam.incident_beam = cat(dim, beam.incident_beam, other_ibeam{:});
    end
    
    function beam = subsrefInternal(beam, subs)
      % Get the subscripted beam
      
      % Ref a/b coefficients (Bsc)
      beam = subsrefInternal@ott.beam.vswf.Bsc(beam, subs);
      
      if numel(subs) > 1
        if subs{1} == 1 || strcmpi(subs{1}, ':')
          subs = subs(2:end);
        end
        assert(numel(subs) == 1, 'Only 1-D indexing supported for now');
      end

      % Ref incident_beam and tmatrix properties
      beam.tmatrix = beam.tmatrix(subs{:});
      beam.incident_beam = beam.incident_beam(subs{:});
    end

    function beam = subsasgnInternal(beam, subs, rem, other)
      % Assign to the subscripted beam
      
      % Ref a/b coefficients (Bsc)
      beam = subsasgnInternal@ott.beam.vswf.Bsc(beam, subs, rem, other);

      if numel(subs) > 1
        if subs(1) == 1
          subs = subs(2:end);
        end
        assert(numel(subs) == 1, 'Only 1-D indexing supported for now');
      end

      assert(isempty(rem), 'Assignment to parts of beams not supported');
      
      % Ref incident_beam and tmatrix properties
      if isempty(other)
        % Delete data
        beam.tmatrix(subs{:}) = ott.scat.vswf.Tmatrix.empty();
        beam.incident_beam(subs{:}) = [];
      else
        beam.tmatrix(subs{:}) = other.tmatrix;
        beam.incident_beam(subs{:}) = other.incident_beam;
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
