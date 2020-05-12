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
