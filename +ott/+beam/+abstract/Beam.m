classdef (Abstract) Beam < ott.beam.Beam ...
    & matlab.mixin.Heterogeneous
% Base class for abstract beam descriptions.
% Inherits from :class:`ott.beam.Beam` and
% :class:`matlab.mixin.Hetrogeneous`.
%
% Abstract beams describe particular properties of beams but do not
% specify the method to calculate the fields.  Most classes declare
% no properties, instead, properties are declared via
% :class:`ott.beam.properties` sub-classes.
%
% These classes define casts to field calculation classes.  Abstract
% classes implement the same interface as full field calculation classes
% but the field calculation methods typically involve a cast to another
% beam type.
%
% Abstract beam can be formed into hetrogeneous arrays.  These arrays
% are assumed to be coherent and can be cast to beam implementations.
%
% Methods
%   - contains      -- Query if a array_type is contained in the array
%
% Supported casts
%   - Beam             -- Default beam cast for arrays, uses Coherent/Bsc
%   - vswf.Bsc         -- (Inherited)
%   - vswf.Pointmatch  -- (Inherited)
%   - vswf.NearfieldPm -- (Inherited)
%   - vswf.FarfieldPm  -- (Inherited)
%   - Array
%   - Coherent
%   - Incoherent

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function beam = Beam(varargin)
      % Construct an abstract beam instance
      %
      % Usage
      %   beam = beam@ott.beam.abstract.Beam(...)
      %
      % For optional parameters, see :class:`ott.beam.properties.Beam`.

      beam = beam@ott.beam.Beam(varargin{:});
    end

    function b = contains(beam, array_type)
      % Query if a array_type is contained in the array.
      %
      % Usage
      %   b = beam.contains(array_type)
      %
      % Parameters
      %   - array_type (enum) -- An array type, must be one of
      %     'array', 'coherent' or 'incoherent'.

      assert(any(strcmpi(array_type, {'array', 'incoherent', 'coherent'})), ...
          'array_type must be ''array'', ''incoherent'' or ''coherent''');

      if numel(beam) > 1
        if strcmpi(array_type, 'coherent')
          b = true;
        else
          b = false;
          for ii = 1:numel(beam)
            b = b | beam(ii).contains(array_type);
          end
        end
      else
        % Other cases handled by ott.beam.abstract.Array etc.
        b = false;
      end
    end

    function beam = ott.beam.Beam(beam, varargin)
      % Cast an array of beams to a coherent array of beams
      %
      % This method is called for hetrogeneous arrays.

      assert(isa(beam, 'ott.beam.abstract.Beam'), ...
          'First argument must be abstract.Beam');

      % TODO: Cast to coherent or array depending on contents
      %   Maybe warn if not coherent?

      if numel(beam) > 1
        beam = ott.beam.Coherent(beam);
      else
        beam = ott.beam.vswf.Bsc(beam);
      end
    end

    function beam = ott.beam.Array(beam, varargin)
      % Construct an array of beam objects

      assert(isa(beam, 'ott.beam.abstract.Beam'), ...
          'First argument must be abstract.Beam');

      beam_array = ott.beam.Array(size(beam));

      for ii = numel(beam)
        beam_array(ii) = ott.beam.Beam(beam(ii));
      end
    end

    function beam_array = ott.beam.Coherent(beam, varargin)
      % Generate incoherent array of beams

      assert(isa(beam, 'ott.beam.abstract.Beam'), ...
          'First argument must be abstract.Beam');

      beam_array = ott.beam.Coherent(size(beam));

      for ii = numel(beam)
        beam_array(ii) = ott.beam.Beam(beam(ii));
      end
    end

    function beam_array = ott.beam.Incoherent(beam, varargin)
      % Generate incoherent array of dipoles

      assert(isa(beam, 'ott.beam.abstract.Beam'), ...
          'First argument must be abstract.Beam');

      beam_array = ott.beam.Incoherent(size(beam));

      for ii = numel(beam)
        beam_array(ii) = ott.beam.Beam(beam(ii));
      end
    end
  end

  methods (Access=protected)
    function beam = castHelper(cast, beam, varargin)
      % Helper for casts
      %
      % This cast converts native arrays to ott.beam.Array

      assert(isa(beam, 'ott.beam.abstract.Beam'), ...
          'First argument must be a abstract.Beam');

      ott.utils.nargoutCheck(beam, nargout);

      if numel(beam) > 1
        oldbeam = beam;
        beam = ott.beam.Array('coherent', size(beam));
        for ii = 1:numel(oldbeam)
          beam(ii) = cast(beam, varargin{:});
        end
      else
        beam = cast(beam, varargin{:});
      end
    end

    function beam = castArrayHelper(cast, beam, varargin)
      % Helper for array casts
      %
      % Checks the argument and applies the cast

      assert(isa(beam, 'ott.beam.abstract.Beam'), ...
          'First argument must be a abstract.Beam');
      ott.utils.nargoutCheck(beam, nargout);

      beam = cast(beam, varargin{:});
    end

    function beam = castAbstractHelper(cast, beam, varargin)
      % Helper for casts to other abstract beams
      %
      % This cast doesn't remove native arrays.

      assert(isa(beam, 'ott.beam.abstract.Beam'), ...
          'First argument must be a abstract.Beam');

      beam = cast(beam, varargin{:});
    end
  end

  methods (Static, Sealed, Access = protected)
    function default_object = getDefaultScalarElement
      % Default object for a Shape array is an empty shape
      default_object = ott.beam.abstract.Empty();
    end
  end
end

