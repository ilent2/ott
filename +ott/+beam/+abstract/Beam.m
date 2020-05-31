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
%   - plus          -- Form coherent array
%   - minus         -- Form coherent array
%   - uminus        -- Negative fields
%
% Supported casts
%   - Beam             -- Default beam cast, uses Bsc
%   - vswf.Bsc         -- (Inherited)
%   - vswf.Pointmatch  -- (Inherited)
%   - vswf.NearfieldPm -- (Inherited)
%   - vswf.FarfieldPm  -- (Inherited)
%   - Array            -- (Sealed) Cast hetrogeneous array to Array

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

    function beam = ott.beam.Beam(varargin)
      % Cast to Bsc
      beam = ott.beam.vswf.Bsc(varargin{:});
    end
  end

  methods (Sealed)
    function beam = plus(a, b)
      % Form coherent array from pair of beams

      if numel(a) > 1 || numel(b) > 1
        beam = ott.beam.Array([1, 2]);
        beam(1) = a;
        beam(2) = b;
      else
        beam = [a, b];
      end
    end

    function beam = minus(a, b)
      % Form coherent array from pair of beams

      beam = a + (-b);
    end

    function beam = uminus(beam)
      % Negate the fields of the beam
      beam = beam.arrayCast(@ott.beam.abstract.Negative);
    end

    function beam = arrayCast(beam, cast, varargin)
      % Cast array of beams to specified type.
      %
      % This function applies the cast to each beam in the array.
      % If the input is not an array, doesn't return an array.
      % If all results are abstract beams, keeps the array as a
      % hetrogeneous array, otherwise tries to create a coherent array
      % using defaultArrayType.
      %
      % If a coherent array can't be created, raises a
      % 'ott:beam:abstract:Beam:not_coherent_array' warning and creates
      % a generic arrays instead.
      %
      % Usage
      %   beam = beam.arrayCast(cast, ...)
      %
      % Example
      %   beam = beam.arrayCast(@ott.beam.Beam)
      %   Casts all elements to Beam instances.

      allAbstract = true;
      allCoherent = true;

      beam_array = cell(size(beam));
      for ii = 1:numel(beam)
        beam_array{ii} = cast(beam(ii), varargin{:});
        allAbstract = allAbstract ...
            & isa(beam_array{ii}, 'ott.beam.abstract.Beam');
        allCoherent = allCoherent & ~beam_array{ii}.contains('incoherent');
      end

      if allAbstract
        beam = reshape([beam_array{:}], size(beam_array));
      elseif numel(beam_array) == 1
        beam = beam_array{1};
      elseif allCoherent
        arrayTypes = cellfun(@(b) b.defaultArrayType('coherent'), ...
            beam_array, 'UniformOutput', false);
        if all(cellfun(@(a) isequal(a, arrayTypes{1}), arrayTypes))
          beam = arrayTypes{1}(beam_array);
        else
          beam = ott.beam.Array('coherent', beam_array);
        end
      else
        warning('ott:beam:abstract:Beam:not_coherent_array', ...
            'Unable to cast to coherent array, casting to generic array');

        arrayTypes = cellfun(@(b) b.defaultArrayType('array'), ...
            beam_array, 'UniformOutput', false);
        if all(cellfun(@(a) isequal(a, arrayTypes{1}), arrayTypes))
          beam = arrayTypes{1}(beam_array);
        else
          beam = ott.beam.Array('array', beam_array);
        end
      end
    end

    function b = arrayContains(beam, array_type)
      % Apply contains method to each beam array element.
      %
      % Usage
      %   b = beams.arrayContains(array_type)
      %
      % See :meth:`contains` for additional information.

      assert(any(strcmpi(array_type, {'array', 'incoherent', 'coherent'})), ...
          'array_type must be ''array'', ''incoherent'' or ''coherent''');

      if strcmpi(array_type, 'coherent')
        b = numel(beam) > 1;
      else
        b = false;
        for ii = 1:numel(beam)
          b = b | beam(ii).contains(array_type);
        end
      end
    end

    function beam = ott.beam.Array(beam, varargin)
      % Construct a new ott.beam.Array with contents of this array
      %
      % Usage
      %   beam = ott.beam.Array(beam, ...)
      %   Additional arguments are passed to ott.beam.Array.

      assert(isa(beam, 'ott.beam.abstract.Beam'), ...
          'First argument must be a abstract.Beam');
      ott.utils.nargoutCheck(beam, nargout);

      allCoherent = ~beam.arrayContains('incoherent') ...
          && ~beam.arrayContains('array');
      if allCoherent
        array_type = 'coherent';
      else
        array_type = 'incoherent';
      end

      args = ott.beam.Array.likeProperties(beam(1), varargin);
      beam = ott.beam.Array('beams', num2cell(beam), ...
          'array_type', array_type, args{:});
    end
  end

  methods (Access=protected)
    function beam = castHelper(cast, beam, varargin)
      % Helper for casts.  See also arrayCast for hetrogeneous arrays.

      assert(isa(beam, 'ott.beam.abstract.Beam'), ...
          'First argument must be a abstract.Beam');

      ott.utils.nargoutCheck(beam, nargout);

      beam = cast(beam, varargin{:});
    end
  end

  methods (Static, Sealed, Access = protected)
    function default_object = getDefaultScalarElement
      % Default object for a Shape array is an empty shape
      default_object = ott.beam.abstract.Empty();
    end

    function cobj = convertObject(~, objToConvert)
      % Apply the objects cast
      cobj = ott.beam.abstract.Beam(objToConvert);
    end
  end
end

