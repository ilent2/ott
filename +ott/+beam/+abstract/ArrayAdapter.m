classdef ArrayAdapter < ott.beam.abstract.Beam
% Adapter for hetrogeneous arrays with ott.beam.Array elements

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    data
  end

  properties (Dependent)
    power
    medium
    omega
    array_type
    beams
  end

  methods
    function beam = ArrayAdapter(data)
      % Construct array adapter
      %
      % Usage
      %   beam = ArrayAdapter(data)

      beam.data = data;
    end

    function b = contains(beam, varargin)
      % Query if a array_type is contained in the array.
      %
      % Defers to the beam.data.contains method.
      %
      % Usage
      %   b = beam.contains(array_type)
      %
      % Parameters
      %   - array_type (enum) -- An array type, must be one of
      %     'array', 'coherent' or 'incoherent'.

      b = beam.data.contains(varargin{:});
    end

    function beam = ott.beam.Beam(beam)
      % Retrieve beam data
      beam = beam.data;
    end
  end

  methods (Hidden)
    function E = efieldInternal(beam, varargin)
      E = beam.data.efieldInternal(beam, varargin{:});
    end
    function E = hfieldInternal(beam, varargin)
      E = beam.data.hfieldInternal(beam, varargin{:});
    end
    function E = efarfieldInternal(beam, varargin)
      E = beam.data.efarfieldInternal(beam, varargin{:});
    end
    function E = hfarfieldInternal(beam, varargin)
      E = beam.data.hfarfieldInternal(beam, varargin{:});
    end
  end

  methods % Getters/setters
    function beam = set.data(beam, val)
      assert(isa(val, 'ott.beam.properties.ArrayType'), ...
          'data must be a Array');
      beam.data = val;
    end

    function v = get.beams(beam)
      v = beam.data.beams;
    end

    function v = get.array_type(beam)
      v = beam.data.array_type;
    end

    function v = get.power(beam)
      v = beam.data.power;
    end
    function v = get.medium(beam)
      v = beam.data.medium;
    end
    function v = get.omega(beam)
      v = beam.data.omega;
    end
  end
end
