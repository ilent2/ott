classdef Negative < ott.beam.abstract.Beam
% Negative beam for coherent addition/subtraction
%
% Negates the fields of the beam.
%
% Properties
%   - data      -- Internal beam data to negate

  properties
    data       % Internal beam data to negate
  end

   properties (Dependent)
    omega     % Beam optical frequency
    medium    % Beam optical medium
    power     % Beam power
  end

  methods
    function beam = Negative(other)
      % Construct a negated beam
      %
      % Usage
      %   beam = Negative(beam)

      beam.data = other;
    end

    function b = contains(beam, array_type)
      % Query if a array_type is contained in the array.
      %
      % Applies contains to the data.
      %
      % Usage
      %   b = beam.contains(array_type)
      %
      % Parameters
      %   - array_type (enum) -- An array type, must be one of
      %     'array', 'coherent' or 'incoherent'.

      b = beam.data.contains(array_type);
    end
  end

  methods (Hidden)
    function E = efieldInternal(beam, varargin)
      % Negate fields
      E = beam.data.efieldInternal(xyz, varargin{:});
      E = -E;
    end

    function H = hfieldInternal(beam, varargin)
      % Negate fields
      H = beam.data.hfieldInternal(xyz, varargin{:});
      H = -H;
    end

    function E = efarfieldInternal(beam, varargin)
      % Negate fields
      E = beam.data.efarfieldInternal(xyz, varargin{:});
      E = -E;
    end

    function H = hfarfieldInternal(beam, varargin)
      % Negate fields
      H = beam.data.hfarfieldInternal(xyz, varargin{:});
      H = -H;
    end
  end

  methods % Getters/setters
    function omega = get.omega(beam)
      omega = beam.data.omega;
    end
    function beam = set.omega(beam, val)
      beam.beams.omega = val;
    end

    function medium = get.medium(beam)
      medium = beam.data.medium;
    end
    function beam = set.medium(beam, val)
      beam.data.medium = val;
    end

    function power = get.power(beam)
      power = beam.data.power;
    end
    function beam = set.power(beam, val)
      beam.data.power = val;
    end
  end
end
