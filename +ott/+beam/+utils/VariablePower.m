classdef (Abstract) VariablePower < ott.beam.Properties
% Adds a variable power property to a Beam.
% Inherits from :class:`ott.beam.Properties`.
%
% Properties (Hidden)
%   - powerInternal       -- The power property
%
% Methods (Hidden)
%   - getBeamPower        -- Get the internal power value
%   - setBeamPower        -- Set the internal power value

  properties (Hidden)
    powerInternal
  end

  methods
    function beam = VariablePower(power, varargin)
      % Construct a beam specifying the power
      %
      % Usage
      %   beam = VariablePower(power, ...)
      %
      % For optional arguments, see :class:`Properties`.

      beam = beam@ott.beam.Properties(varargin{:});
      beam.powerInternal = power;
    end
  end

  methods (Hidden)
    function val = getBeamPower(beam)
      % Get the internal power value
      val = beam.powerInternal;
    end
    function beam = setBeamPower(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'power must be numeric scalar');
      beam.powerInternal = val;
    end
  end
end
