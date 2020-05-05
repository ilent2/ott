classdef (Abstract) Vector2Polarisation < ott.beam.Properties
% Adds a 2-vector polarisation property to a beam.
% Inherits from :class:`ott.beam.Properties`.
%
% Properties
%   - polarisation        -- Polarisation 2-vector [x, y]
%
% Methods
%   - set.polarisation    -- Checks polarisation size is 2x1

  properties
    polarisation          % x and y polarisation 2-vector [x, y]
  end

  methods
    function beam = Vector2Polarisation(pol, varargin)
      % Construct a new Vector2Polarisation layer
      %
      % Usage
      %   beam = Vector2Polarisation(pol, ...)

      beam = beam@ott.beam.Properties(varargin{:});
      beam.polarisation = pol;
    end
  end

  methods % Getters/setters
    function beam = set.polarisation(beam, val)
      assert(isnumeric(val) && numel(val) == 2, ...
        'polarisation must be 2 element vector');
      beam.polarisation = val(:).';
    end
  end
end
