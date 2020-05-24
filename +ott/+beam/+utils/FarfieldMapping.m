classdef FarfieldMapping
% Declares a far-field mapping property for beams
%
% Properties
%   - mapping       -- Far-field mapping (sintheta or tantheta)

  properties
    mapping
  end

  methods (Static)
    function args = likeProperties(other, args)
      % Construct array of like-proeprties
      if isa(other, 'ott.beam.utils.FarfieldMapping')
        args = ott.utils.addDefaultParameter('mapping', other.mapping, args);
      end
    end
  end

  methods
    function beam = FarfieldMapping(varargin)
      % Costruct far-field mapping property
      %
      % Usage
      %   beam = beam@ott.beam.utils.FarfieldMapping(...)
      %
      % Optional parameters
      %   - mapping (enum) -- Initial mapping value.
      %     Default leaves mapping unassigned.

      p = inputParser;
      p.addOptional('mapping', [], ...
          @(x) any(strcmpi(x, {'sintheta', 'tantheta'})));
      p.parse(varargin{:});

      if isempty(p.UsingDefaults)
        beam.mapping = p.Results.mapping;
      end
    end
  end

  methods % Getters/setters
    function beam = set.mapping(beam, val)
      assert(any(strcmpi(val, {'sintheta', 'tantheta'})), ...
          'mapping must be ''sintheta'' or ''tantheta''');
      beam.mapping = val;
    end
  end
end
