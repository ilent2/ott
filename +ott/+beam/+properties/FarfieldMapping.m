classdef FarfieldMapping
% Declares a far-field mapping property for beams
%
% Properties
%   - mapping       -- Far-field mapping (sintheta or tantheta)
%   - hemisphere    -- Far-field hemisphere (pos or neg)

  properties
    mapping
    hemisphere
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
      % Construct far-field mapping property
      %
      % Usage
      %   beam = beam@ott.beam.utils.FarfieldMapping(mapping, hemisphere, ...)
      %
      % Named arguments
      %   - mapping (enum) -- Initial mapping value.
      %   - hemisphere (enum) -- Hemisphere for paraxial field.

      p = inputParser;
      p.addOptional('mapping', [], ...
          @(x) any(strcmpi(x, {'sintheta', 'tantheta'})));
      p.addOptional('hemisphere', [], ...
          @(x) any(strcmpi(x, {'pos', 'neg'})));
      p.parse(varargin{:});

      beam.mapping = p.Results.mapping;
      beam.hemisphere = p.Results.hemisphere;
    end
  end

  methods % Getters/setters
    function beam = set.mapping(beam, val)
      assert(any(strcmpi(val, {'sintheta', 'tantheta'})), ...
          'mapping must be ''sintheta'' or ''tantheta''');
      beam.mapping = val;
    end

    function beam = set.hermisphere(beam, val)
      assert(any(strcmpi(val, {'pos', 'neg'})), ...
          'hemisphere must be pos or neg');
      beam.hemisphere = val;
    end
  end
end
