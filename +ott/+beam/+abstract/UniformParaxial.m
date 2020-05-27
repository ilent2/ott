classdef UniformParaxial < ott.beam.abstract.CastNearfield ...
    & ott.beam.properties.Material ...
    & ott.beam.properties.FarfieldMapping
% A beam with a uniform paraxial far-field.
% Inherits from :class:`CastNearfield`
% and :class:`ott.beam.properties.Material`
% and :class:`ott.beam.properties.FarfieldMapping`.
%
% Properties
%   - polarisation
%   - mapping       -- Mapping for far-field to paraxial conversion
%   - direction     -- Mask direction (positive or negative hemisphere)
%
% All casts inherited from base.

  properties
    polarisation
  end

  properties (Dependent)
    power
  end

  methods (Static)
    function args = likeProperties(other, args)
      % Construct an array of like-properties
      if isa(other, 'ott.beam.abstract.UniformFarfield')
        args = ott.utils.addDefaultParameter(...
            'polarisation', other.polarisation, args);
      end
      args = ott.beam.properties.Material.likeProperties(other, args);
      args = ott.beam.properties.FarfieldMapping.likeProperties(other, args);
    end

    function beam = like(other, varargin)
      % Create a beam like another beam
      %
      % Usage
      %   beam = UniformFarfield.like(other, ...)
      %
      % See constructor for arguments.

      args = ott.beam.abstract.UniformParaxial.likeProperties(other, varargin);
      beam = ott.beam.abstract.UniformParaxial(args{:});
    end
  end

  methods
    function beam = UniformParaxial(varargin)
      % Construct a uniform paraxial far-field instance.
      %
      % Usage
      %   beam = UniformParaxial(polarisation, ...)
      %
      % Parameters
      %   - polarisation (2 numeric) -- [x; y] polarisation of
      %     paraxial far-field.
      %
      % Optional named arguments
      %   - mapping (enum) -- 'sintheta' or 'tantheta'.  Default: 'tantheta'
      %   - hemisphere (enum) -- 'pos' or 'neg'.  Default: 'pos'

      p = inputParser;
      p.addOptional('polarisation', [], @isnumeric);
      p.addParameter('mapping', 'tantheta');
      p.addParameter('hemisphere', 'pos');
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.properties.Material(unmatched{:});
      beam = beam@ott.beam.properties.FarfieldMapping(...
          'mapping', p.Results.mapping, 'hemisphere', p.Results.hemisphere);
      beam.polarisation = p.Results.polarisation;
    end
  end

  methods (Hidden)
    function E = efarfieldInternal(beam, rtp, varargin)
      % Uniform polarisation across paraxial plane

      % TODO: Implement
      error('Not yet implemented');
    end

    function H = hfarfieldInternal(beam, rtp, varargin)
      % Uniform polarisation across sphere

      % TODO: Implement
      error('Not yet implemented');
    end
  end

  methods % Getters/setters
    function beam = set.polarisation(beam, val)
      assert(isnumeric(val) && numel(val) == 2, ...
          'polarisation must be 2 element numeric');
      beam.polarisation = val(:);
    end

    function power = get.power(beam, val)
      % TODO: Implement
      error('Not yet implemented');
    end
  end
end

