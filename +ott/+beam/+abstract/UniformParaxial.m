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
      if isa(other, 'ott.beam.abstract.UniformParaxial')
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
      %   - mapping (enum) -- 'sin' or 'tan' or 'theta'.  Default: 'tan'
      %   - hemisphere (enum) -- 'pos' or 'neg'.  Default: 'pos'

      p = inputParser;
      p.addOptional('polarisation', [], @isnumeric);
      p.addParameter('mapping', 'tan');
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

      [~, dA] = beam.farfield2paraxial(rtp, ...
          'mapping', beam.mapping, 'direction', beam.hemisphere);

      Ex = beam.polarisation(1);
      Ey = beam.polarisation(2);
      Etheta = -Ex .* cos(phi) - Ey .* sin(phi);
      Ephi = -Ex .* sin(phi) + Ey .* cos(phi);

      Etheta = Etheta .* dA;
      Ephi = Ephi .* dA;

      Ertp = [0*Etheta; Etheta; Ephi];

      E = ott.utils.FieldVector(rtp, Ertp, 'spherical');
    end

    function H = hfarfieldInternal(beam, rtp, varargin)
      % Uniform polarisation across sphere

      [~, dA] = beam.farfield2paraxial(rtp, ...
          'mapping', beam.mapping, 'direction', beam.hemisphere);

      % TODO: Unit conversion factor?
      Hy = beam.polarisation(1);
      Hx = beam.polarisation(2);
      phi = rtp(3, :);
      Htheta = -Hx .* cos(phi) - Hy .* sin(phi);
      Hphi = -Hx .* sin(phi) + Hy .* cos(phi);

      Htheta = Htheta .* dA;
      Hphi = Hphi .* dA;

      Hrtp = [0*Htheta; Htheta; Hphi];

      H = ott.utils.FieldVector(rtp, Hrtp, 'spherical');
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

