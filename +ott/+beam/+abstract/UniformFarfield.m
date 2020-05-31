classdef UniformFarfield < ott.beam.abstract.CastNearfield ...
    & ott.beam.properties.Material
% A beam with a uniform far-field.
% Inherits from :class:`Beam` and :class:`ott.beam.properties.Material`.
%
% Properties
%   - polarisation        - Far-field polarisation (theta/phi)
%
% All casts inherited from base.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    polarisation      % [theta; phi] polarisation of far-field
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
    end

    function beam = like(other, varargin)
      % Create a beam like another beam
      %
      % Usage
      %   beam = UniformFarfield.like(other, ...)
      %
      % See constructor for arguments.

      args = ott.beam.abstract.UniformFarfield.likeProperties(other, varargin);
      beam = ott.beam.abstract.UniformFarfield(args{:});
    end
  end

  methods
    function beam = UniformFarfield(varargin)
      % Construct a uniform far-field instance.
      %
      % Usage
      %   beam = UniformFarfield(polarisation, ...)
      %
      % Parameters
      %   - polarisation (2 numeric) -- [theta; phi] polarisation of
      %     far-field.

      p = inputParser;
      p.addOptional('polarisation', [], @isnumeric);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.properties.Material(unmatched{:});
      beam.polarisation = p.Results.polarisation;
    end
  end

  methods (Hidden)
    function E = efarfieldInternal(beam, rtp, varargin)
      % Uniform polarisation across sphere

      if size(rtp, 1) == 2
        rtp = [zeros(1, size(rtp, 2)); rtp];
      end

      E = ott.utils.FieldVector(rtp, ...
          repmat([0; beam.polarisation], 1, size(rtp, 2)), 'spherical');
    end

    function H = hfarfieldInternal(beam, rtp, varargin)
      % Uniform polarisation across sphere

      if size(rtp, 1) == 2
        rtp = [zeros(1, size(rtp, 2)); rtp];
      end

      % TODO: Missing conversion factor?
      H = ott.utils.FieldVector(rtp, ...
          repmat([0; flip(beam.polarisation)], 1, size(rtp, 2)), 'spherical');
    end
  end

  methods % Getters/setters
    function beam = set.polarisation(beam, val)
      assert(isnumeric(val) && numel(val) == 2, ...
          'polarisation must be 2 element numeric');
      beam.polarisation = val(:);
    end

    function power = get.power(beam)
      % TODO: Check this
      power = 4*pi*sum(abs(beam.polarisation).^2);
    end
    function beam = set.power(beam, val)
      % Scale polarisation
      beam.polarisation = beam.polarisation ...
          ./ abs(beam.polarisation) .* sqrt(val ./ (4*pi));
    end
  end
end
