classdef Gaussian < ott.beam.properties.Gaussian ...
    & ott.beam.abstract.Paraxial
% Describes the Beam casts for a Gaussian beam.
% Inherits from :class:`CastBoth`, :class:`ott.beam.properties.Gaussian`,
% and :class:`utils.VariablePower`.
%
% Casts
%   - vswf.Bsc            -- Cast to vswf.Gaussian
%   - paraxial.Paraxial   -- Cast to paraxial.Gaussian
%
% Additional casts inherited from base.
% See also :class:`ott.beam.properties.Gaussian` for properties.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Static)
    function args = likeProperties(other, args)
      % Construct an array of like-properties
      args = ott.beam.properties.VariablePower.likeProperties(other, args);
      args = ott.beam.properties.Gaussian.likeProperties(other, args);
    end

    function beam = like(other, varargin)
      % Create a beam like another beam
      %
      % Usage
      %   beam = Gaussian.like(other, ...)
      %
      % See constructor for arguments.

      args = ott.beam.abstract.Gaussian.likeProperties(other, varargin);
      beam = ott.beam.abstract.Gaussian(args{:});
    end
  end

  methods
    function beam = Gaussian(varargin)
      % Construct a new Abstract Gaussian beam
      %
      % Usage
      %   beam = Gaussian(waist, ...)
      %   Parameters can also be passed as named arguments.
      %
      % Parameters
      %   - waist (numeric) -- Beam waist [L]
      %
      % Optional named arguments
      %   - power (numeric) -- beam power.  Default: ``1.0``.
      %
      % See also :class:`ott.beam.properties.Gaussian` for properties.

      args = ott.utils.addDefaultParameter('power', 1.0, varargin);
      beam = beam@ott.beam.properties.Gaussian(args{:});
    end

    function beam = ott.beam.paraxial.Paraxial(varargin)
      % Construct a paraxial.Gaussian beam
      beam = ott.beam.paraxial.Gaussian(varargin{:});
    end

    function beam = ott.beam.vswf.Bsc(varargin)
      % Cast to VSWF Gaussian
      beam = ott.beam.vswf.Gaussian(varargin{:});
    end
  end
end

