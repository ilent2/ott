classdef HermiteGaussian < ott.beam.properties.HermiteGaussian ...
    & ott.beam.abstract.Beam ...
    & ott.beam.utils.VariablePower
% Abstract representation of a Hermite-Gaussian beam
% Inherits from :class:`Gaussian` and
% :class:`ott.beam.properties.LaguerreGaussian`.
%
% Supported casts
%   - Beam            -- Default Beam cast, uses paraxial.HermiteGaussian
%   - vswf.Bsc        -- Default Bsc cast, uses vswf.HermiteGaussian
%   - vswf.HermiteGaussian
%   - vswf.LaguerreGaussian     -- Only for HG00
%   - vswf.Gaussian             -- Only for HG00
%   - vswf.InceGaussian         -- Only for HG00
%   - paraxial.HermiteGaussian
%   - paraxial.LaguerreGaussian -- Only for HG00
%   - paraxial.Gaussian         -- Only for HG00
%
% See also :class:`ott.beam.properties.HermiteGaussian` for properties.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Static)
    function args = likeProperties(other, args)
      % Construct an array of like-properties
      args = ott.beam.utils.VariablePower.likeProperties(other, args);
      args = ott.beam.properties.HermiteGaussian.likeProperties(other, args);
    end

    function beam = like(other, varargin)
      % Create a beam like another beam
      %
      % Usage
      %   beam = HermiteGaussian.like(other, ...)
      %
      % See constructor for arguments.

      args = ott.beam.abstract.HermiteGaussian.likeProperties(...
          other, varargin);
      beam = ott.beam.abstract.HermiteGaussian(args{:});
    end
  end

  methods
    function beam = HermiteGaussian(varargin)
      % Construct a new Abstract Hermite-Gaussian beam
      %
      % Usage
      %   beam = HermiteGaussian(waist, mmode, nmode, ...)
      %
      % Parameters
      %   - waist (numeric) -- Beam waist
      %   - mmode (integer) -- Mode number
      %   - nmode (integer) -- Mode number
      %
      % For optional parameters, see :class:`Properties`.

      args = ott.utils.addDefaultParameter('power', 1.0, varargin);
      beam = beam@ott.beam.properties.Gaussian(args{:});
    end
  end
end

