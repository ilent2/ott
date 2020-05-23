classdef InceGaussian < ott.beam.properties.InceGaussian ...
    ott.beam.abstract.Gaussian
% Abstract representation of a Ince-Gaussian beam
% Inherits from :class:`Gaussian`.
%
% Supported casts
%   - Beam            -- Default Beam cast, uses vswf.InceGaussian
%   - vswf.Bsc        -- Default Bsc cast, uses vswf.InceGaussian
%   - vswf.InceGaussian
%   - vswf.LaguerreGaussian     -- Only for Gaussian
%   - vswf.Gaussian             -- Only for Gaussian
%   - vswf.HermiteGaussian      -- Only for Gaussian
%   - paraxial.LaguerreGaussian -- Only for Gaussian
%   - paraxial.Gaussian         -- Only for Gaussian
%   - paraxial.HermiteGaussian  -- Only for Gaussian
%
% See also :class:`ott.beam.properties.InceGaussian` for properties.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function beam = InceGaussian(varargin)
      % Construct a new Abstract Ince-Gaussian beam
      %
      % Usage
      %   beam = InceGaussian(waist, lmode, porder, parity, ellipticity, ...)
      %   Parameters can also be passed as named arguments.
      %
      % Parameters
      %   - waist (numeric) -- Beam waist
      %   - lmode (numeric) -- Azimuthal mode number.
      %   - porder (numeric) -- Paraxial mode number.
      %   - parity (enum) -- Either 'even' or 'odd'.
      %   - ellipticity (numeric) -- Ellipticity of coordinates.
      %
      % Optional named arguments
      %   - power (numeric) -- beam power.  Default: ``1.0``.
      %
      % For optional parameters, see :class:`Properties`.

      args = ott.utils.addDefaultParameter('power', 1.0, varargin);
      beam = beam@ott.beam.properties.InceGaussian(args{:});
    end
  end
end

