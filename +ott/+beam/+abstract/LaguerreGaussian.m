classdef LaguerreGaussian < ott.beam.properties.LaguerreGaussian ...
    & ott.beam.abstract.Gaussian
% Abstract representation of a Laguerre-Gaussian beam.
% Inherits from :class:`Gaussian` and
% :class:`ott.beam.properties.LaguerreGaussian`.
%
% Supported casts
%   - Beam            -- Default Beam cast, uses paraxial.LaguerreGaussian
%   - vswf.Bsc        -- Default Bsc cast, uses vswf.LaguerreGaussian
%   - vswf.LaguerreGaussian
%   - vswf.Gaussian             -- Only for LG00
%   - vswf.HermiteGaussian      -- Only for LG00
%   - vswf.InceGaussian         -- Only for LG00
%   - paraxial.LaguerreGaussian
%   - paraxial.Gaussian         -- Only for LG00
%   - paraxial.HermiteGaussian  -- Only for LG00
%
% See also :class:`ott.beam.properties.LaguerreGaussian` for properties.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function beam = LaguerreGaussian(varargin)
      % Construct a new Abstract Laguerre Gaussian beam
      %
      % Usage
      %   beam = LaguerreGaussian(waist, lmode, pmode, ...)
      %   Parameters can also be passed as named arguments.
      %
      % Parameters
      %   - waist (numeric)     -- Beam waist
      %   - lmode (integer)     -- Azimuthal LG mode
      %   - pmode (integer > 0) -- Radial LG mode
      %
      % Optional named arguments
      %   - power (numeric) -- beam power.  Default: ``1.0``.
      %
      % See :class:`ott.beam.properties.LaguerreGaussian`.

      args = ott.utils.addDefaultParameter('power', 1.0, varargin);
      beam = beam@ott.beam.properties.LaguerreGaussian(args{:});
    end
  end
end

