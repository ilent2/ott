classdef InceGaussian < ott.beam.properties.InceGaussian ...
    & ott.beam.abstract.Paraxial
% Abstract representation of a Ince-Gaussian beam
% Inherits from :class:`Paraxial` and
% :class:`ott.beam.properties.InceGaussian`.
%
% Casts
%   - vswf.Bsc        -- Cast to vswf.InceGaussian
%   - vswf.InceGaussian
%
% Additional casts inherited from base.
% See also :class:`ott.beam.properties.InceGaussian` for properties.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Static)
    function args = likeProperties(other, args)
      % Construct an array of like-properties
      args = ott.beam.utils.VariablePower.likeProperties(other, args);
      args = ott.beam.properties.InceGaussian.likeProperties(other, args);
    end

    function beam = like(other, varargin)
      % Create a beam like another beam
      %
      % Usage
      %   beam = InceGaussian.like(other, ...)
      %
      % See constructor for arguments.

      args = ott.beam.abstract.InceGaussian.likeProperties(...
          other, varargin);
      beam = ott.beam.abstract.InceGaussian(args{:});
    end
  end

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

    function beam = ott.beam.vswf.Bsc(varargin)
      % Default Beam cast, uses vswf.InceGaussian
      beam = ott.beam.vswf.InceGaussian(varargin{:});
    end

    function beam = ott.beam.vswf.InceGaussian(beam, varargin)
      % Cast to vswf.InceGaussian

      assert(isa(beam, 'ott.beam.abstract.InceGaussian'), ...
          'first argument must be a abtract.InceGaussian');

      beam = castHelper(@ott.beam.vswf.InceGaussian.like, ...
          beam, 'lmode', beam.lmode, 'porder', beam.porder, ...
          'ellipticity', beam.ellipticity, 'parity', beam.parity, ...
          varargin{:});
    end
  end
end
