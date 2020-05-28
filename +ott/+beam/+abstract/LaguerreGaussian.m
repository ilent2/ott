classdef LaguerreGaussian < ott.beam.properties.LaguerreGaussian ...
    & ott.beam.abstract.Paraxial
% Abstract representation of a Laguerre-Gaussian beam.
% Inherits from :class:`CastBoth` and
% :class:`ott.beam.properties.LaguerreGaussian`.
%
% Casts
%   - vswf.Bsc            -- Cast to vswf.LaguerreGaussian
%   - paraxial.Paraxial   -- Cast to paraxial.LaguerreGaussian
%   - paraxial.LaguerreGaussian
%   - vswf.LaguerreGaussian
%
% Additional casts inherited from base.
% See also :class:`ott.beam.properties.LaguerreGaussian` for properties.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Static)
    function args = likeProperties(other, args)
      % Construct an array of like-properties
      args = ott.beam.utils.VariablePower.likeProperties(other, args);
      args = ott.beam.properties.LaguerreGaussian.likeProperties(other, args);
    end

    function beam = like(other, varargin)
      % Create a beam like another beam
      %
      % Usage
      %   beam = LaguerreGaussian.like(other, ...)
      %
      % See constructor for arguments.

      args = ott.beam.abstract.LaguerreGaussian.likeProperties(...
          other, varargin);
      beam = ott.beam.abstract.LaguerreGaussian(args{:});
    end
  end

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

    function beam = ott.beam.vswf.Bsc(varargin)
      % Default Bsc cast, uses vswf.LaguerreGaussian
      beam = ott.beam.vswf.LaguerreGaussian(varargin{:});
    end

    function beam = ott.beam.paraxial.Paraxial(varargin)
      % Default Bsc cast, uses paraxial.LaguerreGaussian
      beam = ott.beam.paraxial.LaguerreGaussian(varargin{:});
    end

    function beam = ott.beam.vswf.LaguerreGaussian(beam, varargin)
      % Cast to vswf.LaguerreGaussian

      assert(isa(beam, 'ott.beam.abstract.LaguerreGaussian'), ...
          'first argument must be a abtract.LaguerreGaussian');

      beam = castHelper(@ott.beam.vswf.LaguerreGaussian.like, ...
          beam, 'lmode', beam.lmode, 'pmode', beam.pmode, varargin{:});
    end

    function beam = ott.beam.paraxial.LaguerreGaussian(beam, varargin)
      % Cast to paraxial.LaguerreGaussian

      assert(isa(beam, 'ott.beam.abstract.LaguerreGaussian'), ...
          'first argument must be a abtract.LaguerreGaussian');

      beam = castHelper(@ott.beam.paraxial.LaguerreGaussian.like, ...
          beam, 'lmode', beam.lmode, 'pmode', beam.pmode, varargin{:});
    end
  end
end

