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
      % Optional named arguments
      %   - power (numeric) -- Beam power.  Default: ``1.0``.
      %
      % For optional parameters, see :class:`Properties`.

      args = ott.utils.addDefaultParameter('power', 1.0, varargin);
      beam = beam@ott.beam.properties.HermiteGaussian(args{:});
    end

    function beam = ott.beam.Beam(varargin)
      % Default Beam cast, uses paraxial.HermiteGaussian
      beam = ott.beam.paraxial.HermiteGaussian(varargin{:});
    end

    function beam = ott.beam.vswf.Bsc(varargin)
      % Default Bsc cast, uses vswf.HermiteGaussian
      beam = ott.beam.vswf.HermiteGaussian(varargin{:});
    end

    function beam = ott.beam.vswf.Gaussian(beam, varargin)
      % Cast to vswf.Gaussian via abstract.Gaussian
      beam = ott.beam.abstract.Gaussian(beam, varargin{:});
      beam = ott.beam.vswf.Gaussian(beam, varargin{:});
    end

    function beam = ott.beam.vswf.LaguerreGaussian(beam, varargin)
      % Cast to vswf.LaguerreGaussian via abstract.Gaussian
      beam = ott.beam.abstract.Gaussian(beam, varargin{:});
      beam = ott.beam.vswf.LaguerreGaussian(beam, varargin{:});
    end

    function beam = ott.beam.vswf.InceGaussian(beam, varargin)
      % Cast to vswf.InceGaussian via abstract.Gaussian
      beam = ott.beam.abstract.Gaussian(beam, varargin{:});
      beam = ott.beam.vswf.InceGaussian(beam, varargin{:});
    end

    function beam = ott.beam.paraxial.Gaussian(beam, varargin)
      % Cast to paraxial.Gaussian via abstract.Gaussian
      beam = ott.beam.abstract.Gaussian(beam, varargin{:});
      beam = ott.beam.paraxial.Gaussian(beam, varargin{:});
    end

    function beam = ott.beam.paraxial.LaguerreGaussian(beam, varargin)
      % Cast to paraxial.LaguerreGaussian via abstract.Gaussian
      beam = ott.beam.abstract.Gaussian(beam, varargin{:});
      beam = ott.beam.paraxial.LaguerreGaussian(beam, varargin{:});
    end

    function beam = ott.beam.abstract.Gaussian(beam, varargin)
      % Cast the beam to an abstract.Gaussian

      assert(isa(beam, 'ott.beam.abstract.HermiteGaussian'), ...
          'first argument must be a abtract.HermiteGaussian');
      assert(beam.lmode == 0 && beam.pmode == 0, ...
          'Beam must be HG00 for cast to Gaussian');

      beam = castHelper(@ott.beam.abstract.Gaussian.like, ...
          beam, varargin{:});
    end

    function beam = ott.beam.vswf.HermiteGaussian(beam, varargin)
      % Cast to vswf.HermiteGaussian

      assert(isa(beam, 'ott.beam.abstract.HermiteGaussian'), ...
          'first argument must be a abtract.HermiteGaussian');

      beam = castHelper(@ott.beam.vswf.HermiteGaussian.like, ...
          beam, 'mmode', beam.mmode, 'nmode', beam.nmode, varargin{:});
    end

    function beam = ott.beam.paraxial.HermiteGaussian(beam, varargin)
      % Cast to paraxial.HermiteGaussian

      assert(isa(beam, 'ott.beam.abstract.HermiteGaussian'), ...
          'first argument must be a abtract.HermiteGaussian');

      beam = castHelper(@ott.beam.paraxial.HermiteGaussian.like, ...
          beam, 'mmode', beam.mmode, 'nmode', beam.nmode, varargin{:});
    end
  end
end
