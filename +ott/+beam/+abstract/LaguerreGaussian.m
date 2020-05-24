classdef LaguerreGaussian < ott.beam.properties.LaguerreGaussian ...
    & ott.beam.abstract.Beam ...
    & ott.beam.utils.VariablePower
% Abstract representation of a Laguerre-Gaussian beam.
% Inherits from :class:`Gaussian` and
% :class:`ott.beam.properties.LaguerreGaussian`.
%
% Supported casts
%   - Beam            -- Default Beam cast, uses paraxial.LaguerreGaussian
%   - vswf.Bsc        -- Default Bsc cast, uses vswf.LaguerreGaussian
%   - vswf.LaguerreGaussian
%   - vswf.Gaussian             -- Only for LG00, via abstract.Gaussian
%   - vswf.HermiteGaussian      -- Only for LG00, via abstract.Gaussian
%   - vswf.InceGaussian         -- Only for LG00, via abstract.Gaussian
%   - paraxial.LaguerreGaussian
%   - paraxial.Gaussian         -- Only for LG00, via abstract.Gaussian
%   - paraxial.HermiteGaussian  -- Only for LG00, via abstract.Gaussian
%   - abstract.Gaussian
%
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

    function beam = ott.beam.Beam(varargin)
      % Default Beam cast, uses paraxial.LaguerreGaussian
      beam = ott.beam.paraxial.LaguerreGaussian(varargin{:});
    end

    function beam = ott.beam.vswf.Bsc(varargin)
      % Default Bsc cast, uses vswf.LaguerreGaussian
      beam = ott.beam.vswf.LaguerreGaussian(varargin{:});
    end

    function beam = ott.beam.vswf.Gaussian(beam, varargin)
      % Cast to vswf.Gaussian via abstract.Gaussian
      beam = ott.beam.abstract.Gaussian(beam, varargin{:});
      beam = ott.beam.vswf.Gaussian(beam, varargin{:});
    end

    function beam = ott.beam.vswf.HermiteGaussian(beam, varargin)
      % Cast to vswf.HermiteGaussian via abstract.Gaussian
      beam = ott.beam.abstract.Gaussian(beam, varargin{:});
      beam = ott.beam.vswf.HermiteGaussian(beam, varargin{:});
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

    function beam = ott.beam.paraxial.HermiteGaussian(beam, varargin)
      % Cast to paraxial.HermiteGaussian via abstract.Gaussian
      beam = ott.beam.abstract.Gaussian(beam, varargin{:});
      beam = ott.beam.paraxial.HermiteGaussian(beam, varargin{:});
    end

    function beam = ott.beam.abstract.Gaussian(beam, varargin)
      % Cast the beam to an abstract.Gaussian

      assert(isa(beam, 'ott.beam.abstract.LaguerreGaussian'), ...
          'first argument must be a abtract.LaguerreGaussian');

      assert(beam.lmode == 0 && beam.pmode == 0, ...
          'Beam must be LG00 for cast to Gaussian');

      beam = ott.beam.abstract.Gaussian.like(beam, varargin{:});
    end

    function beam = ott.beam.vswf.LaguerreGaussian(beam, varargin)
      % Cast to vswf.LaguerreGaussian
      beam = castHelper(@ott.beam.vswf.LaguerreGaussian.like, ...
          beam, 'lmode', beam.lmode, 'pmode', beam.pmode, varargin{:});
    end

    function beam = ott.beam.paraxial.LaguerreGaussian(beam, varargin)
      % Cast to paraxial.LaguerreGaussian
      beam = castHelper(@ott.beam.paraxial.LaguerreGaussian.like, ...
          beam, 'lmode', beam.lmode, 'pmode', beam.pmode, varargin{:});
    end
  end

  methods (Access=private)
    function beam = castHelper(cast, beam, varargin)
      % Helper for casts

      assert(isa(beam, 'ott.beam.abstract.LaguerreGaussian'), ...
          'First argument must be a abstract.LaguerreGaussian');

      if numel(beam) > 1
        oldbeam = beam;
        beam = ott.beam.Coherent(size(beam));
        for ii = 1:numel(oldbeam)
          beam(ii) = cast(beam, varargin{:});
        end
      else
        beam = cast(beam, varargin{:});
      end
    end
  end
end

