classdef InceGaussian < ott.beam.properties.InceGaussian ...
    & ott.beam.abstract.Beam ...
    & ott.beam.utils.VariablePower
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

    function beam = ott.beam.Beam(varargin)
      % Default Beam cast, uses vswf.InceGaussian
      beam = ott.beam.vswf.InceGaussian(varargin{:});
    end

    function beam = ott.beam.vswf.Bsc(varargin)
      % Default Beam cast, uses vswf.InceGaussian
      beam = ott.beam.vswf.InceGaussian(varargin{:});
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

    function beam = ott.beam.vswf.LaguerreGaussian(beam, varargin)
      % Cast to vswf.LaguerreGaussian via abstract.Gaussian
      beam = ott.beam.abstract.Gaussian(beam, varargin{:});
      beam = ott.beam.vswf.LaguerreGaussian(beam, varargin{:});
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

    function beam = ott.beam.paraxial.LaguerreGaussian(beam, varargin)
      % Cast to paraxial.LaguerreGaussian via abstract.Gaussian
      beam = ott.beam.abstract.Gaussian(beam, varargin{:});
      beam = ott.beam.paraxial.LaguerreGaussian(beam, varargin{:});
    end

    function beam = ott.beam.abstract.Gaussian(beam, varargin)
      % Cast the beam to an abstract.Gaussian

      assert(isa(beam, 'ott.beam.abstract.InceGaussian'), ...
          'first argument must be a abtract.InceGaussian');

      assert(beam.lmode == 0 && beam.porder == 0 && ...
          strcmpi(beam.parity, 'even') && beam.ellipticity == 0.0, ...
          'Beam must be Gaussian for cast to Gaussian');

      beam = castHelper(@ott.beam.abstract.Gaussian.like, ...
          beam, varargin{:});
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
