classdef Gaussian < ott.beam.properties.Gaussian ...
    & ott.beam.abstract.Beam ...
    & ott.beam.utils.VariablePower
% Describes the Beam casts for a Gaussian beam.
% Inherits from :class:`Beam`, :class:`ott.beam.properties.Gaussian`,
% and :class:`utils.VariablePower`.
%
% Supported casts
%   - Beam              -- Default Beam cast, uses GaussianDavis5
%   - vswf.Bsc          -- Default Bsc cast, uses vswf.Gaussian
%   - paraxial.Paraxial -- Default Paraxial, uses paraxial.Gaussian
%   - GaussianDavis5
%   - vswf.Gaussian
%   - vswf.HermiteGaussian
%   - vswf.LaguerreGaussian
%   - vswf.InceGaussian
%   - paraxial.Gaussian
%   - paraxial.LaguerreGaussian
%   - paraxial.HermiteGaussian
%
% See also :class:`ott.beam.properties.Gaussian` for properties.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Static)
    function args = likeProperties(other, args)
      % Construct an array of like-properties
      args = ott.beam.utils.VariablePower.likeProperties(other, args);
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

    function beam = ott.beam.Beam(varargin)
      % Cast the beam to a ott.beam.Beam object
      %
      % The default beam is a ott.beam.GaussianDavis5, since it
      % is a good compromise between speed and accuracy.

      beam = ott.beam.GaussianDavis5(varargin{:});
    end

    function beam = ott.beam.GaussianDavis5(beam, varargin)
      % Cast beam to a GaussianDavis5
      beam = castHelper(@ott.beam.GaussianDavis5.like, ...
          beam, varargin{:});
    end

    function beam = ott.beam.paraxial.Paraxial(varargin)
      % Construct a paraxial.Gaussian beam
      beam = ott.beam.paraxial.Gaussian(varargin{:});
    end

    function beam = ott.beam.paraxial.Gaussian(beam, varargin)
      % Construct a paraxial.Gaussian beam
      beam = castHelper(@ott.beam.paraxial.Gaussian.like, ...
          beam, varargin{:});
    end

    function beam = ott.beam.paraxial.HermiteGaussian(beam, varargin)
      % Construct a paraxial.HermiteGaussian beam
      beam = castHelper(@ott.beam.paraxial.HermiteGaussian.like, ...
          beam, 'mmode', 0, 'nmode', 0, varargin{:});
    end

    function beam = ott.beam.paraxial.LaguerreGaussian(beam, varargin)
      % Construct a paraxial.LaguerreGaussian beam
      beam = castHelper(@ott.beam.paraxial.LaguerreGaussian.like, ...
          beam, 'lmode', 0, 'pmode', 0, varargin{:});
    end

    function beam = ott.beam.vswf.Bsc(varargin)
      % Cast to VSWF Gaussian
      beam = ott.beam.vswf.Gaussian(varargin{:});
    end

    function beam = ott.beam.vswf.Gaussian(beam, varargin)
      % Convert to VSWF Gaussian
      beam = castHelper(@ott.beam.vswf.Gaussian.like, beam, varargin{:});
    end

    function beam = ott.beam.vswf.LaguerreGaussian(beam, varargin)
      % Convert to VSWF Laguerre-Gaussian
      beam = castHelper(@ott.beam.vswf.LaguerreGaussian.like, ...
          beam, 'lmode', 0, 'pmode', 0, varargin{:});
    end

    function beam = ott.beam.vswf.HermiteGaussian(beam, varargin)
      % Convert to VSWF Hermite-Gaussian
      beam = castHelper(@ott.beam.vswf.HermiteGaussian.like, ...
          beam, 'mmode', 0, 'nmode', 0, varargin{:});
    end

    function beam = ott.beam.vswf.InceGaussian(beam, varargin)
      % Convert to VSWF Ince-Gaussian
      beam = castHelper(@ott.beam.vswf.InceGaussian.like, ...
          beam, 'porder', 0, 'lmode', 0, 'parity', 'even', ...
          'ellipticity', 1.0, varargin{:});
    end
  end

  methods (Access=private)
    function beam = castHelper(cast, beam, varargin)
      % Helper for casts

      assert(isa(beam, 'ott.beam.abstract.Gaussian'), ...
          'First argument must be a abstract.Gaussian');

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

