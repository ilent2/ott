classdef Bessel < ott.beam.properties.Bessel ...
    & ott.beam.abstract.Beam ...
    & ott.beam.utils.InfinitePower
% Abstract description of a Bessel beam.
% Inherits from :class:`ott.beam.properties.Bessel` and :class:`Beam`.
%
% Supported casts
%   Beam              -- Default Beam cast, uses vswf.Bessel
%   Bsc               -- Default Bsc cast, uses vswf.Bessel
%   vswf.Bessel
%   vswf.BesselBasis

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function beam = Bessel(varargin)
      % Construct a new abstract Bessel beam
      %
      % Usage
      %   beam = Bessel(theta, ...)
      %   Parameters can also be passed as named arguments.
      %
      % Parameters
      %   - theta (numeric) -- Angle of the Bessel beam from beam axis.

      beam = beam@ott.beam.properties.Bessel(varargin{:});
    end

    function beam = ott.beam.Beam(varargin)
      % Cast to vswf.Bessel beam
      beam = ott.beam.vswf.Bessel(varargin{:});
    end

    function beam = ott.beam.vswf.Bsc(varargin)
      % Cast to vswf.Bessel beam
      beam = ott.beam.vswf.Bessel(varargin{:});
    end

    function beam = ott.beam.vswf.Bessel(beam, varargin)
      % Cast to vswf.Bessel beam

      assert(isa(beam, 'ott.beam.abstract.Bessel'), ...
          'First argument must be a abstract.Bessel');

      % TODO: Other beam properties
      beam = ott.beam.vswf.Bessel(beam.theta, varargin{:});
    end

    function beam = ott.beam.vswf.BesselBasis(beam, varargin)
      % Cast to vswf.Bessel beam

      assert(isa(beam, 'ott.beam.abstract.Bessel'), ...
          'First argument must be a abstract.Bessel');

      % TODO: Other beam properties
      beam = ott.beam.vswf.BesselBasis(beam.theta, varargin{:});
    end
  end
end
