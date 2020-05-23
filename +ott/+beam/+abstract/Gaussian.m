classdef Gaussian < ott.beam.properties.Gaussian ...
    & ott.beam.abstract.Beam ...
    & ott.beam.utils.VariablePower
% Describes the Beam casts for a Gaussian beam.
% Inherits from :class:`Beam`, :class:`ott.beam.properties.Gaussian`,
% and :class:`utils.VariablePower`.
%
% Supported casts
%   - Beam            -- Default Beam cast, uses GaussianDavis5
%   - vswf.Bsc        -- Default Bsc cast, uses vswf.Gaussian
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

  methods
    function beam = Gaussian(waist, varargin)
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

    function beam = ott.beam.Beam(oldbeam, varargin)
      % Cast the beam to a ott.beam.Beam object
      %
      % The default beam is a ott.beam.GaussianDavis5, since it
      % is a good compromise between speed and accuracy.

      beam = ott.beam.GaussianDavis5(oldbeam, varargin{:});
    end

    function beam = ott.beam.GaussianDavis5(oldbeam, varargin)
      % Cast beam to a GaussianDavis5

      beam = ott.beam.GaussianDavis5(oldbeam.waist, ...
        'omega', oldbeam.omega, 'medium', oldbeam.medium, ...
        'power', oldbeam.power, 'position', oldbeam.position, ...
        'rotation', oldbeam.rotation, varargin{:});
    end

    function beam = ott.beam.vswf.Bsc(oldbeam, varargin)
      % Convert to a Bsc object
      %
      % Usage
      %   bsc = ott.beam.vswf.Bsc(abstract_beam, ...)
      %
      % Optional named arguments
      %   - suggestedNmax (numeric) -- Suggested Nmax for the Bsc.
      %     This parameter is ignored.

      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('suggestedNmax', []);
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      NA = oldbeam.wavelength0 ./ (pi * oldbeam.waist);

      beam = ott.beam.vswf.PmGauss('lg', [0, 0], ...
          'omega', oldbeam.omega, 'medium', oldbeam.medium, ...
          'NA', NA, 'power', oldbeam.power, 'position', oldbeam.position, ...
          'rotation', oldbeam.rotation, unmatched{:});
    end
  end
end

