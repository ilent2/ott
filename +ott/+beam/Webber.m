classdef Webber < ott.beam.BscWBessel & ott.beam.properties.Webber
% Construct a VSWF representation of a Bessel beam.
% Inherits from :class:`+ott.+beam.BscInfinite` and
% :class:`+ott.+beam.+properties.Webber`.
%
% Properties
%   - theta       -- Annular angle [radians]
%   - alpha       -- Parameter describe Webber beam
%   - parity      -- Parity of beam (either 'even' or 'odd')
%   - data        -- Internal BSC instance describing beam

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function beam = Webber(varargin)
      % Construct a new VSWF Webber beam representation.
      %
      % Usage
      %   beam = Webber(theta, alpha, parity, ...)
      %
      % Optional named parameters
      %   - theta (numeric) -- Annular beam angle from -z direction (radians).
      %     Default: ``pi/4``.
      %
      %   - alpha (numeric) -- Webber beam parameter.  Default: ``1``.
      %
      %   - parity (enum) -- Parity of beam.  Either 'even' or 'odd'
      %     Default: ``'even'``.
      %
      %   - Nmax (numeric) -- Initial beam Nmax.  Default: ``0``.
      %     This parameter automatically grows when the beam is used.
      %     See also :meth:`recalculate` and :meth:`getData`.

      p = inputParser;
      p.addOptional('theta', pi/4, @isnumeric);
      p.addOptional('alpha', 1, @isnumeric);
      p.addOptional('parity', 'even', ...
          @(x) sum(strcmpi(x, {'even', 'odd'})) == 1);
      p.addParameter('Nmax', 0, @isnumeric);
      p.parse(varargin{:});

      beam.theta = p.Results.theta;
      beam.alpha = p.Results.alpha;
      beam.parity = p.Results.parity;
      beam = beam.recalculate(p.Results.Nmax);
    end

    function beam = recalculate(beam, Nmax)
      % Re-calculate BSC data for specified Nmax.

      % TODO: Support for different polarisations?

      [beam.data, beam.besselWeights] = ott.bsc.Annular.FromWebber(...
          beam.theta, beam.alpha, beam.parity, 'Nmax', Nmax);
    end
  end
end

