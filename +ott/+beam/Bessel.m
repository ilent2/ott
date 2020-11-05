classdef Bessel < ott.beam.BscWBessel & ott.beam.properties.Bessel
% Construct a VSWF representation of a Bessel beam.
% Inherits from :class:`ott.beam.BscInfinite` and
% :class:`ott.beam.properties.Bessel`.
%
% See base classes for list of properties/methods.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function beam = Bessel(varargin)
      % Construct a new VSWF Bessel beam representation.
      %
      % Usage
      %   beam = Bessel(theta, polfield, lmode, ...)
      %
      % Optional named arguments
      %   - theta (numeric) -- Annular beam angle from -z direction (radians).
      %     Default: ``pi/4``.
      %
      %   - polfield (2 numeric) -- Field in the [theta, phi] or [x, y]
      %     directions.  Default: ``[1, 1i]``.
      %
      %   - polbasis (enum) -- Polarisation basis.  Can be either
      %     'polar' or 'cartesian'.  Default: ``'polar'``.
      %
      %   - lmode (numeric) -- Azimuthal mode number.  Default: ``0``.
      %
      %   - Nmax (numeric) -- Initial beam Nmax.  Default: ``0``.
      %     This parameter automatically grows when the beam is used.
      %     See also :meth:`recalculate` and :meth:`getData`.

      p = inputParser;
      p.addOptional('theta', pi/4);
      p.addOptional('polfield', [1, 1i]);
      p.addOptional('lmode', 0);
      p.addParameter('polbasis', 'polar');
      p.addParameter('Nmax', 0);
      p.parse(varargin{:});

      beam.theta = p.Results.theta;
      beam.polfield = p.Results.polfield;
      beam.polbasis = p.Results.polbasis;
      beam.lmode = p.Results.lmode;
      beam = beam.recalculate(p.Results.Nmax);
    end

    function beam = recalculate(beam, Nmax)
      % Re-calculate BSC data for specified Nmax.

      % TODO: Support for different polarisations

      beam.data = ott.bsc.Annular.FromBessel(Nmax, beam.theta, ...
          beam.polfield(:), beam.lmode);
    end
  end
end
