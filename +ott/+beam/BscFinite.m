classdef BscFinite < ott.beam.BscBeam
% A beam represented by a finite VSWF expansion.
% Inherits from :class:`BscBeam`.
%
% This class describes beams which can be represented using a finite VSWF
% expansion. The class stores a :class:`Bsc` instance internally.
% BSC coefficients at any other location can be found by applying a
% translation to the beam data.
%
% As with the :class:`BscBeam` class, this class assumes a regular beam.
%
% Supported casts
%   - ott.bsc.Bsc   -- Get the BSC data after applying transformations
%
% Additional properties/methods inherited from :class:`BscBeam`.

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  methods
    function beam = BscFinite(varargin)
      % Construct a new Finite Bsc instance
      %
      % Usage
      %   beam = BscFinite(data, ...)
      %
      % Optional named parameters
      %   - data (ott.bsc.Bsc) -- Initial data for beam.
      %     Default: ``ott.bsc.Bsc.empty()``.
      %
      %   - index_medium (numeric) -- Refractive index of the medium.
      %     Default: ``1.0``.
      %
      %   - omega (numeric) -- Optical angular frequency [Hz].
      %     Default: ``3e8/1064e-9*2*pi`` (i.e., default vacuum
      %     wavelength is 1064 nm).

      beam = beam@ott.beam.BscBeam(varargin{:});
    end

    function bsc = ott.bsc.Bsc(beam, Nmax)
      % Get the beam data for the specified Nmax.
      %
      % If the current data is empty, calls :meth:`recalculate`.
      %
      % If the position and rotation have been set, rotations and translations
      % are applied to the beam (rotations first, translations second).
      % The translated beam will have the specified Nmax.
      %
      % Usage
      %   bsc = ott.bsc.Bsc(beam, Nmax)
      %
      % Parameters
      %   - beam (ott.beam.BscFinite) -- The beam to get data from.
      %
      %   - Nmax (numeric) -- Desired Nmax.  This parameter is ignored
      %     except when a translation is applied, in which case the
      %     translated beam has this Nmax.  Default: current Nmax or 0.

      if nargin < 2
        Nmax = max([0, beam.data.Nmax]);
      end

      if isempty(beam.data)
        beam = beam.recalcualte(Nmax);
      end

      bsc = beam.data;

      % Rotate the beam
      bsc = bsc.rotate(beam.rotation);

      % Translate the beam
      % Assumes beam is a regular beam
      bsc = bsc.translateXyz(beam.position ./ beam.wavelength, 'Nmax', Nmax);
    end
  end
end
