classdef BscWBessel < ott.beam.BscInfinite
% Base class for weighted BscBessel beams.
% Inherits from :class:`BscInfinite`.
%
% Properties
%   - besselWeights    -- Apply weights at end of cast (default: [])
%
% Static methods
%   - ott.bsc.Bsc       -- Cast to a bsc instance (smartly)

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    besselWeights = []
  end

  methods
    function bsc = ott.bsc.Bsc(beam, Nmax)
      % Get the BSC data for a specific Nmax
      %
      % If the current data is empty or the current Nmax is less than
      % the desired Nmax, calls :meth:`recalculate` and then returns the
      % internal data.
      %
      % Rotations are applied without changing Nmax.
      % Translations parallel to the beam axis are applied without
      % changing Nmax.
      %
      % Usage
      %   bsc = ott.bsc.Bsc(beam, Nmax)
      %
      % Parameters
      %   - beam (ott.beam.BscBeam) -- The beam to get data from.
      %
      %   - Nmax (numeric) -- Desired Nmax.  The default is either 0
      %     or the current beam Nmax minus the Nmax required for any
      %     radial translations, whichever is greater.

      % Get parallel and perpendicular component of translation
      axial = dot(beam.position, beam.rotation(:, 3));
      radial = beam.position - axial * beam.rotation(:, 3);

      % Calculate Nmax needed for radial part
      rNmax = ott.utils.ka2nmax(vecnorm(radial) * beam.wavenumber);

      if nargin < 2
        Nmax = max(0, max([0, beam.Nmax]) - rNmax);
      end

      if (Nmax+rNmax) > beam.Nmax
        beam = beam.recalculate(Nmax + rNmax);
      end

      % Get data and apply rotation
      bsc = beam.data.rotate(beam.rotation);

      % Apply Nmax preserving translation
      bsc = bsc.translateZ(axial);

      % Apply Nmax non-preserving translation
      bsc = bsc.translateXyz(radial, 'Nmax', Nmax);

      % Apply weights (if available)
      if ~isempty(beam.besselWeights)
        bsc = bsc * beam.besselWeights(:);
      end

      % Apply scale
      bsc = bsc * beam.scale;
    end
  end
end

