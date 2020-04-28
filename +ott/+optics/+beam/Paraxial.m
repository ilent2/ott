classdef (Abstract) Paraxial < ott.optics.beam.Beam
% A abstract beam class for Paraxial beams.
% Inherits from :class:`Beam`.
%
% Assumes the beam is paraxial, i.e. :math:`E_z = H_z = 0` and::
%
%    H = \hat{z} \times \frac{E}{Z}
%
% where :math:`Z` is the medium impedance.
%
% Abstract methods:
%   - efieldInternal    -- Called by efield
%   - efarfieldInternal -- Called by efarfield
%   - getBeamPower      -- get method called by dependent property power

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Hidden)
    function H = hfieldInternal(beam, xyz)
      % Calculate h field from e field
      %
      % Usage
      %   H = beam.hfield(xyz)
      %   Calculates the fields at the specified locations (3xN matrix).
      %   Returns a :class:`ott.utils.FieldVector`.

      E = beam.efieldInternal(xyz);

      % Construct H-field (
      H = zeros(size(E));
      H([2, 1], :) = E.vxyz([1, 2], :) ./ beam.impedance;

      % Package output
      H = ott.utils.FieldVector(xyz, H, 'cartesian');
    end

    function H = hfarfieldInternal(beam, rtp)
      error('Not yet implemented');
    end
  end
end

