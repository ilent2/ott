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
% Methods
%   - efarfieldInternal -- Approximates far-field using efieldInternal
%   - hfieldInternal    -- Calculates hfield from efieldInternal
%
% Abstract methods:
%   - efieldInternal    -- Called by efield
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

    function E = efarfieldInternal(beam, rtp, varargin)
      % Approximates the e far-field from the e near-field
      %
      % Calculates the fields at a distance of 10 wavelengths from
      % the coordinate origin.  This function should be overloaded
      % if a better approximation exists for your beam.
      %
      % Usage
      %   E = beam.efarfieldInternal(rtp)
      %
      % Optional named arguments
      %   - mapping (enum) -- Either 'theta' or 'rho' for field direction
      %     to keep.  Radial component for 'theat' is scaled by cos(theta).
      %     Default: ``'theta'``.
      %
      %   - far_distance (numeric) -- A far distance to calculate fields.
      %     Default: ``10*beam.wavelength``.

      p = inputParser;
      p.addParameter('mapping', 'theta');
      p.addParameter('far_distance', 10*max([beam.wavelength]));
      p.parse(varargin{:});

      % Choose a distance
      r = p.Results.far_distance;
      if size(rtp, 1) == 2
        rtp = [p.Results.far_distance; rtp];
      else
        rtp(1, :) = p.Results.far_distance;
      end

      % Calculate the fields at this distance
      xyz = ott.utils.rtp2xyz(rtp);
      E = beam.efieldInternal(xyz);

      onesrow = ones(size(rtp(1, :)));
      zerosrow = zeros(size(rtp(1, :)));

      phi_hat = ott.utils.rtpv2xyzv([zerosrow; zerosrow; onesrow], rtp);
      theta_hat = ott.utils.rtpv2xyzv([zerosrow; onesrow; zerosrow], rtp);

      switch p.Results.mapping
        case 'theta'
          % Nothing to do
        case 'rho'
          theta_hat(3, :) = 0;
          theta_hat = theta_hat ./ vecnorm(theta_hat);
        otherwise
          error('Mapping must be ''theta'' or ''rhow''');
      end

      % Calculate fields in specified direction
      Ephi = dot(phi_hat, E.vxyz);
      Etheta = dot(theta_hat, E.vxyz);

      % Package output
      Ertp = [ones(size(Ephi)); Etheta; Ephi];
      E = ott.utils.FieldVector(rtp, Ertp, 'spherical');
    end

    function H = hfarfieldInternal(beam, rtp)
      % Calculate h field from e far-field
      %
      % Usage
      %   H = beam.hfarfield(rtp)
      %   Calculates the fields at the specified locations (3xN matrix).
      %   Returns a :class:`ott.utils.FieldVector`.

      E = beam.efarfieldInternal(rtp);

      % Construct H-field (might have 2 or three rows)
      H = zeros(size(E));
      swprows = [size(H, 1)-1, size(H, 1)];
      H(flip(swprows), :) = E.vrtp(swprows, :) ./ beam.impedance;

      % Package output
      H = ott.utils.FieldVector(rtp, H, 'spherical');
    end
  end
end

