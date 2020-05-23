classdef Dda < ott.scat.utils.BeamForce
% Method to simulate scattering using the discrete dipole approximation.
% Inherits from :class:`ott.scat.utils.BeamForce`
%
% This class stores the dipole locations and interactions of many discrete
% dipoles and calculates the scattered fields by solving::
%
%    A p = Ei
%
% and::
%
%   Es = F p
%
% where :math:`A` is the dipole interaction matrix, :math:`p` is the
% dipole polarisability induced by the incident field :math:`Ei`,
% :math:`Es` is the scattered field and :math:`F` is the transfer
% matrix relating the polarisabilities to the field locations.
%
% Properties
%   - dipole_xyz    -- Dipole locations
%   - interaction   -- Dipole interaction matrix

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    dipole_xyz      % Dipole locations
    interaction     % Dipole interaction matrix
  end

  methods
    function sbeam = scatter(dda, ibeam, varargin)
      % Calculates the scattered beam
      %
      % Solves the equation::
      %
      %   A p = Ei
      %
      % for the polarizabilities and returns a new dipole beam object.

      % Get the incident field at each dipole
      Ei = ibeam.ehfield(dda.dipole_xyz);

      % Solve the system for p
      p = A \ Ei;

      % Construct a new beam
      sbeam = ott.beam.Scattered('scattered', ...
          ott.beam.Dipole(dipole_xyz, p), ...
          'incident_beam', ibeam, ...
          'particle', dda);
    end
  end
end

