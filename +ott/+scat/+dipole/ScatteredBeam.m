classdef ScatteredBeam < ott.optics.beam.Beam
% Describes a beam that has been scattered by a dipole particle.
% Inherits from :class:`ott.optics.beam.Beam`
%
% Methods
%   - efield      -- Calculate the electric field
%   - hfield      -- Calculate the magnetic field
%   - ehfield     -- Calculate electric and magnetic fields
%   - visualise   -- Generate a visualisation of the fields
%
% Properties
%   - position    -- The dipole position in the original beam
%   - dipole      -- The dipole object that cause the scattering
%   - beam        -- The original beam
%   - type        -- This beam type (total or scattered)
%
% Dependent properties
%   - force_scat  -- The scattered force on the dipole (3x1 numeric)
%   - force_grad  -- The gradient force on the dipole (3x1 numeric)
%   - force       -- The total force on the dipole (3x1 numeric)

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

% TODO: Compare all these with T-matrix limit
% TODO: This has/should moved to dipole beam type???

  properties
    position    % The dipole position in the original beam
    dipole      % The dipole object that cause the scattering
    beam        % The original beam
  end

  properties (Dependent)
    force       % The scattered force on the dipole (3x1 numeric)
    force_scat  % The gradient force on the dipole (3x1 numeric)
    force_grad  % The total force on the dipole (3x1 numeric)
  end

  methods (Hidden)
    function E = efieldInternal(sbeam, xyz)
      % Calculate the total scattered field
      %
      % Uses equation 2.06 from
      %
      %   B. T. Draine, The Discrete-Dipole Approximation and Its
      %   Application to Interstellar Graphite Grains.
      %   The Astrophysical Journal, 333:848-872,1988 October 15

      % TODO: Do we use the polarizability or the polarizability
      % with the correction factor

      % TODO: Implement me
      error('Not yet implemented');
    end

    function H = hfieldInternal(sbeam, xyz)
      % TODO: Implement
      error('Not yet implemented');
    end
  end

  methods
    function sbeam = ScatteredBeam(beam, dipole, position)
      % Construct a new dipole scattered beam
      %
      % Usage
      %   sbeam = ScatteredBeam(beam, dipole, position)
      %
      % Parameters
      %   - beam -- A :class:`ott.optics.beam.Beam`
      %   - dipole -- A :class:`Dipole`
      %   - position (3x1 numeric) -- Location of the scattering event

      sbeam = sbeam@ott.optics.beam.Beam();

      sbeam.beam = beam;
      sbeam.dipole = dipole;
      sbeam.position = position;
    end
  end

  methods
    function sbeam = set.beam(sbeam, val)
      assert(isa(val, 'ott.optics.beam.Beam'), ...
        'beam must be a ott.optics.beam.Beam objects');
      sbeam.beam = val;
    end
    function sbeam = set.dipole(sbeam, val)
      assert(isa(val, 'ott.optics.dipole.Dipole'), ...
        'dipole must be a ott.optics.dipole.Dipole');
      sbeam.dipole = val;
    end
    function sbeam = set.position(sbeam, val)
      assert(isnumeric(val) && all(size(val) == [3, 1]), ...
        'position must be 3x1 numeric');
      sbeam.position = val;
    end

    function val = get.force(sbeam)
      % Call the dipole force method
      force = sbeam.dipole.force(sbeam.beam, sbeam.position);
    end
    function val = get.force_scat(sbeam)
      % Call the dipole force method
      val = sbeam.dipole.force_scat(sbeam.beam, sbeam.position);
    end
    function val = get.force_grad(sbeam)
      % Call the dipole force method
      val = sbeam.dipole.force_grad(sbeam.beam, sbeam.position);
    end
  end
end

