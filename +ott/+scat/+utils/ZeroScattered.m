classdef (Abstract) ZeroScattered < ott.scat.utils.Particle
% Base class for particles which produce Zero-scattered beams.
% Inherits from :class:`ott.scat.utils.Particle`.
% This class is incompatible with :class:`ott.scat.utils.BeamForce`.
%
% Abstract properties
%   - position    -- Position of the particle [3xN]
%   - rotation    -- Orientation of the particle [3x3]
%
% Methods
%   - scatter     -- Constructs a :class:`ott.beam.ZeroScattered` beam.
%   - force       -- Calculate force on particle in beam
%   - torque      -- Calculate torque on particle in beam
%   - forcetorque -- Calculate force and torque on particle in beam
%
% Abstract methods
%   - rotate      -- Apply rotation to particle
%   - forceInternal       -- Method called by `force`
%   - torqueInternal      -- Method called by `torque`
%   - forcetorqueInternal -- Method called by `forcetorque`

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  methods (Hidden)
    function sbeam = scatterInternal(particle, ibeam)
      % Construct a new :class:`ott.beam.ZeroScattered` beam.
      %
      % Usage
      %   sbeam = particle.scatter(beam)

      sbeam = ott.beam.ZeroScattered('total', ibeam, particle);
    end
  end
end
