classdef BeamForce < ott.scat.utils.Particle
% Provides force calculation methods using scatter and beam calculation.
% Inherits from :class:`ott.scat.utils.Particle`.
%
% This class is useful for particles which calculate how a beam is
% scattered and where the force should be calculated from the change
% in momentum of the two beam objects.
%
% Zero-order (non-scattering) methods should not inherit from this
% class since the scattered beam is equal to the incident beam.
%
% Hidden methods
%   - forceInternal       -- Method called by `force`
%   - torqueInternal      -- Method called by `torque`
%   - forcetorqueInternal -- Method called by `forcetorque`
%
% Abstract methods
%   - rotate      -- Apply rotation to particle
%   - scatter     -- Calculate scattered beam
%
% For full methods and properties, see :class:`ott.scat.utils.Particle`.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  methods (Hidden)
    function force = forceInternal(particle, ibeam, varargin)
      % Calculate the force using the beam method

      scattered = particle.scatter(ibeam);
      force = scattered.force();
    end

    function torque = torqueInternal(particle, ibeam, varargin)
      % Calculate the torque using the beam method

      scattered = particle.scatter(ibeam);
      force = scattered.torque();
    end

    function varargout = forcetorqueInternal(particle, ibeam, varargin)
      % Calculate the force-torque using the beam method

      scattered = particle.scatter(ibeam);
      [varargout{1:nargout}] = scattered.forcetorque();
    end
  end
end

