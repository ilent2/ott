classdef BeamForce < ott.scat.Particle
% Provides force calculation methods using scatter and beam calculation.
% Inherits from :class:`ott.scat.Particle`.
%
% This class is useful for particles which calculate how a beam is
% scattered and where the force should be calculated from the change
% in momentum of the two beam objects.
%
% Zero-order (non-scattering) methods should not inherit from this
% class since the scattered beam is equal to the incident beam.
%
% Redefined methods
%   - force        -- Calculate scattered then calculate force
%   - torque       -- Calculate scattered then calculate torque
%   - forcetorque  -- Calculate scattered then calculate force/torque
%
% Hidden methods
%   - forceInternal       -- Unused method (same effect as `force`)
%   - torqueInternal      -- Unused method (same effect as `torque`)
%   - forcetorqueInternal -- Unused method (same effect as `forcetorque`)
%
% Abstract methods
%   - rotate      -- Apply rotation to particle
%   - scatter     -- Calculate scattered beam
%
% For full methods and properties, see :class:`ott.scat.Particle`.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  methods
    function varargout = force(particle, beam, varargin)
      % Calculate the force from a beam on a particle.
      %
      % Uses the scatter method and the scattered beams force method.
      %
      % Usage
      %   force = particle.force(beam, ...)
      %
      % For details on usage/arguments see :meth:`forcetorque`.
      
      % In most cases, it is more optimal to scatter all beams first
      scattered = particle.scatter(beam, varargin{:});
      [varargout{1:nargout}] = scattered.force();
    end

    function varargout = torque(particle, beam, varargin)
      % Calculate the torque from a beam on a particle.
      %
      % Uses the scatter method and the scattered beams torque method.
      %
      % Usage
      %   torque = particle.torque(beam, ...)
      %
      % For details on usage/arguments see :meth:`forcetorque`.
      
      % In most cases, it is more optimal to scatter all beams first
      scattered = particle.scatter(beam, varargin{:});
      [varargout{1:nargout}] = scattered.torque();
    end

    function varargout = forcetorque(particle, beam, varargin)
      % Calculate force and torque from a beam on a particle.
      %
      % Uses the scatter method and the scattered beams forcetorque method.
      %
      % Usage
      %   [force, torque] = particle.forcetorque(beam, ...)
      %
      % Optional named arguments
      %   - position (3xN numeric) -- Positions of the particle to
      %     calculate properties for.
      %     Default: ``[]``.
      %
      %   - rotation (3x3N numeric) -- Orientations of the particle
      %     to calculate properties for.
      %     Default: ``[]``.
      
      % In most cases, it is more optimal to scatter all beams first
      scattered = particle.scatter(beam, varargin{:});
      [varargout{1:nargout}] = scattered.forcetorque();
    end
  end

  methods (Hidden)
    function force = forceInternal(particle, ibeam, varargin)
      % Calculate the force using the beam method
      %
      % This method is not used by :meth:`force`.  For most beams it
      % is better to scatter all beams first.  This method allows
      % scattering individual beams.

      scattered = particle.scatter(ibeam);
      force = scattered.force();
    end

    function torque = torqueInternal(particle, ibeam, varargin)
      % Calculate the torque using the beam method
      %
      % This method is not used by :meth:`torque`.  For most beams it
      % is better to scatter all beams first.  This method allows
      % scattering individual beams.

      scattered = particle.scatter(ibeam);
      torque = scattered.torque();
    end

    function varargout = forcetorqueInternal(particle, ibeam, varargin)
      % Calculate the force-torque using the beam method
      %
      % This method is not used by :meth:`forcetorque`.  For most beams it
      % is better to scatter all beams first.  This method allows
      % scattering individual beams.

      scattered = particle.scatter(ibeam);
      [varargout{1:nargout}] = scattered.forcetorque();
    end
  end
end

