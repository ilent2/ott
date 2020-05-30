classdef (Abstract) Particle < ott.utils.RotationPositionEntity
% Describes the abstract interface for particle.
% Inherits from :class:`ott.utils.RotationPositionEntity`.
%
% Particle-centric scattering methods should inherit from this class.
% Particles which describe geometric shapes can choose to inherit from
% this class and the geometric shape class, or use the
% :class:`ShapeProperty` to declare a `shape` property in the class.
%
% Abstract properties
%   - positionInternal    -- Position of the particle [3xN]
%   - rotationInternal    -- Orientation of the particle [3x3]
%
% Methods
%   - force       -- Calculate force on particle in beam
%   - torque      -- Calculate torque on particle in beam
%   - forcetorque -- Calculate force and torque on particle in beam
%   - scatter     -- Calculate scattered beam
%
% Hidden methods
%   - forcetorqueInternal -- Method called by `forcetorque`
%
% Abstract methods
%   - rotate              -- Apply rotation to particle
%   - forceInternal       -- Method called by `force`
%   - torqueInternal      -- Method called by `torque`
%   - scatterInternal     -- Method called by `scatter`

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  methods (Abstract, Hidden)
    scatterInternal      % Calculate scattered beam
    forceInternal        % Calculate force on particle in beam
    torqueInternal       % Calculate torque on particle in beam
  end

  methods
    function varargout = force(particle, beam, varargin)
      % Calculate the force from a beam on a particle.
      %
      % Usage
      %   force = particle.force(beam, ...)
      %
      % For details on usage/arguments see :meth:`forcetorque`.

      [varargout{1:nargout}] = ...
          ott.utils.RotationPositionProp.rotationPositionHelper(...
          @particle.forceInternal, beam, 'prxcat', 2, varargin{:});
    end

    function varargout = torque(particle, beam, varargin)
      % Calculate the torque from a beam on a particle.
      %
      % Usage
      %   torque = particle.torque(beam, ...)
      %
      % For details on usage/arguments see :meth:`forcetorque`.

      [varargout{1:nargout}] = ...
          ott.utils.RotationPositionProp.rotationPositionHelper(...
          @particle.torqueInternal, beam, 'prxcat', 2, varargin{:});
    end

    function varargout = forcetorque(particle, beam, varargin)
      % Calculate force and torque from a beam on a particle.
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

      [varargout{1:nargout}] = ...
          ott.utils.RotationPositionProp.rotationPositionHelper(...
          @particle.forcetorqueInternal, beam, 'prxcat', 2, varargin{:});
    end

    function varargout = scatter(particle, beam, varargin)
      % Calculate beams scattered by a particle
      %
      % Usage
      %   sbeam = particle.scatter(ibeam, ...)
      %
      % Parameters and outputs
      %   - ibeam (ott.beam.abstract.Beam) -- Incident beam object.
      %     A copy of the incident beam is typically contained inside
      %     the scattered beam.  The copy may be modified from the
      %     original, for example, when translations and rotations
      %     are applied to the beam.
      %
      %   - sbeam (ott.beam.abstract.Scattered) -- The scattered beam
      %     or an array of scattered beams.
      %
      % Named parameters
      %   - position (3xN numeric) -- Position for the particle.
      %     Applied by translating the beam in the opposite direction.
      %     Default: ``[]``.
      %
      %   - rotation (3x3N numeric) -- Rotation for the particle.
      %     Applied by rotating the beam in the opposite direction.
      %     Inverse rotation is applied after scattering, effectively
      %     rotating the particle in the beam.  Default: ``[]``.
      %
      %   - array_type (enum) -- Array type for scattered beam with respect
      %     to rotation and position arguments.  Default: ``'array'``.
      %     Only used when position/rotation are arrays.

      p = inputParser;
      p.addParameter('array_type', 'array');
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Scatter the beam after applying position/rotation
      [varargout{1:nargout}] = ...
          ott.utils.RotationPositionProp.rotationPositionHelper(...
          @particle.scatterInternal, beam, 'prxcat', 2, unmatched{:});

      % Set the output array type
      if isa(varargout{1}, 'ott.beam.utils.ArrayType')
        varargout{1}.array_type = p.Results.array_type;
      end
    end
  end

  methods (Hidden)
    function [force, torque] = forcetorqueInternal(particle, varargin)
      % Calculate force and torque on particle in beam
      %
      % The default implementation simply calls ``force`` and ``torque``.

      force = particle.force(varargin{:});
      torque = particle.torque(varargin{:});
    end
  end
end

