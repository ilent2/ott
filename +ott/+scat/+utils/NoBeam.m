classdef NoBeam < ott.scat.utils.ZeroScattered
% Provides force and torque methods which don't require a beam.
% Inherits from :class:`ZeroScattered`.
%
% This class is useful for methods that produce the same result
% regardless of if a beam is supplied (i.e. interpolation methods).
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

  methods
    function varargout = force(particle, varargin)
      % Calculate the force from a beam on a particle.
      %
      % Usage
      %   force = particle.force(beam, ...)
      %
      %   force = particle.force(...)
      %   Same as above, the beam argument is ignored.
      %
      % For details on usage/arguments see :meth:`forcetorque`.

      p = inputParser;
      p.addOptional('beam', ott.beam.abstract.Empty, ...
          @(x) isa(x, 'ott.beam.abstract.Beam'));
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      [varargout{1:nargout}] = force@ott.scat.utils.ZeroScattered(...
          particle, p.Results.beam, unmatched{:});
    end

    function varargout = torque(particle, varargin)
      % Calculate the torque from a beam on a particle.
      %
      % Usage
      %   torque = particle.torque(beam, ...)
      %
      %   torque = particle.torque(...)
      %   Same as above, an empty beam is used as the default.
      %
      % For details on usage/arguments see :meth:`forcetorque`.

      p = inputParser;
      p.addOptional('beam', ott.beam.abstract.Empty, ...
          @(x) isa(x, 'ott.beam.abstract.Beam'));
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      [varargout{1:nargout}] = torque@ott.scat.utils.ZeroScattered(...
          particle, p.Results.beam, unmatched{:});
    end

    function varargout = forcetorque(particle, varargin)
      % Calculate the force and torque from a beam on a particle.
      %
      % Usage
      %   [force, torque] = particle.forcetorque(beam, ...)
      %
      %   [force, torque] = particle.forcetorque(...)
      %   Same as above, an empty beam is used as the default.
      %
      % Optional named arguments
      %   - position (3xN numeric) -- Positions of the particle to
      %     calculate properties for.
      %     Default: ``[]``.
      %
      %   - rotation (3x3N numeric) -- Orientations of the particle
      %     to calculate properties for.
      %     Default: ``[]``.

      p = inputParser;
      p.addOptional('beam', ott.beam.abstract.Empty, ...
          @(x) isa(x, 'ott.beam.abstract.Beam'));
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      [varargout{1:nargout}] = forcetorque@ott.scat.utils.ZeroScattered(...
          particle, p.Results.beam, unmatched{:});
    end
  end
end
