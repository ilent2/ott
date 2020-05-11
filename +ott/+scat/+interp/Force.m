classdef Force < ott.scat.utils.NoBeam & ott.utils.RotationPositionProp
% Method to estimate forces from interpolated force data.
% Inherits from :class:`ott.scat.utils.ZeroScattered`.
%
% Interpolated particles do not scatter optical beams, hence they
% inherit from :class:`ZeroScattered`.  This class uses local
% interpolation to estimate the force around a set of positions.
%
% Properties
%   - interpolants     -- Underlying interpolate objects

% TODO: Support for torques and scattering
%   Should this all be in one class or separate classes?

% TODO: Support for casting from Particle to scattered interpolant
%   This should be implemented in ott.scat.utils.Particle (probably)

  properties
    interpolants
  end

  methods
    function particle = Force(position, force, varargin)
      % Construct a new force interpolation particle
      %
      % Usage
      %   particle = Force(position, force, ...)
      %
      % Parameters
      %   - position (3xN numeric)
      %   - force    (3xN numeric)

      fx = scatteredInterpolant(position.', force(1, :).');
      fy = scatteredInterpolant(position.', force(2, :).');
      fz = scatteredInterpolant(position.', force(3, :).');
      particle.interpolants = {fx, fy, fz};

    end
  end

  methods (Hidden)
    function f = forceInternal(particle, beam, varargin)
      % Evaluate interpolation

      f = [particle.interpolants{1}(beam.position.');
           particle.interpolants{2}(beam.position.');
           particle.interpolants{3}(beam.position.')];
    end

    function t = torqueInternal(particle, beam, varargin)
      % Constant zero torque
      t = zeros(3, numel(beam));
    end

    function [f, t] = forcetorqueInternal(particle, beam, varargin)
      % Calculate force and torque

      f = particle.forceInternal(beam, varargin{:});
      t = particle.torqueInternal(beam, varargin{:});
    end
  end
end
