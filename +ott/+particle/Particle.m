classdef Particle < matlab.mixin.Heterogeneous ...
    & ott.utils.RotationPositionProp
% Base class for particles in optical tweezers simulations.
% Inherits from :class:`ott.utils.RotationPositionProp` and
% `matlab.mixin.Hetrogeneous`.
%
% This class combined the optical scattering methods, drag calculation
% methods and other properties required to simulate the dynamics of a
% particle in an optical tweezers simulation.  In future version of OTT,
% this class may change to support multiple scattering methods.
%
% This is an abstract class.  For instances of this class see
% :class:`Variable` and :class:Fixed`.
%
% Abstract properties
%   - shape       -- Geometric shape representing the particle
%   - tmatrix     -- Describes the scattering (particle-beam interaction)
%   - drag        -- Describes the drag (particle-fluid interaction)
%   - tinternal   -- T-matrix for internal scattered field (optional)
%   - mass        -- Particle mass [kg] (optional)
%
% Properties
%   - position    -- Particle position [m]
%   - rotation    -- Particle orientation (3x3 rotation matrix)
%
% Methods
%   - surf        -- Uses the shape surf method to visualise the particle
%   - setMassFromDensity -- Calcualte mass from homogeneous density

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (Abstract)
    shape        % Geometric shape representing the particle
    tmatrix      % Describes the scattering (particle-beam interaction)
    drag         % Describes the drag (particle-fluid interaction)
    tinternal    % T-matrix for internal scattered field (optional)
  end

  properties
    mass         % Particle mass [kg] (optional)
  end

  methods
    function varargout = surf(particle, varargin)
      % Generate visualisation of the shape using the shape.surf method.
      %
      % Applies the particle rotation and translation to the shape
      % before visualisation.
      %
      % Usage
      %   [...] = particle.surf(...)
      %
      % For full details and usage, see :meth:`+ott.+shape.Shape.surf`.

      oshape = particle.shape;
      oshape = oshape.rotate(particle.rotation);
      oshape.position = oshape.position + particle.position;
      [varargout{1:nargout}] = oshape.surf(varargin{:});
    end

    function particle = setMassFromDensity(particle, density)
      % Uses the particle shape and a homogeneous density to calculate mass
      %
      % Particle shape must have a valid volume method.
      %
      % Usage
      %   particle = particle.setMassFromDensity(density)
      %
      % Parameters
      %   - density (numeric) -- Homogeneous density [kg/m^3]

      assert(~isempty(particle.shape), ...
          'particle shape must be set');
      assert(isnumeric(density) && isscalar(density) && density >= 0, ...
          'density must be positive numeric scalar');

      particle.mass = density * particle.shape.volume;
    end
  end

  methods % Getters/setters
    function particle = particle.mass(particle, val)
      if isempty(val)
        particle.mass = [];
      else
        assert(isnumeric(val) && isscalar(val) && val >= 0, ...
            'mass must be positive numeric scalar or empty');
        particle.mass = val;
      end
    end
  end
end


