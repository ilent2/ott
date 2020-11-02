classdef Particle < matlab.mixin.Heterogeneous ...
    & ott.utils.RotationPositionProp
% Base class describing particles in optical tweezers simulations.
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
%
% Properties
%   - position    -- Particle position [m]
%   - rotation    -- Particle orientation (3x3 rotation matrix)
%
% Methods
%   - surf        -- Uses the shape surf method to visualise the particle

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (Abstract)
    shape        % Geometric shape representing the particle
    tmatrix      % Describes the scattering (particle-beam interaction)
    drag         % Describes the drag (particle-fluid interaction)
    tinternal    % T-matrix for internal scattered field (optional)
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
  end
end

