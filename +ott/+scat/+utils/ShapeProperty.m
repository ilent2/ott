classdef ShapeProperty < ott.scat.utils.Particle
% Declares a shape property and dependent properties for shapes.
% Inherits from :class:`ott.scat.utils.Particle`.
%
% Use this class when the scattering method (particle) can describe
% multiple geometric shapes.  If the particle only describes one type of
% shape it may be better to simply inherit from the Shape class.
%
% Properties
%   - shape         -- The encapsulated shape object.
%
% Dependent properties
%   - position      -- Uses the shapes location.
%   - rotation      -- Uses the shapes rotation.
%
% Dependent methods
%   - rotate        -- Uses the shapes rotate method.
%
% Methods
%   - force       -- Calculate force on particle in beam
%   - torque      -- Calculate torque on particle in beam
%   - forcetorque -- Calculate force and torque on particle in beam
%   - validateShape   -- Method called to validate the shape
%
% Abstract methods
%   - scatter     -- Calculate scattered beam
%   - forceInternal       -- Method called by `force`
%   - torqueInternal      -- Method called by `torque`
%   - forcetorqueInternal -- Method called by `forcetorque`

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  properties
    shape         % The encapsulated shape object.
  end

  properties (Dependent)
    position      % Uses the shapes location.
    rotation      % Uses the shapes rotation.
  end

  methods
    function particle = ShapeProperty(shape)
      % Construct and initialise the shape property
      %
      % Usage
      %   particle = ShapeProperty(shape)

      particle.shape = shape;
    end

    function particle = rotate(particle, varargin)
      % Apply rotation to the particle (and shape)
      %
      % See also :class:`ott.shapes.Shape/rotate`.

      particle.shape = particle.shape.rotate(varargin{:});
    end
  end

  methods (Hidden)
    function shape = validateShape(particle, shape)
      % Method called to validate the shape.
      %
      % If your scattering method requires a specific shape,
      % overload this method with the additional conditions or casts.

      assert(isa(shape, 'ott.shapes.Shape'), ...
          'shape must be of type ''ott.shapes.Shape''');
    end
  end

  methods % Getters/setters
    function particle = set.shape(particle, val)
      particle.shape = particle.validateShape(val);
    end

    function position = get.position(particle)
      position = particle.position;
    end
    function particle = set.position(particle, val)
      particle.shape.position = val;
    end

    function rotation = get.rotation(particle)
      rotation = particle.rotation;
    end
    function particle = set.rotation(particle, val)
      particle.shape.rotation = val;
    end
  end
end

