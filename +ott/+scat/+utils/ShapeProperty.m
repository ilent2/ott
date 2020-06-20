classdef ShapeProperty
% Declares a shape property and dependent properties for shapes.
%
% Use this class when the scattering method (particle) can describe
% multiple geometric shapes.  If the particle only describes one type of
% shape it may be better to simply inherit from the Shape class.
%
% Does not implement a class constructor.  The shape property is
% not initialised and must be set in the sub-class constructor.
%
% Properties
%   - shape         -- The encapsulated shape object.
%
% Dependent properties
%   - positionInternal      -- Uses the shapes location.
%   - rotationInternal      -- Uses the shapes rotation.
%
% Methods
%   - validateShape   -- Method called to validate the shape

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  properties
    shape         % The encapsulated shape object.
  end

  properties (Dependent, Hidden)
    positionInternal      % Uses the shapes location.
    rotationInternal      % Uses the shapes rotation.
  end

  methods (Hidden)
    function shape = validateShape(~, shape)
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

    function position = get.positionInternal(particle)
      position = particle.shape.position;
    end
    function particle = set.positionInternal(particle, val)
      particle.shape.position = val;
    end

    function rotation = get.rotationInternal(particle)
      rotation = particle.shape.rotation;
    end
    function particle = set.rotationInternal(particle, val)
      particle.shape.rotation = val;
    end
  end
end

