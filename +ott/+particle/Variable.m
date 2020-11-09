classdef Variable < ott.particle.Particle
% A particle whose drag/tmatrix are automatically recomputed.
% Inherits from :class:`ott.particle.Particle`.
%
% Changing the shape or refractive index of this particle causes the
% T-matrix and drag data to be re-calculated.  This is useful for modelling
% particles with time-varying properties.
%
% Properties
%   - shape             -- Shape describing the object (changes tmatrix/drag)
%   - relative_index    -- Particle refractive index (changes tmatrix)
%   - tmatrix_method    -- Method for T-matrix calculation
%   - tinternal_method  -- Internal T-matrix calculation methdo
%   - drag_method       -- Method for drag calculation
%
% Dependent properties
%   - drag              -- Description of drag properties
%   - tmatrix           -- Description of optical scattering properties
%   - tinternal         -- Internal T-matrix (optional)
%
% Methods
%   - setProperties     -- Set shape and relative_index simultaneously
%
% Static methods
%   - FromShape         -- Create instance using FromShape of Tmatirx/Stokes
%   - Sphere            -- Construct instance using Sphere approximation
%   - StarShaped        -- Construct instance for star shaped particles

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    tmatrix_method     % Method for T-matrix calculation
    tinternal_method   % Internal T-matrix calculation methdo
    drag_method        % Method for drag calculation
  end

  properties (Hidden, SetAccess=protected)
    dragInternal        % Internal data for drag property
    tmatrixInternal     % Internal data for tmatrix property
    tinternalInternal   % Internal data for tinternal property
    shapeInternal       % Internal shape data
    riInternal          % Internal refractive index data
  end

  properties (Dependent)
    shape              % Shape describing the object (changes tmatrix/drag)
    relative_index     % Particle refractive index (changes tmatrix)
    drag               % Description of drag properties
    tmatrix            % Description of optical scattering properties
    tinternal          % Internal T-matrix (optional)
  end

  methods (Static)
    function particle = FromShape(varargin)
      % Construct a new Variable particle using the FromShape method.
      %
      % This function sets up the ``tmatrix_method``, ``tinternal_method``
      % and ``drag_method`` functions to use the ``FromShape`` methods
      % from :class:`ott.tmatrix.Tmatrix` and :class:`ott.drag.Stokes`.
      %
      % Usage
      %   particle = ott.particle.Variable.FromShape(...)
      %
      % Optional named arguments
      %   - wavelength_medium (numeric) -- Medium wavelength.  [m]
      %     Used for converting shape units [m] to wavelength units.
      %     Default: ``1064e-9`` (a common IR trapping wavelength).
      %
      %   - viscosity (numeric) -- Viscosity of medium. [Ns/m2]
      %     Default: ``8.9e-4`` (approximate viscosity of water).
      %
      %   - internal (logical) -- If the internal T-matrix calculation
      %     method should also be set.  Default: ``false``.
      %
      % Unmatched parameters passed to class constructor.

      p = inputParser;
      p.addParameter('wavelength_medium', 1064e-9);
      p.addParameter('viscosity', 8.9e-4);
      p.addParameter('internal', false);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      wavelength = p.Results.wavelength_medium;

      tmatrix_method = @(shape, ri, ~) ott.tmatrix.Tmatrix.FromShape(...
          shape ./ wavelength, 'relative_index', ri);

      if p.Results.internal
        tinternal_method = @(shape, ri, ~) ott.tmatrix.Tmatrix.FromShape(...
            shape ./ wavelength, 'relative_index', ri, 'internal', true);
      else
        tinternal_method = [];
      end

      drag_method = @(shape, ~) ott.drag.Stokes.FromShape(shape, ...
          'viscosity', p.Results.viscosity);

      particle = ott.particle.Variable(...
          'tmatrix_method', tmatrix_method, ...
          'tinternal_method', tinternal_method, ...
          'drag_method', drag_method, unmatched{:});
    end

    function particle = Sphere(varargin)
      % Construct a variable particle for a spherical geometry
      %
      % Uses the ``FromShape`` methods from :class:`ott.tmatrix.Mie` and
      % :class:`ott.drag.StokesSphere`.
      %
      % Usage
      %   particle = ott.particle.Variable.Sphere(...)
      %
      % Optional named arguments
      %   - wavelength_medium (numeric) -- Medium wavelength.  [m]
      %     Used for converting shape units [m] to wavelength units.
      %     Default: ``1064e-9`` (a common IR trapping wavelength).
      %
      %   - viscosity (numeric) -- Viscosity of medium. [Ns/m2]
      %     Default: ``8.9e-4`` (approximate viscosity of water).
      %
      %   - internal (logical) -- If the internal T-matrix calculation
      %     method should also be set.  Default: ``false``.
      %
      % Unmatched parameters passed to class constructor.

      p = inputParser;
      p.addParameter('wavelength_medium', 1064e-9);
      p.addParameter('viscosity', 8.9e-4);
      p.addParameter('internal', false);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      wavelength = p.Results.wavelength_medium;

      tmatrix_method = @(shape, ri, ~) ott.tmatrix.Mie.FromShape(...
          shape ./ wavelength, 'relative_index', ri);

      if p.Results.internal
        tinternal_method = @(shape, ri, ~) ott.tmatrix.Mie.FromShape(...
            shape ./ wavelength, 'relative_index', ri, 'internal', true);
      else
        tinternal_method = [];
      end

      drag_method = @(shape, ~) ott.drag.StokesSphere.FromShape(shape, ...
          'viscosity', p.Results.viscosity);

      particle = ott.particle.Variable(...
          'tmatrix_method', tmatrix_method, ...
          'tinternal_method', tinternal_method, ...
          'drag_method', drag_method, unmatched{:});
    end

    function particle = StarShaped(varargin)
      % Construct a variable particle for a star shaped geometry.
      %
      % Uses the ``FromShape`` methods from :class:`ott.tmatrix.Tmatrix` and
      % :class:`ott.drag.StokesStarShaped`.
      %
      % Usage
      %   particle = ott.particle.Variable.StarShaped(...)
      %
      % Optional named arguments
      %   - wavelength_medium (numeric) -- Medium wavelength.  [m]
      %     Used for converting shape units [m] to wavelength units.
      %     Default: ``1064e-9`` (a common IR trapping wavelength).
      %
      %   - viscosity (numeric) -- Viscosity of medium. [Ns/m2]
      %     Default: ``8.9e-4`` (approximate viscosity of water).
      %
      %   - internal (logical) -- If the internal T-matrix calculation
      %     method should also be set.  Default: ``false``.
      %
      % Unmatched parameters passed to class constructor.

      p = inputParser;
      p.addParameter('wavelength_medium', 1064e-9);
      p.addParameter('viscosity', 8.9e-4);
      p.addParameter('internal', false);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      wavelength = p.Results.wavelength_medium;

      tmatrix_method = @(shape, ri, ~) ott.tmatrix.Tmatrix.FromShape(...
          shape ./ wavelength, 'relative_index', ri);

      if p.Results.internal
        tinternal_method = @(shape, ri, ~) ott.tmatrix.Tmatrix.FromShape(...
            shape ./ wavelength, 'relative_index', ri, 'internal', true);
      else
        tinternal_method = [];
      end

      drag_method = @(shape, ~) ott.drag.StokesStarShaped.FromShape(shape, ...
          'viscosity', p.Results.viscosity);

      particle = ott.particle.Variable(...
          'tmatrix_method', tmatrix_method, ...
          'tinternal_method', tinternal_method, ...
          'drag_method', drag_method, unmatched{:});
    end
  end

  methods
    function particle = Variable(varargin)
      % Construct a new Variable particle instance.
      %
      % Usage
      %   particle = Variable(...)
      %
      % Optional named arguments
      %   - tmatrix_method ([] | function_handle) -- T-matrix calculation
      %     method. Signature: ``@(shape, ri, old_tmatrix)``.
      %     Default: ``[]``.
      %
      %   - drag_method ([] | function_handle) -- Drag calculation method.
      %     Signature: ``@(shape, old_drag)``.
      %     Default: ``[]``.
      %
      %   - tinternal_method ([] | function_handle) -- internal T-matrix
      %     method. Signature: ``@(shape, ri, old_tmatrix)``.
      %     Default: ``[]``.
      %
      %   - initial_shape ([] | ott.shape.Shape) -- Initial particle
      %     geometry.  Default: ``[]``.
      %
      %   - initial_relative_index (numeric) -- Initial particle relative
      %     refractive index.  Default: ``[]``.

      p = inputParser;
      p.addParameter('drag_method', []);
      p.addParameter('tinternal_method', []);
      p.addParameter('tmatrix_method', []);
      p.addParameter('initial_shape', []);
      p.addParameter('initial_relative_index', []);
      p.parse(varargin{:});

      particle.drag_method = p.Results.drag_method;
      particle.tmatrix_method = p.Results.tmatrix_method;
      particle.tinternal_method = p.Results.tinternal_method;
      particle = particle.setProperties(...
          p.Results.initial_shape, ...
          p.Results.initial_relative_index);
    end

    function particle = setProperties(particle, shape, relative_index)
      % Set particle shape and relative index simultaneously.
      %
      % The shape and relative index properties can be changed directly
      % using the corresponding class properties, however this causes
      % two updates of the T-matrix/drag data.  This method provided a
      % option for adjusting both properties simultaniously.
      %
      % Usage
      %   particle = particle.setProperties(shape, relative_index)
      %
      % Parameters
      %   - shape ([] | ott.shape.Shape) -- Particle geometry.
      %
      %   - relative_index ([] | numeric) -- Relative refractive index.

      assert(isnumeric(relative_index) || isempty(relative_index), ...
          'relative_index must be numeric or empty');
      assert(isa(shape, 'ott.shape.Shape') || isempty(shape), ...
          'shape must be empty or ot.shape.Stape');

      particle.riInternal = relative_index;
      particle.shapeInternal = shape;

      if ~isempty(particle.relative_index) && ~isempty(particle.shape)
        % Calculate new T-matrix and drag data

        if ~isempty(particle.tmatrix_method)
          particle.tmatrixInternal = particle.tmatrix_method(...
              particle.shape, particle.relative_index, ...
              particle.tmatrixInternal);
        end

        if ~isempty(particle.drag_method)
          particle.dragInternal = particle.drag_method(...
              particle.shape, particle.dragInternal);
        end

        if ~isempty(particle.tinternal_method)
          particle.tinernalInternal = particle.tinternal_method(...
              particle.shape, particle.relative_index, ...
              particle.tinternalInternal);
        end
      end
    end
  end

  methods % Getters/setters
    function particle = set.drag(particle, ~) %#ok<INUSD>
      error('ott:particle:variable:set_drag', ...
        'Cannot set drag properties of particle.Variable class');
    end
    function drag = get.drag(particle)
      drag = particle.dragInternal;
    end

    function particle = set.tmatrix(particle, ~) %#ok<INUSD>
      error('ott:particle:variable:set_tmatrix', ...
        'Cannot set tmatrix property of particle.Variable class');
    end
    function tmatrix = get.tmatrix(particle)
      tmatrix = particle.tmatrixInternal;
    end

    function particle = set.tinternal(particle, ~) %#ok<INUSD>
      error('ott:particle:variable:set_tinternal', ...
        'cannot set tinternal property of particle.Variable class');
    end
    function tinternal = get.tinternal(particle)
      tinternal = particle.tinternalInternal;
    end

    function particle = set.shape(particle, val)
      particle = particle.setProperties(val, particle.relative_index);
    end
    function shape = get.shape(particle)
      shape = particle.shapeInternal;
    end

    function particle = set.relative_index(particle, val)
      particle = particle.setProperties(particle.shape, val);
    end
    function ri = get.relative_index(particle)
      ri = particle.riInternal;
    end

    function particle = set.drag_method(particle, val)
      assert(isempty(val) || isa(val, 'function_handle'), ...
          'drag_method must be empty or a function handle');
      particle.drag_method = val;
    end

    function particle = set.tmatrix_method(particle, val)
      assert(isempty(val) || isa(val, 'function_handle'), ...
          'tmatrix_method must be empty or a function handle');
      particle.tmatrix_method = val;
    end

    function particle = set.tinternal_method(particle, val)
      assert(isempty(val) || isa(val, 'function_handle'), ...
          'tinternal_method must be empty or a function handle');
      particle.tinternal_method = val;
    end
  end
end


