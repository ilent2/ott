classdef Fixed < ott.particle.Particle
% A particle with stored drag/tmatrix properties.
% Inherits from :class:`ott.particle.Particle`.
%
% Properties
%   - drag        -- Description of drag properties
%   - tmatrix     -- Description of optical scattering properties
%   - shape       -- Description of the geometry [m]
%   - tinternal   -- Internal T-matrix (optional)
%
% Static methods
%   - FromShape   -- Construct a particle from a shape description

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    drag         % Description of drag properties
    tmatrix      % Description of optical scattering properties
    shape        % Description of the geometry [m]
    tinternal    % Internal T-matrix (optional)
  end

  methods (Static)
    function particle = FromShape(shape, varargin)
      % Construct a particle from a shape description.
      %
      % Usage
      %   particle = Fixed.FromShape(shape, ...)
      %
      % Parameters
      %   - shape (ott.shape.Shape) -- Particle geometry. [m]
      %
      % Named parameters
      %   - index_relative (numeric) -- Relative refractive index of particle.
      %     Default: ``[]``.  Must be supplied or computable from
      %     `index_particle` and `index_medium`.
      %
      %   - index_particle (numeric) -- Refractive index in particle.
      %     Default: ``[]``.
      %
      %   - index_medium (numeric) -- Refractive index in medium.
      %     Default: ``1.0`` unless both `index_relative` and `index_particle`
      %     are supplied.  In which case `index_medium` is ignored.
      %
      %   - wavelength0 (numeric) -- Vacuum wavelength.  [m]
      %     Default: ``1064e-9`` (a common IR trapping wavelength).
      %
      %   - viscosity (numeric) -- Viscosity of medium. [Ns/m2]
      %     Default: ``8.9e-4`` (approximate viscosity of water).
      %
      %   - internal (logical) -- If the internal T-matrix should also
      %     be computed.  Default: ``false``.
      %
      %   - mass (numeric) -- Particle mass.  Default: ``[]``.

      p = inputParser;
      p.addParameter('index_relative', []);
      p.addParameter('index_particle', []);
      p.addParameter('index_medium', 1.0);
      p.addParameter('viscosity', 8.9e-4);
      p.addParameter('wavelength0', 1064e-9);
      p.addParameter('internal', false);
      p.addParameter('mass', []);
      p.parse(varargin{:});

      % Get index_relative
      index_relative = p.Results.index_relative;
      index_particle = p.Results.index_particle;
      assert(~(isempty(index_relative) && isempty(index_particle)), ...
          'Either or both index_relative and index_particle must be supplied');
      if isempty(index_relative)
        index_relative = index_particle ./ p.Results.index_medium;
      end

      % Get wavelength_medium
      if isempty(index_particle)
        index_medium = p.Results.index_medium;
      else
        index_medium = index_particle ./ index_relative;
      end
      wavelength_medium = p.Results.wavelength0 ./ index_medium;

      % Calculate internal T-matrix
      tmatrix = ott.tmatrix.Tmatrix.FromShape(...
          shape ./ wavelength_medium, 'relative_index', index_relative);

      % Calculate drag
      drag = ott.drag.Stokes.FromShape(shape, ...
          'viscosity', p.Results.viscosity);

      tinternal = [];
      if p.Results.internal
        tinternal = ott.tmatrix.Tmatrix.FromShape(...
            shape ./ wavelength_medium, ...
            'relative_index', index_relative, ...
            'internal', true, 'mass', mass);
      end

      % Construct particle
      particle = ott.particle.Fixed(shape, drag, tmatrix, ...
          'tinternal', tinternal);
    end
  end

  methods
    function particle = Fixed(varargin)
      % Construct a new particle instance with fixed properties
      %
      % Usage
      %   particle = Fixed(shape, drag, tmatrix, ...)
      %
      % Optional named arguments
      %   - shape (ott.shape.Shape) -- Geometry description.  Units: m.
      %
      %   - drag (ott.drag.Stokes) -- Drag tensor description.
      %
      %   - tmatrix (ott.tmatrix.Tmatrix) -- Scattering description.
      %
      %   - tinternal (ott.tmatrix.Tmatrix) -- Internal T-matrix.
      %
      %   - mass (numeric) -- Particle mass.  Default: ``[]``.

      p = inputParser;
      p.addOptional('shape', [], @(x) isa(x, 'ott.shape.Shape'));
      p.addOptional('drag', [], @(x) isa(x, 'ott.drag.Stokes'));
      p.addOptional('tmatrix', [], @(x) isa(x, 'ott.tmatrix.Tmatrix'));
      p.addParameter('tinternal', []);
      p.addParameter('mass', []);
      p.parse(varargin{:});

      particle.drag = p.Results.drag;
      particle.shape = p.Results.shape;
      particle.tmatrix = p.Results.tmatrix;
      particle.tinternal = p.Results.tinternal;
      particle.mass = p.Results.mass;
    end
  end

  methods % Getters/setters
    function particle = set.shape(particle, val)
      if isempty(val)
        particle.shape = val;
      else
        assert(isa(val, 'ott.shape.Shape'), ...
          'shape must be a ott.shape.Shape instance');
        particle.shape = val;
      end
    end

    function particle = set.drag(particle, val)
      if isempty(val)
        particle.drag = [];
      else
        assert(isa(val, 'ott.drag.Stokes'), ...
          'drag must be a ott.drag.Stokes instance');
        particle.drag = val;
      end
    end

    function particle = set.tmatrix(particle, val)
      if isempty(val)
        particle.tmatrix = [];
      else
        assert(isa(val, 'ott.tmatrix.Tmatrix'), ...
          'tmatrix must be a ott.tmatrix.Tmatrix instance');
        particle.tmatrix = val;
      end
    end

    function particle = set.tinternal(particle, val)
      if isempty(val)
        particle.tinternal = [];
      else
        assert(isa(val, 'ott.tmatrix.Tmatrix'), ...
            'tinternal must be empty or ott.tmatrix.Tmatrix instance');
        particle.tinternal = val;
      end
    end
  end
end


