classdef StokesStarShaped < ott.drag.Stokes ...
    & ott.drag.mixin.VarViscosity ...
    & ott.drag.mixin.CalcInvDrag
% Abstract base class for star shaped particle methods in an unbounded medium.
% Inherits from :class:`Stokes`.
%
% Properties
%   - inverse     -- Calculated from forward
%   - viscosity   -- Viscosity of medium
%   - shape       -- A star shaped particle describing the geometry
%
% Abstract properties
%   - forwardInternal    -- Drag tensor calculated by method
%
% Static methods
%   FromShape     -- Defers to StokesLambNn (may change in future)
%
% Supported casts
%   - ott.shape.Shape -- Retrieves the geometry (see `shape`)

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    shape       % Geometry for drag calculation
  end

  methods (Static)
    function drag = FromShape(shape, varargin)
      % Construct drag tensor for star shaped particle
      %
      % Usage
      %   StokesStarShaped.FromShape(shape, ...)
      %
      % Passes all arguments to :class:`StokesLambNn`.

      drag = ott.drag.StokesLambNn(shape, varargin{:});
    end
  end

  methods
    function drag = StokesStarShaped(varargin)
      % Construct base class for star-shaped particle methods
      %
      % Usage
      %   drag = drag@ott.drag.StokesStarShaped(shape, viscosity, ...)
      %
      % Parameters
      %   - shape (ott.shape.Shape) -- A shape object with a starShaped
      %     property with the `true` value.
      %   - viscosity (numeric) -- Viscosity of medium (default: 1.0)
      %
      % Parameters can also be passed as named arguments.
      % Unmatched parameters are passed to :class:`Stokes`.

      p = inputParser;
      p.addOptional('shape', [], @(x) isa(x, 'ott.shape.Shape'));
      p.addOptional('viscosity', 1.0, @isnumeric);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      drag = drag@ott.drag.Stokes(unmatched{:});
      drag.shape = p.Results.shape;
      drag.viscosity = p.Results.viscosity;
    end

    function shape = ott.shape.Shape(drag)
      % Retrieve the internal geometry
      shape = drag.shape;
      shape.position = [0;0;0];
    end
  end

  methods % Getters/setters
    function drag = set.shape(drag, val)
      assert(isa(val, 'ott.shape.Shape') && val.starShaped, ...
          'shape must be a star-shaped ott.shape.Shape instance');
      drag.shape = val;
    end
  end
end

