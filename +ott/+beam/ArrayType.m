classdef ArrayType < ott.beam.Beam
% Provides methods for beams with array data
% Inherits from :class:`Beam`.
%
% This class declares an ``arrayType`` property and modified visualisation
% functions to handle combining the result of incoherent beams.
% Not all visualisations are compatible with incoherent array types.
%
% Properties
%   - arrayType     -- Type of array (either coherent/incoherent/array)
%
% Methods
%   - containsIncoherent    -- Returns true if the has incoherent elements

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    arrayType     % Type of array (coherent/incoherent/array)
  end

  methods
    function beam = ArrayType(varargin)
      % Construct a new array type.
      %
      % This class is abstract, call this method from a sub-class.
      %
      % Parameters
      %   - arrayType (enum) -- Either coherent/incoherent/array
      %     Default: ``'array'``.
      %
      % All other parameters passed to base.

      p = inputParser;
      p.addParameter('arrayType', 'array');
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.Beam(unmatched{:});
      beam.arrayType = p.Results.arrayType;
    end

    function b = containsIncoherent(beam)
      % Check if the beam contains any incoherent components
      %
      % Usage
      %   b = baem.containsIncoherent()

      b = strcmpi(beam.arrayType, 'incoherent');
    end

    function [imout, XY, vswfData] = visNearfield(beam, varargin)
      % Create a visualisation of the beam
      %
      % If the beam is incoherent, combines the third dimension of the
      % data after applying :meth:`VisualisationData`.
      %
      % Usage
      %   beam.visNearfield(...) -- display an image of the beam in
      %   the current figure window.
      %
      %   [im, XY, data] = beam.visualise(...) -- returns a image of the beam.
      %   If the beam object is an array, returns an image for each beam.
      %   Also returns a :class:`ott.utils.VswfData` structure for fast
      %   repeated calculations.
      %
      % Optional named arguments
      %   - size (2 numeric) -- Number of rows and columns in image.
      %     Default: ``[80, 80]``.
      %
      %   - field (enum) -- Type of visualisation to generate, see
      %     :meth:`VisualisationData` for valid options.
      %     Default: ``'irradiance'``.
      %
      %   - axis (enum|cell) -- Describes the slice to visualise.
      %     Can either be 'x', 'y' or 'z' for a plane perpendicular to
      %     the x, y and z axes respectively; or a cell array with 2
      %     or 3 unit vectors (3x1) for x, y, [z] directions.
      %     Default: ``'z'``.
      %
      %   - offset (numeric) -- Plane offset along axis (default: 0.0)
      %
      %   - range (numeric|cell) -- Range of points to visualise.
      %     Can either be a cell array { x, y }, two scalars for
      %     range [-x, x], [-y, y] or 4 scalars [ x0, x1, y0, y1 ].
      %     Default: ``[1, 1].*beam.wavelength``.
      %
      %   - mask (function_handle) Describes regions to remove from the
      %     generated image.  Function should take one argument for the
      %     3xN field xyz field locations and return a logical array mask.
      %     Values with ``false`` are removed.  Default: ``[]``.
      %
      %   - plot_axes (axes handle) -- Axes to place the visualisation in.
      %     If empty, no visualisation is generated.
      %     Default: ``gca()`` if ``nargout == 0`` otherwise ``[]``.
      %
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.
      %
      % Unmatched parameters are passed to the field calculation function.
      
      assert(numel(beam) == 1, ...
          'Expected only a single beam.  Use ott.beam.Array for arrays');

      [imout, XY, vswfData] = visNearfield@ott.beam.Beam(...
          beam, varargin{:});

      if strcmpi(beam.arrayType, 'incoherent')
        imout = sum(imout, 3);
      end

      % Generate visualisation
      beam.visShowPlot(imout, XY, nargout, varargin{:});

      if nargout == 0
        clear imout XY vswfData
      end
    end

    function [imout, XY, vswfData] = visFarfield(beam, varargin)
      % Create a visualisation of the beam by projecting the far-field
      % onto a plane.
      %
      % If the beam is incoherent, combines the third dimension of the
      % data after applying :meth:`VisualisationData`.
      %
      % Usage
      %   beam.visualiseFarfield(...) displays an image of the beam
      %   in the current axes.
      %
      %   [im, XY, data] = beam.visualise(...) returns a image of the beam.
      %   If the beam object is an array, returns an image for each beam.
      %
      % Optional named arguments
      %   - size (2 numeric) -- Number of rows and columns in image.
      %     Default: ``[80, 80]``.
      %
      %   - field (enum) -- Type of visualisation to generate, see
      %     :meth:`VisualisationData` for valid options.
      %     Default: ``'irradiance'``.
      %
      %   - range (numeric|cell) -- Range of points to visualise.
      %     Can either be a cell array { x, y }, two scalars for
      %     range [-x, x], [-y, y] or 4 scalars [ x0, x1, y0, y1 ].
      %     Default: ``[1, 1]``.
      %
      %   - direction (enum|2 numeric|3x3 numeric) -- Direction to visualise.
      %     Either 'pos', 'neg' for the positive and negative z hemispheres,
      %     [roty, rotz] for the rotation angles (in degrees), or
      %     a 3x3 rotation matrix to apply to the coordinates.
      %     Default: ``'pos'``.
      %
      %   - mapping (enum) -- Mapping from theta-phi to far-field.
      %     Must be one of 'sin', 'tan' or 'theta'.
      %     Default: ``'sin'``.
      %
      %   - plot_axes (axes handle) -- Axes to place the visualisation in.
      %     If empty, no visualisation is generated.
      %     Default: ``gca()`` if ``nargout == 0`` otherwise ``[]``.
      %
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.
      %
      %   - basis (enum) -- Vector spherical wave function basis to
      %     visualise.  Must be one of 'incoming' or 'outgoing'.
      %     Default: ``'incoming'``.

      [imout, XY, vswfData] = visFarfield@ott.beam.Beam(...
          beam, varargin{:});

      if strcmpi(beam.arrayType, 'incoherent')
        imout = sum(imout, 3);
      end

      % Generate visualisation
      beam.visShowPlot(imout, XY, nargout, varargin{:});

      if nargout == 0
        clear imout XY vswfData
      end
    end

    function [imout, XYZ, vswfData] = visFarfieldSphere(beam, varargin)
      % Generate a spherical surface visualisation of the far-field
      %
      % If the beam is incoherent, combines the third dimension of the
      % data after applying :meth:`VisualisationData`.
      %
      % Usage
      %   beam.visualiseFarfieldSphere(...)
      %   Generate a visualisation of the far-field in the current axes.
      %
      %   [I, XYZ, data] = beam.visualiseFarfieldSphere(...)
      %   Outputs the field value and three coordinate matrices that
      %   can be passed to ``surf(XYZ{1}, XYZ{2}, XYZ{3})``
      %
      % Optional named arguments
      %   - field (enum) -- Type of visualisation to generate, see
      %     :meth:`VisualisationData` for valid options.
      %     Default: ``'irradiance'``.
      %
      %   - normalise (logical) -- If the field value should be normalized.
      %     Default: ``false``.
      %
      %   - type (enum) -- Type of surface visualisation.
      %     Can be either 'sphere' or '3dpolar'.
      %     Default: ``'sphere'``.
      %
      %   - npts (numeric) -- Number of points for sphere surface.
      %     Passed to ``sphere(npts)``.  Default: ``100``.
      %
      %   - plot_axes (axes handle) -- Axes to place the visualisation in.
      %     If empty, no visualisation is generated.
      %     Default: ``gca()`` if ``nargout == 0`` otherwise ``[]``.
      %
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.
      %
      %   - basis (enum) -- Vector spherical wave function basis to
      %     visualise.  Must be one of 'incoming' or 'outgoing'.
      %     Default: ``'incoming'``.

      [imout, XYZ, vswfData] = visFarfieldSphere@ott.beam.Beam(...
          beam, varargin{:});

      if strcmpi(beam.arrayType, 'incoherent')
        imout = sum(imout, 3);
      end

      % Generate visualisation
      beam.visFarfieldSpherePlot(imout, XYZ, nargout, varargin{:});

      if nargout == 0
        clear imout XYZ vswfData
      end
    end

    function [imout, ptheta, vswfData] = visFarfieldSlice(beam, varargin)
      % Generate a 2-D slice through the far-field
      %
      % If the beam is incoherent, combines the third dimension of the
      % data after applying :meth:`VisualisationData`.
      %
      % Usage
      %   beam.visualiseFarfieldSlice(phi, ...)
      %   Generates a 2-D slice at angle phi around the z-axis.
      %   If `phi` is not prsent, uses default of ``[0, pi/2]``.
      %   Plots into the current axes.
      %
      %   [im, theta, data] = beam.visualiseFarfieldSlice(...)
      %   Outputs the calculated values and corresponding angles.
      %
      % Optional named arguments
      %   - field (enum) -- Type of visualisation to generate, see
      %     :meth:`VisualisationData` for valid options.
      %     Default: ``'irradiance'``.
      %
      %   - normalise (logical) -- If the field value should be normalized.
      %     Default: ``false``.
      %
      %   - npts (numeric) -- Number of points for sphere surface.
      %     Passed to ``sphere(npts)``.  Default: ``100``.
      %
      %   - plot_axes (axes handle) -- Axes to place the visualisation in.
      %     If empty, no visualisation is generated.
      %     Default: ``gca()`` if ``nargout == 0`` otherwise ``[]``.
      %
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.
      %
      %   - basis (enum) -- Vector spherical wave function basis to
      %     visualise.  Must be one of 'incoming' or 'outgoing'.
      %     Default: ``'incoming'``.

      [imout, ptheta, vswfData] = visFarfieldSlice@ott.beam.Beam(...
          beam, varargin{:});

      if strcmpi(beam.arrayType, 'incoherent')
        imout = sum(imout, 3);
      end

      % Generate visualisation
      beam.visFarfieldSlicePlot(imout, ptheta, nargout, varargin{:});

      if nargout == 0
        clear imout ptheta vswfData
      end
    end
  end

  methods % Getters/setters
    function beam = set.arrayType(beam, val)
      assert(sum(strcmpi(val, {'coherent', 'incoherent', 'array'})) == 1, ...
          'arrayType must be coherent/incoherent/array');
      beam.arrayType = val;
    end
  end
end

