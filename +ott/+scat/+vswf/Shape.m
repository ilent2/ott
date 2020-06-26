classdef Shape < ott.scat.vswf.Tmatrix
% Constructs a T-matrix from a geometric shape.
% Inherits from :class:`ott.scat.vswf.Tmatrix`.
%
% Unlike other T-matrix classes, this class has a variable shape.
% Changing the shape causes the T-matrix to be reconstructed using
% either a user selected method or a default method.
%
% Properties
%   - shape     -- Shape to construct T-matrix
%   - relativeMedium   -- Material to use for the shape
%   - tmatrixMethod    -- Method to use for T-matrix calculation
%
% Static methods
%   - SmartCylinder    -- Smart method choice for cylindrical particles
%
% Additional properties inherited from :class:`Tmatrix`.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Hidden)
    shapeInternal           % Internal shape property
    relativeMediumInternal  % Internal relativeMedium property
    tmatrixMethodInternal   % Internal tmatrixMethod property
  end

  properties (Dependent)
    shape                   % Geometry (causes update)
    relativeMedium          % Optical properties (causes  update)
    tmatrixMethod           % Method for calculation (causes update)
  end

  methods (Static)
    function tmatrix = SmartCylinder(varargin)
      % Constructs a T-matrix for a cylinder using smart method selection
      %
      % Either uses DDA, EBCM or Pointmatch depending on the cylinder's
      % aspect ratio and material.  Pointmatching is the preferred method,
      % otherwise defaults to EBCM or DDA if error tolerance is not met.
      %
      % Uses the results from
      %
      %   Qi et al.,
      %   Optics Letters Vol. 39, Issue 16, pp. 4827-4830 (2014)
      %   https://doi.org/10.1364/OL.39.004827
      %
      % Usage
      %   tmatrix = SmartCylinder(shape, relativeMedium, ...)
      %
      % Optional named arguments
      %   - tolerance (enum) -- Error tolerance, can either be
      %     'one' or 'ten' for approximately 1% and 10% contours
      %     from the paper. Default: ``'ten'``.
      %
      % Unmatched arguments are passed to the appropriate constructor.

      p = inputParser;
      p.addOptional('shape', [], @(x) isa(x, 'ott.shapes.Shape'));
      p.addOptional('relativeMedium', [], ...
          @(x) isa(x, 'ott.beam.medium.Relative'));
      p.addParameter('wavelength', 1.0, @(x) isnumeric(x) & isscalar(x));
      p.addParameter('tolerance', 'ten');
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      args = ott.utils.unmatchedArgs(p);

      shape = p.Results.shape;
      if ~isa(shape, 'ott.shapes.Cylinder')
        shape = ott.shapes.Cylinder(shape);
      end

      assert(numel(shape) == 1, 'shape must be a single shape');
      shape = shape ./ p.Results.wavelength;

      relativeMedium = p.Results.relativeMedium;

      % Construct argument list
      args = [{'shape', shape, 'relativeMedium', relativeMedium}, args];

      % EBCM 1% data
      ebcm1 = {};
      ebcm1.x = [73, 169, 198, 228, 261, 391, 586, 718, 718, ...
          657, 523, 457, 262, 73];
      ebcm1.y = [409, 406, 418, 423, 397, 412, 400, 375, 223, ...
          193, 195, 165, 204, 390];
      ebcm1.x = (ebcm1.x - 73) * 2.0 / (718 - 73);
      ebcm1.y = -(ebcm1.y - 438) * 6.0 / (438 - 9);

      % PM 1% data
      pm1 = {};
      pm1.x = [297, 355, 394, 718, 718, 591, 525, 391, 361, 297];
      pm1.y = [943, 933, 946, 894, 868, 846, 874, 864, 913, 913];
      pm1.x = (pm1.x - 73) * 2.0 / (718 - 73);
      pm1.y = -(pm1.y - 985) * 6.0 / (985 - 555);

      % EBCM 10% data
      ebcm10 = {};
      ebcm10.x = [73, 193, 718, 718, 525, 328, 229, 160, 73];
      ebcm10.y = [430, 426, 381, 37, 94, 177, 214, 274, 375];
      ebcm10.x = (ebcm10.x - 73) * 2.0 / (718 - 73);
      ebcm10.y = -(ebcm10.y - 438) * 6.0 / (438 - 9);

      % PM 10% data
      pm10 = {};
      pm10.x = [130, 160, 328, 397, 462, 522, 589, 718, 718, ...
          654, 589, 522, 328, 265, 130];
      pm10.y = [961, 970, 967, 951, 946, 946, 925, 912, 753, ...
          784, 798, 798, 865, 874, 948];
      pm10.x = (pm10.x - 73) * 2.0 / (718 - 73);
      pm10.y = -(pm10.y - 985) * 6.0 / (985 - 555);

      % Paper parameters were in physical units
      lambda = 1.064;
      diameter = lambda * shape.radius * 2;
      len = shape.height;

      switch p.Results.tolerance
        case 'ten'
          if inpolygon(diameter, len, pm10.x, pm10.y)
            method = 'pm';
          elseif inpolygon(diameter, len, ebcm10.x, ebcm10.y)
              method = 'ebcm';
          else
              method = 'dda';
          end
        case 'one'
          if inpolygon(diameter, len, pm1.x, pm1.y)
            method = 'pm';
          elseif inpolygon(diameter, len, ebcm1.x, ebcm1.y)
            method = 'ebcm';
          else
            method = 'dda';
          end
        otherwise
          error('tolerance must be ''one'' or ''ten''');
      end

      % If shape is anisotropic use DDA
      if ~relativeMedium.isIsotropic
        method = 'dda';
      end

      assert(~relativeMedium.isMagnetic, ...
        'Magnetic particles not supported, permeability must be 1');

      switch method
        case 'dda'
          tmatrix = ott.scat.vswf.Dda.FromShape(args{:});
        case 'ebcm'
          tmatrix = ott.scat.vswf.Ebcm.FromAxisymInterpShape(args{:});
        case 'pm'
          tmatrix = ott.scat.vswf.Pointmatch.FromStarShape(args{:});
        otherwise
          error('Internal error');
      end
    end
  end

  methods
    function tmatrix = Shape(varargin)
      % Construct a T-matrix shape wrapper object
      %
      % This object has variable properties for the shape, material and
      % method which, when changed, cause the T-matrix data to be updated.
      %
      % Usage
      %   tmatrix = Shape(shape, relativeMedium, tmatrixMethod, ...)
      %
      % Parameters
      %   - shape (ott.shapes.Shape) -- Description of particle geometry.
      %
      %   - relativeMedium (ott.beam.medium.Relative) -- Description of
      %     particle optical properties.
      %
      %   - tmatrixMethod (function_handle) -- Method to use for
      %     calculating the T-matrix data.  Must take ``shape`` and
      %     ``relativeMedium`` as the first two arguments and return
      %     a valid T-matrix (i.e., it supports class constructors,
      %     user defined functions and ``FromShape``-like static methods).
      %     Default: ``@ott.scat.vswf.TmatrixShape.SphereMaxRadius``.
      %
      % Optional named arguments
      %   - wavelength (numeric) -- Used to convert `shape` input to
      %     relative units, i.e. `radius_rel = radius ./ wavelength`.
      %     This parameter not used for setting the T-matrix material.
      %     Default: ``1.0`` (i.e., `shape` is already in relative units).

      p = inputParser;
      p.addOptional('shape', [], @(x) isa(x, 'ott.shapes.Shape'));
      p.addOptional('relativeMedium', [], ...
          @(x) isa(x, 'ott.beam.medium.Relative'));
      p.addOptional('tmatrixMethod', ...
          @ott.scat.vswf.TmatrixShape.SphereMaxRadius, ...
          @(x) isa(x, 'function_handle'));
      p.addParameter('wavelength', 1.0);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      tmatrix = tmatrix@ott.scat.vswf.Tmatrix(unmatched{:});
      tmatrix = tmatrix.changeProperties(...
          p.Results.shape, p.Results.relativeMedium, ...
          p.Results.tmatrixMethod, 'wavelength', p.Results.wavelength);
    end

    function tmatrix = changeProperties(tmatrix, varargin)
      % Change the T-matrix properties and re-calculate the data
      %
      % This method is more efficient than changing the properties
      % individually as it only re-calculate the T-matrix once.
      %
      % Usage
      %   tmatrix = tmatrix.changeProperties(shape, medium, method)
      %
      % For parameters, see class constructor (:meth:`Shape`).

      p = inputParser;
      p.addOptional('shape', [], @(x) isa(x, 'ott.shapes.Shape'));
      p.addOptional('relativeMedium', [], ...
          @(x) isa(x, 'ott.beam.medium.Relative'));
      p.addOptional('tmatrixMethod', ...
          @ott.scat.vswf.TmatrixShape.SphereMaxRadius, ...
          @(x) isa(x, 'function_handle'));
      p.addParameter('wavelength', 1.0);
      p.parse(varargin{:});

      % Check output arguments
      tmatrix.nargoutCheck(nargout);

      S = warning('off', 'ott:scat:vswf:Shape:property_change');
      tmatrix.shapeInternal = p.Results.shape ./ p.Results.wavelength;
      tmatrix.relativeMediumInternal = p.Results.relativeMedium;
      tmatrix.tmatrixMethodInternal = p.Results.tmatrixMethod;
      warning(S);

      % Update the T-matrix
      tmatrix = tmatrix.calculateData();
    end

    function tmatrix = calculateData(tmatrix)
      % Calculate data for the T-matrix
      %
      % This function is called when visible properties are changed.

      % Check output arguments
      tmatrix.nargoutCheck(nargout);

      % Construct the new T-matrix data
      newTmatrix = tmatrix.tmatrixMethod(tmatrix.shape, ...
          tmatrix.relativeMaterial);
      assert(isa(newTmatrix, 'ott.scat.vswf.Tmatrix') ...
          && numel(newTmatrix) == 1, ...
          'tmatrixMethod must return a single T-matrix');

      % Assign to the Shape object
      tmatrix.data = newTmatrix.data;
      tmatrix = tmatrix.setType(newTmatrix.type);
    end
  end

  methods (Hidden)
    function shape = getGeometry(tmatrix, wavelength)
      shape = tmatrix.shape .* wavelength;
    end
  end

  methods % Getters/setters
    function val = get.shape(tmatrix)
      val = tmatrix.shapeInternal;
    end
    function tmatrix = set.shape(tmatrix, val)
      S = warning('off', 'ott:scat:vswf:Shape:property_change');
      tmatrix.shapeInternal = val;
      warning(S);
      tmatrix = tmatrix.calculateData();
    end

    function val = get.relativeMedium(tmatrix)
      val = tmatrix.relativeMediumInternal;
    end
    function tmatrix = set.relativeMedium(tmatrix, val)
      S = warning('off', 'ott:scat:vswf:Shape:property_change');
      tmatrix.relativeMediumInternal = val;
      warning(S);
      tmatrix = tmatrix.calculateData();
    end

    function val = get.tmatrixMethod(tmatrix)
      val = tmatrix.tmatrixMethodInternal;
    end
    function tmatrix = set.tmatrixMethod(tmatrix, val)
      S = warning('off', 'ott:scat:vswf:Shape:property_change');
      tmatrix.tmatrixMethodInternal = val;
      warning(S);
      tmatrix = tmatrix.calculateData();
    end

    function tmatrix = set.shapeInternal(tmatrix, val)
      warning('ott:scat:vswf:Shape:property_change', ...
          'Changing this Hidden property doesn''t update the T-matrix data');
      assert(isa(val, 'ott.shapes.Shape'), ...
          'shape must be a ott.shapes.Shape');
      tmatrix.shapeInternal = val;
    end

    function tmatrix = set.relativeMediumInternal(tmatrix, val)
      warning('ott:scat:vswf:Shape:property_change', ...
          'Changing this Hidden property doesn''t update the T-matrix data');
      assert(isa(val, 'ott.beam.medium.Relative'), ...
          'relatveMedium must be a ott.beam.medium.Relative');
      tmatrix.relativeMediumInternal = val;
    end

    function tmatrix = set.tmatrixMethodInternal(tmatrix, val)
      warning('ott:scat:vswf:Shape:property_change', ...
          'Changing this Hidden property doesn''t update the T-matrix data');
      assert(isa(val, 'function_handle'), ...
          'tmatrixMethod must be a function_handle');
      tmatrix.tmatrixMethodInternal = val;
    end
  end
end

