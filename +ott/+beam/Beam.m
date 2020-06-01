classdef (Abstract) Beam < ott.beam.properties.Beam
% Abstract base class for beam approximations.
% Inherits from :class:`ott.beam.utils.BeamInterface`.
%
% Field calculation methods
%   - efield     -- Calculate electric field around the origin
%   - hfield     -- Calculate magnetic field around the origin
%   - ehfield    -- Calculate electric and magnetic fields around the origin
%   - egradient  -- Gradient of electric field
%   - hgradient  -- Gradient of magnetic field
%   - ehgradient -- Gradient of electric and magnetic fields
%   - poynting   -- Calculate the Poynting vector field
%   - efarfield  -- Calculate electric fields in the far-field
%   - hfarfield  -- Calculate magnetic fields in the far-field
%   - ehfarfield -- Calculate electric and magnetic fields in the far-field
%   - poyntingFarfield -- Calculate the Poynting vector in the far-field
%   - eparaxial  -- Calculate electric fields in the paraxial far-field
%   - hparaxial  -- Calculate magnetic fields in the paraxial far-field
%   - ehparaxial -- Calculate electric and magnetic paraxial far-fields
%
% Force and torque related methods
%   - intensityMoment -- Calculate moment of beam intensity in the far-field
%   - force           -- Calculate the change in momentum between two beams
%   - torque          -- Calculate change in angular momentum between beams
%   - forcetorque     -- Calculate the force and the torque between beams
%
% Field visualisation methods
%   - visualise -- Generate a visualisation around the origin
%   - visualiseFarfield -- Generate a visualisation at the far-field
%   - visualiseFarfieldSlice -- Visualise the field on a angular slice
%   - visualiseFarfieldSphere -- Visualise the filed on a sphere
%
% Properties
%   - wavelength    -- Wavelength of beam in medium (default: 1.0)
%
% Dependent properties
%   - wavenumber    -- Wave-number of beam in medium
%
% Abstract properties
%   - power       -- The power of the beam (may be infinite)
%
% Abstract methods
%   - efieldInternal    -- Called by efield
%   - hfieldInternal    -- Called by hfield
%   - efarfieldInternal -- Called by efarfield
%   - hfarfieldInternal -- Called by hfarfield
%   - getBeamPower      -- get method called by dependent property power
%
% Hidden methods
%   - egradientInternal   -- Called by egradient (estimates numerically)
%   - hgradientInternal   -- Called by hgradient (estimates numerically)
%   - ehgradientInternal  -- Called by ehgradient (estimates numerically)
%   - forceInternal       -- Force calculation method (to be overridden)
%   - torqueInternal      -- Torque calculation method (to be overridden)
%   - forcetorqueInternal -- Force/Torque calculation method (to be overridden)
%
% Supported casts
%   - Ray                 -- Constructs a Ray set from the far-field
%   - vswf.Bsc            -- Default Bsc cast, uses vswf.FarfieldPm
%   - vswf.Pointmatch     -- Default Bsc.Pm cast, uses vswf.FarfieldPm
%   - vswf.FarfieldPm
%   - vswf.NearfieldPm
%   - vswf.BesselBasis
%   - vswf.PlaneBasis
%   - vswf.LgParaxialBasis
%   - abstract.InterpFarfield
%   - abstract.InterpNearfield2d
%   - abstract.InterpNearfield3d

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

% TODO: Add support for re-usable visualisation data

  methods (Abstract)
    efieldInternal(obj)      % Method called by efield
    hfieldInternal(obj)      % Method called by hfield
    efarfieldInternal(obj)   % Called by efarfield
    hfarfieldInternal(obj)   % Called by hfarfield
  end

  methods (Static)
    function data = VisualisationData(field_type, field)
      % Helper to generate the visualisation data output.
      % This function is not intended to be called directly, instead
      % see :meth:`visualise` or :meth:`visualiseFarfield`.
      %
      % Usage
      %   VisualisationData(field_type, fieldvector)
      %   Takes a field_type string and a :class:`ott.utils.FieldVector`
      %   with the coordinates and field vectors.
      %
      % Parameters
      %   - fieldvector (3xN FieldVector) -- Field vector values
      %     and coordinate locations.
      %
      %   - field_type -- (enum) Type of field to calculate.
      %     Supported types include:
      %
      %     - irradiance  -- :math:`\sqrt{Ex^2 + Ey^2 + Ez^2}`
      %     - E2          -- :math:`Ex^2 + Ey^2 + Ez^2`
      %     - Sum(Abs(E)) -- :math:`|Ex| + |Ey| + |Ez|`
      %
      %     - Re(Er), Re(Et), Re(Ep), Re(Ex), Re(Ey), Re(Ez)
      %     - Abs(Er), Abs(Et), Abs(Ep), Abs(Ex), Abs(Ey), Abs(Ez)
      %     - Arg(Er), Arg(Et), Arg(Ep), Arg(Ex), Arg(Ey), Arg(Ez)

      assert(size(field, 1) == 3, 'fieldvector must be 3xN matrix');
      assert(isa(field, 'ott.utils.FieldVector'), ...
          'fieldvector must be a instance of FieldVector');

      % Generate the requested field
      switch field_type
        case 'irradiance'
          data = sqrt(sum(abs(field.vxyz).^2, 1));
        case 'E2'
          data = sum(abs(field.vxyz).^2, 1);
        case 'sum(Abs(E))'
          data = sum(abs(field.vxyz), 1);

        case 'Re(Er)'
          data = real(field.vrtp(1, :));
        case 'Re(Et)'
          data = real(field.vrtp(2, :));
        case 'Re(Ep)'
          data = real(field.vrtp(3, :));

        case 'Re(Ex)'
          data = real(field.vxyz(1, :));
        case 'Re(Ey)'
          data = real(field.vxyz(2, :));
        case 'Re(Ez)'
          data = real(field.vxyz(3, :));

        case 'Abs(Er)'
          data = abs(field.vrtp(1, :));
        case 'Abs(Et)'
          data = abs(field.vrtp(2, :));
        case 'Abs(Ep)'
          data = abs(field.vrtp(3, :));

        case 'Abs(Ex)'
          data = abs(field.vxyz(1, :));
        case 'Abs(Ey)'
          data = abs(field.vxyz(2, :));
        case 'Abs(Ez)'
          data = abs(field.vxyz(3, :));

        case 'Arg(Er)'
          data = angle(field.vrtp(1, :));
        case 'Arg(Et)'
          data = angle(field.vrtp(2, :));
        case 'Arg(Ep)'
          data = angle(field.vrtp(3, :));

        case 'Arg(Ex)'
          data = angle(field.vxyz(1, :));
        case 'Arg(Ey)'
          data = angle(field.vxyz(2, :));
        case 'Arg(Ez)'
          data = angle(field.vxyz(3, :));

        case 'Er'
          data = field.vrtp(1, :);
        case 'Et'
          data = field.vrtp(2, :);
        case 'Ep'
          data = field.vrtp(3, :);

        case 'Ex'
          data = field.vxyz(1, :);
        case 'Ey'
          data = field.vxyz(2, :);
        case 'Ez'
          data = field.vxyz(3, :);

        otherwise
          error('OTT:BSC:GetVisualisationData:unknown_field_type', ...
              'Unknown value for field_type.');
      end
    end
  end

  methods
    function beam = Beam(varargin)
      % Pass-through constructor for beam interface
      %
      % Usage
      %   beam = beam@ott.beam.utils.BeamInterface(...)
      %
      % All properties passed to ott.beam.properties.Beam

      beam = beam@ott.beam.properties.Beam(varargin{:});
    end

    function beam = ott.beam.Ray(beam, varargin)

      % TODO: Multiple beams and type check (look at abstract.Bessel)
      % TODO: Should FocussedTopHat have a specialisation of this?

      theta = linspace(0, beam.angle, 50);
      phi = linspace(0, 2*pi, 100);
      radius = 100;
      [R, T, P] = meshgrid(radius, theta, phi);
      rtp = [R(:), T(:), P(:)].';
      vrtp = [0*R(:), ones(size(T(:))), 0*P(:)].';
      [vxyz, xyz] = ott.utils.rtp2xyz(vrtp, rtp);

      Etp = beam.efarfield(rtp).vrtp(2:3, :);

      beam = ott.beam.Ray('origin', xyz, 'dierction', -xyz, ...
          'polarisation1', vxyz, 'field', Etp, varargin{:});
    end

    function beam = ott.beam.vswf.Pointmatch(varargin)
      % Cast to vswf.FarfieldPm
      beam = ott.beam.vswf.FarfieldPm(varargin{:});
    end

    function beam = ott.beam.vswf.Bsc(varargin)
      % Cast to vswf.Pointmatch
      beam = ott.beam.vswf.Pointmatch(varargin{:});
    end

    function beam = ott.beam.vswf.FarfieldPm(beam, varargin)
      % Construct a far-field point matched beam

      Nmax = 20;
      ntheta = 2*(Nmax+1);
      nphi = 2*(Nmax+1);

      % Sample point in far-field
      theta = linspace(0, pi, ntheta);
      phi = linspace(0, 2*pi, nphi);
      [T, P] = meshgrid(theta, phi);
      tp = [T(:), P(:)].';

      % For abstract beams, this invokes the Beam cast.
      E = beam.efarfield(tp);
      Etp = [E.vrtp(2, :).'; E.vrtp(3, :).'];

      ci = ott.utils.combined_index(Nmax, Nmax);
      [nn, mm] = ott.utils.combined_index(1:ci);
      beam = ott.beam.vswf.FarfieldPm(nn, mm, tp, Etp);

    end

    function beam = ott.beam.abstract.InterpFarfield(beam, varargin)
      % Construct a far-field interpolated beam

      % TODO: Multiple beams and type check

      % Calculate far-field pattern
      % TODO: Might be smarter choice
      theta = linspace(0, pi, 100);
      phi = linspace(0, 2*pi, 100);
      [T, P] = meshgrid(theta, phi);
      tp = [T(:), P(:)].';

      % For abstract beams, this invokes the Beam cast.
      [Etp, Htp] = beam.ehfarfield(tp);

      beam = ott.beam.abstract.InterpFarfield(tp, Etp, Htp, varargin{:});
    end

    function E = efield(beam, xyz, varargin)
      % Calculate E and H field
      %
      % Usage
      %   E = beam.efield(xyz)
      %   Calculates the fields at the specified locations (3xN matrix).
      %   Returns a :class:`ott.utils.FieldVector`.
      %
      % Optional named arguments
      %   - origin (enum) -- Specify coordinate origin.  Can be
      %     'local' or 'world'.  Default: ``'world'``.
      %
      % Other arguments are passed to efieldInternal.

      % Transform coordinates
      [xyz, unmatched] = beam.transformWorldToLocal(xyz, varargin{:});

      E = beam.efieldInternal(xyz, unmatched{:});
    end

    function H = hfield(beam, xyz, varargin)
      % Calculate E and H field
      %
      % Usage
      %   H = beam.hfield(xyz)
      %   Calculates the fields at the specified locations (3xN matrix).
      %   Returns a :class:`ott.utils.FieldVector`.
      %
      % Optional named arguments
      %   - origin (enum) -- Specify coordinate origin.  Can be
      %     'local' or 'world'.  Default: ``'world'``.
      %
      % Other arguments are passed to hfieldInternal.

      % Transform coordinates
      [xyz, unmatched] = beam.transformWorldToLocal(xyz, varargin{:});

      H = beam.hfieldInternal(xyz, unmatched{:});
    end

    function [E, H] = ehfield(beam, xyz, varargin)
      % Calculate E and H field
      %
      % Usage
      %   [E, H] = beam.ehfield(xyz)
      %   Calculates the fields at the specified locations (3xN matrix).
      %   Returns a :class:`ott.utils.FieldVector`.
      %
      % Optional named arguments
      %   - origin (enum) -- Specify coordinate origin.  Can be
      %     'local' or 'world'.  Default: ``'world'``.
      %
      % Other arguments are passed to ehfieldInternal.

      % Transform coordinates
      [xyz, unmatched] = beam.transformWorldToLocal(xyz, varargin{:});

      [E, H] = beam.ehfieldInternal(xyz, unmatched{:});
    end

    function gradE = egradient(beam, xyz, varargin)
      % Calculate the gradient of the E-field at the specified location.
      %
      % Usage
      %   gradE = beam.egradient(xyz)
      %   Calculates the fields at the specified locations (3xN matrix).
      %   Returns a 3x3xN matrix of gradient tensors.
      %
      % Optional named arguments
      %   - origin (enum) -- Specify coordinate origin.  Can be
      %     'local' or 'world'.  Default: ``'world'``.
      %
      % Other arguments are passed to egradientInternal.

      % Transform coordinates
      [xyz, unmatched] = beam.transformWorldToLocal(xyz, varargin{:});

      gradE = beam.egradientInternal(xyz, unmatched{:});
    end

    function gradH = hgradient(beam, xyz, varargin)
      % Calculate the gradient of the H-field at the specified location.
      %
      % Usage
      %   gradE = beam.hgradient(xyz)
      %   Calculates the fields at the specified locations (3xN matrix).
      %   Returns a 3x3xN matrix of gradient tensors.
      %
      % Optional named arguments
      %   - origin (enum) -- Specify coordinate origin.  Can be
      %     'local' or 'world'.  Default: ``'world'``.
      %
      % Other arguments are passed to hgradientInternal.

      % Transform coordinates
      [xyz, unmatched] = beam.transformWorldToLocal(xyz, varargin{:});

      gradH = beam.hgradientInternal(xyz, unmatched{:});
    end

    function [gradE, gradH] = ehgradient(beam, xyz, varargin)
      % Calculate the gradient of the E-field at the specified location.
      %
      % Usage
      %   gradE = beam.ehgradient(xyz)
      %   Calculates the fields at the specified locations (3xN matrix).
      %   Returns a 3x3xN matrix of gradient tensors.
      %
      % Optional named arguments
      %   - origin (enum) -- Specify coordinate origin.  Can be
      %     'local' or 'world'.  Default: ``'world'``.
      %
      % Other arguments are passed to ehgradientInternal.

      % Transform coordinates
      [xyz, unmatched] = beam.transformWorldToLocal(xyz, varargin{:});

      [gradE, gradH] = beam.ehgradientInternal(xyz, unmatched{:});
    end

    function S = poynting(beam, xyz, varargin)
      % Calculate the Poynting vector
      %
      % This function calculates::
      %
      %     S = E \times H
      %
      % Usage
      %   S = beam.poynting(xyz)
      %   Calculates the Poynting vector at the locations (3xN matrix)
      %   Returns a :class:`ott.utils.FieldVector`.

      [E, H] = beam.ehfield(xyz, varargin{:});
      S = ott.utils.FieldVector(xyz, cross(E.vxyz, H.vxyz), 'cartesian');
    end

    function S = poyntingFarfield(beam, rtp, varargin)
      % Calculate the Poynting vector in the far-field
      % See :meth:`poynting` for further information.
      %
      % Usage
      %   S = beam.poyntingFarfield(rtp)
      %   Returns a :class:`ott.utils.FieldVector`.

      [E, H] = beam.ehfarfield(rtp, varargin{:});
      xyz = ott.utils.rtp2xyz(E.locations);
      S = ott.utils.FieldVector(xyz, cross(E.vxyz, H.vxyz), 'cartesian');
    end

    function E = efarfield(beam, rtp, varargin)
      % Calculate E and H field
      %
      % Usage
      %   E = beam.efarfield(rtp)
      %   Calculates the fields at the specified locations.
      %   Returns a :class:`ott.utils.FieldVector`.
      %
      % Parameters
      %   - rtp (2xN | 3xN numeric) -- Far-field locations.
      %     radial coordinate is ignored by most methods.

      E = beam.efarfieldInternal(rtp, varargin{:});
    end

    function H = hfarfield(beam, rtp, varargin)
      % Calculate E and H field
      %
      % Usage
      %   H = beam.hfarfield(rtp)
      %   Calculates the fields at the specified locations.
      %   Returns a :class:`ott.utils.FieldVector`.
      %
      % See :meth:`efarfield` for parameters.

      H = beam.hfarfieldInternal(rtp, varargin{:});
    end

    function [E, H] = ehfarfield(beam, rtp, varargin)
      % Calculate E and H far-fields
      %
      % Usage
      %   [E, H] = beam.ehfarfield(rtp)
      %   Calculates the fields at the specified locations.
      %   Returns a :class:`ott.utils.FieldVector`.
      %
      % See :meth:`efarfield` for parameters.

      % Default implementation just calls each part, but for some
      % beams it may be more efficient to calculate these together,
      % in which case this method should be over-written.
      E = beam.efarfield(rtp, varargin{:});
      H = beam.hfarfield(rtp, varargin{:});
    end

    function E = eparaxial(beam, nxyz, varargin)
      % Calculates the far-fields with a paraxial projection
      %
      % Usage
      %   E = beam.eparaxial(nxyz, ...)
      %   Calculates the fields at the normalized coordinates.
      %   Returns a :class:`ott.utils.FieldVector`.
      %
      % Parameters
      %   - nxyz (2xN|3xN numeric) -- Normalized locations of the
      %     far-field points.  z value is ignored, x and y values
      %     should be between [-1, 1].
      %
      % Optional named arguments
      %   - mapping (enum) -- Mapping from theta-phi to far-field.
      %     Must be one of 'sin', 'tan' or 'theta'.
      %     Default: ``'sin'``.
      %
      %   - direction (enum) -- Beam hemisphere to calculate.
      %     Must be 'neg' or 'pos'.  Default: ``pos``.

      % Convert normalized coordinates to rtp
      rtp = ott.beam.properties.FarfieldMapping.paraxial2farfield(...
          nxyz, varargin{:});

      % Calculate fields
      E = beam.efarfield(rtp, varargin{:});
    end

    function H = hparaxial(beam, nxyz, varargin)
      % Calculates the far-fields with a paraxial projection
      %
      % Usage
      %   H = beam.hparaxial(nxyz, ...)
      %   Calculates the fields at the normalized coordinates.
      %   Returns a :class:`ott.utils.FieldVector`.
      %
      % See :meth:`eparaxial` for parameters.

      % Convert normalized coordinates to rtp
      rtp = ott.beam.properties.FarfieldMapping.paraxial2farfield(...
          nxyz, varargin{:});

      % Calculate fields
      H = beam.hfarfield(rtp, varargin{:});
    end

    function [E, H] = ehparaxial(beam, nxyz, varargin)
      % Calculates the far-fields with a paraxial projection
      %
      % Usage
      %   [E, H] = beam.ehparaxial(nxyz, ...)
      %   Calculates the fields at the normalized coordinates.
      %   Returns a :class:`ott.utils.FieldVector`.
      %
      % See :meth:`eparaxial` for parameters.

      % Convert normalized coordinates to rtp
      rtp = ott.beam.properties.FarfieldMapping.paraxial2farfield(...
          nxyz, varargin{:});

      % Calculate fields
      E = beam.efarfield(rtp, varargin{:});
      H = beam.hfarfield(rtp, varargin{:});
    end

    function varargout = visualise(beam, varargin)
      % Create a visualisation of the beam
      %
      % Usage
      %   beam.visualise(...) displays an image of the beam in the current
      %   figure window.
      %
      %   im = beam.visualise(...) returns a image of the beam.
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
      %     Default: ``[]``.
      %
      %   - plot_axes (axes handle) -- Axes to place the visualisation in.
      %     If empty, no visualisation is generated.
      %     Default: ``gca()`` if ``nargout == 0`` otherwise ``[]``.
      %
      % Unmatched arguments are passed to the field function.

      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('field', 'irradiance');
      p.addParameter('size', []);
      p.addParameter('axis', 'z');
      p.addParameter('offset', 0.0);
      p.addParameter('range', [1, 1].*max([beam(:).wavelength]));
      p.addParameter('mask', []);
      p.addParameter('plot_axes', []);
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Get range and size of data
      default_sz = [80, 80];
      [xrange, yrange, sz] = beam.visualiseGetRange(p.Results.range, ...
          p.Results.size, default_sz);

      % Generate grid of positions
      assert(isscalar(p.Results.offset) && isnumeric(p.Results.offset), ...
          'offset must be numeric scalar');
      [xx, yy, zz] = meshgrid(xrange, yrange, p.Results.offset);

      % Change the grid to the user requested orientation
      [xyz, labels] = beam.visualiseGetXyz(xx, yy, zz, p.Results.axis);

      % Calculate mask
      if isempty(p.Results.mask)
        mask = true(sz);
      else
        mask = p.Results.mask(xyz);
      end
      
      % Calculate which points are needed
      our_xyz = xyz(:, mask(:));
     
      % Calculate visualisation data
      imout = beam.visualiseImoutXyzMasked(sz, our_xyz, ...
          p.Results.field, mask, unmatched{:});

      % Display the visualisation
      ott.beam.Beam.visualiseShowPlot(...
          nargout, p.Results.plot_axes, imout, ...
          {xrange, yrange}, labels);

      % Assign outputs if requested
      if nargout == 1
        varargout{1} = imout;
      end
    end

    function varargout = visualiseFarfield(beam, varargin)
      % Create a visualisation of the beam by projecting the far-field
      % onto a plane.
      %
      % Usage
      %   beam.visualiseFarfield(...) displays an image of the beam
      %   in the current axes.
      %
      %   im = beam.visualise(...) returns a image of the beam.
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
      % Unmatched arguments are passed to the field function.

      % Parse arguments
      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('field', 'irradiance');
      p.addParameter('size', []);
      p.addParameter('direction', 'pos');
      p.addParameter('range', [1, 1]);
      p.addParameter('mapping', 'sin');
      p.addParameter('plot_axes', []);
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Generate grid of coordinates
      default_sz = [80, 80];
      [xrange, yrange, sz] = beam.visualiseGetRange(p.Results.range, ...
          p.Results.size, default_sz);
      [xx, yy] = meshgrid(xrange, yrange);
      nxyz = [xx(:), yy(:), zeros(size(xx(:)))].';

      % Convert coordinates to spherical
      rtp = ott.beam.properties.FarfieldMapping.paraxial2farfield(...
          nxyz, 'mapping', p.Results.mapping);
      rtp(imag(rtp) ~= 0) = nan;

      % Apply direction coordinate transformation
      if ischar(p.Results.direction)
        if strcmpi(p.Results.direction, 'pos')
          % Nothing to do
        elseif strcmpi(p.Results.direction, 'neg')
          rtp(2, :) = pi - rtp(2, :);
        else
          error('Unknonw direction char string');
        end
      elseif isnumeric(p.Results.direction)
        if all(size(p.Results.direction) == [3, 3])
          R = p.Results.direction;
        elseif numel(p.Results.direction) == 2
          R = ott.utils.roty(p.Results.direction(1));
          R = ott.utils.roty(p.Results.direction(2)) * R;
        else
          error('direction must be char, 3x3 numeric or 2-numeric')
        end

        xyz = R * ott.utils.rtp2xyz(rtp);
        rtp = ott.utils.xyz2rtp(xyz);
      end

      % Calculate field data
      imout = beam.visualiseImoutRtp(sz, rtp, p.Results.field, unmatched{:});

      % Display the visualisation
      ott.beam.Beam.visualiseShowPlot(...
          nargout, p.Results.plot_axes, imout, ...
          {xrange, yrange}, {'Direction 1', 'Direction 2'});

      % Assign outputs if requested
      if nargout == 1
        varargout{1} = imout;
      end
    end

    function varargout = visualiseFarfieldSphere(beam, varargin)
      % Generate a spherical surface visualisation of the far-field
      %
      % Usage
      %   beam.visualiseFarfieldSphere(...)
      %   Generate a visualisation of the far-field in the current axes.
      %
      %   [I, X, Y, Z] = beam.visualiseFarfieldSphere(...)
      %   Outputs the field value and three coordinate matrices that
      %   can be passed to ``surf(X, Y, Z)``.
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
      % Unmatched arguments are passed to the field function.

      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('field', 'irradiance');
      p.addParameter('normalise', false);
      p.addParameter('type', 'sphere');
      p.addParameter('npts', 100);
      p.addParameter('plot_axes', []);
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Generate grid
      [X,Y,Z] = sphere(p.Results.npts);

      % Calculate field data
      rtp = ott.utils.xyz2rtp(X, Y, Z);
      imout = beam.visualiseImoutRtp(size(X), rtp, ...
          p.Results.field, unmatched{:});

      % Generate visualisation
      if nargout == 0 || ~isempty(p.Results.plot_axes)

        % Get the default axes
        our_axes = p.Results.plot_axes;
        if isempty(our_axes)
          our_axes = gca();
        end

        % Check that we only have one visualisation to show
        if ~ismatrix(imout)
          warning('ott:beam:Beam:non_matrix_plot', ...
              'Multiple beams generated, only showing first');
        end

        I = imout(:, :, 1);
        if p.Results.normalise
          I = I ./ max(abs(I(:)));
        end

        switch p.Results.type
          case 'sphere'
            surf(our_axes, X, Y, Z, I, ...
                'facecolor','interp','edgecolor','none');
          case '3dpolar'
            surf(our_axes, abs(I).*X,abs(I).*Y,abs(I).*Z,I,...
              'facecolor','interp','edgecolor','none');
          otherwise
            error('Unknown visualisation type');
        end

        zlabel(our_axes, 'Z');
        xlabel(our_axes, 'X');
        ylabel(our_axes, 'Y');
        view(our_axes, 50, 20);
        axis(our_axes, 'equal');
      end

      % Assign output data
      if nargout >= 1
        varargout{1} = imout;
        if nargout == 4
          varargout{2:4} = {X, Y, Z};
        end
      end
    end

    function varargout = visualiseFarfieldSlice(beam, phi, varargin)
      % Generate a 2-D slice through the far-field
      %
      % Usage
      %   beam.visualiseFarfieldSlice(phi, ...)
      %   Generates a 2-D slice at angle phi around the z-axis.
      %   Plots into the current axes.
      %
      %   im = beam.visualiseFarfieldSlice(...)
      %   Outputs the calculated values.
      %
      %   [theta, im] = beam.visualiseFarfieldSlice(...)
      %   Also outputs the corresponding angles.
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
      % Unmatched arguments are passed to the field function.

      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('field', 'irradiance');
      p.addParameter('normalise', false);
      p.addParameter('npts', 100);
      p.addParameter('plot_axes', []);
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Generate grid
      ptheta = linspace(0, 2*pi, p.Results.npts);
      [r, theta, phi] = ott.utils.matchsize(0, ptheta(:), phi);
      rtp = [r(:), theta(:), phi(:)].';

      % Calculate fields
      imout = beam.visualiseImoutRtp(p.Results.npts, rtp, ...
          p.Results.field, unmatched{:});

      % Generate visualisation
      if nargout == 0 || ~isempty(p.Results.plot_axes)

        % Get the default axes
        our_axes = p.Results.plot_axes;
        if isempty(our_axes)
          our_axes = gca();
        end

        % Check that we only have one visualisation to show
        if ~ismatrix(imout)
          warning('ott:beam:Beam:non_matrix_plot', ...
              'Multiple beams generated, only showing first');
        end

        I = imout(:, :, 1);
        if p.Results.normalise
          I = I ./ max(abs(I(:)));
        end

        % Generate plot
        polaraxes(our_axes);
        polarplot(ptheta, I);
      end

      % Assign outputs
      if nargout == 2
        varargout{1:2} = {theta, imout};
      elseif nargout == 1
        varargout{1} = imout;
      end
    end

    function [moments, ints] = intensityMoment(beam, varargin)
      % Calculate moment of the beam intensity in the far-field
      %
      % Usage
      %   [moment, int] = beam.intensityMoment(...)
      %
      % Optional named arguments
      %   - theta_range (2 numeric) -- Range of angles to integrate over.
      %     Default: ``[0, pi]``.
      %
      %   - ntheta (numeric) -- Number of theta points.  (Default: 100)
      %   - nphi (numeric) -- Number of phi points.  (Default: 100)
      %
      % Unmatched arguments are passed to the field function.

      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('theta_range', [0, pi]);
      p.addParameter('ntheta', 100);
      p.addParameter('nphi', 100);
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Setup grid
      [theta, phi] = ott.utils.angulargrid(p.Results.ntheta, p.Results.nphi);
      dtheta = diff(theta(1:2));
      dphi = diff(phi(1:2));

      % Truncate the theta range
      keep = theta > p.Results.theta_range(1) & theta < p.Results.theta_range(2);
      theta = theta(keep);
      phi = phi(keep);

      rtp = [ones(numel(theta), 1), theta(:), phi(:)].';

      % So integrals match sign convention used in ott.forcetorque
      uxyz = ott.utils.rtp2xyz(rtp);
      uxyz(3, :) = -uxyz(3, :);

      imout = beam.intensityMomentImout(4, rtp, uxyz, theta, phi, dtheta, dphi, unmatched{:});

      moments = imout(1:3, :);
      ints = imout(4, :);
    end

    function varargout = force(beam, other, varargin)
      % Calculate change in linear momentum between beams.
      %
      % Usage
      %   force = beam.force(other_beam, ...)
      %
      %   force = beam.force(particle, ...)
      %   Use the particle's force method with flipped position/rotation
      %   and flipped output.
      %
      % For details on usage/arguments see :meth:`forcetorque`.

      if isa(other, 'ott.scat.utils.Particle')
        [varargout{1:nargout}] = other.force(beam, varargin{:});
      elseif isa(other, 'ott.beam.abstract.Beam')
        [varargout{1:nargout}] = beam.callParticleMethod(...
            @beam.forceInternal, other, varargin{:});
      end
    end

    function varargout = torque(beam, other, varargin)
      % Calculate change in angular momentum between beams.
      %
      % Usage
      %   torque = beam.torque(other_beam, ...)
      %
      %   torque = beam.torque(particle, ...)
      %   Use the particle's torque method with flipped position/rotation
      %   and flipped output.
      %
      % For details on usage/arguments see :meth:`forcetorque`.

      if isa(other, 'ott.scat.utils.Particle')
        [varargout{1:nargout}] = other.torque(beam, varargin{:});
      elseif isa(other, 'ott.beam.abstract.Beam')
        [varargout{1:nargout}] = beam.callParticleMethod(...
            @beam.torqueInternal, other, varargin{:});
      end
    end

    function varargout = forcetorque(beam, other, varargin)
      % Calculate change in momentum between beams.
      %
      % Usage
      %   [force, torque] = beam.forcetorque(other_beam, ...)
      %   Returns 3xN matrices for the force and torque in Cartesian
      %   coordinates.
      %
      %   [force, torque] = beam.forcetorque(particle, ...)
      %   Use the particle's `forcetorque` method.
      %   Use the particle's `forcetorque` method with flipped
      %   position/rotation and flipped output.
      %
      % Parameters
      %   - other (Beam|scat.utils.Particle) -- A beam to compare the force
      %     with or a particle with a forcetorque method.
      %
      %   - position (3xN numeric) -- Distance to translate beam before
      %     calculating the scattered beam using the T-matrix.
      %     Default: ``[]``.
      %
      %   - rotation (3x3N numeric) -- Angle to rotate beam before
      %     calculating the scattered beam using the T-matrix.
      %     Inverse rotation is applied to scattered beam, effectively
      %     rotating the particle.
      %     Default: ``[]``.

      if isa(other, 'ott.scat.utils.Particle')
        [varargout{1:nargout}] = other.forcetorque(beam, varargin{:});
      elseif isa(other, 'ott.beam.abstract.Beam')
        [varargout{1:nargout}] = beam.callParticleMethod(...
            @beam.forcetorqueInternal, other, varargin{:});
      end
    end

    function varargout = scatter(beam, particle, varargin)
      % Calculate the beam scattered by a particle
      %
      % Usage
      %   sbeam = ibeam.scatter(particle, ...)
      %
      % Uses the particle's scatter method.  For details and function
      % arguments, see the particle scatter method's documentation.

      [position, rotation, unmatched] = beam.flipPositionRotation(varargin{:});

      [varargout{1:nargout}] = particle.scatter(beam, ...
          'position', position, 'rotation', rotation, unmatched{:});
    end
  end

  methods (Hidden)
    function varargout = callParticleMethod(beam, method, other, varargin)
      % Helper for calling particle methods with position/rotation inputs
      %
      % Forces should be flipped when calculated with the beam as the
      % primary (force on the beam vs force on the particle).
      %
      % Position and rotation should be inverted.

      % Flip position/rotation
      [position, rotation, unmatched] = beam.flipPositionRotation(varargin{:});

      % Call particle method
      [varargout{1:nargout}] = ...
        ott.utils.RotationPositionProp.rotationPositionHelper(method, ...
        other, 'position', position, 'rotation', rotation, ...
        'prxcat', 2, unmatched{:});

      % Flip outputs
      for ii = 1:nargout
        varargout{ii} = -varargout{ii};
      end
    end

    function [position, rotation, unmatched] = flipPositionRotation(...
        beam, varargin)
      % Flip position and rotation variable arguments

      p = inputParser;
      p.addParameter('position', []);
      p.addParameter('rotation', []);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Flip position/rotation
      position = -p.Results.position;
      rotation = p.Results.rotation;
      for ii = 1:size(rotation, 2)./3
        rotation(:, (1:3) + (ii-1)*3) = rotation(:, (1:3) + (ii-1)*3).';
      end
    end

    function [xyz, unmatched] = transformWorldToLocal(beam, xyz, varargin)
      % Helper for world to local transforms with arguments
      %
      % Usage
      %   xyz = transformWorldToLocal(beam, xyz, ...)
      %   Parses inputs for origin parameter.  Error on unmatched.
      %
      %   [xyz, unmatched] = transformWorldToLocal(beam, xyz, ...)
      %   Keeps unmatched arguments

      p = inputParser;
      p.KeepUnmatched = nargout == 2;
      p.addParameter('origin', 'world');
      p.parse(varargin{:});

      if nargout == 2
        unmatched = ott.utils.unmatchedArgs(p);
      end

      % Translate locations to beam coordinates
      if strcmpi(p.Results.origin, 'world')
        xyz = beam.global2local(xyz);
      end
    end

    function div = numericalDivergence(beam, xyz, field, varargin)
      % Estimate divergence of field numerically (forward difference)
      %
      % Usage
      %   gray = numericalGradient(beam, xyz, field, ...)
      %
      % Parameters
      %   - field (function_handle) -- Function returning a FieldVector field.

      % TODO: Array support

      p = inputParser;
      p.addParameter('dx', 1.0e-3*beam.wavelength);
      p.parse(varargin{:});

      dx = p.Results.dx;

      E0 = field(xyz);

      Ex = field(xyz + [dx;0;0]);
      Ey = field(xyz + [0;dx;0]);
      Ez = field(xyz + [0;0;dx]);
      Eforward = [Ex.vxyz(1, :); Ey.vxyz(2, :); Ez.vxyz(3, :)];

      div = (Eforward - E0.vxyz) ./ dx;

    end

    function grad = numericalGradient(beam, xyz, field, varargin)
      % Calculate gradient numerically
      %
      % Usage
      %   gray = numericalGradient(beam, xyz, field, ...)
      %
      % Parameters
      %   - field (function_handle) -- Function returning a scalar or
      %     `ott.utils.FieldVector` field.  If the type is scalar,
      %     returns a 3xN matrix, otherwise returns a 3x3xN matrix.

      % TODO: Array support

      p = inputParser;
      p.addParameter('dx', 1.0e-3*beam.wavelength);
      p.parse(varargin{:});

      dx = p.Results.dx;

      E0 = field(xyz);

      Ex = field(xyz + [dx;0;0]);
      Ey = field(xyz + [0;dx;0]);
      Ez = field(xyz + [0;0;dx]);

      if isa(E0, 'ott.utils.FieldVector')

        % Vector gradient

        E0 = reshape(E0.vxyz, 3, 1, []);
        Ex = reshape(Ex.vxyz, 3, 1, []);
        Ey = reshape(Ey.vxyz, 3, 1, []);
        Ez = reshape(Ez.vxyz, 3, 1, []);

        Eforward = [Ex, Ey, Ez];
        grad = (Eforward - E0) ./ dx;

      else

        % Scalar gradient
        Eforward = [Ex; Ey; Ez];
        grad = (Eforward - E0) ./ dx;

      end

    end

    function gradE = egradientInternal(beam, xyz, varargin)
      % Default implementation estimates gradient numerically
      % If your beam has a nicer method, overload this function.

      gradE = beam.numericalGradient(xyz, @beam.efield, varargin{:});
    end

    function gradH = hgradientInternal(beam, xyz, varargin)
      % Default implementation estimates gradient numerically
      % If your beam has a nicer method, overload this function.

      gradH = beam.numericalGradient(xyz, @beam.hfield, varargin{:});
    end

    function [gradE, gradH] = ehgradientInternal(beam, xyz, varargin)
      % Default implementation estimates gradient numerically
      % If your beam has a nicer method, overload this function.

      gradE = beam.egradientInternal(xyz, varargin{:});
      gradH = beam.hgradientInternal(xyz, varargin{:});
    end

    function [E, H] = ehfieldInternal(beam, xyz, varargin)
      % Default implementation just calls each part, but for some
      % beams it may be more efficient to calculate these together,
      % in which case this method should be over-written.

      E = beam.efieldInternal(xyz, varargin{:});
      H = beam.hfieldInternal(xyz, varargin{:});
    end

    function data = calculateVisualisationData(beam, field_type, data)
      % Apply the VisualisationData function to a data array
      %
      % This function is overloaded by Array types in order to
      % implement incoherent combination.
      
      data = beam.arrayApply(...
        @(x) beam.VisualisationData(field_type, x), data);
      
      % Combine incoherent components
      if isa(beam, 'ott.beam.properties.ArrayType')
        data = beam.combineIncoherentArray(data, 3);
      end
    end

    function imout = intensityMomentImout(beam, ~, rtp, uxyz, ...
        theta, ~, dtheta, dphi, varargin)
      % Calculate the imout data for the intensityMoment function

      % Calculate fields
      % TODO: Should this use poyntingFarfield?
      E = beam.efarfield(rtp, varargin{:});
      Eirr = beam.calculateVisualisationData('E2', E);

      ints = sum(Eirr .* sin(theta.') .* dtheta .* dphi, 2);

      % Calculate moment in Cartesian coordinates
      Eirr_xyz = uxyz .* Eirr;
      moments = sum(Eirr_xyz .* sin(theta.') .* dtheta .* dphi, 2);

      imout = [moments(:); ints];
    end
    
    function arr = applyMask(~, arr, mask)
      arr(~mask) = nan;
    end

    function imout = visualiseImoutXyzMasked(beam, sz, xyz, field_type, mask, varargin)
      % Calculate the imout data for the visualise function

      % Calculate field and visualisation data
      E = beam.efield(xyz, varargin{:});
      imout = beam.calculateVisualisationData(field_type, E);
      
      % Reshape output to match required size
      imout = beam.arrayApply(@(x) reshape(x, sz), imout);

      % Create layer from masked data
      if ~all(mask(:))
        imout = beam.arrayApply(@(x) beam.applyMask(x, mask), imout);
      end
    end

    function visdata = visualiseImoutRtp(beam, sz, rtp, field_type, varargin)
      % Calculate the imout data for farfield visualisation functions

      E = beam.efarfield(rtp, varargin{:});
      visdata = beam.calculateVisualisationData(field_type, E);
      
      % Reshape output to match required size
      visdata = beam.arrayApply(@(x) reshape(x, [sz, 1]), visdata);
    end

    function force = forceInternal(beam, other, varargin)
      % Calculate the force using the intensityMoment method
      %
      % This function should be overridden in sub-classes which
      % support other (possibly more efficient) force calculation methods.

      fbeam = beam.intensityMoment();
      obeam = other.intensityMoment();

      force = fbeam - obeam;
    end

    function torque = torqueInternal(beam, other, varargin)
      % Calculate the torque
      %
      % This function should be overridden in sub-classes which
      % support other (possibly more efficient) torque calculation methods.

      % TODO: Implement this
      warning('Not yet implemented');
      torque = [0;0;0];
    end

    function [force, torque] = forcetorqueInternal(beam, other, varargin)
      % Calculate the force and the torque.
      %
      % This function uses the forceInternal and torqueInternal methods.
      %
      % This function should be overridden in sub-classes which
      % support other (possibly more efficient) torque calculation methods.

      force = beam.forceInternal(other, varargin{:});
      torque = beam.torqueInternal(other, varargin{:});
    end
  end

  methods (Static)

    function [xrange, yrange, sz] = visualiseGetRange(range, sz, default_sz)
      % Get the xrange and yrange dimensions of the plot
      %
      % Usage
      %   [xrange, yrange, sz] = visualiseGetRange(range, sz, default_sz)

      if iscell(range)

        if ~isempty(sz)
          warning('Ignoring size parameter');
        end

        xrange = range{1};
        yrange = range{2};
        sz = [length(yrange), length(xrange)];

      elseif length(range) == 2

        if isempty(sz)
          sz = default_sz;
        end
        assert(isnumeric(sz) && numel(sz) == 2, ...
            'size must be 2 element numeric vector');

        xrange = linspace(-1, 1, sz(2))*range(1);
        yrange = linspace(-1, 1, sz(1))*range(2);

      elseif length(range) == 4

        if isempty(sz)
          sz = default_sz;
        end
        assert(isnumeric(sz) && numel(sz) == 2, ...
            'size must be 2 element numeric vector');

        xrange = linspace(range(1), range(2), sz(2));
        yrange = linspace(range(3), range(4), sz(1));

      else
        error('ott:Bsc:visualise:range_error', ...
            'Incorrect type or size of range arguments');
      end
    end

    function [xyz, labels] = visualiseGetXyz(xx, yy, zz, our_axis)
      % Helper function for :meth:`visualise`
      %
      % Usage
      %   [xyz, labels] = beam.visualiseGetXyz(xx, yy, zz, axis)
      %
      %   [~, labels] = beam.visualiseGetXyz([], [], [], axis)

      % TODO: Should we split this function?

      % Change the grid to the user requested orientation
      if ischar(our_axis)
        switch our_axis
          case 'x'
            xyz = [zz(:), yy(:), xx(:)].';
            labels = {'Z', 'Y'};
          case 'y'
            xyz = [yy(:), zz(:), xx(:)].';
            labels = {'Z', 'X'};
          case 'z'
            xyz = [xx(:), yy(:), zz(:)].';
            labels = {'X', 'Y'};
          otherwise
            error('Unknown axis name specified');
        end
      elseif iscell(our_axis) && numel(our_axis) >= 2
        dir1 = our_axis{1}(:);
        dir2 = our_axis{2}(:);
        if numel(our_axis) == 3
          dir3 = our_axis{3}(:);
        else
          dir3 = cross(dir1(:), dir2(:));
        end

        assert(isnumeric(dir1) && numel(dir1) == 3 ...
            && isnumeric(dir2) && numel(dir2) == 3 ...
            && isnumeric(dir3) && numel(dir3) == 3, ...
            'direction vectors must be 3 element numeric');

        labels = {'Direction 1', 'Direction 2'};

        if ~isempty(xx) && ~isempty(yy) && ~isempty(zz)
          xyz = dir1.*xx(:).' + dir2.*yy(:).' + dir3.*zz(:).';
        else
          xyz = [];
        end
      else
        error('axis must be character or a 2 or 3 element cell array');
      end
    end

    function visualiseShowPlot(vis_nargout, our_axes, imout, ...
        ranges, labels)
      % Helper function for generating the visualisation
      % Called by :meth:`visualise`.
      %
      % Usage
      %   beam.visualiseShowPlot(vis_nargout, our_axes, imout, ranges, labels)

      if vis_nargout == 0 || ~isempty(our_axes)

        % Get the default axes
        if isempty(our_axes)
          our_axes = gca();
        end

        % Check that we only have one visualisation to show
        if ~ismatrix(imout)
          warning('ott:beam:Beam:non_matrix_plot', ...
              'Multiple beams generated, only showing first');
        end

        % Generate the plot
        imagesc(our_axes, ranges{1}, ranges{2}, imout(:, :, 1), ...
            'AlphaData', ~isnan(imout(:, :, 1)));
        xlabel(our_axes, labels{1});
        ylabel(our_axes, labels{2});
        axis(our_axes, 'image');
      end
    end
  end
end

