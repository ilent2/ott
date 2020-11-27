classdef Beam < matlab.mixin.Heterogeneous ...
    & ott.utils.RotationPositionProp
% Provides a high-level view of a optical tweezers toolbox beam.
% Inherits from :class:`ott.utils.RotationPositionProp` and
% `matlab.mixin.Hetrogeneous`.
%
% This class is the base class for OTT beams.  :class:`Beam` and its
% sub-classes provide a description of beams but do not implement any
% of the field or scattering calculation methods.
% These classes separate properties describing the beam (such as position,
% power and waist) from properties specific to the numerical calculation
% (VSWF data, Nmax, apparent power).
% Internally, :class:`Beam` uses :class:`ott.bsc.Bsc` for field calculations.
% Depending on the implementation, the :class:`ott.bsc.Bsc` data is
% either stored or calculated when required.
%
% The other major difference from :class:`ott.bsc.Bsc` is the units.
% :class:`Beam` uses SI units for all quantities, making integration
% with dynamics simulations easer.
%
% In OTTv2, this interface will likely be extended to other types of
% beam/scattering methods (such as paraxial beams or other approximations).
%
% Properties
%   - position        -- (3x1 numeric) Location of the beam
%   - rotation        -- (3x3 numeric) Orientation of the beam
%   - wavelength      -- Wavelength in medium [m]
%   - wavenumber      -- Wavenumber in medium [1/m]
%   - speed           -- Speed of light in the medium [m/s]
%   - speed0          -- Speed of light in vacuum [m/s]
%   - scale           -- Scaling parameter for beam fields (default: 1)
%
% Abstract properties
%   - omega           -- Optical angular frequency of light [1/s]
%   - index_medium    -- Refractive index of the medium
%
% Methods
%   - setWavelength -- Set wavelength property
%   - setWavenumber -- Set wavenumber property
%   - scatter       -- Calculate how a beam is scattered
%
% Field calculation methods
%   - ehfield    -- Calculate electric and magnetic fields around the origin
%   - ehfieldRtp -- Calculate electric and magnetic fields around the origin
%   - ehfarfield -- Calculate electric and magnetic fields in the far-field
%   - eparaxial  -- Calculate electric fields in the paraxial far-field
%   - hparaxial  -- Calculate magnetic fields in the paraxial far-field
%   - ehparaxial -- Calculate electric and magnetic paraxial far-fields
%
% Force and torque related methods
%   - forcetorque     -- Calculate the force and the torque between beams
%   - intensityMoment -- Calculate moment of beam intensity in far-field
%   - force      -- Calculate the change in momentum between two beams
%   - torque     -- Calculate change in angular momentum between beams
%   - spin       -- Calculate change in spin momentum between beams
%
% Field visualisation methods
%   - visNearfield      -- Generate a visualisation around the origin
%   - visFarfield       -- Generate a visualisation at the far-field
%   - visFarfieldSlice  -- Visualise the field on a angular slice
%   - visFarfieldSphere -- Visualise the filed on a sphere
%
% Mathematical operations
%   - times,mtimes     -- Scalar multiplication of beam fields
%   - rdivide,mrdivide -- Scalar division of beam fields
%   - uminus           -- Flip beam field intensity
%   - plus, minus      -- Combine beams coherently
%   - or               -- Combine beams incoherently
%
% Abstract methods
%   - efield     -- Calculate electric field around the origin
%   - hfield     -- Calculate magnetic field around the origin
%   - efieldRtp  -- Calculate electric field around the origin (sph. coords.)
%   - hfieldRtp  -- Calculate magnetic field around the origin (sph. coords.)
%   - efarfield  -- Calculate electric fields in the far-field
%   - hfarfield  -- Calculate magnetic fields in the far-field
%   - scatterInternal -- Called to calculate the scattered beam

% Copyright 2018-2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    scale          % Scaling parameter for beam fields
  end

  properties (Abstract)
    index_medium   % Refractive index of the medium
    omega          % Optical angular frequency of light [1/s]
  end

  properties (Dependent)
    speed          % Speed of light in the medium [m/s]
    wavelength     % Wavelength in medium [m]
    wavelength0    % Wavelength in vacuum [m]
    wavenumber     % Wavenumber in medium [1/m]
    impedance      % Impedance of the medium
  end

  properties (Constant)
    speed0 = 3e8   % Speed of light in vacuum [m/s]
  end

  methods (Abstract)
    efield(~)      % Calculate electric field around the origin
    hfield(~)      % Calculate magnetic field around the origin
    efieldRtp(~)   % Calculate electric field around the origin (sph. coords.)
    hfieldRtp(~)   % Calculate magnetic field around the origin (sph. coords.)
    efarfield(~)   % Calculate electric fields in the far-field
    hfarfield(~)   % Calculate magnetic fields in the far-field
    scatterInternal(~)     % Called to calculate the scattered beam
  end

  methods (Static)
    function data = VisualisationData(field_type, field)
      % Helper to generate the visualisation data output.
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
      %     - irradiance  -- :math:`\sqrt{|Ex|^2 + |Ey|^2 + |Ez|^2}`
      %     - E2          -- :math:`|Ex|^2 + |Ey|^2 + |Ez|^2`
      %     - Sum(Abs(E)) -- :math:`|Ex| + |Ey| + |Ez|`
      %
      %     - Re(Er), Re(Et), Re(Ep), Re(Ex), Re(Ey), Re(Ez)
      %     - Abs(Er), Abs(Et), Abs(Ep), Abs(Ex), Abs(Ey), Abs(Ez)
      %     - Arg(Er), Arg(Et), Arg(Ep), Arg(Ex), Arg(Ey), Arg(Ez)

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
          error('ott:bsc:Bsc:GetVisualisationData:unknown_field_type', ...
              'Unknown value for field_type.');
      end
    end
  end

  methods (Access=protected)
    function beam = Beam(varargin)
      % Construct a new Beam instance (should only be called by sub-classes)
      %
      % Usage
      %   beam = Beam(...)
      %
      % Optional named arguments
      %   - scale (numeric) -- Initial scale.  Default: ``1``.
      %
      %   - position (3x1 numeric) -- Initial position [m].
      %     Default: ``[0;0;0]``.
      %
      %   - rotation (3x3 numeric) -- Initial rotation matrix.
      %     Default: ``eye(3)``.

      p = inputParser;
      p.addParameter('scale', 1.0);
      p.addParameter('position', [0;0;0]);
      p.addParameter('rotation', eye(3));
      p.parse(varargin{:});

      beam.position = p.Results.position;
      beam.rotation = p.Results.rotation;
      beam.scale = p.Results.scale;
    end
  end

  methods

    function beam = setWavelength(beam, val, mode)
      % Set the beam wavelength property
      %
      % Usage
      %   beam = beam.setWavelength(value, mode)
      %
      % Parameters
      %   - value (numeric) -- Wavelength values [meters]
      %   - mode (enum) -- Either 'fixedFrequency' or 'fixedSpeed'
      %     to set the wavelength by chaining ``speed`` or ``omega``.

      ott.utils.nargoutCheck(beam, nargout);

      if strcmpi(mode, 'fixedSpeed')
        beam.omega = 2*pi*beam.speed0./val;
      elseif strcmpi(mode, 'fixedFrequency')
        beam.speed = beam.omega./(2*pi).*val;
      else
        error('Unknown mode value');
      end
    end

    function beam = setWavenumber(beam, val, mode)
      % Set the beam wavenumber property
      %
      % Usage
      %   beam = beam.setWavenumber(value, mode)
      %
      % Parameters
      %   - value (numeric) -- Wavenumber values [1/meters]
      %   - mode (enum) -- Either 'fixedFrequency' or 'fixedSpeed'
      %     to set the wavelength by chaining ``speed`` or ``omega``.

      ott.utils.nargoutCheck(beam, nargout);
      beam = beam.setWavelength(2*pi./val, mode);
    end

    %
    % Field calculation functions
    %

    function varargout = eparaxial(beam, xy, varargin)
      % Calculate E paraxial far-field
      %
      % Usage
      %   [E, data] = beam.eparaxial(xy, ...)
      %   E is of type :class:`ott.utils.FieldVector`.
      %
      % Parameters
      %   - xy (2xN numeric) -- Paraxial coordinates for field
      %     calculation. Packaged [x; y].
      %
      % Optional named parameters
      %   - mapping (enum) -- Mapping for paraxial projection.
      %     See :func:`ott.utils.paraxial2rtp` for details.
      %     Default: ``'sin'``.
      %
      %   - direction (enum | 2 numeric) -- Mapping direction.
      %     See :func:`ott.utils.paraxial2rtp` for details.
      %     Default: ``'pos'``.
      %
      % Additional parameters passed to :meth:`efarfield`.

      p = inputParser;
      p.addParameter('mapping', 'sin');
      p.addParameter('direction', 'pos');
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      rtp = ott.utils.paraxial2rtp(xy, ...
          p.Results.mapping, p.Results.direction);
      [varargout{1:nargout}] = beam.efarfield(rtp, unmatched{:});
    end

    function varargout = hparaxial(beam, xy, varargin)
      % Calculate H paraxial far-field
      %
      % Usage
      %   [H, data] = beam.hparaxial(xy, ...)
      %   H is of type :class:`ott.utils.FieldVector`.
      %
      % Parameters
      %   - xy (2xN numeric) -- Paraxial coordinates for field
      %     calculation. Packaged [x; y].
      %
      % Optional named parameters
      %   - mapping (enum) -- Mapping for paraxial projection.
      %     See :func:`ott.utils.paraxial2rtp` for details.
      %     Default: ``'sin'``.
      %
      %   - direction (enum | 2 numeric) -- Mapping direction.
      %     See :func:`ott.utils.paraxial2rtp` for details.
      %     Default: ``'pos'``.
      %
      % Additional parameters passed to :meth:`hfarfield`.

      p = inputParser;
      p.addParameter('mapping', 'sin');
      p.addParameter('direction', 'pos');
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      rtp = ott.utils.paraxial2rtp(xy, ...
          p.Results.mapping, p.Results.direction);
      [varargout{1:nargout}] = beam.hfarfield(rtp, unmatched{:});
    end

    function varargout = ehparaxial(beam, xy, varargin)
      % Calculate E and H paraxial far-fields
      %
      % Usage
      %   [E, H, data] = beam.ehparaxial(xy, ...)
      %   E and H are of type :class:`ott.utils.FieldVector`.
      %
      % Parameters
      %   - xy (2xN numeric) -- Paraxial coordinates for field
      %     calculation. Packaged [x; y].
      %
      % Optional named parameters
      %   - mapping (enum) -- Mapping for paraxial projection.
      %     See :func:`ott.utils.paraxial2rtp` for details.
      %     Default: ``'sin'``.
      %
      %   - direction (enum | 2 numeric) -- Mapping direction.
      %     See :func:`ott.utils.paraxial2rtp` for details.
      %     Default: ``'pos'``.
      %
      % Additional parameters passed to :meth:`ehfarfield`.

      p = inputParser;
      p.addParameter('mapping', 'sin');
      p.addParameter('direction', 'pos');
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      rtp = ott.utils.paraxial2rtp(xy, ...
          p.Results.mapping, p.Results.direction);
      [varargout{1:nargout}] = beam.ehfarfield(rtp, unmatched{:});
    end

    function [E, H, vswfData] = ehfarfield(beam, rtp, varargin)
      % Calculate E and H fields in Cartesian coordinates. (SI units)
      %
      % Usage
      %   [E, H, data] = beam.ehfarfield(rtp, ...)
      %   E and H are of type :class:`ott.utils.FieldVector`.
      %
      % Parameters
      %   - rtp (3xN numeric) -- Spherical coordinates for field calculation.
      %     Units of meters/radians.  Packaged [r; theta; phi]
      %
      % Optional named parameters
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.

      [E, vswfData] = beam.efarfield(rtp, varargin{:});
      [H, vswfData] = beam.hfarfield(rtp, 'data', vswfData);
    end

    function [E, H, vswfData] = ehfield(beam, xyz, varargin)
      % Calculate E and H fields in Cartesian coordinates. (SI units)
      %
      % Usage
      %   [E, H, data] = beam.ehfield(xyz, ...)
      %   E and H are of type :class:`ott.utils.FieldVector`.
      %
      % Parameters
      %   - rtp (3xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of meters.  Packaged [x; y; z]
      %
      % Optional named parameters
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.

      [E, vswfData] = beam.efield(xyz, varargin{:});
      [H, vswfData] = beam.hfield(xyz, 'data', vswfData);
    end

    function [E, H, vswfData] = ehfieldRtp(beam, rtp, varargin)
      % Calculate E and H fields in spherical coordinates. (SI units)
      %
      % Usage
      %   [E, H, data] = beam.ehfieldRtp(rtp, ...)
      %   E and H are of type :class:`ott.utils.FieldVector`.
      %
      % Parameters
      %   - rtp (3xN numeric) -- Spherical coordinates for field calculation.
      %     Units of meters/radians.  Packaged [r; theta; phi].
      %
      % Optional named parameters
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.

      [E, vswfData] = beam.efieldRtp(rtp, varargin{:});
      [H, vswfData] = beam.hfieldRtp(rtp, 'data', vswfData);
    end
    
    function [Efv, psi] = translateFarfields(beam, Efv, direction)
      % Applies a phase shift to the far-fields similar to a translation
      %
      % Usage
      %   [Efv, psi] = beam.translateFarfields(Efv, direction)
      %
      % Returns
      %   - Efv (ott.utils.FieldVector) -- Scaled fields.
      %   - psi (N numeric) -- Phase factors.
      %
      % Parameters
      %   - Efv (ott.utils.FieldVector) -- Far-fields to shift.
      %
      %   - direction (3x1 numeric) -- Direction and distance to shift.
      %     Cartesian coordinates.
      
      if all(direction == [0;0;0])
        % No work to do
        sz = size(Efv);
        psi = ones(1, sz(2:end));
      else
        % Convert to cartesian field vector
        Efv = ott.utils.FieldVectorCart(Efv);

        % Get coordinates and normlise
        uxyz = double(Efv(4:6, :));
        uxyz = uxyz ./ vecnorm(uxyz);
        
        psi = sum(uxyz .* direction, 1) * beam.wavenumber;
        Efv = Efv .* exp(1i*psi);
      end
    end

    %
    % Force and torque related methods
    %

    function f = force(beam, other, varargin)
      % Calculate the difference in momentum between beams.
      %
      % The default implementation uses :meth:`intensityMoment`.
      % This method should be overloaded if there is a more optimal
      % method for your time of beam.
      %
      % Usage
      %   force = beam.force(other_beam)
      %
      %   force = beam.force(particle, ...)
      %
      % Returns
      %   - 3x1 numeric vector.
      
      assert(isa(beam, 'ott.beam.Beam'), ...
        'first argument must be Beam instance');
      
      if ~isa(other, 'ott.beam.Beam')
        % TODO: This is a bit of duplication between BscBeam and in the
        % future (once implemented) Array.  Can we improve it?
        
        p = inputParser;
        p.addParameter('position', other.position);
        p.addParameter('rotation', other.rotation);
        p.parse(varargin{:});
        
        rotation = p.Results.rotation;
        position = p.Results.position;
        
        Nrot = size(rotation, 2)/3;
        Npos = size(position, 2);
        assert(Nrot == Npos || Nrot == 1 || Npos == 1, ...
          'number of rotations/positions must mach or be singular');
        if Npos == 1, position = repmat(position, 1, Nrot); end
        if Nrot == 1, rotation = repmat(rotation, 1, Npos); end
        
        f = zeros(3, size(position, 2));
        for ii = 1:size(position, 2)
          sbeam = beam.scatter(other, 'position', position(:, ii), ...
            'rotation', rotation(:, (1:3) + (ii-1)*3));
          f(:, ii) = sbeam.force();
        end
        return;
      end

      fbeam = beam.intensityMoment() ./ beam.speed;
      obeam = other.intensityMoment() ./ other.speed;

      f = obeam - fbeam;
    end

    function t = torque(~, varargin)
      % Calculate the angular momentum difference between beams.
      %
      % There is no default implementation of this method.  This function
      % simply returns nans.

      t = nan(3, 1);
    end

    function t = spin(~, varargin)
      % Calculate the spin angular momentum between beams.
      %
      % There is no default implementation of this method.  This function
      % simply returns nans.

      t = nan(3, 1);
    end

    function varargout = forcetorque(ibeam, sbeam, varargin)
      % Calculate force, torque and spin (in SI units)
      %
      % Usage
      %   [f, t, s] = incident_beam.forcetorque(scattered_beam) --
      %   Calculates the force/torque/spin between two beams.
      %
      %   [f, t, s] = incident_beam.forcetorque(particle, ...) --
      %   First calculates the scattering, then calculates the
      %   force/torque/spin between resulting beams.
      %   Additional parameters are passed to :meth:`scatter`.

      % Calculate scattering if required
      if ~isa(sbeam, 'ott.beam.Beam')
        sbeam = ibeam.scatter(sbeam, varargin{:});
        [varargout{1:nargout}] = sbeam.forcetorque();
        return;
      end

      varargout{1} = ibeam.force(sbeam);
      if nargout > 1
        varargout{2} = ibeam.torque(sbeam);
        if nargout > 2
          varargout{3} = ibeam.spin(sbeam);
        end
      end
    end

    function varargout = intensityMoment(beam, varargin)
      % Calculate moments of the beam intensity in the far-field.
      %
      % Calculates the moment of far-field energy flux::
      %
      %   S = |E|^2
      %
      % Using this quantities, it should be possible to get a reasonable
      % estimate for the force.
      %
      % Usage
      %   [S, I, ...] = beam.intensityMoment(...)
      %
      % Returns
      %   - S (3x1 numeric) -- Moment of energy density.
      %   - I (1 numeric) -- Total outgoing energy flux [W].
      %
      % Optional named arguments
      %   - theta_range (2 numeric) -- Range of angles to integrate over.
      %     Default: ``[0, pi]``.
      %
      %   - ntheta (numeric) -- Number of theta points.  (Default: 100)
      %
      %   - nphi (numeric) -- Number of phi points.  (Default: 100)
      %
      % Additional parameters passed to far-field calculation functions.
      % Additional outputs from far-field calculation functions are
      % returned after first three arguments.

      p = inputParser;
      p.addParameter('theta_range', [0, pi]);
      p.addParameter('ntheta', 100);
      p.addParameter('nphi', 100);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Setup grid
      [theta, phi] = ott.utils.angulargrid(p.Results.ntheta, p.Results.nphi);
      dtheta = diff(theta(1:2));
      dphi = diff(phi([1, p.Results.ntheta+1]));

      % Truncate the theta range
      keep = theta > p.Results.theta_range(1) ...
          & theta < p.Results.theta_range(2);
      theta = theta(keep);
      phi = phi(keep);

      rtp = [ones(numel(theta), 1), theta(:), phi(:)].';

      % Calculate Cartesian coordinate unit vectors
      uxyz = ott.utils.rtp2xyz(rtp);

      % Calculate field and E2
      [E, varargout{3:nargout}] = beam.efarfield(rtp, ...
          'basis', 'outgoing', unmatched{:});
      Eirr = beam.VisualisationData('E2', E);

      % Make sure units are correct (S = E x H = 1/Z*E^2)
      Eirr = Eirr ./ beam.impedance;

      % Calculate intensity
      varargout{2} = sum(Eirr .* sin(theta.') .* dtheta .* dphi, 2);

      % Calculate moment in Cartesian coordinates
      Eirr_xyz = uxyz .* Eirr;
      varargout{1} = sum(Eirr_xyz .* sin(theta.') .* dtheta .* dphi, 2);
    end

    %
    % Visualisation functions
    %

    function varargout = visNearfield(beam, varargin)
      % Create a visualisation of the beam
      %
      % Usage
      %   beam.visNearfield(...) -- display an image of the beam in
      %   the current figure window.
      %
      %   [im, XY, data] = beam.visualise(...) -- returns a image of the
      %   beam and a cell array of the coordinates used in the calculation.
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

      p = inputParser;
      p.addParameter('field', 'irradiance');
      p.addParameter('size', []);
      p.addParameter('axis', 'z');
      p.addParameter('offset', 0.0);
      p.addParameter('range', beam.defaultVisRange);
      p.addParameter('mask', []);
      p.addParameter('plot_axes', []);   % Not used (just filtered)
      p.addParameter('data', ott.utils.VswfData(), ...
          @(x) isa(x, 'ott.utils.VswfData'));
      p.KeepUnmatched = true;
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
      [xyz, ~] = beam.visualiseGetXyz(xx, yy, zz, p.Results.axis);

      % Calculate mask
      if isempty(p.Results.mask)
        mask = true(sz);
      else
        mask = p.Results.mask(xyz);
      end

      % Calculate which points are needed
      our_xyz = xyz(:, mask(:));

      % Calculate fields
      [E, data] = beam.efield(our_xyz, 'data', p.Results.data, ...
          unmatched{:});

      % Unpack masked image (use nans for other values)
      Ef = nan([size(xyz), size(E, 3)]);
      Ef(:, mask(:), :) = E.vxyz;
      xyz = repmat(xyz, [1, 1, size(Ef, 3)]);
      Ef = ott.utils.FieldVectorCart(Ef, xyz);

      % Generate visualisation data
      imout = beam.VisualisationData(p.Results.field, Ef);
      imout = reshape(imout, [sz, size(E, 3)]);

      % Display the visualisation
      beam.visShowPlot(imout, {xrange, yrange}, nargout, varargin{:});

      % Assign output
      if nargout >= 1
        varargout{1} = imout;
        if nargout >= 2
          varargout{2} = {xrange, yrange, xyz};
          if nargout >= 3
            varargout{3} = data;
          end
        end
      end
    end

    function varargout = visFarfield(beam, varargin)
      % Create a visualisation of the beam by projecting the far-field
      % onto a plane.
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

      % Parse arguments
      p = inputParser;
      p.addParameter('field', 'irradiance');
      p.addParameter('size', []);
      p.addParameter('direction', 'pos');
      p.addParameter('range', [1, 1]);
      p.addParameter('mapping', 'sin');
      p.addParameter('plot_axes', []);
      p.addParameter('data', ott.utils.VswfData(), ...
          @(x) isa(x, 'ott.utils.VswfData'));
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Generate grid of coordinates
      default_sz = [80, 80];
      [xrange, yrange, sz] = beam.visualiseGetRange(p.Results.range, ...
          p.Results.size, default_sz);
      [xx, yy] = meshgrid(xrange, yrange);
      nxy = [xx(:), yy(:)].';

      % Calculate fields
      [E, data] = beam.eparaxial(nxy, 'data', p.Results.data, ...
          unmatched{:}, 'mapping', p.Results.mapping, ...
          'direction', p.Results.direction);

      % Generate visualisation data
      imout = beam.VisualisationData(p.Results.field, E);
      imout = reshape(imout, [sz, size(E, 3)]);

      % Display the visualisation
      beam.visShowPlot(imout, {xrange, yrange}, nargout, varargin{:});

      % Assign output
      if nargout >= 1
        varargout{1} = imout;
        if nargout >= 2
          varargout{2} = {xrange, yrange};
          if nargout >= 3
            varargout{3} = data;
          end
        end
      end
    end

    function varargout = visFarfieldSphere(beam, varargin)
      % Generate a spherical surface visualisation of the far-field
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

      p = inputParser;
      p.addParameter('field', 'irradiance');
      p.addParameter('normalise', false);
      p.addParameter('type', 'sphere');
      p.addParameter('npts', 100);
      p.addParameter('plot_axes', []);
      p.addParameter('data', ott.utils.VswfData(), ...
          @(x) isa(x, 'ott.utils.VswfData'));
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Generate grid
      [X,Y,Z] = sphere(p.Results.npts);

      % Calculate fields
      rtp = ott.utils.xyz2rtp(X, Y, Z);
      [E, data] = beam.efarfield(rtp, 'data', p.Results.data, ...
          unmatched{:});

      % Generate visualisation data
      imout = beam.VisualisationData(p.Results.field, E);
      imout = reshape(imout, size(X));

      % Generate visualisation
      beam.visFarfieldSpherePlot(imout, {X, Y, Z}, nargout, varargin{:});

      % Assign output data
      if nargout >= 1
        varargout{1} = imout;
        if nargout >= 2
          varargout{2} = {X, Y, Z};
          if nargout >= 3
            varargout{3} = data;
          end
        end
      end
    end

    function [imout, ptheta, vswfData] = visFarfieldSlice(beam, varargin)
      % Generate a 2-D slice through the far-field
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

      p = inputParser;
      p.addOptional('phi', [0, pi/2], @isnumeric);
      p.addParameter('field', 'irradiance');
      p.addParameter('basis', 'incoming');
      p.addParameter('normalise', false);
      p.addParameter('npts', 100);
      p.addParameter('plot_axes', []);
      p.addParameter('data', ott.utils.VswfData(), ...
          @(x) isa(x, 'ott.utils.VswfData'));
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Generate grid
      ptheta = linspace(0, 2*pi, p.Results.npts);
      [r, ptheta, phi] = ott.utils.matchsize(0, ptheta(:), p.Results.phi);
      rtp = [r(:), ptheta(:), phi(:)].';
      rtp = ott.utils.sanitiseRtp(rtp);

      % Calculate fields
      [E, vswfData] = beam.efarfield(rtp, 'data', p.Results.data, ...
          unmatched{:});

      % Generate visualisation data
      imout = beam.VisualisationData(p.Results.field, E);

      % Generate visualisation
      beam.visFarfieldSlicePlot(imout, ptheta, nargout, varargin{:});

      % Assign outputs
      if nargout == 0
        clear imout ptheta vswfData;
      end
    end
  end

  methods (Sealed)
    
    function sbeam = scatter(ibeam, particle, varargin)
      % Calculate how a particle scatters the beam
      %
      % Usage
      %   sbeam = scatter(ibeam, particle, ...)
      %
      % Returns
      %   - sbeam (ott.beam.Scattered) -- Scattered beam encapsulating
      %     the particle, incident beam, scattered beams(s).  For a
      %     method which doesn't create a scattered beam, see
      %     :meth:`scatterBsc`.
      %
      % Parameters
      %   - particle (ott.particle.Particle) -- Particle with
      %     T-matrix properties (possibly internal and external).
      %
      % Optional named arguments
      %   - position (3x1 numeric) -- Particle position.
      %     Default: ``particle.position``.
      %
      %   - rotation (3x3 numeric) -- Particle rotation.
      %     Default: ``particle.rotation``.
      
      if isa(particle, 'ott.tmatrix.Tmatrix')
        particle = ott.particle.Fixed('tmatrix', particle);
      end
      
      p = inputParser;
      p.addParameter('position', particle.position);
      p.addParameter('rotation', particle.rotation);
      p.parse(varargin{:});
      
      particle.position = p.Results.position;
      particle.rotation = p.Results.rotation;
      
      sbeam(numel(ibeam)) = ott.beam.Empty();
      for ii = 1:numel(ibeam)
        sbeam(ii) = ibeam(ii).scatterInternal(particle);
      end
    end

    %
    % Mathematical operators
    %

    function beam = times(beam, other)
      % Scale the beam fields (only supports scalar multiplication)
      %
      % Usage
      %   beam = beam .* scale

      % Defer for now, may change in future
      beam = beam * other;
    end

    function beam = mtimes(beam, other)
      % Scale the beam fields (only supports scalar multiplication)
      %
      % Usage
      %   beam = beam * scalar

      % Swap inputs if required
      if isnumeric(beam)
        [beam, other] = deal(other, beam);
      end

      assert(~isa(other, 'ott.beam.Beam'), ...
          'At least one input must be numeric');

      for ii = 1:numel(beam)
        beam(ii).scale = beam(ii).scale * other;
      end
    end

    function beam = rdivide(beam, other)
      % Scalar division of beam fields
      %
      % Usage
      %   beam = beam ./ value

      assert(isnumeric(other), 'second argument must be numeric');

      for ii = 1:numel(beam)
        beam(ii).scale = beam(ii).scale ./ other;
      end
    end

    function beam = mrdivide(beam, other)
      % Scalar division of beam fields
      %
      % Usage
      %   beam = beam / value

      % Defer
      beam = beam ./ other;
    end

    function beam = uminus(beam)
      % Unitary minus: flips beam field intensity
      %
      % Usage
      %   beam = -beam;

      for ii = 1:numel(beam)
        beam(ii).scale = -beam(ii).scale;
      end
    end

    function beam = plus(a, b)
      % Combine beams coherently
      %
      % Usage
      %   beam = beam1 + beam2

      assert(isscalar(a) && isscalar(b), 'only works for scalar inputs');

      beam = ott.beam.Array([a, b], 'arrayType', 'coherent');
    end

    function beam = minus(a, b)
      % Combine beams coherently
      %
      % Usage
      %   beam = beam1 - beam2

      assert(isscalar(a) && isscalar(b), 'only works for scalar inputs');

      beam = ott.beam.Array([a, -b], 'arrayType', 'coherent');
    end

    function beam = or(a, b)
      % Combine beams incoherently
      %
      % Usage
      %   beam = beam1 | beam2

      assert(isscalar(a) && isscalar(b), 'only works for scalar inputs');

      beam = ott.beam.Array([a, -b], 'arrayType', 'incoherent');
    end
  end

  methods (Hidden)
    function rng = defaultVisRangeInternal(~)
      rng = [1,1];
    end
  end
  
  methods (Sealed, Hidden)
    function rng = defaultVisRange(beam)
      % Get the default visualisation range
      
      rng = beam(1).defaultVisRangeInternal();
      for ii = 2:numel(beam)
        rng = max(rng, beam(ii).defaultVisRangeInternal());
      end
    end
  end

  methods (Static, Hidden)
    function visShowPlot(imout, ranges, our_nargout, varargin)
      % Generate plots for visNearfield and visFarfield

      p = inputParser;
      p.addParameter('plot_axes', []);
      p.KeepUnmatched = true;
      p.parse(varargin{:});

      our_axes = p.Results.plot_axes;

      if our_nargout == 0 || ~isempty(our_axes)

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
        axis(our_axes, 'image');
      end
    end

    function visFarfieldSpherePlot(imout, XYZ, our_nargout, varargin)
      % Generate the plot for far-field sphere

      p = inputParser;
      p.addParameter('normalise', false);
      p.addParameter('type', 'sphere');
      p.addParameter('plot_axes', []);
      p.KeepUnmatched = true;
      p.parse(varargin{:});

      % Generate visualisation
      if our_nargout == 0 || ~isempty(p.Results.plot_axes)

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
            surf(our_axes, XYZ{1}, XYZ{2}, XYZ{3}, I, ...
                'facecolor','interp','edgecolor','none');
          case '3dpolar'
            surf(our_axes, abs(I).*XYZ{1},abs(I).*XYZ{2},abs(I).*XYZ{3},I,...
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
    end

    function visFarfieldSlicePlot(imout, ptheta, our_nargout, varargin)
      % Generate the plot for far-field slice

      p = inputParser;
      p.addOptional('phi', [0, pi/2], @isnumeric);
      p.addParameter('normalise', false);
      p.addParameter('plot_axes', []);
      p.KeepUnmatched = true;
      p.parse(varargin{:});

      % Generate visualisation
      if our_nargout == 0 || ~isempty(p.Results.plot_axes)

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
    end

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
  end
  
  methods (Static, Sealed, Access = protected)
    function default_object = getDefaultScalarElement
      default_object = ott.beam.Empty();
    end
  end

  methods % Getters/setters
    function speed = get.speed(beam)
      % Get the speed in the medium
      speed = beam.speed0 ./ beam.index_medium;
    end

    function beam = set.speed(beam, val)
      % Set the speed in the medium
      assert(isnumeric(val) && isscalar(val) && val > 0, ...
          'speed must be positive numeric scalar');
      beam.index_medium = beam.speed0 ./ val;
    end

    function wavenumber = get.wavenumber(beam)
      wavenumber = beam.omega ./ beam.speed;
    end
    function beam = set.wavenumber(beam, ~) %#ok<INUSD>
      error('ott:beam:Beam:set_wavenumber', ...
        'Cannot set wavenumber, set speed, omega or use setWavenumber');
    end

    function wavelength = get.wavelength(beam)
      wavelength = 2*pi ./ beam.wavenumber;
    end
    function beam = set.wavelength(beam, ~) %#ok<INUSD>
      error('ott:beam:Beam:set_wavelength', ...
        'Cannot set wavelength, set speed, omega or use setWavelength');
    end
    
    function wavelength = get.wavelength0(beam)
      wavelength = beam.speed0 ./ beam.omega .* 2*pi;
    end
    function beam = set.wavelength0(beam, val)
      assert(isnumeric(val) && isscalar(val) && val > 0, ...
          'wavelength0 must be positive numeric scalar');
      beam.omega = 3e8/val*2*pi;
    end

    function val = get.impedance(beam)
      % Assuming mu_r = 1 and non-conductive medium
      val = 376.73 ./ beam.index_medium;
    end
  end
end

