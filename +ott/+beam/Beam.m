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
%   - index_medium    -- Refractive index of the medium
%   - wavelength      -- Wavelength in medium [m]
%   - wavenumber      -- Wavenumber in medium [1/m]
%   - omega           -- Optical angular frequency of light [1/s]
%   - speed           -- Speed of light in the medium [m/s]
%   - speed0          -- Speed of light in vacuum [m/s]
%
% Abstract properties
%   - data            -- Bsc data describing beam
%
% Force and torque related methods
%   - force           -- Calculate the change in momentum between two beams
%   - torque          -- Calculate change in angular momentum between beams
%   - spin            -- Calculate change in spin momentum between beams
%   - forcetorque     -- Calculate the force and the torque between beams
%
% Field visualisation methods
%   - visNearfield      -- Generate a visualisation around the origin
%   - visFarfield       -- Generate a visualisation at the far-field
%   - visFarfieldSlice  -- Visualise the field on a angular slice
%   - visFarfieldSphere -- Visualise the filed on a sphere
%
% Methods
%   - fields etc...
%   - containsIncoherent    -- Returns true if the has incoherent elements

% Copyright 2018-2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    index_medium   % Refractive index of the medium
    omega          % Optical angular frequency of light [1/s]
  end

  properties (Dependent)
    speed          % Speed of light in the medium [m/s]
    wavelength     % Wavelength in medium [m]
    wavenumber     % Wavenumber in medium [1/m]
  end

  properties (Constant)
    speed0 = 3e8    % Speed of light in vacuum [m/s]
  end

  methods

    %
    % Force and torque related methods
    %

    function varargout = force(ibeam, sbeam)
      % Calculate change in linear momentum between beams.
      %
      % Uses the internal bsc data of the two beams to calculate the
      % force (in Newtons).
      %
      % Usage
      %   force = ibeam.force(sbeam)
      %
      %   [fx, fy, fz] = ibeam.force(sbeam)

      % Calculate force using internal methods
      [varargout{1:nargout}] = ibeam.data.force(sbeam.data);

      % Convert units of output to SI
      for ii = 1:nargout
        varargout{ii} = varargout{ii} .* beam.power ./ beam.speed;
      end
    end

    function varargout = torque(ibeam, sbeam)
      % Calculate change in angular momentum between beams.
      %
      % Uses the internal bsc data of the two beams to calculate the
      % torque (in Newton meters).
      %
      % Usage
      %   torque = ibeam.torque(sbeam)
      %
      %   [tx, ty, tz] = ibeam.torque(sbeam)

      % Calculate force using internal methods
      [varargout{1:nargout}] = ibeam.data.torque(sbeam.data);

      % Convert units of output to SI
      for ii = 1:nargout
        varargout{ii} = varargout{ii} .* beam.power ./ beam.omega;
      end
    end

    function varargout = spin(ibeam, sbeam)
      % Calculate change in spin angular momentum between beams.
      %
      % Uses the internal bsc data of the two beams to calculate the
      % spin torque (in Newton meters).
      %
      % Usage
      %   spin = ibeam.spin(sbeam)
      %
      %   [sx, sy, sz] = ibeam.spin(sbeam)

      % Calculate force using internal methods
      [varargout{1:nargout}] = ibeam.data.spin(sbeam.data);

      % Convert units of output to SI
      for ii = 1:nargout
        varargout{ii} = varargout{ii} .* beam.power ./ beam.omega;
      end
    end

    function [force, torque, spin] = forcetorque(ibeam, sbeam)
      % Calculate change in momentum between beams
      %
      % Usage
      %   [f, t, s] = ibeam.forcetorque(sbeam) calculates the force,
      %   torque and spin between the incident beam ``ibeam`` and
      %   scattered beam ``sbeam``.
      %   Outputs 3x[N...] matrix depending on the number and shape of beams.

      % Calculate force using internal methods
      [varargout{1:nargout}] = ibeam.data.forcetorque(sbeam.data);

      % Convert to SI units
      if nargout >= 1
        varargout{1} = varargout{1} .* beam.power ./ beam.speed;
        for ii = 2:nargout
          varargout{ii} = varargout{ii} .* beam.power ./ beam.omega;
        end
      end
    end

    %
    % Field visualisation methods
    %

    function varargout = visNearfield(beam, varargin)
      % Create a visualisation of the beam
      %
      % Usage
      %   beam.visualise(...) displays an image of the beam in the current
      %   figure window.
      %
      %   [im, data] = beam.visualise(...) returns a image of the beam.
      %   If the beam object is an array, returns an image for each beam.
      %   Also returns a :class:`ott.utils.VswfData` structure for fast
      %   repeated calculations.
      %
      % See :meth:`+ott.+bsc.Bsc.visNearfield` for parameters and details.

      [varargout{1:nargout}] = beam.data.visNearfield(varargin{:});
    end

    function varargout = visFarfield(beam, varargin)
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
      % See :meth:`+ott.+bsc.Bsc.visFarfield` for parameters and details.

      [varargout{1:nargout}] = beam.data.visFarfield(varargin{:});
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
      % See :meth:`+ott.+bsc.Bsc.visFarfieldSphere` for parameters and details.

      [varargout{1:nargout}] = beam.data.visFarfieldSphere(varargin{:});
    end

    function varargout = visFarfieldSlice(beam, varargin)
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
      % See :meth:`+ott.+bsc.Bsc.visFarfieldSlice` for parameters and details.

      [varargout{1:nargout}] = beam.data.visFarfieldSlice(varargin{:});
    end

    %
    % Generic methods
    %

    function b = containsIncoherent(beam)
      % TODO
      error('not yet implemented');
    end
  end

  methods % Getters/setters
    function speed = get.speed(beam)
      % Get the speed in the medium
      speed = beam.speed0 ./ index_medium;
    end

    function beam = set.speed(beam, val)
      % Set the speed in the medium
      assert(isnumeric(val) && isscalar(val) && val > 0, ...
          'speed must be positive numeric scalar');
      beam.index_medium = beam.speed0 ./ val;
    end

    function beam = set.index_medium(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
          'index_medium must be numeric scalar');
      beam.index_medium = val;
    end

    function wavenumber = get.wavenumber(beam)
      wavenumber = beam.omega ./ beam.speed;
    end
    function beam = set.wavenumber(beam, val)
      error('Cannot set wavenumber, set speed or omega instead');
    end

    function wavelength = get.wavelength(beam)
      wavelength = 2*pi ./ beam.wavenumber;
    end
    function beam = set.wavelength(beam, val)
      error('Cannot set wavelength, set speed or omega instead');
    end

    function beam = set.omega(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
          'omega must be numeric scalar');
      beam.omega = val;
    end
  end
end

