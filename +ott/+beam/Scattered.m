classdef Scattered < ott.beam.Beam
% Describes the beam scattered from a particle.
% Inherits from :class:`ott.beam.Beam`.
%
% When the internal field and the particle are supplied, the visualisation
% and field calculation method calculate either the external or internal
% fields depending on the requested points location.  The visualisation
% methods also provide an option for showing the particle outline.
%
% Properties
%   - incident    -- Beam data incident on the particle
%   - scattered   -- Beam data scattered by the particle
%   - outgoing    -- Outgoing modes radiation from the particle
%   - incoming    -- Incoming modes scattered by the particle
%   - internal    -- Internal beam data (optional)
%   - particle    -- Particle responsible for scattering (optional)
%
% Field calculation methods
%   - efield     -- Calculate electric field around the origin
%   - hfield     -- Calculate magnetic field around the origin
%   - ehfield    -- Calculate electric and magnetic fields around the origin
%   - efieldRtp  -- Calculate electric field around the origin (sph. coords.)
%   - hfieldRtp  -- Calculate magnetic field around the origin (sph. coords.)
%   - ehfieldRtp -- Calculate electric and magnetic fields around the origin
%   - efarfield  -- Calculate electric fields in the far-field
%   - hfarfield  -- Calculate magnetic fields in the far-field
%   - ehfarfield -- Calculate electric and magnetic fields in the far-field
%   - eparaxial  -- Calculate electric fields in the paraxial far-field
%   - hparaxial  -- Calculate magnetic fields in the paraxial far-field
%   - ehparaxial -- Calculate electric and magnetic paraxial far-fields
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

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    incident    % Incident beam
    scattered   % Scattered beam (optional)
    internal    % Internal beam (optional)
    particle    % Particle instance (optional)
  end

  properties (Dependent)
    outgoing     % Outgoing modes radiation from the particle
    incoming     % Incoming modes scattered by the particle
    shape        % Either particle or the shape property of particle
    index_medium % Refractive index of incident/scattered medium
    omega        % Optical frequency
  end

  methods
    function bm = Scattered(varargin)
      % Construct a new scattered beam representation
      %
      % Usage
      %   beam = Scattered(scattered, incident, particle, internal, ...)
      %
      % Optional named arguments
      %   - scattered (ott.beam.Beam) -- The scattered beam data.  This will
      %     typically be an outgoing scattered beam (such as that generated
      %     using most T-matrix classes).  Default: ``[]``.
      %
      %   - incident (ott.beam.Beam) -- The incident beam data.  This is
      %     required for total fields to look correct.
      %     Default: ``ott.beam.Empty()``.
      %
      %   - particle (ott.shape.Shape | ott.particle.Particle | []) -- An
      %     object with a geometry, used for deciding if points are inside
      %     the particle.  If omitted, all points are assumed valid.
      %     Default: ``[]``.
      %
      %   - internal (ott.beam.Beam | []) -- Fields to display for points
      %     that are inside the shape.  If omitted, NaN is used for internal
      %     points.  Default: ``[]``.

      p = inputParser;
      p.addOptional('scattered', [], ...
          @(x) isempty(x) || isa(x, 'ott.beam.Beam'));
      p.addOptional('incident', ott.beam.Empty(), ...
          @(x) isa(x, 'ott.beam.Beam'));
      p.addOptional('particle', [], ...
          @(x) isempty(x) || isa(x, 'ott.particle.Particle') ...
          || isa(x, 'ott.shape.Shape'));
      p.addOptional('internal', [], ...
          @(x) isempty(x) || isa(x, 'ott.beam.Beam'));
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      bm = bm@ott.beam.Beam(unmatched{:});
      bm.scattered = p.Results.scattered;
      bm.incident = p.Results.incident;
      bm.particle = p.Results.particle;
      bm.internal = p.Results.internal;
    end

    %
    % Field calculation methods
    %

    function [E, vswfData] = efield(beam, xyz, varargin)
      % Calculate the electric field (in SI units)
      %
      % Usage
      %   [E, vswfData] = beam.efield(xyz, ...)
      %
      % Parameters
      %   - xyz (3xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of meters.  Packaged [x;y;z].
      %
      % Optional named parameters
      %   - type (enum) -- Type of field to calculate.  Can be either
      %     'incident', 'scattered', or 'total'.  When only the
      %     ``scattered`` property is set, 'scattered' and 'total' behave
      %     the same.  Default: ``'total'``.
      %
      % Unmatched parameters are passed to the corresponding method
      % of the incident/scattered :class:`Beam` class.

      % Get external data and unmatched arguments
      [odata, unmatched] = beam.getAndParseType(varargin{:});

      % Find internal points
      binside = beam.insideXyz(xyz);

      % Calculate external fields
      E = zeros(size(xyz));
      [Ef, vswfData] = odata.efield(xyz(:, ~binside), unmatched{:});
      E(:, ~binside) = Ef.vxyz;

      % Calculate internal fields
      if ~isempty(beam.internal)
        xyzI = beam.global2local(xyz);
        [Ef, vswfData] = beam.internal.efield(...
            xyzI(:, binside), unmatched{:}, 'data', vswfData);
        E(:, binside) = Ef.vxyz;
      else
        E(:, binside) = nan;
      end
      
      % Package output
      E = ott.utils.FieldVectorCart(E, xyz);

      if nargout == 1
        clear vswfData;
      end
    end

    function [H, vswfData] = hfield(beam, xyz, varargin)
      % Calculate the magnetic field (in SI units)
      %
      % Usage
      %   [H, vswfData] = beam.hfield(xyz, ...)
      %
      % Parameters
      %   - xyz (3xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of meters.  Packaged [x;y;z].
      %
      % Optional named parameters
      %   - type (enum) -- Type of field to calculate.  Can be either
      %     'incident', 'scattered', or 'total'.  When only the
      %     ``scattered`` property is set, 'scattered' and 'total' behave
      %     the same.  Default: ``'total'``.
      %
      % Unmatched parameters are passed to the corresponding method
      % of the incident/scattered :class:`Beam` class.

      % Get external data and unmatched arguments
      [odata, unmatched] = beam.getAndParseType(varargin{:});

      % Find internal points
      binside = beam.insideXyz(xyz);

      % Calculate external fields
      H = zeros(size(xyz));
      [Hf, vswfData] = odata.hfield(xyz(:, ~binside), unmatched{:});
      H(:, ~binside) = Hf.vxyz;

      % Calculate internal fields
      if ~isempty(beam.internal)
        xyzI = beam.global2local(xyz);
        [Hf, vswfData] = beam.internal.hfield(...
            xyzI(:, binside), unmatched{:}, 'data', vswfData);
        H(:, binside) = Hf.vxyz;
      else
        H(:, binside) = nan;
      end
      
      % Package output
      H = ott.utils.FieldVectorCart(H, xyz);

      if nargout == 1
        clear vswfData;
      end
    end

    function [E, H, vswfData] = ehfield(beam, xyz, varargin)
      % Calculate the electric and magnetic field (in SI units)
      %
      % Usage
      %   [E, H, vswfData] = beam.ehfield(xyz, ...)
      %
      % Parameters
      %   - xyz (3xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of meters.  Packaged [x;y;z].
      %
      % Optional named parameters
      %   - type (enum) -- Type of field to calculate.  Can be either
      %     'incident', 'scattered', or 'total'.  When only the
      %     ``scattered`` property is set, 'scattered' and 'total' behave
      %     the same.  Default: ``'total'``.
      %
      % Unmatched parameters are passed to the corresponding method
      % of the incident/scattered :class:`Beam` class.

      % Get external data and unmatched arguments
      [odata, unmatched] = beam.getAndParseType(varargin{:});

      % Find internal points
      binside = beam.insideXyz(xyz);

      % Calculate external fields
      E = zeros(size(xyz));
      H = E;
      [Ef, Hf, vswfData] = odata.ehfield(xyz(:, ~binside), unmatched{:});
      E(:, ~binside) = Ef.vxyz;
      H(:, ~binside) = Hf.vxyz;

      % Calculate internal fields
      if ~isempty(beam.internal)
        xyzI = beam.global2local(xyz);
        [Ef, Hf, vswfData] = beam.internal.ehfield(...
            xyzI(:, binside), unmatched{:}, 'data', vswfData);
        E(:, binside) = Ef.vxyz;
        H(:, binside) = Hf.vxyz;
      else
        E(:, binside) = nan;
        H(:, binside) = nan;
      end
      
      % Package output
      E = ott.utils.FieldVectorCart(E, xyz);
      H = ott.utils.FieldVectorCart(H, xyz);

      if nargout < 3
        clear vswfData;
      end
    end

    function [E, vswfData] = efieldRtp(beam, rtp, varargin)
      % Calculate the electric field (in SI units)
      %
      % Usage
      %   [E, vswfData] = beam.efieldRtp(rtp, ...)
      %
      % Parameters
      %   - rtp (3xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of meters/radians.  Packaged [r;theta;phi].
      %
      % Optional named parameters
      %   - type (enum) -- Type of field to calculate.  Can be either
      %     'incident', 'scattered', or 'total'.  When only the
      %     ``scattered`` property is set, 'scattered' and 'total' behave
      %     the same.  Default: ``'total'``.
      %
      % Unmatched parameters are passed to the corresponding method
      % of the incident/scattered :class:`Beam` class.

      % Get external data and unmatched arguments
      [odata, unmatched] = beam.getAndParseType(varargin{:});

      % Find internal points
      binside = beam.insideRtp(rtp);

      % Calculate external fields
      E = zeros(size(rtp));
      [Ef, vswfData] = odata.efieldRtp(rtp(:, ~binside), unmatched{:});
      E(:, ~binside) = Ef.vrtp;

      % Calculate internal fields
      if ~isempty(beam.internal)
        rtpI = beam.global2localRtp(rtp);
        [Ef, vswfData] = beam.internal.efieldRtp(...
            rtpI(:, binside), unmatched{:}, 'data', vswfData);
        E(:, binside) = Ef.vrtp;
      else
        E(:, binside) = nan;
      end
      
      % Package output
      E = ott.utils.FieldVectorSph(E, rtp);

      if nargout == 1
        clear vswfData;
      end
    end

    function [H, vswfData] = hfieldRtp(beam, rtp, varargin)
      % Calculate the magnetic field (in SI units)
      %
      % Usage
      %   [E, vswfData] = beam.hfieldRtp(rtp, ...)
      %
      % Parameters
      %   - rtp (3xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of meters/radians.  Packaged [r;theta;phi].
      %
      % Optional named parameters
      %   - type (enum) -- Type of field to calculate.  Can be either
      %     'incident', 'scattered', or 'total'.  When only the
      %     ``scattered`` property is set, 'scattered' and 'total' behave
      %     the same.  Default: ``'total'``.
      %
      % Unmatched parameters are passed to the corresponding method
      % of the incident/scattered :class:`Beam` class.

      % Get external data and unmatched arguments
      [odata, unmatched] = beam.getAndParseType(varargin{:});

      % Find internal points
      binside = beam.insideRtp(rtp);

      % Calculate external fields
      H = zeros(size(rtp));
      [Hf, vswfData] = odata.hfieldRtp(rtp(:, ~binside), unmatched{:});
      H(:, ~binside) = Hf.vrtp;

      % Calculate internal fields
      if ~isempty(beam.internal)
        rtpI = beam.global2localRtp(rtp);
        [Hf, vswfData] = beam.internal.hfieldRtp(...
            rtpI(:, binside), unmatched{:}, 'data', vswfData);
        H(:, binside) = Hf.vrtp;
      else
        H(:, binside) = nan;
      end
      
      % Package output
      H = ott.utils.FieldVectorSph(H, rtp);
      
      if nargout == 1
        clear vswfData;
      end
    end

    function [E, H, vswfData] = ehfieldRtp(beam, rtp, varargin)
      % Calculate the electric and magnetic field (in SI units)
      %
      % Usage
      %   [E, H, vswfData] = beam.ehfieldRtp(rtp, ...)
      %
      % Parameters
      %   - rtp (3xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of meters/radians.  Packaged [r;theta;phi].
      %
      % Optional named parameters
      %   - type (enum) -- Type of field to calculate.  Can be either
      %     'incident', 'scattered', or 'total'.  When only the
      %     ``scattered`` property is set, 'scattered' and 'total' behave
      %     the same.  Default: ``'total'``.
      %
      % Unmatched parameters are passed to the corresponding method
      % of the incident/scattered :class:`Beam` class.

      % Get external data and unmatched arguments
      [odata, unmatched] = beam.getAndParseType(varargin{:});

      % Find internal points
      binside = beam.insideRtp(rtp);

      % Calculate external fields
      E = zeros(size(rtp));
      H = E;
      [Ef, Hf, vswfData] = odata.ehfieldRtp(rtp(:, ~binside), unmatched{:});
      E(:, ~binside) = Ef.vrtp;
      H(:, ~binside) = Hf.vrtp;

      % Calculate internal fields
      if ~isempty(beam.internal)
        rtpI = beam.global2localRtp(rtp);
        [Ef, Hf, vswfData] = beam.internal.ehfieldRtp(...
            rtpI(:, binside), unmatched{:}, 'data', vswfData);
        E(:, binside) = Ef.vrtp;
        H(:, binside) = Hf.vrtp;
      else
        E(:, binside) = nan;
        H(:, binside) = nan;
      end
      
      % Package output
      E = ott.utils.FieldVectorSph(E, rtp);
      H = ott.utils.FieldVectorSph(H, rtp);

      if nargout < 3
        clear vswfData;
      end
    end

    function varargout = efarfield(beam, rtp, varargin)
      % Calculate the electric field (in SI units)
      %
      % Calculates the fields in the local reference frame and then
      % applies phase/rotations to give the global far-field.
      %
      % Usage
      %   [E, vswfData] = beam.efarfield(rtp, ...)
      %
      % Parameters
      %   - rtp (3xN | 2xN numeric) -- Spherical coordinates for field
      %     calculation. Packaged [r; theta; phi] or [theta; phi].
      %
      % Optional named parameters
      %   - type (enum) -- Type of field to calculate.  Can be either
      %     'incident', 'scattered', or 'total'.  When only the
      %     ``scattered`` property is set, 'scattered' and 'total' behave
      %     the same.  Default: ``'total'``.
      %
      % Unmatched parameters are passed to the corresponding method
      % of the incident/scattered :class:`Beam` class.
      
      % TODO: Rotations
      
      scat_position = beam.position;
      beam.position = [0;0;0];

      % Get external data and unmatched arguments
      [odata, unmatched] = beam.getAndParseType(varargin{:});

      % Calculate fields (in local reference frame)
      [varargout{1:nargout}] = odata.efarfield(rtp, unmatched{:});
      
      if nargout >= 1
        varargout{1} = beam.translateFarfields(varargout{1}, -scat_position);
      end
    end

    function varargout = hfarfield(beam, rtp, varargin)
      % Calculate the magnetic field (in SI units)
      %
      % Calculates the fields in the local reference frame and then
      % applies phase/rotations to give the global far-field.
      %
      % Usage
      %   [H, vswfData] = beam.hfarfield(rtp, ...)
      %
      % Parameters
      %   - rtp (3xN | 2xN numeric) -- Spherical coordinates for field
      %     calculation. Packaged [r; theta; phi] or [theta; phi].
      %
      % Optional named parameters
      %   - type (enum) -- Type of field to calculate.  Can be either
      %     'incident', 'scattered', or 'total'.  When only the
      %     ``scattered`` property is set, 'scattered' and 'total' behave
      %     the same.  Default: ``'total'``.
      %
      % Unmatched parameters are passed to the corresponding method
      % of the incident/scattered :class:`Beam` class.
      
      % TODO: Rotations
      
      scat_position = beam.position;
      beam.position = [0;0;0];

      % Get external data and unmatched arguments
      [odata, unmatched] = beam.getAndParseType(varargin{:});

      % Calculate fields
      [varargout{1:nargout}] = odata.hfarfield(rtp, unmatched{:});
      
      if nargout >= 1
        varargout{1} = beam.translateFarfields(varargout{1}, -scat_position);
      end
    end

    function varargout = ehfarfield(beam, rtp, varargin)
      % Calculate the electric and magnetic field (in SI units)
      %
      % Calculates the fields in the local reference frame and then
      % applies phase/rotations to give the global far-field.
      %
      % Usage
      %   [E, H, vswfData] = beam.ehfarfield(rtp, ...)
      %
      % Parameters
      %   - rtp (3xN | 2xN numeric) -- Spherical coordinates for field
      %     calculation. Packaged [r; theta; phi] or [theta; phi].
      %
      % Optional named parameters
      %   - type (enum) -- Type of field to calculate.  Can be either
      %     'incident', 'scattered', or 'total'.  When only the
      %     ``scattered`` property is set, 'scattered' and 'total' behave
      %     the same.  Default: ``'total'``.
      %
      % Unmatched parameters are passed to the corresponding method
      % of the incident/scattered :class:`Beam` class.
      
      % TODO: Rotations
      
      scat_position = beam.position;
      beam.position = [0;0;0];

      % Get external data and unmatched arguments
      [odata, unmatched] = beam.getAndParseType(varargin{:});

      % Calculate fields
      [varargout{1:nargout}] = odata.ehfarfield(rtp, unmatched{:});
      
      if nargout >= 1
        [varargout{1}, psi] = beam.translateFarfields(varargout{1}, -scat_position);
        if nargout >= 2
          varargout{2} = varargout{2} .* psi;
        end
      end
    end
    
    %
    % Field visualisation methods
    %
    
    function [im, XY, data] = visNearfield(beam, varargin)
      % Create a visualisation of the beam.
      %
      % If the shape property of the beam is set, also draws a contour
      % around the edge of the shape.  Does this recursively for all
      % incident scattered beams.
      %
      % Usage
      %   beam.visNearfield(...) -- display an image of the beam in
      %   the current figure window.
      %
      %   [im, XY, data] = beam.visualise(...) -- returns a image of the
      %   beam and a cell array of the coordinates used in the calculation.
      %   If the beam object is an array, returns an image for each beam.
      %
      % Optional named parameters
      %   - drawContour (logical) -- If the contour should be drawn.
      %     Default: ``true``.  Only has effect if plot shown.
      %
      %   - contourSmooth (numeric) -- Smothing factor (pixels) for
      %     Gaussian blur of contour.  Good for smooth shapes, not good
      %     for cubes.  Default: ``2``.
      %
      % For additional parameters/usage information, see :class:`Beam`.
      
      p = inputParser;
      p.addParameter('plot_axes', []);
      p.addParameter('drawContour', true);
      p.addParameter('contourSmooth', 2);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);
      
      % Get/setup the axes before calling parent
      our_axes = p.Results.plot_axes;
      if (nargout == 0 || ~isempty(our_axes))

        % Get the default axes
        if isempty(our_axes)
          our_axes = gca();
        end
      end
      
      % Call base for most of the work
      [im, XY, data] = visNearfield@ott.beam.Beam(beam, ...
          'plot_axes', our_axes, unmatched{:});

      % Do we need to draw a contour
      if ~isempty(our_axes) && p.Results.drawContour

        % Get the default axes
        if isempty(our_axes)
          our_axes = gca();
        end
      
        % Get data for inside/outside
        inside = beam.insideXyz(XY{3}, 'recursive', true);
        inside = reshape(inside, numel(XY{2}), numel(XY{1}));
        
        if ~all(inside(1) == inside(:))

          % Smooth contour data (makes things look nicer)
          if p.Results.contourSmooth ~= 0
            inside = ott.utils.gaussfilt(double(inside), p.Results.contourSmooth);
          end

          % Draw contour on plot
          hold on;
          [~, ch] = contour(our_axes, XY{1}, XY{2}, inside, [0.5, 0.5], 'w--');
          ch.LineWidth = 2;
          hold off;
        end
      end
      
      if nargout == 0
        clear im XY data
      end
    end

    %
    % Force/torque related methods
    %

    function varargout = force(beam, obeam, varargin)
      % Calculate change in linear momentum between beams.
      %
      % Uses the internal bsc data of the two beams to calculate the
      % force (in Newtons).
      %
      % Usage
      %   force = beam.force() -- Calculate force for this scattered
      %   beam (uses ``incident`` and ``scattered``).
      %
      %   force = ibeam.force(sbeam) -- Calculate force between this
      %   beam and another beam.
      %
      %   [fx, fy, fz] = ibeam.force(...) -- As above, but unpacking
      %   the outputs for X/Y/Z.

      if nargin == 1
        assert(~isempty(beam.scattered) && ~isempty(beam.incident), ...
          'scattered and incident must be set for force without arguments');
        [odata, ~] = beam.getAndParseType('type', 'outgoing');
        [idata, ~] = beam.getAndParseType('type', 'incoming');
        [varargout{1:nargout}] = idata.force(odata);
      else
        [odata, ~] = beam.getAndParseType('type', 'total');
        [varargout{1:nargout}] = odata.force(obeam, varargin{:});
      end
    end

    function varargout = torque(beam, obeam, varargin)
      % Calculate change in angular momentum between beams.
      %
      % Uses the internal bsc data of the two beams to calculate the
      % torque (in Newton meters).
      %
      % Usage
      %   torque = beam.torque() -- Calculate torque for this scattered
      %   beam (uses ``incident`` and ``scattered``).
      %
      %   torque = ibeam.torque(sbeam) -- Calculate torque between this
      %   beam and another beam.
      %
      %   [tx, ty, tz] = ibeam.torque(...) -- As above, but unpacking
      %   the outputs for X/Y/Z.

      if nargin == 1
        assert(~isempty(beam.scattered) && ~isempty(beam.incident), ...
          'scattered and incident must be set for torque without arguments');
        [odata, ~] = beam.getAndParseType('type', 'outgoing');
        [idata, ~] = beam.getAndParseType('type', 'incoming');
        [varargout{1:nargout}] = idata.torque(odata);
      else
        [odata, ~] = beam.getAndParseType('type', 'total');
        [varargout{1:nargout}] = odata.torque(obeam, varargin{:});
      end
    end

    function varargout = spin(beam, obeam, varargin)
      % Calculate change in spin angular momentum between beams.
      %
      % Uses the internal bsc data of the two beams to calculate the
      % spin torque (in Newton meters).
      %
      % Usage
      %   spin = beam.spin() -- Calculate spin for this scattered
      %   beam (uses ``incident`` and ``scattered``).
      %
      %   spin = ibeam.spin(sbeam) -- Calculate spin between this
      %   beam and another beam.
      %
      %   [sx, sy, sz] = ibeam.spin(...) -- As above, but unpacking
      %   the outputs for X/Y/Z.

      if nargin == 1
        assert(~isempty(beam.scattered) && ~isempty(beam.incident), ...
          'scattered and incident must be set for spin without arguments');
        [odata, ~] = beam.getAndParseType('type', 'outgoing');
        [idata, ~] = beam.getAndParseType('type', 'incoming');
        [varargout{1:nargout}] = idata.spin(odata);
      else
        [odata, ~] = beam.getAndParseType('type', 'total');
        [varargout{1:nargout}] = odata.spin(obeam, varargin{:});
      end
    end

    function varargout = forcetorque(beam, obeam, varargin)
      % Calculate change in momentum between beams
      %
      % Usage
      %   [f, t, s] = ibeam.forcetorque(sbeam) -- calculates the force,
      %   torque and spin for this scattered beam (uses ``incident``
      %   and ``scattered``).
      %
      %   [f, t, s] = ibeam.forcetorque(sbeam) -- calculates the force,
      %   torque and spin between the incident beam ``ibeam`` and
      %   scattered beam ``sbeam``.
      %
      %   Outputs 3x[N...] matrix depending on the number and shape of beams.

      if nargin == 1
        assert(~isempty(beam.scattered) && ~isempty(beam.incident), ...
          'scattered and incident must be set for forcetorque without arguments');
        [odata, ~] = beam.getAndParseType('type', 'outgoing');
        [idata, ~] = beam.getAndParseType('type', 'incoming');
        [varargout{1:nargout}] = idata.forcetorque(odata);
      else
        [odata, ~] = beam.getAndParseType('type', 'total');
        [varargout{1:nargout}] = odata.forcetorque(obeam, varargin{:});
      end
    end

    %
    % Generic methods
    %
    
    function sbeam = scatterInternal(beam, particle, varargin)
      % Calculate how a particle scatters the beam.
      %
      % This does not include optical interaction between particles.
      %
      % Usage
      %   sbeam = scatter(ibeam, particle)
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
      
      if beam.insideXyz(particle.position)
        % Internal scattering
        [odata, ~] = beam.getAndParseType('type', 'internal', varargin{:});
      else
        % External scattering
        [odata, ~] = beam.getAndParseType('type', 'total', varargin{:});
      end
      
      % TODO: Check if particle overlaps boundary
      
      sbeam = odata.scatterInternal(particle);
      sbeam.incident = beam.translateXyz(-sbeam.position) ...
          .rotate(sbeam.rotation.');
    end

    function binside = insideRtp(beam, rtp)
      % Calculate if points are inside the particle (Spherical coordi.)
      %
      % Uses the :meth:`+ott.+shape.Shape.insideRtp` method of particle.
      % If no particle is provided, returns a vector of false.

      rtp = beam.global2local(rtp);
      
      if ~isempty(beam.shape)
        rtp = beam.particle.global2localRtp(rtp);
        binside = beam.shape.insideRtp(rtp);
      else
        if isempty(beam.scattered) && ~isempty(beam.internal)
          binside = true(1, size(rtp, 2));
        else % Only have incident and/or scattered
          binside = false(1, size(rtp, 2));
        end
      end
    end

    function binside = insideXyz(beam, xyz, varargin)
      % Calculate if points are inside the particle.
      %
      % Uses the :meth:`+ott.+shape.Shape.insideXyz` method of particle.
      % If no particle is provided, returns a vector of false.
      %
      % Usage
      %   b = beam.insideXyz(xyz, ...)
      %
      % Optional named arguments
      %   - recursive (logical) -- If true, and the incident beam is
      %     a scattered beam, also calls the incident beam internal method.
      %     Default: ``false``
      
      p = inputParser;
      p.addParameter('recursive', false);
      p.parse(varargin{:});
      
      xyz = beam.global2local(xyz);
      
      if ~isempty(beam.shape)
        xyz = beam.particle.global2local(xyz);
        binside = beam.shape.insideXyz(xyz);
      else
        if isempty(beam.scattered) && ~isempty(beam.internal)
          binside = true(1, size(xyz, 2));
        else % Only have incident and/or scattered
          binside = false(1, size(xyz, 2));
        end
      end
      
      if p.Results.recursive && isa(beam.incident, 'ott.beam.Scattered')
        binside = binside | beam.incident.insideXyz(xyz, 'recursive', true);
      end
    end
  end

  methods (Hidden)
    function [odata, unmatched] = getAndParseType(beam, varargin)
      % Get and parse Type parameter

      p = inputParser;
      p.addParameter('type', 'total');
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Get beam data for exterior
      switch p.Results.type
        case 'total'
          odata = ott.beam.Array([beam.incident, 2*beam.scattered], ...
              'arrayType', 'coherent');
        case 'incident'
          odata = beam.incident;
        case 'scattered'
          odata = beam.scattered;
        case 'outgoing'
          odata = beam.outgoing;
        case 'incoming'
          odata = beam.incoming;
        otherwise
          error('ott:beam:Scattered:unknown_type', ...
              'Unknown type specified, must be incident, scattered or total');
      end
      
      odata = odata.rotate(beam.rotation).translateXyz(beam.position);
    end

    function val = defaultVisRangeInternal(beam)
      [odata, ~] = beam.getAndParseType('type', 'total');
      val = odata.defaultVisRange;
    end
  end

  methods % Getters/setters
    function val = get.omega(beam)
      if ~isempty(beam.scattered)
        val = beam.scattered.omega;
      elseif ~isempty(beam.internal)
        val = beam.internal.omega;
      elseif ~isempty(beam.incident)
        val = beam.incident.omega;
      else
        error('Unable to get beam optical frequency');
      end
    end
    
    function val = get.index_medium(beam)
      if ~isempty(beam.scattered)
        val = beam.scattered.index_medium;
      elseif ~isempty(beam.incident)
        val = beam.incident.index_medium;
      else
        % TODO: We could also use the particle relative index and internal?
        error('unable to get beam refractive index');
      end
    end
    
    function beam = set.incident(beam, val)
      assert(isa(val, 'ott.beam.Beam') && isscalar(val), ...
          'incident must be a single ott.beam.Beam instance');
      beam.incident = val;
    end

    function beam = set.scattered(beam, val)
      if isempty(val)
        beam.scattered = [];
      else
        assert(isa(val, 'ott.beam.Beam'), ...
            'scattered must be [] | ott.beam.Beam');
        assert(isscalar(val), 'scattered must be a scalar');
        beam.scattered = val;
      end
    end
    
    function beam = get.outgoing(beam)
      sbsc = ott.bsc.Bsc(beam.scattered);
      ibsc = ott.bsc.Bsc(beam.incident, [sbsc.Nmax]);
    
      arrayType = 'coherent';
      if isa(beam.incident, 'ott.beam.ArrayType')
        arrayType = beam.incident.arrayType;
      end

      beam = ott.beam.BscOutgoing(ibsc + 2*sbsc, ...
        'omega', beam.omega, 'index_medium', beam.index_medium, ...
        'arrayType', arrayType);
    end
    
    function beam = get.incoming(beam)
      
      % Getting a sbsc instance just for Nmax seems wasteful
      sbsc = ott.bsc.Bsc(beam.scattered);
      ibsc = ott.bsc.Bsc(beam.incident, [sbsc.Nmax]);
    
      arrayType = 'coherent';
      if isa(beam.incident, 'ott.beam.ArrayType')
        arrayType = beam.incident.arrayType;
      end
      
      beam = ott.beam.BscFinite(ibsc, ...
        'omega', beam.omega, 'index_medium', beam.index_medium, ...
        'arrayType', arrayType);
    end
      
    function beam = set.particle(beam, val)
      if isempty(val)
        beam.particle = [];
      else
        assert(isa(val, 'ott.particle.Particle') ...
            || isa(val, 'ott.shape.Shape'), ...
            'particle must be [] | ott.shape.Shape | ott.particle.Particle');
        assert(isscalar(val), 'particle must be a scalar');
        beam.particle = val;
      end
    end

    function beam = set.internal(beam, val)
      if isempty(val)
        beam.internal = [];
      else
        assert(isa(val, 'ott.beam.Beam'), ...
            'internal must be [] | ott.beam.Beam');
        assert(isscalar(val), 'internal must be a scalar');
        beam.internal = val;
      end
    end

    function val = get.shape(beam)
      if isempty(beam.particle) || isa(beam.particle, 'ott.shape.Shape')
        val = beam.particle;
      else
        val = beam.particle.shape;
      end
    end
    function beam = set.shape(beam, val)
      if isempty(val)
        beam.particle = [];
      else
        assert(isa(val, 'ott.shape.Shape') && isscalar(val), ...
            'shape must be singular shape object');
        beam.particle = val;
      end
    end
  end
end
