classdef PlaneWave < ott.beam.properties.PlaneWaveArray ...
    & ott.beam.properties.AnyArrayType ...
    & ott.beam.utils.InfinitePower ...
    & ott.beam.Beam
% Description of a plane wave beam.
% Inherits from :class:`abstract.PlaneWave`.
%
% The class has separate vectors for the field and the polarisation.
% This allows the field to be zero without loosing information about
% the polarisation of the plane wave.
%
% Units of the fields depend on units used for the parameters.
%
% Static methods
%   - FromFieldVectors -- Construct a spectrum from EH field vectors
%   - FromNearfield -- Discrete plane wave spectrum from near-field slice
%   - FromFarfield  -- Discrete plane wave spectrum from a beam's far-field
%   - FromParaxial  -- Discrete plane wave spectrum from paraxial far-field
%
% Methods
%   - rotate      -- Rotate the direction and polarisation
%   - rotate*     -- Rotate the particle around the X,Y,Z axis
%
% Field visualisation methods
%   - visualise -- Generate a visualisation around the origin
%   - visualiseRays -- Generate a quiver-style visualisation of directions
%   - visualiseFarfield -- Generate a visualisation at the far-field
%   - visualiseFarfieldSlice -- Visualise the field on a angular slice
%   - visualiseFarfieldSphere -- Visualise the filed on a sphere
%
% Properties
%   - origin        -- Ray origins, 3xN array (default [0;0;0])
%   - direction     -- Direction of propagation (3xN Cartesian)
%   - field         -- Field parallel and perpendicular to polarisation
%   - polarisation  -- Primary polarisation direction
%
% Properties inherited from Beam
%   - power         -- The power of the beam (infinite)
%   - wavelength    -- Wavelength of beam in medium (default: 1.0)
%   - permittivity  -- Material relative permittivity (default: 1.0)
%   - permeability  -- Material relative permeability (default: 1.0)
%   - wavenumber    -- Wave-number of beam in medium
%   - impedance     -- Impedance of the medium
%
% Supported casts
%   - vswf.Bsc      -- Default Bsc cast, uses vswf.PlaneWave
%   - vswf.PlaneWave

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Static)
    function args = likeProperties(other, args)
      % Construct an array of like-properties
      args = ott.beam.properties.AnyArrayType.likeProperties(other, args);
      args = ott.beam.properties.PlaneWaveArray.likeProperties(other, args);
      args = ott.beam.Beam.likeProperties(other, args);
    end

    function beam = like(other, varargin)
      % Create a beam like another beam
      %
      % Usage
      %   beam = Gaussian.like(other, ...)
      %
      % See constructor for arguments.

      args = ott.beam.PlaneWave.likeProperties(other, varargin);
      beam = ott.beam.PlaneWave(args{:});
    end

    function beam = empty(varargin)
      % Construct an emtpy beam array
      %
      % Usage
      %   beam = ott.beam.PlaneWave.empty(...)
      %
      % Additional parameters are passed to the constructor.
      
      empt = zeros(3, 0);
      beam = ott.beam.PlaneWave('direction', empt, ...
        'polarisation', empt, ...
        'field', empt(1, :), 'origin', empt, varargin{:});
    end
    
    function beam = FromFieldVectors(E, H)
      % Construct a new Ray from E and H field vectors
      %
      % Usage
      %   beam = FromFieldVectors(E, H)

      assert(isa(E, 'ott.utils.FieldVector'), 'E must be a FieldVector');
      assert(isa(H, 'ott.utils.FieldVector'), 'H must be a FieldVector');

      % Get coordinates
      if strcmpi(E.basis, 'spherical')
        xyz = ott.utils.rtp2xyz(E.locations);
      elseif strcmpi(E.basis, 'cartesian')
        xyz = E.locations;
      else
        error('Unsupported field vector type');
      end

      Exyz = E.vxyz;
      Hxyz = H.vxyz;

      % Calculate Poynting vector
      S = 0.5.*real(cross(Exyz, conj(Hxyz)));
      Snorm = vecnorm(S, 2, 1);

      % Filter points with no intensity
      mask = Snorm == 0;
      xyz(:, mask) = [];
      S(:, mask) = [];
      Snorm(:, mask) = [];

      % Normalise S afterwards (to avoid nan)
      S = S ./ Snorm;

      % Calculate orthogonal vectors (if present/needed)
      Ex = real(Exyz(:, ~mask));
      Ey = imag(Exyz(:, ~mask));
      Exnorm = vecnorm(Ex, 2, 1) ./ vecnorm(abs(Exyz(:, ~mask)), 2, 1);
      Eynorm = vecnorm(Ey, 2, 1) ./ vecnorm(abs(Exyz(:, ~mask)), 2, 1);
      Snorm = [Snorm, Snorm] .* [Exnorm, Eynorm];

      % Filter properties (and duplicate as required)
      xyz = [xyz(:, Exnorm ~= 0), xyz(:, Eynorm ~= 0)];
      S = [S(:, Exnorm ~= 0), S(:, Eynorm ~= 0)];
      Snorm = Snorm(:, [Exnorm ~= 0, Eynorm ~= 0]);
      pol = [Ex(:, Exnorm ~= 0), Ey(:, Eynorm ~= 0)];

      % Construct beam object
      % TODO: Should use use a 2xN field instead of a 1xN field?
      %       This was leftover from when it was in the Ray class
      beam = ott.beam.PlaneWave('origin', xyz, ...
        'direction', S, 'field', Snorm, 'polarisation', pol);
    end

    function beam = FromFarfield(oldbeam)
      % Constructs a new Ray instance from the beam far-field
      %
      % Uses the ``ehfarfield()`` methods from a :class:`ott.beam.Beam`
      % to sample a slice of the near-field.  The beam direction and intensity
      % is set from the cross product of the E and H fields.
      %
      % Each sample location produces zeros, one or two rays, corresponding
      % to the real and imaginary parts of E.  If either the real or
      % imaginary part is zero, only one ray is produced.
      %
      % Usage
      %   ray = FromFarfield(beam, ...)

      % Calculate fields on a grid
      [xx, yy, zz] = sphere(40);
      xyz = [xx(:), yy(:), zz(:)].';
      rtp = ott.utils.xyz2rtp(xyz);
      [E, H] = oldbeam.ehfarfield(rtp);

      % Generate from field vectors
      beam = ott.beam.PlaneWave.FromFieldVectors(E, H);
    end

    function beam = FromNearfield(oldbeam)
      % Constructs a new Ray instance from a near-field beam slice
      %
      % Uses the ``ehfield()`` methods from a :class:`ott.beam.Beam`
      % to sample a slice of the near-field.  The ray direction and intensity
      % is set from the cross product of the E and H fields.
      %
      % Each sample location produces zero, one or two rays, corresponding
      % to the real and imaginary parts of E.  If either the real or
      % imaginary part is zero, only one ray is produced.
      %
      % Usage
      %   ray = FromNearfield(beam, ...)

      % Calculate fields on a grid
      [xx, yy, zz] = meshgrid(linspace(-1, 1), linspace(-1, 1), 0);
      xyz = [xx(:), yy(:), zz(:)].';
      [E, H] = oldbeam.ehfield(xyz);

      % Generate from field vectors
      beam = ott.beam.PlaneWave.FromFieldVectors(E, H);
    end

    function beam = FromParaxial(oldbeam)
      % Constructs a new Ray instance from the beam paraxial far-field
      %
      % Uses the ``ehparaxial()`` methods from a :class:`ott.beam.Beam`
      % to sample a slice of the near-field.  The ray direction and intensity
      % is set from the cross product of the E and H fields.
      %
      % Each sample location produces zero, one or two rays, corresponding
      % to the real and imaginary parts of E.  If either the real or
      % imaginary part is zero, only one ray is produced.
      %
      % Usage
      %   ray = FromNearfield(beam, ...)

      % Calculate fields on a grid
      [xx, yy, zz] = meshgrid(linspace(-1, 1), linspace(-1, 1), 0);
      xyz = [xx(:), yy(:), zz(:)].';
      [E, H] = oldbeam.ehparaxial(xyz);

      % Generate from field vectors
      beam = ott.beam.PlaneWave.FromFieldVectors(E, H);
    end
  end

  methods
    function bm = PlaneWave(varargin)
      % Construct plane wave array properties
      %
      % Usage
      %   beam = PlaneWave(origin, direction, polarisation1, ...)
      %
      %   beam = PlaneWave(origin, 'directionSet', directionSet, ...)
      %
      % Parameters
      %   - origin (3xN numeric) -- Plane wave origins.
      %   - direction (3xN numeric) -- Plane wave direction.
      %   - polarisation1 (3xN numeric) -- Primary polarisation direction.
      %
      %   - directionSet (3x3N numeric) -- Array formed by combining
      %     direction/polarisation vectors into rotation matrices.  The
      %     direction vector should be the last column of the matrix.
      %     Default: ``eye(3)``.
      %
      % Optional named arguments
      %   - field (2xN numeric) -- Field parallel and perpendicular to
      %     plane wave polarisation direction.
      %     Default: ``[1, 1i]``.
      %
      %   - array_type (enum) -- Beam array type.  Can be
      %     'coherent', 'incoherent' or 'array'.  Default: ``'coherent'``.
      
      [origin, directionSet, field, unmatched] = ...
          ott.beam.properties.PlaneWaveArray.parseArgs(varargin{:});

      p = inputParser;
      p.addParameter('array_type', 'coherent');
      p.KeepUnmatched = true;
      p.parse(unmatched{:});
      unmatched = ott.utils.unmatchedArgs(p);

      bm = bm@ott.beam.properties.AnyArrayType(p.Results.array_type);
      bm = bm@ott.beam.properties.PlaneWaveArray(...
          origin, 'directionSet', directionSet, 'field', field, ...
          unmatched{:});
    end
  end

  methods (Hidden)
    function E = efieldInternal(beam, xyz, varargin)
      [E, ~] = beam.ehfieldInternal(xyz, varargin{:});
    end

    function H = hfieldInternal(beam, xyz, varargin)
      [~, H] = beam.ehfieldInternal(xyz, varargin{:});
    end

    function [E, H] = ehfieldInternal(beam, xyz, varargin)
      % Calculate the E-field
      %
      % Usage
      %   E = beam.efieldInternal(xyz, ...)
      %
      % Optional named arguments
      %   - method (enum) -- Method to use when drawing rays
      %     Can be one of
      %     - inf -- Infinite extend plane waves
      %     - ray_invr -- Power scaled with 1/distance from centre
      %     - invr -- Power scaled with 1/distance from origin
      %     Default: ``'inf'``.

      p = inputParser;
      p.addParameter('method', 'inf');
      p.parse(varargin{:});
      
      assert(ismatrix(xyz) && size(xyz, 1) == 3 && isnumeric(xyz), ...
        'xyz must be 3xN numeric matrix');
      
      E = zeros([3, numel(beam), size(xyz, 2)]);
      H = zeros([3, numel(beam), size(xyz, 2)]);
      
      % Get dependent properties from beam
      Edir = beam.polarisation1;
      Hdir = beam.polarisation2;
      wavenumber = beam.wavenumber;
      origin = beam.origin;
      direction = beam.direction;
      
      % Get method from properties outside loop (opt R2018a)
      method = p.Results.method;

      for ii = 1:size(xyz, 2)

        % Calculate distance along ray from origin
        rlocal = xyz(:, ii) - origin;
        dist = dot(direction, rlocal);

        % Calculate normalisation factor
        switch method
          case 'inf'
            scale = 1.0;

          case 'invr'
            scale = 1.0 ./ vecnorm(rlocal);

          case 'ray_invr'
            scale = 1.0 ./ vecnorm(rlocal - dist.*beam.direction);

          otherwise
            error('Unknown method specified');
        end

        % Calculate the polarisation (possibly in two directions)
        P0 = Edir .* beam.field(1, :) + Hdir .* beam.field(2, :);
        H0 = Hdir .* beam.field(1, :) + Edir .* beam.field(2, :);

        % Calculate the field at the location
        E(:, :, ii) = scale .* P0 .* exp(1i.*wavenumber.*dist);
        H(:, :, ii) = scale .* H0 .* exp(1i.*wavenumber.*dist);
      end
      
      if strcmpi(beam.array_type, 'coherent') || numel(beam) == 1
        E = squeeze(sum(E, 2));
        H = squeeze(sum(H, 2));
        
        % Package output
        E = ott.utils.FieldVector(xyz, E, 'cartesian');
        H = ott.utils.FieldVector(xyz, H, 'cartesian');
      else
        E = mat2cell(E, 3, ones(1, numel(beam)), size(xyz, 2));
        H = mat2cell(H, 3, ones(1, numel(beam)), size(xyz, 2));
        
        % Remove single dimension
        E = cellfun(@(x) squeeze(x), E, 'UniformOutput', false);
        H = cellfun(@(x) squeeze(x), H, 'UniformOutput', false);
        
        % Package output
        E = cellfun(@(x) ott.utils.FieldVector(xyz, x, 'cartesian'), ...
            E, 'UniformOutput', false);
        H = cellfun(@(x) ott.utils.FieldVector(xyz, x, 'cartesian'), ...
            H, 'UniformOutput', false);
      end
    end

    function E = efarfieldInternal(beam, rtp, varargin)
      [E, ~] = beam.ehfarfieldInternal(rtp, varargin{:});
    end

    function H = hfarfieldInternal(beam, rtp, varargin)
      [~, H] = beam.ehfarfieldInternal(rtp, varargin{:});
    end
    
    function force = forceInternal(beam, incident, varargin)
      % Calculate the from the change in poynting vectors
      
      [E1, H1] = beam.ehfieldInternal([0;0;0]);
      S1 = beam.arrayApply(@(E, H) real(cross(E.vxyz, conj(H.vxyz))), E1, H1);
      S1 = beam.combineIncoherentArray(S1, 2);
      if iscell(S1)
        S1 = cell2mat(S1);
      end
      
      [E2, H2] = incident.ehfieldInternal([0;0;0]);
      S2 = incident.arrayApply(@(E, H) real(cross(E.vxyz, conj(H.vxyz))), E2, H2);
      S2 = incident.combineIncoherentArray(S2, 2);
      if iscell(S2)
        S2 = cell2mat(S2);
      end

      force = S2 - S1;
    end

    function [E, H] = ehfarfieldInternal(beam, rtp, varargin)
      % Calculate the E-field in the far-field
      %
      % Usage
      %   E = efarfieldInternal(beam, rtp, ...)
      %
      % Optional named arguments
      %   - method (enum) -- Method to use when drawing rays
      %     Can be one of
      %     - delta -- Only draw the ray if rtp is a ray direction
      %     - dot   -- Draw dot(rtp, ray direction)
      %     Default: ``''dot'``.

      p = inputParser;
      p.addParameter('method', 'dot');
      p.parse(varargin{:});

      % It's easier to work in Cartesian coordinates
      rtp(1, :) = 1.0;
      xyz = ott.utils.rtp2xyz(rtp);

      Edir = beam.polarisation;
      Hdir = cross(beam.direction, beam.polarisation);

      % Calculate the polarisation (possibly in two directions)
      E0 = Edir .* beam.field(1, :);
      H0 = Hdir .* beam.field(1, :);
      if size(beam.field, 1) == 2
        E0 = E0 + Hdir .* beam.field(2, :);
        H0 = H0 + Edir .* beam.field(2, :);
      end
      
      E = zeros(size(xyz));
      H = zeros(size(xyz));

      for ii = 1:size(xyz, 2)

        % Calculate component of ray in specified direction
        intensities = sum(beam.direction .* xyz(:, ii), 1)./vecnorm(beam.direction);

        switch p.Results.method
          case 'dot'
            % Nothing to do
          case 'delta'
            tol = 1.0e-6;
            intensities = double(abs(1 - intensities) < tol);
          otherwise
            error('Unknown method specified');
        end

        % Calculate the field at the location
        E(:, ii) = sum(intensities .* E0, 2);
        H(:, ii) = sum(intensities .* H0, 2);
      end

      % Return result in Cartesian coordinates
      E = ott.utils.FieldVector(xyz, E, 'cartesian');
      H = ott.utils.FieldVector(xyz, H, 'cartesian');
    end
  end
end
