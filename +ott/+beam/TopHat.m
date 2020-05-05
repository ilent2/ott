classdef TopHat < ott.beam.Beam & ott.beam.paraxial.Paraxial ...
    & ott.utils.Vector
% Representation of a diffraction-free beam with a Top-Hat profile.
% Inherits from :class:`Beam`, :class:`Paraxial`
% and :class:`ott.utils.Vector`.
%
% A top-hat beam has a sharp edge intensity profile, for example, in
% polar coordinates::
%
%   E(r, z) = \left\{\begin{array}{ll} 1    & r \leq 1 \\
%                                   0    & \text{otherwise}\end{array}\right.
%
% This class has two properties which control the shape: ``profile`` and
% ``field``.  Profile is a function handle which returns logical
% values for inside/outside the top-hat; Field can either be a scalar
% or a function handle for the polarisation distribution across the top-hat.
%
% This function is similar to :class:`PlaneWave` and can be cast to
% :class:`PlaneWave` when the field profile is uniform across the beam.
%
% Static methods
%   - ProfileSquare     -- Square top-hat profile
%   - ProfileRectangle  -- Rectangle top-hat profile
%   - ProfileCircle     -- Circle top-hat profile
%   - ProfileEllipse    -- Ellipse top-hat profile
%
% Properties
%   - profile       -- Function handle for top-hat profile
%   - field         -- Field distribution function, scalar or 2-vector
%   - origin        -- Ray origins, 3xN array (default [0;0;0])
%   - direction     -- Direction of propagation (3xN Cartesian)
%   - polarisation  -- Primary polarisation direction
%   - powerInternal -- Internal power value (numeric | empty)
%
% Methods
%   - rotate      -- Rotate the direction and polarisation
%   - rotate*     -- Rotate the particle around the X,Y,Z axis
%
% Field visualisation methods
%   - visualise -- Generate a visualisation around the origin
%   - visualiseFarfield -- Generate a visualisation at the far-field
%   - visualiseFarfieldSlice -- Visualise the field on a angular slice
%   - visualiseFarfieldSphere -- Visualise the filed on a sphere
%
% Properties inherited from Beam
%   - power         -- The power of the beam (infinite)
%   - wavelength    -- Wavelength of beam in medium (default: 1.0)
%   - permittivity  -- Material relative permittivity (default: 1.0)
%   - permeability  -- Material relative permeability (default: 1.0)
%   - wavenumber    -- Wave-number of beam in medium
%   - impedance     -- Impedance of the medium

  properties
    field           % Field parallel and perpendicular to polarisation
    polarisation    % Primary polarisation direction
    profile         % Function handle for beam profile
  end

  properties (Hidden)
    powerInternal   % Internal power value (numeric | empty)
  end

  methods (Static)
    function fcn = ProfileSquare(width, direction)
      % Creates a square top-hat profile function handle
      %
      % Usage
      %   function_handle = ProfileSquare(width, direction)
      %
      % Parameters
      %   - width (numeric) -- Side length of the square
      %
      %   - direction (2x1 numeric) -- Orientation of the square
      %     with respect to the beam primary polarisation direction.
      %     This vector should point from the centre to the middle of
      %     one of the square's side.
      %     It is recommended to use ``[1; 0]`` and change the beam
      %     polarisation if you are using automatic power estimation.

      assert(isnumeric(width) && isscalar(width), ...
          'width must be numeric scalar');
      assert(isnumeric(direction) && all(size(direction) == [2, 1]), ...
          'direction must be numeric 2x1 vector');

      % Ensure direction is normalised
      direction = direction ./ vecnorm(direction);

      dir1 = direction;
      dir2 = [-dir1(2); dir1(1)];

      % Expression that can be converted to a symbol
      fcn = @(x, y) (dir1.' * [x;y]).^2 <= (width./2).^2 ...
          & (dir2.' * [x;y]).^2 <= (width./2).^2;
    end

    function fcn = ProfileCircle(radius)
      % Creates a circle top-hat profile function handle
      %
      % Usage
      %   function_handle = ProfileCircle(radius)

      assert(isnumeric(radius) && isscalar(radius), ...
          'radius must be numeric scalar');

      % Expression that can be converted to a symbol
      fcn = @(x,y) x.^2 + y.^2 <= radius.^2;
    end

    function fcn = ProfileEllipse(radii, direction)
      % Creates a elliptical top-hat profile function handle
      %
      % Usage
      %   function_handle = ProfileEllipse(radii, direction)
      %
      % Parameters
      %   - widths (2-numeric) -- Width of the rectangle sides.
      %
      %   - direction (2x1 numeric) -- Orientation of the square
      %     with respect to the beam primary polarisation direction.
      %     This vector points in the direction with length ``widths(1)``.
      %     It is recommended to use ``[1; 0]`` and change the beam
      %     polarisation if you are using automatic power estimation.

      assert(isnumeric(radii) && numel(radii) == 2, ...
          'radii must be 2-numeric');
      assert(isnumeric(direction) && all(size(direction) == [2, 1]), ...
          'direction must be numeric 2x1 vector');

      % Ensure direction is normalised
      direction = direction ./ vecnorm(direction);

      rot = [direction, flip(direction).*[-1;1]];

      % Expression that can be converted to a symbol
      fcn = @(x, y) sum(((rot * [x;y])./radii(:)).^2, 1) <= 1.0;
    end

    function fcn = ProfileRectangle(widths, direction)
      % Creates a rectangle top-hat profile function handle
      %
      % Usage
      %   function_handle = ProfileRectangle(widths, direction)
      %
      % Parameters
      %   - widths (2-numeric) -- Width of the rectangle sides.
      %
      %   - direction (2x1 numeric) -- Orientation of the square
      %     with respect to the beam primary polarisation direction.
      %     This vector points in the direction with length ``widths(1)``.
      %     It is recommended to use ``[1; 0]`` and change the beam
      %     polarisation if you are using automatic power estimation.

      assert(isnumeric(widths) && numel(widths) == 2, ...
          'width must be 2-numeric');
      assert(isnumeric(direction) && all(size(direction) == [2, 1]), ...
          'direction must be numeric 2x1 vector');

      % Ensure direction is normalised
      direction = direction ./ vecnorm(direction);

      dir1 = direction;
      dir2 = [-dir1(2); dir1(1)];

      % Expression that can be converted to a symbol
      fcn = @(x, y) (dir2.' * [x;y]).^2 <= (widths(1)./2).^2 ...
          & (dir1.' * [x;y]).^2 <= (widths(2)./2).^2;
    end
  end

  methods
    function bm = TopHat(varargin)
      % Construct a new top hat beam or top hat beam array
      %
      % Usage
      %   beam = TopHat(...)
      %
      % Optional named arguments
      %   - profile (function_handle) -- Describes the shape of the top-hat.
      %     See also :meth:`ProfileCircle` and similar functions.
      %     Default: ``ProfileCircle(1.0)``.
      %
      %   - direction (3xN numeric) -- direction vectors (Cartesian)
      %     Default: ``[0;0;1]``.
      %
      %   - polarisation (3xN numeric) -- polarisation vectors (Cartesian)
      %     Default: ``[1;0;0]``.
      %
      %   - field (1xN|2xN numeric|function_handle)) -- Function handle
      %     or field vectors parallel and (optionally) perpendicular to
      %     the polarisation direction.
      %     Allows for 0 intensity with finite polarisation direction.
      %     Default: ``1``.
      %
      %   - origin (3xN numeric) -- Origin of plane waves.
      %     Default: ``[0;0;0]``.
      %
      %   - power (numeric | enum) -- Power of the beam.  Must be either
      %     a numeric value or 'symbolic' to use the symbolic toolbox
      %     or 'numeric' to use ``integral2`` to estimate the power
      %     numerically.
      %     Default: ``'numeric'``.

      % Parse inputs
      p = inputParser;
      p.addParameter('profile', ott.beam.TopHat.ProfileCircle(1.0));
      p.addParameter('direction', [0;0;1]);
      p.addParameter('polarisation', [1;0;0]);
      p.addParameter('field', 1.0);
      p.addParameter('origin', [0;0;0]);
      p.addParameter('power', 'numeric');
      p.parse(varargin{:});

      % Call Vector base class
      % Using 'bm' instead of 'beam', possible bug in R2018a.
      bm = bm@ott.utils.Vector(p.Results.origin, p.Results.direction);

      % Store remaining parameters
      bm.polarisation = p.Results.polarisation;
      bm.field = p.Results.field;
      bm.profile = p.Results.profile;
      bm.powerInternal = p.Results.power;
    end

    function beam = ott.beam.PlaneWave(beam)
      % Cast the object to a PlaneWave instance
      %
      % This requires the polarisation profile to be uniform.

      assert(isnumeric(beam.field), 'Field profile must be numeric');

      beam = ott.beam.PlaneWave(...
          'direction', beam.direction, ...
          'polarisation', beam.polarisation, ...
          'origin', beam.origin, ...
          'field', beam.field);
    end
  end

  methods (Hidden)

    function p = getBeamPower(beam)
      % Gets the internal beam power or attempts to calculate it

      if ischar(beam.powerInternal)
        % Attempt to calculate power

        if strcmpi(beam.powerInternal, 'symbolic')
          % Use the symbolic toolbox and the profile function
          func = piecewise(sym(beam.profile), 1, 0);

          % Add field intensity profile
          if isa(beam.field, 'function_handle')
            func = abs(beam.field(x, y)) .* func;
          else
            func = abs(beam.field) .* func;
          end

          p = double(int(int(func, -Inf, Inf), -Inf, Inf));

        else
          % Evaluate numerically

          % Get function and add field
          if isa(beam.field, 'function_handle')
            func = @(x, y) abs(beam.field(x, y)) .* beam.profile(x, y);
          else
            func = @(x, y) double(beam.profile(x, y)) .* abs(beam.field);
          end

          p = integral2(func, -Inf, Inf, -Inf, Inf);

        end

      else
        % Get the internally stored result
        p = beam.powerInternal;
      end
    end

    function E = efieldInternal(beam, xyz)
      % Calculate the electric field of the Top-hat beam

      E = zeros(size(xyz));
      
      for ii = 1:size(xyz, 2)

        % Calculate distance along ray from origin
        rlocal = xyz(:, ii) - beam.origin;
        dist = dot(beam, rlocal);

        % Calculate xy coordinates for profile evaluation
        xydir = rlocal - dist.*beam.direction;
        xdir = dot(xydir, beam.polarisation);
        ydir = vecnorm(xydir - xdir);

        % Evaluate beam profile
        profile = double(beam.profile(xdir, ydir));

        % Add field components
        if isa(beam.field, 'function_handle')
          profile = profile .* beam.field(xdir, ydir);
        else
          profile = profile .* beam.field;
        end

        % Calculate the field at the location
        E(:, ii) = profile .* exp(1i.*beam.wavenumber.*dist);
      end

      % Package output
      E = ott.utils.FieldVector(xyz, E, 'cartesian');
    end
  end

  methods % Getters/setters
    % Properties
    %   - field         -- Field parallel and perpendicular to polarisation
    %   - polarisation  -- Primary polarisation direction
    %   - profile       -- Function handle for profile

    function beam = set.profile(beam, val)

      % Check type
      assert(isa(val, 'function_handle'), ...
          'profile must be a function handle');

      % Check length
      assert(size(val, 2) == size(beam.direction, 2), ...
          'polarisation must have same length as direction');

      beam.profile = val;
    end

    function beam = set.field(beam, val)

      % Check type and rows
      assert(isa(val, 'function_handle') ...
        || (isnumeric(val) && any(size(val, 1) == [1, 2])), ...
        'field must be 1xN, 2xN numeric or function handle');

      % Check length
      assert(size(val, 2) == size(beam.direction, 2), ...
          'polarisation must have same length as direction');

      beam.field = val;
    end

    function beam = set.polarisation(beam, val)

      % Check type and rows
      assert(isnumeric(val) && size(val, 1) == 3, ...
          'polarisation must be numeric 3xN matrix');

      % Check length
      assert(size(val, 2) == size(beam.direction, 2), ...
          'polarisation must have same length as direction');

      beam.polarisation = val;
    end

    function beam = set.powerInternal(beam, val)
      assert((ischar(val) && any(strcmpi(val, {'numeric', 'symbolic'}))) ...
          || (isnumeric(val) && isscalar(val)), ...
          'Power must be numeric value or ''numeric'' or ''symbolic''');
      beam.powerInternal = val;
    end
  end
end

