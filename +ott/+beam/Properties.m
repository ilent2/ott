classdef (Abstract) Properties
% A base class for Beam and abstract.Beam representations.
%
% This class defines the common properties and methods to these
% two classes.
%
% Properties
%   - power         -- The power of the beam (may be infinite)
%   - wavelength    -- Wavelength of beam in medium (default: 1.0)
%   - permittivity  -- Material relative permittivity (default: 1.0)
%   - permeability  -- Material relative permeability (default: 1.0)
%   - combine_rules -- Combination rules for beam array dimensions
%
% Dependent properties
%   - wavenumber    -- Wave-number of beam in medium
%   - impedance     -- Impedance of the medium
%   - coherent_dims -- Coherent dimensions of beam array (logical mask)
%   - incoherent_dims -- Incoherent dimensions of beam array (logical mask)
%
% Methods
%   - size          -- Size of the beam array
%   - cat           -- Concatenation method for beam arrays
%   - combineArray  -- Combine an array using combination rules
%   - combinedSize  -- Size of the beam after applying combination rules
%
% Abstract methods
%   - getBeamPower      -- get method called by dependent property power

% TODO: Some inconsistencies with relative and absolute permittivity

  properties
    wavelength     % Wavelength of beam in medium (default: 1.0)
    permittivity   % Material relative permittivity (default: 1.0)
    permeability   % Material relative permeability (default: 1.0)
    speed0         % Speed of light in vacuum (default: 1.0)
    combine_rules  % Combination rules for beam array dimensions (cell array)
  end

  properties (Dependent)
    power           % The power of the beam (may be infinite)
    wavenumber      % Wave-number of beam in medium
    impedance       % Impedance of the medium
    omega          % Optical frequency of light
    index_medium   % Refractive index in medium
    speed          % Speed of light in medium
    wavelength0    % Vacuum wavelength of beam

    coherent_dims   % Coherent dimensions of beam array (logical mask)
    incoherent_dims % Incoherent dimensions of beam array (logical mask)
  end

  methods (Abstract)
    getBeamPower    % get method called by dependent property power
  end

  methods (Static)
    function [permittivity, wavelength, speed0] = parseInputs(varargin)
      % Helper to parse the inputs to the constructor

      % Pre-parse (for like argument)
      pre = inputParser;
      pre.KeepUnmatched = true;
      pre.addParameter('like', []);
      pre.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(pre);

      % Get defaults from like (if available)
      def_permittivity = 1.0;
      def_wavelength = 1.0;
      def_speed0 = 1.0;
      if ~isempty(pre.Results.like)
        def_permittivity = pre.Results.like.permittivity;
        def_wavelength = pre.Results.like.wavelength;
        def_speed0 = pre.Results.like.speed0;
      end

      % All values should be numeric scalars
      type_check = @(x) isnumeric(x) & isscalar(x);

      % Setup input arguments
      p = ott.utils.RelatedArgumentParser;
      p.addRequired('permittivity', def_permittivity, type_check);
      p.addRequired('wavelength', def_wavelength, type_check);
      p.addRequired('speed0', def_speed0, type_check);
      p.addOptional('omega', [], type_check);
      p.addOptional('index_medium', [], type_check);
      p.addOptional('wavenumber', [], type_check);
      p.addOptional('speed', [], type_check);
      p.addOptional('wavelength0', [], type_check);

      % Add rules relating required to optional
      p.addRule('speed = speed0 ./ index_medium');
      p.addRule('wavenumber = 2*pi / wavelength');
      p.addRule('permittivity = index_medium.^2');
      p.addRule('speed0 = omega / (2*pi) * wavelength0');
      p.addRule('speed = omega / (2*pi) * wavelength');

      % Parse inputs
      p.parse(unmatched{:});

      % Assign variables to output
      speed0 = p.RequiredResults.speed0;
      permittivity = p.RequiredResults.permittivity;
      wavelength = p.RequiredResults.wavelength;
    end
  end

  methods
    function beam = Properties(varargin)
      % Initialize properties to defaults
      %
      % Optional named arguments
      %   - permittivity (numeric) -- Relative permittivity of medium
      %   - wavelength (numeric) -- Wavelength in medium [L]
      %   - speed0 (numeric) -- Speed of light in vacuum [L/T]
      %
      %   - omega (numeric) -- Optical frequency [2*pi/T]
      %   - index_medium (numeric) -- Refractive index in medium
      %   - wavenumber (numeric) -- Wave-number in medium [2*pi/L]
      %   - speed (numeric) -- Speed of light in medium [L/T]
      %   - wavelength0 (numeric) -- Wavelength in medium [L]
      %
      %   - like (Beam) -- Uses this beam for default properties.

      beam.wavelength = 1.0;
      beam.permittivity = 1.0;
      beam.permeability = 1.0;

      % Parse optional inputs
      [beam.permittivity, beam.wavelength, beam.speed0] = ...
          beam.parseInputs(varargin{:});
    end

    function arr = combineArray(beam, arr, base_sz, combine_type)
      % Apply combination rules to a arbitrary array.
      %
      % Usage
      %   arr = beam.combineArray(arr, base_sz, combine_type)
      %
      % Parameters
      %   arr -- The array to apply the rules too.
      %   The size of arr (excluding leading dimensions) must match the
      %   size of the beam array or the beam array after applying
      %   the incoherent or coherent rules first.
      %
      %   base_sz (numeric) -- The number of leading dimensions to
      %   ignore in arr.
      %
      %   combine_type (enum) -- Combination rules to apply.
      %   Specify the combination type, can either be 'both', 'coherent'
      %   or 'incoherent'.  Default: 'both'.

      if nargin < 4
        combine_type = 'both';
      end

      assert(isnumeric(base_sz) && isscalar(base_sz) && base_sz >= 0 ...
          && round(base_sz) == base_sz, ...
          'base_sz must be positive integer');

      % Get array of dimensions to combine
      switch combine_type
        case 'both'
          mask = beam.coherent_dims | beam.incoherent_dims;
          assert(base_sz + ndims(arr) == numel(mask), ...
              'Array has too few dimensions');
        case 'coherent'
          mask = beam.coherent_dims;

          if base_sz + ndims(arr) < numel(mask)
            % Try discarding incoherent dimensions
            mask(beam.incoherent_dims) = [];
            assert(base_sz + ndims(arr) == numel(mask), ...
                'Array has too few dimensions');
          end

        case 'incoherent'
          mask = beam.incoherent_dims;

          if base_sz + ndims(arr) < numel(mask)
            % Try discarding coherent dimensions
            mask(beam.coherent_dims) = [];
            assert(base_sz + ndims(arr) == numel(mask), ...
                'Array has too few dimensions');
          end

        otherwise
          error('Unknown combine type specified');
      end
      sum_dims = base_sz + find(mask);

      % Sum dimensions in reverse order
      for dim = flip(sum_dims)
        arr = sum(arr, dim);
      end
    end

    function sz = combinedSize(beam, combine_type, dim)
      % Calculate the size of the beam after applying combination rules
      %
      % Usage
      %   sz = beam.combinedSize()
      %
      %   sz = beam.combinedSize(combine_type)
      %   Specify the combination type, can either be 'both', 'coherent'
      %   or 'incoherent'.  Default: 'both'.
      %
      %   sz = beam.combineSize(combine_type, dim)
      %   Only return the size of the specified dimension.
      %
      % The minimum size is [1, 1].

      % Get the un-combined size
      sz = size(beam);

      if nargin < 2
        combine_type = 'both';
      end

      switch combine_type
        case 'both'
          mask = beam.coherent_dims | beam.incoherent_dims;
        case 'coherent'
          mask = beam.coherent_dims;
        case 'incoherent'
          mask = beam.incoherent_dims;
        otherwise
          error('Unknown combine type specified');
      end

      % Remove combined dimensions
      sz(mask) = [];
      if numel(sz) == 0
        sz = [1, 1];
      elseif numel(sz) == 1
        sz = [1, sz];
      end

      % Get specified dimension
      if nargin == 3
        sz = sz(dim);
      end
    end
  end

  methods (Hidden)
    function beam = setBeamPower(beam, val)
      % Function to set the beam power (if supported)
      % Override this function if your beam supports this feature
      error('Setting beam power not supported');
    end
  end

  methods % Getters/setters
    function val = get.power(beam)
      val = beam.getBeamPower();
    end
    function beam = set.power(beam, val)
      beam = beam.setBeamPower(val);
    end

    function beam = set.wavelength(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'wavelength must be numeric scalar');
      beam.wavelength = val;
    end

    function beam = set.wavenumber(beam, val)
      % Set the wavelength
      assert(isnumeric(val), 'wavenumber must be numeric');
      beam.wavelength = 2*pi./val;
    end
    function val = get.wavenumber(beam)
      val = 2*pi/beam.wavelength;
    end

    function beam = set.permittivity(beam, val)
      assert(isscalar(val) && isnumeric(val), ...
          'permittivity must be numeric scalar');
      beam.permittivity = val;
    end

    function beam = set.permeability(beam, val)
      assert(isscalar(val) && isnumeric(val), ...
          'permeability must be numeric scalar');
      beam.permeability = val;
    end

    function val = get.impedance(beam)
      % TODO: What about tensor properties
      val = sqrt(beam.permeability ./ beam.permittivity);
    end

    function beam = set.speed0(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'speed0 must be numeric scalar');
      beam.speed0 = val;
    end

    function val = get.omega(beam)
      wavelength0 = beam.wavelength .* beam.index_medium;
      val = 2*pi*beam.speed0./wavelength0;
    end
    function val = get.index_medium(beam)
      val = sqrt(beam.permittivity);
    end
    function val = get.speed(beam)
      val = beam.speed0 ./ beam.index_medium;
    end
    function val = get.wavelength0(beam)
      val = beam.wavelength ./ beam.index_medium;
    end

    function val = get.coherent_dims(beam)
      val = strcmpi('coherent', beam.combine_rules);
    end
    function val = get.incoherent_dims(beam)
      val = strcmpi('incoherent', beam.combine_rules);
    end
  end
end
