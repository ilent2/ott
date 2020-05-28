classdef DipoleArray < ott.beam.properties.Dipole
% Base class for dipole array objects
%
% Properties
%   - position      -- Position of the dipole collection
%   - rotation      -- Orientation of the dipole collection
%
% Methods
%   - setDipoles    -- Set dipole positions/polarizations
%   - size          -- Size of the dipole array
%
% Properties
%   - location      -- Locations of individual dipoles (3xN)
%   - polarization  -- Polarization of the dipoles (3NxM)

  properties (SetAccess=protected)
    location          % Dipole locations
    polarization      % Dipole polarization
  end

  methods
    function varargout = size(beam, varargin)
      % Get the number of beams contained in this object
      %
      % Usage
      %   sz = size(beam)   or    sz = beam.size()
      %   For help on arguments, see builtin ``size``.
      %
      % The leading dimension is always 1.  May change in future.

      sz = size(beam.location);
      sz(1) = 1;

      [varargout{1:nargout}] = ott.utils.size_helper(sz, varargin{:});
    end
  end

  methods (Hidden)
    function beam = catInternal(dim, beam, varargin)
      % Concatenate beams

      assert(dim == 2, 'Only 1xN arrays supported (may change in future)');

      other_loc = {};
      other_pol = {};
      for ii = 1:length(varargin)
        other_loc{ii} = varargin{ii}.location;
        other_pol{ii} = varargin{ii}.polarization;
      end

      beam.location = cat(dim, beam.location, other_loc{:});
      beam.polarization = cat(dim, beam.polarization, other_pol{:});
    end

    function beam = plusInternal(beam1, beam2)
      % Concatenate two coherent beams together

      beam = cat(2, beam1, beam2);
    end

    function beam = subsrefInternal(beam, subs)
      % Get the subscripted beam

      if numel(subs) > ndims(beam.location)
        if subs(1) == 1
          subs = subs(2:end);
        end
        assert(numel(subs) > ndims(beam.location), ...
            'Too many subscript indices');
      end

      beam.location = beam.location(:, subs{:});

      idx = (1:3).' + 3*(subs{1}-1);
      beam.polarization = beam.polarization(idx(:), :);
    end

    function beam = subsasgnInternal(beam, subs, rem, other)
      % Assign to the subscripted beam

      if numel(subs) > ndims(beam.location)
        if subs(1) == 1
          subs = subs(2:end);
        end
        assert(numel(subs) > ndims(beam.location), ...
            'Too many subscript indices');
      end

      assert(isempty(rem), 'Assignment to parts of beams not supported');

      idx = (1:3).' + 3*(subs{1}-1);

      if isempty(other)
        % Delete data
        beam.location(:, subs{:}) = other;
        beam.polarization(idx(:), :) = other;

      else
        % Ensure we have a dipole
        assert(isa(other, 'ott.beam.properties.DipoleArray'), ...
            'Only DipoleArray beams supported for now');

        % Must set data first!
        beam.location(:, subs{:}) = other.location;
        beam.polarization(:, subs{:}) = other.polarization;
      end
    end
  end

  methods % Getters/setters
    function beam = set.polarization(beam, val)
      assert(isnumeric(val) && ismatrix(val) && mod(size(val, 1), 3) == 0, ...
          'polarization must be a 3NxM numeric matrix');
      beam.polarization = val;
    end

    function beam = set.location(beam, val)
      assert(isnumeric(val) && ismatrix(val) && size(val, 1) == 3, ...
          'location must be 3xN numeric matrix');
      beam.location = val;
    end
  end
end

