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
%
% Dependent properties
%   - ndipoles      -- Number of dipoles in the array N
%   - nbeams        -- Number of beams in the array M

  properties (SetAccess=protected)
    location          % Dipole locations
    polarization      % Dipole polarization
  end

  properties (Dependent)
    nbeams
    ndipoles
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

      sz = [beam.ndipoles, beam.nbeams];

      [varargout{1:nargout}] = ott.utils.size_helper(sz, varargin{:});
    end
  end

  methods (Hidden)
    function beam = catInternal(dim, beam, varargin)
      % Concatenate beams

      assert(dim == 1 || dim == 2, 'Dimension must be 1 or 2');

      other_loc = {};
      other_pol = {};
      for ii = 1:length(varargin)
        other_loc{ii} = varargin{ii}.location;
        other_pol{ii} = varargin{ii}.polarization;
      end

      if dim == 1
        % Add new dipoles
        beam.location = cat(2, beam.location, other_loc{:});
        beam.polarization = cat(1, beam.polarization, other_pol{:});

      else
        % Check dipole locations are equal
        locEqual = cellfun(@(x) all(x(:) == beam.location(:)), other_loc);
        assert(locEqual, 'All locations must match');

        % Add new polarizations
        beam.polarization = cat(2, beam.polarization, other_pol{:});
      end
    end

    function beam = plusInternal(beam1, beam2)
      % Concatenate two coherent beams together

      beam = cat(2, beam1, beam2);
    end

    function beam = subsrefInternal(beam, subs)
      % Get the subscripted beam

      assert(numel(subs) <= 2, 'Maximum 2 subscripts');

      if numel(subs) == 1
        if beam.nbeams == 1
          % Index dipoles
          beam.location = beam.location(:, subs{:});
          idx = (1:3).' + 3*(subs{1}-1);
          beam.polarization = beam.polarization(idx(:), :);
        else
          % Warn when ambiguous
          if beam.ndipoles ~= 1
            warning('ott:beam:properties:DipoleArray:subsref_single_index', ...
                'Index is ambiguous, indexing beam polarizations');
          end

          % Index beams
          beam.polarization = beam.polarization(:, subs{1});
        end
      else
        % Index beams and dipoles
        beam.location = beam.location(:, subs{1});
        idx = (1:3).' + 3*(subs{1}-1);
        beam.polarization = beam.polarization(idx(:), subs{2});
      end
    end

    function beam = subsasgnInternal(beam, subs, rem, other)
      % Assign to the subscripted beam

      assert(isempty(rem), 'Assignment to parts of beams not supported');
      assert(numel(subs) <= 2, 'Maximum 2 subscripts');

      % TODO: Update this for multiple subscripts
      error('Not yet implemented');

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

    function n = get.ndipoles(beam)
      n = size(beam.location, 2);
    end

    function n = get.nbeams(beam)
      n = size(beam.polarization, 2);
    end
  end
end

