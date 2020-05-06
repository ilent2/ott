classdef Array < ott.beam.abstract.Beam & ott.beam.utils.ArrayType
% Abstract array of beams.
% Inherits from :class:`ott.beam.utils.ArrayType`.
%
% Properties
%   - beams       -- Internal array of beam objects
%   - array_type  -- Type of array ('coherent', 'array' or 'incoherent')
%
% Methods
%   - size          -- Size of the beam array
%   - plus          -- Provides addition of coherent beams
%   - cat           -- Concatenation of beams and arrays
%   - vertcat       -- Vertical concatenation of beams and arrays
%   - horzcat       -- Horizontal concatenation of beams and arrays
%   - subsref       -- For direct indexing of the beams array
%   - combineIncoherentArray  -- Combine cell array of beam data

% TODO: Should arrays have beam properties?

  properties
    beams         % Beams contained in the array
  end

  properties (Dependent)
    contains_incoherent     % True if the beam contains incoherent beams
  end

  methods (Hidden)
    function beam = plusInternal(b1, b2)
      % Add beams coherently

      isArr1 = isa(b1, 'ott.beam.abstract.Array');
      isArr2 = isa(b2, 'ott.beam.abstract.Array');

      if isArr1 && isArr2
        b1.beams = [b1.beams, b2.beams];
        beam = b1;
      elseif isArr1
        b1.beams = [b1.beams, {b2}];
        beam = b1;
      elseif isArr2
        b2.beams = [{b1}, b2.beams];
        beam = b2;
      end
    end

    function beam = catInternal(dim, beam, varargin)
      % Concatenate arrays

      other_beams = {};
      for ii = 1:length(varargin)
        other_beams{ii} = varargin{ii}.beams;
      end

      beam.beams = cat(dim, beam.beams, other_beams{:});
    end

    function beams = subsrefInternal(beam, subs)
      % Get the subscripted beam

      beams = beam.beams{subs{:}};
    end
  end

  methods
    function beam = Array(array_type, sz, varargin)
      % Construct a new abstract beam array
      %
      % Usage
      %   beam = Array(array_type, dim, beam1, beam2, ...)
      %
      % Parameters
      %   - array_type (enum) -- Type of beam array.  Either
      %     'array', 'coherent' or 'incoherent'.
      %
      %   - sz (N numeric) -- Size of the beam array.
      %
      %   - beam1, beam2, ... -- Beams to include in the array.
      %     If empty, creates an empty cell array internally.

      beam = beam@ott.beam.utils.ArrayType('array_type', array_type);

      % Store beam array and reshape
      if length(varargin) == 0
        beam.beams = {};
      else
        beam.beams = reshape(varargin, sz);
      end
    end

    function sz = size(beam, varargin)
      % Get the size of the beam array
      %
      % Usage
      %   sz = size(beam, ...)  or beam.size(...)
      %   For arguments, see help on Matlab's built-in ``size``.

      sz = size(beam.beams, varargin{:});
    end

    function arr = combineIncoherentArray(beam, arr, dim)
      % Combine incoherent layers as required

      % Check if we have work to do
      if strcmpi(beam.array_type, 'coherent')
        return;
      end

      assert(numel(arr) == numel(beam), ...
          'Number of elements in beam must match number of elements in array');

      % Combine child beams
      for ii = 1:numel(beam)
        arr{ii} = beam.beams{ii}.combineIncoherentArray(arr{ii}, dim);
      end
      
      % Combine this array
      arr = combineIncoherentArray@ott.beam.utils.ArrayType(beam, arr, dim);
    end
    
    function data = arrayApply(beam, data, func)
      % Apply function to each array in the beam array output.
      %
      % This function is overloaded by Array types in order to
      % implement incoherent combination.
      
      % Apply visualisatio funtion to sub-beams
      for ii = 1:numel(beam)
        data{ii} = beam.beams{ii}.arrayApply(data{ii}, func);
      end
    end
  end

  methods % Getters/setters
    function val = get.contains_incoherent(beam)

      val = false;

      for ii = 1:numel(beam)
        if isa(beam.beams{ii}, 'ott.beam.utils.ArrayType')
          val = val | strcmpi(beam.beams{ii}, 'incoherent');
          if isa(beam.beams{ii}, 'ott.beam.abstract.Array')
            val = val | beam.beams{ii}.contains_incoherent;
          end
        end
      end
    end

    function beam = set.beams(beam, val)

      assert(iscell(val), 'beams property be cell array');

      % Cann't have coherent arrays of incoherent beams
      if strcmpi(beam.array_type, 'coherent')
        for ii = 1:numel(val)
          if isa(val{ii}, 'ott.beam.utils.ArrayType')

            assert(~strcmpi(val{ii}.array_type, 'incoherent'), ...
              'ott:beam:utils:ArrayType:coherent_with_incoherent', ...
              'Cannot have coherent array of incoherent beams');

            if isa(val{ii}, 'ott.beam.abstract.Array')
              assert(~val{ii}.contains_incoherent, ...
                'ott:beam:utils:ArrayType:coherent_with_incoherent', ...
                'Cannot have coherent array of incoherent beams');
            end
          end
        end
      end

      beam.beams = val;
    end
  end
end
