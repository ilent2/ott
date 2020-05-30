classdef Array < ott.beam.Beam ...
    & ott.beam.properties.AnyArrayType
% A class representing arrays of beams.
% Inherits from :class:`Beam` and :class:`ott.beam.utils.AnyArrayType`.
%
% This class declares the methods needed for :class:`Beam`.
% Output of the :class:`Beam` methods is packaged into arrays of
% the appropriate size.
%
% The :class:`Coherent` and :class:`Incoherent` classes are pseudonyms
% for this class class.  There is no generic array pseudonyms: i.e.,
% casting directly to this class creates a generic array by default.
% The other consequence is the `isa` method returns true for coherent
% and incoherent arrays when applied to instance of :class:`Array`.
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
%
% Static methods
%   - CombineCoherent   -- Combine coherent data from cell arrays
%   - empty         -- Create an empty array

  properties
    beams           % Internal cell array of beam objects
  end

  properties (Dependent)
    omega
    medium
    power
  end

  methods (Static)
    function beam = empty(varargin)
      % Construct an empty beam array with 'array' type.
      %
      % Usage
      %   beam = ott.beam.Array.empty()

      sz = [0, 0];
      beam = ott.beam.Array('array', sz);
    end

    function D = CombineCoherent(D)
      % Combine data coherently (by recursion)
      %
      % Usage
      %   D = CombineCoherent(D)
      %
      % Parameters
      %   D -- Cell array or other data type to be combined.
      %   Recursively combines the cells of D until D is no longer
      %   a cell array, then returns D.

      % End of recursion
      if ~iscell(D)
        return;
      end

      % For each cell in D, call CombineCoherent
      oldD = D;
      D = CombineCoherent(oldD);
      for ii = 2:numel(oldD)
        D = D + CombineCoherent(oldD{ii});
      end
    end
  end

  methods
    function beam = Array(array_type, arg)
      % Construct a new abstract beam array
      %
      % Usage
      %   beam = Array(array_type, dim)
      %
      %   beam = Array(array_type, beams)
      %
      % Parameters
      %   - array_type (enum) -- Type of beam array.  Either
      %     'array', 'coherent' or 'incoherent'.
      %
      %   - sz (N numeric) -- Size of the beam array.
      %
      %   - beams (cell) -- Cell array of beams to form beam array.

      beam = beam@ott.beam.properties.AnyArrayType(array_type);

      % Store beams or create empty
      if iscell(arg)
        % Argument is beams
        beam.beams = arg;
      else
        % Argument is size
        beam.beams = repmat({ott.beam.abstract.Empty}, arg);
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

    function data = arrayApply(beam, func, varargin)
      % Apply function to each array in the beam array output.
      %
      % Usage
      %   data = beam.arrayApply(func, ...)
      %   Additional parameters are passed to the function.
      %   All inputs must be cell arrays of the same size as the beam.
      %
      % This function is overloaded by Array types in order to
      % implement incoherent combination.

      % Apply visualisatio funtion to sub-beams
      if numel(beam) > 1

        areCells = cellfun(@iscell, varargin);
        assert(all(areCells), 'All inputs must be cell arrays');

        data = cell(size(beam));
        for ii = 1:numel(beam)
          sub_data = cellfun(@(x) x{ii}, varargin, 'UniformOutput', false);
          data{ii} =  beam.beams{ii}.arrayApply(func, sub_data{:});
        end
      else
        data = beam.beams{1}.arrayApply(func, varargin{:});
      end
    end
  end

  methods (Hidden)


    function E = deferWithCoherent(beam, func)
      % Helper for defer to beam with combine coherent check

      E = {};

      % Evaluate each beam
      for ii = 1:numel(beam)
        E{ii} = func(beam.beams{ii});
      end

      % Combine if requested
      if strcmpi(beam.array_type, 'coherent')
        E = beam.CombineCoherent(E);
      end
    end

    function E = efieldInternal(beam, varargin)
      E = beam.deferWithCoherent(@(b) b.efieldInternal(varargin{:}));
    end
    function E = hfieldInternal(beam, varargin)
      E = beam.deferWithCoherent(@(b) b.hfieldInternal(varargin{:}));
    end
    function E = efarfieldInternal(beam, varargin)
      E = beam.deferWithCoherent(@(b) b.efarfieldInternal(varargin{:}));
    end
    function E = hfarfieldInternal(beam, varargin)
      E = beam.deferWithCoherent(@(b) b.hfarfieldInternal(varargin{:}));
    end

    function validateArrayType(newType)
      % Check if beam contains incoherent beams
      if strcmpi(newType, 'coherent')
        assert(~beam.contains_incoherent, ...
          'ott:beam:Array:coherent_with_incoherent', ...
          'Cannot have coherent array of incoherent beams');
      end
    end

    function beam = subsasgnInternal(beam, subs, ~, other)
      % Assign to the subsripted beam

      % TODO: Should this support varargin
      % TODO: Should this have remaining subscripts?

      if isa(other, 'ott.beam.Array')
        beam.beams(subs{:}) = other.beams;
      else
        [beam.beams{subs{:}}] = deal(other);
      end
    end

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

      % Check if single index
      all_single = true;
      for ii = 1:numel(subs)
        all_single = all_single & numel(subs{ii}) == 1;
      end

      % TODO: Should we have () and {} indexing?  Would remove need for
      % all_single and make things a bit more 'natural'.
      if all_single
        % Return the individual beam (this would be {})
        beams = beam.beams{subs{:}};
      else
        % Form an array with the selection of new beams (behaviour of ())
        beams = [beam.beams{subs{:}}];

        % Ensure output type is not double (is this a good idea?)
        if isequaln(beams, [])
          beams = ott.beam.Array(beam.array_type, [0, 0]);
        end
      end
    end
  end

  methods % Getters/setters
    function beam = set.beams(beam, val)
      assert(iscell(val), 'beams property be cell array');

      % Cann't have coherent arrays of incoherent beams
      % TODO: Can we have incoherent/coherent array of generic arrays?
      if strcmpi(beam.array_type, 'coherent')
        for ii = 1:numel(val)
          if isa(val{ii}, 'ott.beam.properties.ArrayType')
            assert(~val{ii}.contains('incoherent'), ...
              'ott:beam:utils:ArrayType:coherent_with_incoherent', ...
              'Cannot have coherent array of incoherent beams');
          end
        end
      end

      beam.beams = val;
    end

    function power = get.power(beam)
      if strcmpi(beam.type, 'incoherent')
        power = sum(cellfun(@(x) x.power, beam.beams));
      else
        error('Power not supported for coherent/array yet');
      end
    end

    function omega = get.omega(beam)
      omega = beam.beams{1}.omega;
    end
    function beam = set.omega(beam, val)
      for ii = 1:numel(beam.beams)
        beam.beams{ii}.omega = val;
      end
    end

    function medium = get.medium(beam)
      medium = beam.beams{1}.medium;
    end
    function beam = set.medium(beam, val)
      for ii = 1:numel(beam.beams)
        beam.beams{ii}.medium = val;
      end
    end
  end
end

