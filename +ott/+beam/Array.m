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

      beam = beam@ott.beam.abstract.Array(array_type, sz, varargin{:});
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
end

