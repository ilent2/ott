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

% TODO: Should arrays have beam properties?

  properties
    beams         % Beams contained in the array
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
  end

  methods % Getters/setters
  end
end
