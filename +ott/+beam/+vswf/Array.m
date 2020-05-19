classdef Array < ott.beam.Array
% Specialisation of beam array for arrays of VSWF coefficients.
%
% Adds additional methods for applying transformations to arrays of Bscs.
%
% Properties
%   - beams       -- Internal array of beam objects
%   - array_type  -- Type of array ('coherent', 'array' or 'incoherent')
%
% Methods
%   - applyTransformation
%   - applyTranslation
%   - applyRotation
%
% Type conversion methods
%   - Bsc         -- Combine beam array into single BSC

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function beam = Array(array_type, sz, varargin)
      % Construct a new beam array
      %
      % Usage
      %   beam = Array(bsc)
      %   Convert an existing bsc object to an array of separate bscs.
      %
      %   beam = Array(array_type, sz, beam1, beam2, ...)
      %   Construct a new array containing the specified BSC beams.
      %
      % Parameters
      %   - array_type (enum) -- Type of beam array.  Either
      %     'array', 'coherent' or 'incoherent'.
      %
      %   - sz (N numeric) -- Size of the beam array.
      %
      %   - beam1, beam2, ... -- Beams to include in the array.
      %     If empty, creates an empty cell array internally.

      if ischar(array_type)
        beam_data = varargin;

        for ii = 1:numel(beam_data)
          assert(isa(beam_data{ii}, 'ott.beam.vswf.Bsc'), ...
              'vswf.Array must contain only Bsc objects');
        end
      else
        assert(nargin == 1, 'Incorrect number of arguments');

        bsc = array_type;
        array_type = bsc.array_type;
        sz = size(bsc);

        beam_data = cell(1, numel(bsc));
        for ii = 1:numel(bsc)
          beam_data{ii} = bsc(ii);
        end
      end

      beam = beam@ott.beam.Array(array_type, sz, beam_data{:});
    end

    function bsc = ott.beam.vswf.Bsc(array)
      % Merge the beam shape coefficients into a single BSC object
      %
      % This action discards individual beam information such as
      % position, rotation and information about beam specialisations.

      a = array.beams{1}.a;
      b = array.beams{1}.b;

      % TODO: Grow Nmax if required

      for ii = 2:numel(array)
        a = [a, array.beams{2}.a];
        b = [b, array.beams{2}.b];
      end

      bsc = ott.beam.vswf.Bsc(a, b, 'array_type', array.array_type);
    end

    function bsc = applyTransformation(bsc, varargin)
      % Apply transformation to each BSC in the beam array
      %
      % For parameters and usage, see :meth:`Bsc.applyTransformation`.

      % TODO: Similar translations can be optimised
      %     Output should be a single Bsc with position/rotation merged

      ott.utils.nargoutCheck(bsc, nargout);

      for ii = 1:numel(bsc.beams)
        bsc.beams{ii} = bsc.beams{ii}.applyTransformation(varargin{:});
      end

      % Merge beams into one bsc
      bsc = ott.beam.vswf.Bsc(bsc);
    end

    function bsc = applyRotation(bsc, varargin)
      % Apply rotation transformation to each BSC in beam array.
      %
      % For parameters and usage, see :meth:`Bsc.applyRotation`.

      ott.utils.nargoutCheck(bsc, nargout);

      for ii = 1:numel(bsc.beams)
        bsc.beams{ii} = bsc.beams{ii}.applyRotation(varargin{:});
      end

      % Check if positions are all equal
      rots = cellfun(@(x) x.position(:), bsc.beams, 'UniformOutput', false);
      urots = unique(cell2mat(rots).', 'rows');
      if size(urots, 1) == 1
        bsc = ott.beam.vswf.Bsc(bsc);
      end
    end

    function bsc = applyTranslation(bsc, varargin)
      % Apply translation transformation to each BSC in beam array.
      %
      % For parameters and usage, see :meth:`Bsc.applyTranslation`.

      ott.utils.nargoutCheck(bsc, nargout);

      for ii = 1:numel(bsc.beams)
        bsc.beams{ii} = bsc.beams{ii}.applyTranslation(varargin{:});
      end

      % Check if rotations are all equal
      rots = cellfun(@(x) x.rotation(:), bsc.beams, 'UniformOutput', false);
      urots = unique(cell2mat(rots).', 'rows');
      if size(urots, 1) == 1
        bsc = ott.beam.vswf.Bsc(bsc);
      end
    end
  end
end
