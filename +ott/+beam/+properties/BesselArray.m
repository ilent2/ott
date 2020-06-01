classdef BesselArray < ott.beam.properties.Bessel
% Properties of Bessel beam arrays
%
% Properties
%   - angle       -- Far-field angle of Bessel beam (radians)
%   - field       -- Field in theta and phi directions
%   - lmode       -- Azimuthal angular momentum number

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (SetAccess=protected)
    angle       % Far-field angle of Bessel beam (radians)
    field       % Field in theta and phi directions
    lmode       % Azimuthal angular momentum number
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

      sz = size(beam.angle);
      [varargout{1:nargout}] = ott.utils.size_helper(sz, varargin{:});
    end
  end

  methods (Hidden)
    function beam = catInternal(dim, beam, varargin)
      % Concatenate beams

      assert(dim == 2, 'Only 1xN arrays supported (may change in future)');

      other_angle = {};
      other_field = {};
      other_lmode = {};
      for ii = 1:length(varargin)
        other_angle{ii} = varargin{ii}.angle;
        other_field{ii} = varargin{ii}.field;
        other_lmode{ii} = varargin{ii}.lmode;
      end

      beam.angle = cat(dim, beam.angle, other_angle{:});
      beam.field = cat(dim, beam.field, other_field{:});
      beam.lmode = cat(dim, beam.lmode, other_lmode{:});
    end

    function beam = plusInternal(beam1, beam2)
      % Concatenate two coherent beams together

      beam = cat(2, beam1, beam2);
    end

    function beam = subsrefInternal(beam, subs)
      % Get the subscripted beam

      if numel(subs) > ndims(beam.angle)
        if subs(1) == 1
          subs = subs(2:end);
        end
        assert(numel(subs) > ndims(beam.angle), ...
            'Too many subscript indices');
      end

      beam.angle = beam.angle(:, subs{:});
      beam.field = beam.field(:, subs{:});
      beam.lmode = beam.lmode(:, subs{:});
    end

    function beam = subsasgnInternal(beam, subs, rem, other)
      % Assign to the subscripted beam

      if numel(subs) > ndims(beam.angle)
        if subs(1) == 1
          subs = subs(2:end);
        end
        assert(numel(subs) > ndims(beam.angle), ...
            'Too many subscript indices');
      end

      assert(isempty(rem), 'Assignment to parts of beams not supported');

      if isempty(other)
        % Delete data
        beam.angle(:, subs{:}) = [];
        beam.field(:, subs{:}) = [];
        beam.lmode(:, subs{:}) = [];

      else
        assert(isa(other, 'ott.beam.properties.Bessel'), ...
            'Only PlaneWave beams supported for now');

        beam.angle(:, subs{:}) = other.angle;
        beam.field(:, subs{:}) = other.field;
        beam.lmode(:, subs{:}) = other.lmode;
      end
    end
  end
end
