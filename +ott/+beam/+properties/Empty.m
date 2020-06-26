classdef Empty < ott.beam.properties.Beam ...
    & ott.beam.properties.ZeroPower ...
    & ott.beam.properties.VariableMedium
% Properties of scattered beams
%
% Properties
%   power         -- Constant, 0 (empty beams have no power)

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Static)
    function args = likeProperties(other, args)
      args = ott.beam.properties.Beam.likeProperties(other, args);
      args = ott.beam.properties.VariableMedium.likeProperties(other, args);
    end
  end

  methods
    function beam = Empty(varargin)
      beam = beam@ott.beam.properties.Beam(varargin{:});
    end
  end

  methods (Hidden)
    function E = efieldInternal(beam, xyz, varargin)
      % Returns zeros
      E = ott.utils.FieldVector(xyz, 0*xyz, 'cartesian');
    end

    function H = hfieldInternal(beam, xyz, varargin)
      % Returns zeros
      H = ott.utils.FieldVector(xyz, 0*xyz, 'cartesian');
    end

    function E = efarfieldInternal(beam, rtp, varargin)
      % Returns zeros
      E = ott.utils.FieldVector(rtp, 0*rtp, 'spherical');
    end

    function H = hfarfieldInternal(beam, rtp, varargin)
      % Returns zeros
      H = ott.utils.FieldVector(rtp, 0*rtp, 'spherical');
    end
  end
end
