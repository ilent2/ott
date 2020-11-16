classdef FieldVectorSph < ott.utils.FieldVector
% Spherical field vector instance
% Inherits from :class:`FieldVector`.
%
% Supported casts
%   - FieldVectorCart    -- Convert to Cartesian coordinates
%
% See base class for supported methods.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function fv = FieldVectorSph(varargin)
      % Construct a new field vector instance.
      %
      % If coordiantes are provided, verifies that they satisfy
      % ``r >= 0; 0 <= t <= pi; 0 < p <= 2*pi`` otherwise issues a warning.
      %
      % Usage
      %   fv = FieldVectorSph(rtpv, rtp)
      %
      %   fv = FieldVectorSph([rtpv; rtp])
      %
      % Parameters
      %   - rtpv (3xN numeric) -- Field vectors
      %
      %   - rtp (3xN numeric) -- Coordinate locations (optional).
      %     Default: ``zeros(size(rtpv))``.

      fv = fv@ott.utils.FieldVector(varargin{:});
      
      % Check for valid range of values
      if size(fv, 1) == 6 && (any(fv(4, :)) < 0 ...
          || any(fv(5, :) < 0 | fv(5, :) > pi) ...
          || any(fv(6, :) < 0 | fv(6, :) > 2*pi))
        warning('ott:utils:FieldVectorSph:rtp_outside_range', ...
          'rtp values outside expected range, results may be unexpected');
      end
    end

    function sfv = ott.utils.FieldVectorCart(fv)
      % Cast to a Cartesian field vector

      % Get data
      data = double(fv);
      if size(data, 1) == 3
        data = [data; zeros(size(data))];
      end

      % Change coordinates
      [xyzv, xyz] = ott.utils.rtpv2xyzv(data(1:3, :), data(4:6, :));

      % Package new data
      sfv = ott.utils.FieldVectorCart(xyzv, xyz);
      sfv = reshape(sfv, size(data));
    end
  end
end
