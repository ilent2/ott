classdef FieldVectorCart < ott.utils.FieldVector
% Cartesian field vector instance
% Inherits from :class:`FieldVector`.
%
% Supported casts
%   - FieldVectorSph    -- Convert to spherical coordiantes
%
% See base class for supported methods.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function fv = FieldVectorCart(varargin)
      % Construct a new field vector instance
      %
      % Usage
      %   fv = FieldVectorCart(xyzv, xyz)
      %
      %   fv = FieldVectorCart([xyzv; xyz])
      %
      % Parameters
      %   - xyzv (3xN numeric) -- Field vectors
      %
      %   - xyz (3xN numeric) -- Coordinate locations (optional).
      %     Default: ``zeros(size(xyzv))``.

      fv = fv@ott.utils.FieldVector(varargin{:});
    end

    function sfv = ott.utils.FieldVectorSph(fv)
      % Cast the Cartesian field vector to a spherical field vector

      % Get data
      data = double(fv);
      if size(data, 1) == 3
        data = [data; zeros(size(data))];
      end

      % Change coordinates
      [rtpv, rtp] = ott.utils.xyzv2rtpv(data(1:3, :), data(4:6, :));

      % Package new data
      sfv = ott.utils.FieldVectorSph(rtpv, rtp);
      sfv = reshape(sfv, size(data));
    end
  end
end
