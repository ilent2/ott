classdef AxisymShape
% Base class for axis-symmetric shapes
%
% Abstract properties
%   - perimeter       -- Perimeter from in axis plane

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (Abstract)
    perimeter         % Perimeter from in axis plane
  end

  properties (Dependent)
    zRotSymmetry      % Constant: 0
  end

  methods (Abstract)
    normalsRzInternal(obj)    % normals Cylindrical coordinates
    normalsRtInternal(obj)    % normals Polar coordinates
    insideRzInternal(obj)     % inside Cylindrical coordinates
    insideRtInternal(obj)     % inside Polar coordinates
  end

  methods (Hidden)
    function b = insideXyzInternal(shape, xyz, varargin)
      % Convert xy to radial
      rz = [vecnorm(xyz(1:2, :)); xyz(3, :)];
      b = shape.insideRzInternal(rz, varargin{:});
    end

    function b = insideRtpInternal(shape, rtp, varargin)
      % Drop phi coordinate
      b = shape.insideRtInternal(rtp(1:2, :), varargin{:});
    end

    function nxyz = normalsXyzInternal(shape, xyz, varargin)

      % Convert xy to radial
      rz = [vecnorm(xyz(1:2, :)); xyz(3, :)];

      nxz = shape.normalsRzInternal(rz, varargin{:});

      % Convert to Xyz coordinates
      phi = atan2(xyz(2, :), xyz(1, :));
      nxyz = [nxz(1, :).*cos(phi); nxz(1, :).*sin(phi); nxz(2, :)];
    end

    function nxyz = normalsRtpInternal(shape, rtp, varargin)

      % Drop phi coordinate
      nxz = shape.normalsRtInternal(rtp(1:2, :), varargin{:});

      % Convert to Xyz coordinates
      phi = rtp(3, :);
      nxyz = [nxz(1, :).*cos(phi); nxz(1, :).*sin(phi); nxz(2, :)];
    end
  end

  methods % Getters/setters
    function z = get.zRotSymmetry(~)
      z = 0;
    end
  end
end
