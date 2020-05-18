classdef (Abstract) TranslateHelper
% A helper class providing translation methods.
%
% This class provides methods for global translations which can be
% applied irrespective of the object orientation.
%
% Methods
%   - translateXyz      -- Apply translation in Cartesian coordinates
%   - translateRtp      -- Apply translation in Spherical coordinates
%
% Abstract methods
%   - translateXyzInternal -- Method called by translate*

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Abstract)
    translateXyzInternal   % Method called by translate*
  end

  methods
    function varargout = translateXyz(obj, xyz, varargin)
      % Translate the object by a distance specified in Cartesian coordinates
      %
      % Usage
      %   [obj, ...] = obj.translateXyz(xyz, ...)

      ott.utils.nargoutCheck(obj, nargout);
      assert(isnumeric(xyz) && numel(xyz) == 3, ...
          'xyz must be 3 element numeric vector');

      [varargout{1:nargout}] = obj.translateXyzInternal(xyz, varargin{:});
    end

    function varargout = translateRtp(obj, rtp, varargin)
      % Translate the object by a distance specified in Spherical coordinates
      %
      % Usage
      %   [obj, ...] = obj.translateRtp(rtp, ...)
      %   Angles are specified in radians.

      assert(isnumeric(rtp) && numel(rtp) == 3, ...
          'rtp must be 3 element numeric vector');

      xyz = ott.utils.rtp2xyz(rtp(:));
      [varargout{1:nargout}] = obj.translateXyzInternal(xyz, varargin{:});
    end
  end

end
