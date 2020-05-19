classdef (Abstract) TranslateHelper
% A helper class providing translation methods.
%
% This class provides methods for global translations which can be
% applied irrespective of the object orientation.
%
% Supports array objects which implement `repelem`.
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
      %
      % Parameters
      %   - obj -- An instance of the class or array instance.
      %   - xyz (3xN numeric) -- Positions to translate to.
      %
      % If both `obj` and `xyz` are arrays, they must have the same
      % length, otherwise one must be a scalar.

      ott.utils.nargoutCheck(obj, nargout);

      % Check inputs
      assert(isnumeric(xyz) && ismatrix(xyz) && size(xyz, 1) == 3, ...
          'xyz must be 3xN numeric matrix');

      % Get input sizes and check
      Nxyz = size(xyz, 2);
      Nobj = numel(obj);
      assert(Nxyz == 1 || Nobj == 1 || Nobj == Nxyz, ...
          'Number of objects and positions must be same length or singular');

      Nwork = max([Nxyz, Nobj]);

      % Only replicate objects, let internal method choose how to
      % handle single or multiple positions
      if Nobj == 1, obj = repelem(obj, 1, Nwork); end;

      [varargout{1:nargout}] = obj.translateXyzInternal(xyz, varargin{:});
    end

    function varargout = translateRtp(obj, rtp, varargin)
      % Translate the object by a distance specified in Spherical coordinates
      %
      % Usage
      %   [obj, ...] = obj.translateRtp(rtp, ...)
      %   Angles are specified in radians.
      %
      % Parameters
      %   - obj -- An instance of the class or array instance.
      %   - rtp (3xN numeric) -- Positions to translate to.
      %
      % If both `obj` and `rtp` are arrays, they must have the same
      % length, otherwise one must be a scalar.

      assert(isnumeric(rtp) && ismatrix(rtp) && size(rtp, 1) == 3, ...
          'rtp must be 3xN numeric matrix');

      xyz = ott.utils.rtp2xyz(rtp(:));
      [varargout{1:nargout}] = obj.translateXyz(xyz, varargin{:});
    end
  end

end
