classdef RotationPositionProp < ott.utils.RotationPositionEntity
% Adds rotation and position properties to a class.
% Inherits from :class:`ott.utils.RotationPositionEntity`.
%
% Properties
%   - position (3x1 numeric) -- Global position of the particle
%   - rotation (3x3 numeric) -- Global orientation of the particle
%
% Methods
%   - rotate*     -- Functions for rotating the entity
%   - translate*  -- Functions for translating the entity
%   - local2global -- Transform local coordinate to global space
%   - global2local -- Transform global coordinate to local space
%
% Static methods
%   - rotationPositionHelper -- Helper for processing object arrays

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Hidden)
    positionInternal    % Internal position property
    rotationInternal    % Internal rotation property
  end

  methods
    function obj = RotationPositionProp(varargin)
      % Construct a new entity with position/rotation properties.
      %
      % Usage
      %   obj = ott.utils.RotationPositionEntity(...)
      %
      %   obj = obj@ott.utils.RotationPositionEntity(varargin{:})
      %
      % Named arguments
      %   - rotation (3x3 numeric) -- Entity global orientation.
      %     Default: ``eye(3)``.
      %
      %   - position (3x1 numeric) -- Entity global position.
      %     Default: ``[0;0;0]``.

      obj = obj@ott.utils.RotationPositionEntity(varargin{:});
    end
  end
end
