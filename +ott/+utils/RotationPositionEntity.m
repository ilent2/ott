classdef (Abstract) RotationPositionEntity < ...
    ott.utils.RotateHelper & ott.utils.TranslateHelper
% Describes the rotation/position interface for OTT entities.
% Inherits from :class:`RotateHelper` and :class:`TranslateHelper`.
%
% This class describes a common interface for OTT entities including
% beams, shapes and particles.  In addition to implementing the
% rotate/translate interfaces, this class also adds support for
% local translation/rotations.
%
% The simplest implementation of this class is to define position
% and rotation properties which store the position/rotation.
%
% Dependent properties
%   - position (3x1 numeric) -- Global position of the particle
%   - rotation (3x3 numeric) -- Global orientation of the particle
%
% Methods
%   - rotate*     -- Functions for rotating the entity
%   - translate*  -- Functions for translating the entity
%   - local2global -- Transform local coordinate to global space
%   - global2local -- Transform global coordinate to local space
%
% Abstract properties
%   - positionInternal -- Class specific implementation of position
%   - rotationInternal -- Class specific implementation of rotation

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Abstract, Hidden)
    positionInternal  % Class specific implementation of position
    rotationInternal  % Class specific implementation of rotation
  end

  properties (Dependent)
    position % (3x1 numeric) -- Global position of the particle
    rotation % (3x3 numeric) -- Global orientation of the particle
  end

  methods
    function obj = RotationPositionEntity(varargin)
      % Construct a new entity with position/rotation properties.
      %
      % Usage
      %   obj = obj@ott.utils.RotationPositionEntity(varargin{:})
      %
      % Named arguments
      %   - rotation (3x3 numeric) -- Entity global orientation.
      %     Default: ``eye(3)``.
      %
      %   - position (3x1 numeric) -- Entity global position.
      %     Default: ``[0;0;0]``.

      p = inputParser;
      p.addParameter('position', [0;0;0]);
      p.addParameter('rotation', eye(3));
      p.parse(varargin{:});

      obj.position = p.Results.position;
      obj.rotation = p.Results.rotation;
    end

    function xyz = local2global(obj, xyz)
      % Apply position/rotation transformations to the coordinates
      %
      % Usage
      %   xyz = obj.local2global(xyz)

      assert(isnumeric(xyz) && ismatrix(xyz) && size(xyz, 1) == 3, ...
          'xyz must be 3xN numeric matrix');

      xyz = obj.rotation * xyz + obj.position;
    end

    function xyz = global2local(obj, xyz)
      % Apply position/rotation transformations to the coordinates.
      %
      % Usage
      %   xyz= obj.global2local(xyz)

      assert(isnumeric(xyz) && ismatrix(xyz) && size(xyz, 1) == 3, ...
          'xyz must be 3xN numeric matrix');

      xyz = obj.rotation.' * (xyz - obj.position);
    end
  end

  methods (Hidden)
    function varargout = rotateInternal(obj, mat, varargin)
      % Implement global and local rotations

      p = inputParser;
      p.addParameter('origin', 'global');
      p.parse(varargin{:});

      switch p.Results.origin
        case 'global'
          obj.rotation = mat * obj.rotation;
        case 'local'
          obj.rotation = obj.rotation * mat;
        otherwise
          error('Unknown origin for rotation');
      end

      if nargout > 0
        varargout{1} = obj;
      end
    end

    function varargout = translateXyzInternal(obj, xyz, varargin)
      % Implement global and local translations

      p = inputParser;
      p.addParameter('origin', 'global');
      p.parse(varargin{:});
      
      if size(xyz, 2) == 1 && numel(obj) > 1
        xyz = repmat(xyz, 1, numel(obj));
      end

      switch p.Results.origin
        case 'global'
          for ii = 1:numel(obj)
            obj(ii).position = obj(ii).position + xyz(:, ii);
          end
        case 'local'
          for ii = 1:numel(obj)
            obj(ii).position = obj(ii).position + obj(ii).rotation*xyz(:, ii);
          end
        otherwise
          error('Unknown origin for translation');
      end

      if nargout > 0
        varargout{1} = obj;
      end
    end
  end

  methods % Getters/setters
    % position <-> positionInternal
    % rotation <-> rotationInternal

    function obj = set.position(obj, val)
      assert(isnumeric(val) && numel(val) == 3, ...
          'position must be 3 element numeric vector');
      obj.positionInternal = val;
    end
    function val = get.position(obj, val)
      val = obj.positionInternal;
    end

    function obj = set.rotation(obj, val)
      assert(isnumeric(val) && ismatrix(val) && all(size(val) == [3, 3]), ...
          'position must be 3x3 numeric matrix');
      obj.rotationInternal = val;
    end
    function val = get.rotation(obj, val)
      val = obj.rotationInternal;
    end
  end
end
