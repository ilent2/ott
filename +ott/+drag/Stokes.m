classdef Stokes
% Base class for 6-vector force/torque drag tensors.
% Inherits from :class:`Stokes`.
%
% This class is the base class for drag tensors which can be described
% by a 3x3 translational, rotational, and cross-term matrices in
% Cartesian coordinates.
%
% Properties
%  - forward      -- Rank 2 forward tensor
%  - inverse      -- Rank 2 inverse tensor
%  - allowchanges -- If false, prevents further changes to dependent proeprties
%  - viscosity    -- Viscosity of surrounding medium (default: 1.0)
%  - orientation  -- Current orientation matrix for tensor (default: eye)
%
% Dependent properties
%  - translation      -- Translational component of tensor
%  - rotation         -- Rotational component of tensor
%  - crossterms       -- Cross-terms component of tensor (UD)
%  - itranslation     -- Inverse translational component of tensor
%  - irotation        -- Inverse rotational component of tensor
%  - icrossterms      -- Inverse cross-terms component of tensor (UD)
%
% Methods
%  - inv         -- Return the inverse or calculate the inverse
%  - mtimes      -- Multiply the tensor or calculate the forward and mul
%  - rotate      -- Rotate the tensor specifying a rotation matrix
%  - rotate*     -- Rotate the tensor around the X,Y,Z axis
%  - setOrientation -- Clear the current orientation and apply a new one
%
% Static methods
%  - defaultMethod -- Determine the appropriate methods for a particular shape.
%  - simple        -- Calculate the Stokes drag for a simple shape.
%
% See also StokesSphere and StokesLambNn

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (Dependent)
    translation       % Translational component of tensor
    rotation          % Rotational component of tensor
    crossterms        % Cross-terms component of tensor (UD)

    itranslation      % Inverse translational component of tensor
    irotation         % Inverse rotational component of tensor
    icrossterms       % Inverse cross-terms component of tensor (UD)
  end

  properties
    allowchanges      % Prevent further changes after finalize
    forward           % Rank 2 tensor for the drag
    inverse           % Rank 2 inverse tensor for the drag
  end

  properties (SetAccess=protected)
    viscosity         % Viscosity of surrounding medium (default: 1.0)
    orientation       % Current orientation of the tensor
  end

  methods (Static)
    function method = defaultMethod(shape, varargin)
      % Determine the appropriate methods for a particular shape.
      % Returns one of 'sphere', 'lamb-nn', 'lamb-pm'
      % TODO: Documentation

      % TODO
    end

    function drag = simple(shape, varargin)
      % Calculate the Stokes drag for a simple shape
      % TODO: Documentation

      % Parse inputs
      p = inputParser;
      p.KeepUnmatched = true;
      p.addOptional('parameters', []);
      p.parse(varargin{:});

      % Get a shape object from the inputs
      if ischar(shape) && ~isempty(p.Results.parameters)
        shape = ott.shapes.Shape.simple(shape, p.Results.parameters);
        varargin = varargin(2:end);
      elseif ~isa(shape, 'ott.shapes.Shape') || ~isempty(p.Results.parameters)
        error('Must input either Shape object or string and parameters');
      end

      % Call the appropriate class to do the work
      method = ott.drag.Stokes.defaultMethod(shape);
      switch method
        case 'sphere'
          drag = ott.drag.StokesSphere.simple(shape, varargin{:});
        case 'lamb-pm'
          drag = ott.drag.StokesLambPm.simple(shape, varargin{:});
        case 'lamb-nn'
          drag = ott.drag.StokesLambNn.simple(shape, varargin{:});
        otherwise
          error('Internal error: unsupported method string');
      end
    end
  end

  methods (Access=protected)
    function obj = rotate_data(obj, rot)
      % Apply a rotation to the tensor data

      old_allowchanges = obj.allowchanges;
      obj.allowchanges = true;

      if ~isempty(obj.forward)
        obj.translation = rot.' * obj.translation * rot;
        obj.rotation = rot.' * obj.rotation * rot;
        obj.crossterms = rot.' * obj.crossterms * rot;
      end

      if ~isempty(obj.inverse)
        obj.itranslation = rot.' * obj.itranslation * rot;
        obj.irotation = rot.' * obj.irotation * rot;
        obj.icrossterms = rot.' * obj.icrossterms * rot;
      end

      obj.allowchanges = old_allowchanges;
    end
  end

  methods
    function obj = Stokes(varargin)
      % Construct a new drag tensor instance.
      %
      % Optional named parameters:
      %   - forward (6x6 numeric)   -- Full forward tensor
      %   - inverse (6x6 numeric)   -- Full inverse tensor
      %
      %   - translation (3x3 numeric)
      %   - itranslation (3x3 numeric)
      %   - rotation (3x3 numeric)
      %   - irotation (3x3 numeric)
      %   - crossterms (3x3 numeric)
      %   - icrossterms (3x3 numeric)
      %
      %   - finalize (logical)  -- If inverse/forward should be calculated
      %     when no values are given for the opposite part.
      %
      %   - viscosity (numeric) -- Viscosity of surrounding medium.
      %     Not used for internal calculations.  Default: 1.0.

      p = inputParser();
      p.addParameter('forward', []);
      p.addParameter('inverse', []);
      p.addParameter('translation', []);
      p.addParameter('rotation', []);
      p.addParameter('crossterms', []);
      p.addParameter('itranslation', []);
      p.addParameter('irotation', []);
      p.addParameter('icrossterms', []);
      p.addParameter('finalize', true);
      p.addParameter('viscosity', 1.0);
      p.addParameter('orientation', eye(3));
      p.parse(varargin{:});

      assert(~(~isempty(p.Results.forward) & (...
        ~isempty(p.Results.translation) | ~isempty(p.Results.rotation) ...
          | ~isempty(p.Results.crossterms))), ...
        ['Only forward or any of (forward|rotation|crossterms) ', ...
          'must be set, not both']);

      assert(~(~isempty(p.Results.inverse) & ( ...
        ~isempty(p.Results.itranslation) | ~isempty(p.Results.irotation) ...
          | ~isempty(p.Results.icrossterms))), ...
        ['Only inverse or any of (iforward|irotation|icrossterms) ', ...
          'must be set, not both']);

      % Allow changes to properties (default value)
      obj.allowchanges = true;

      % Store properties
      obj.viscosity = p.Results.viscosity;
      obj.orientation = p.Results.orientation;

      % Set forward matrix if possible
      if ~isempty(p.Results.forward)
        obj.forward = p.Results.forward;
      else
        if ~isempty(p.Results.translation)
          obj.translation = p.Results.translation;
        end
        if ~isempty(p.Results.rotation)
          obj.rotation = p.Results.rotation;
        end
        if ~isempty(p.Results.crossterms)
          obj.crossterms = p.Results.crossterms;
        end
      end

      % Set inverse matrix if possible
      if ~isempty(p.Results.inverse)
        obj.inverse = p.Results.inverse;
      else
        if ~isempty(p.Results.itranslation)
          obj.itranslation = p.Results.itranslation;
        end
        if ~isempty(p.Results.irotation)
          obj.irotation = p.Results.irotation;
        end
        if ~isempty(p.Results.icrossterms)
          obj.icrossterms = p.Results.icrossterms;
        end
      end

      % If no inverse/forward, calculate them
      if p.Results.finalize
        obj = obj.finalize();
      end
    end

    function val = get.translation(obj)
      val = obj.forward(1:3, 1:3);
    end
    function obj = set.translation(obj, val)
      assert(obj.allowchanges, 'Change will invalidate inverse');
      assert(ismatrix(val) && all(size(val) == [3, 3]), 'Must be 3x3 matrix');
      if isempty(obj.forward)
        obj.forward = zeros(6, 6);
      end
      obj.forward(1:3, 1:3) = val;
    end

    function val = get.rotation(obj)
      val = obj.forward(4:6, 4:6);
    end
    function obj = set.rotation(obj, val)
      assert(obj.allowchanges, 'Change will invalidate inverse');
      assert(ismatrix(val) && all(size(val) == [3, 3]), 'Must be 3x3 matrix');
      if isempty(obj.forward)
        obj.forward = zeros(6, 6);
      end
      obj.forward(4:6, 4:6) = val;
    end

    function val = get.crossterms(obj)
      val = obj.forward(1:3, 4:6);
    end
    function obj = set.crossterms(obj, val)
      assert(obj.allowchanges, 'Change will invalidate inverse');
      assert(ismatrix(val) && all(size(val) == [3, 3]), 'Must be 3x3 matrix');
      if isempty(obj.forward)
        obj.forward = zeros(6, 6);
      end
      obj.forward(4:6, 1:3) = val;
      obj.forward(1:3, 4:6) = val.';
    end

    function val = get.itranslation(obj)
      val = obj.inverse(1:3, 1:3);
    end
    function obj = set.itranslation(obj, val)
      assert(obj.allowchanges, 'Change will invalidate forward');
      assert(ismatrix(val) && all(size(val) == [3, 3]), 'Must be 3x3 matrix');
      if isempty(obj.inverse)
        obj.inverse = zeros(6, 6);
      end
      obj.inverse(1:3, 1:3) = val;
    end

    function val = get.irotation(obj)
      val = obj.inverse(4:6, 4:6);
    end
    function obj = set.irotation(obj, val)
      assert(obj.allowchanges, 'Change will invalidate forward');
      assert(ismatrix(val) && all(size(val) == [3, 3]), 'Must be 3x3 matrix');
      if isempty(obj.inverse)
        obj.inverse = zeros(6, 6);
      end
      obj.inverse(4:6, 4:6) = val;
    end

    function val = get.icrossterms(obj)
      val = obj.inverse(1:3, 4:6);
    end
    function obj = set.icrossterms(obj, val)
      assert(obj.allowchanges, 'Change will invalidate forward');
      assert(ismatrix(val) && all(size(val) == [3, 3]), 'Must be 3x3 matrix');
      if isempty(obj.inverse)
        obj.inverse = zeros(6, 6);
      end
      obj.inverse(1:3, 4:6) = val;
      obj.inverse(4:6, 1:3) = val.';
    end

    function iTensor = inv(obj)
      % Return the inverse tensor if set or the inverse of the forward tensor
      %
      % Usage:
      %    inverse_tensor = inv(tensor) returns a 2-D matrix representing
      %    the tensor inverse.
      %
      % If tensor.inverse is unset, attempts to calculate the inverse from
      % the forward tensor.

      if ~isempty(obj.inverse)
        iTensor = obj.inverse;
      elseif ~isempty(obj.forward)
        warning('ott:drag:Stokes:inv_calculation', ...
          ['Calculating forward tensor inverse.  ', ...
           'May not work in some coordinates, may be slow for ', ...
           'repeated calculation']);
        iTensor = inv(obj.forward);
      else
        error('ott:drag:Stokes:no_inverse', 'Cannot calculate inverse');
      end
    end

    function drag = mtimes(obj, vec)
      % Calculate the drag using the tensor and a velocity vector
      %
      % Usage:
      %   drag = tensor * vec; where vec is the velocity vector.
      %
      % If tensor.forward is not set, attempts to calculate the drag
      % from the inverse of the inverse drag tensor.

      if ~isempty(obj.forward)
        tensor = obj.forward;
      elseif ~isempty(obj.inverse)
        warning('ott:drag:Stokes:forward_calculation', ...
          ['Calculating inverse tensor inverse.  ', ...
           'May not work in some coordinates, may be slow for ', ...
           'repeated calculation']);
        tensor = inv(obj.inverse);
      else
        error('ott:drag:Stokes:no_forward', 'Cannot calculate drag');
      end

      drag = tensor * vec;
    end

    function obj = set.forward(obj, val)
      assert((ismatrix(val) & all(size(val) == [6, 6])) | isempty(val), ...
        'Must be 6x6 matrix or empty');
      obj.forward = val;
    end

    function obj = set.inverse(obj, val)
      assert((ismatrix(val) & all(size(val) == [6, 6])) | isempty(val), ...
        'Must be 6x6 matrix or empty');
      obj.inverse = val;
    end
    
    function obj = finalize(obj)
      % If no inverse/forward, calculate them
      %
      % Sets the allowchanges value to false when calculation done.
      
      if isempty(obj.forward) && ~isempty(obj.inverse)
        obj.forward = inv(obj.inverse);
        obj.allowchanges = false;
      elseif isempty(obj.inverse) && ~isempty(obj.forward)
        obj.inverse = inv(obj.forward);
        obj.allowchanges = false;
      end

    end

    function mat = diag(obj)
      % Return the diagonal of the forward drag tensor
      mat = diag(obj.forward);
    end

    function obj = set.orientation(obj, val)
      % Update the orientation
      assert(ismatrix(val) && all(size(val) == [3, 3]), ...
        'Orientation must be 3x3 matrix');
      obj.orientation = val;
    end

    function obj = setOrientation(obj, val)
      % Clear the current orientation and set a new orientation

      rot = eye(3);

      % Clear the current orientation
      if any(any(obj.orientation ~= eye(3)))
        rot = obj.orientation.' * rot;
      end

      % Apply the new rotation
      if any(any(val ~= eye(3)))
        rot = val * rot;
      end

      % Apply all at once and store results
      obj = obj.rotate_data(rot);
      obj.orientation = val;
    end

    function obj = rotateX(obj, angle_rad)
      % Rotate the particle about the x-axis an angle in radians
      import ott.utils.*;
      obj = obj.rotate(rotx(180*angle_rad/pi));
    end

    function obj = rotateY(obj, angle_rad)
      % Rotate the particle about the y-axis an angle in radians
      import ott.utils.*;
      obj = obj.rotate(roty(180*angle_rad/pi));
    end

    function obj = rotateZ(obj, angle_rad)
      % Rotate the particle about the z-axis an angle in radians
      import ott.utils.*;
      obj = obj.rotate(rotz(180*angle_rad/pi));
    end

    function obj = rotateXy(obj, anglex, angley)
      % Rotate the particle about the x then y axis (units radians)
      import ott.utils.*;
      obj = obj.rotate(roty(180*angley/pi) * rotx(180*anglex/pi));
    end

    function obj = rotateXz(obj, anglex, anglez)
      % Rotate the particle about the x then z axis (units radians)
      import ott.utils.*;
      obj = obj.rotate(rotz(180*anglez/pi) * rotx(180*anglex/pi));
    end

    function obj = rotateYz(obj, angley, anglez)
      % Rotate the particle about the y then z axis (units radians)
      import ott.utils.*;
      obj = obj.rotate(rotz(180*anglez/pi) * roty(180*angley/pi));
    end

    function obj = rotateXyz(obj, anglex, angley, anglez)
      % Rotate the particle about the x, y then z axis (units radians)
      import ott.utils.*;
      obj = obj.rotate(rotz(180*anglez/pi) * ...
          roty(180*angley/pi) * rotx(180*anglex/pi));
    end

    function obj = rotate(obj, rot)
      % Applies a 3x3 rotation matrix to the current orientation
      obj = obj.rotate_data(rot);
      obj.orientation = rot * obj.orientation;
    end

    function obj = clearRotation(obj)
      % Clear the current rotation orientation
      obj = obj.rotate(obj.orientation.');
      obj.orientation = eye(3);
    end

    function mat = double(obj)
      mat = obj.forward;
    end

  end % methods

end
