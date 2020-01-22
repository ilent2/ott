classdef Stokes
% Base class for Stokes drag tensors
%
% This class contains a matrix for the forward and inverse drag tensors.
% Depending on the instance, either or both of these may be set.
% The only conditions on these tensors is that they be rank 2, the size,
% and coordinate system depends on the instance.
%
% To calculate the drag::
%
%   tensor = Stokes(eye(3));
%   drag = tensor * [1;2;3];
%
% To calculate the inverse drag::
%
%   tensor = Stokes([], eye(3));
%   idrag = inv(tensor) * [1,2,3];
%
% Properties:
%  - forward     -- Rank 2 forward tensor
%  - inverse     -- Rank 2 inverse tensor
%
% Methods:
%  - inv         -- Return the inverse or calculate the inverse
%  - mtimes      -- Multiply the tensor or calculate the forward and mul
%
% See also ott.drag.Sphere

  properties
    forward         % Rank 2 tensor for the drag
    inverse         % Rank 2 inverse tensor for the drag
  end

  methods
    function obj = Stokes(forward, inverse)
      % Construct a new tensor instance
      %
      % Usage:
      %   tensor = Stokes(forward)
      %
      %   tensor = Stokes(forward, inverse)
      %
      % Parameters:
      %   - forward  -- Forward tensor (optional)
      %   - inverse  -- Inverse tensor (optional)

      assert(nargin >= 1 && nargin <= 2, 'Input must be 1 or 2 elements');
      if nargin < 2
        inverse = [];
        if nargin < 1
          forward = [];
        end
      end

      obj.forward = forward;
      obj.inverse = inverse;
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
          ['Calculating forward tensor inverse.  ',
           'May not work in some coordinates, may be slow for ',
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
          ['Calculating inverse tensor inverse.  ',
           'May not work in some coordinates, may be slow for ',
           'repeated calculation']);
        tensor = inv(obj.inverse);
      else
        error('ott:drag:Stokes:no_forward', 'Cannot calculate drag');
      end

      drag = tensor * vec;
    end

    function obj = set.forward(obj, val)
      assert(ismatrix(val) | isempty(val), 'Must be matrix');
      assert(diff(size(val)) == 0, 'Must be square');
      obj.forward = val;
    end

    function obj = set.inverse(obj, val)
      assert(ismatrix(val) | isempty(val), 'Must be matrix');
      assert(diff(size(val)) == 0, 'Must be square');
      obj.inverse = val;
    end

  end % methods

end

