classdef Stokes6 < ott.drag.Stokes
% Base class for 6-vector force/torque drag tensors.
% Inherits from :class:`Stokes`.
%
% This class is the base class for drag tensors which can be described
% by a 3x3 translational, rotational, and cross-term matrices in
% Cartesian coordinates.
%
% See also Sphere

  properties (Dependent)
    translation       % Translational component of tensor
    rotation          % Rotational component of tensor
    crossterms        % Cross-terms component of tensor (UD)

    itranslation      % Inverse translational component of tensor
    irotation         % Inverse rotational component of tensor
    icrossterms       % Inverse cross-terms component of tensor (UD)
  end

  properties
    allowchanges         % Prevent further changes after finalize
  end

  methods
    function obj = Stokes6(varargin)
      % Construct a new Stokes6 drag tensor instance
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
      p.parse(varargin{:});

      assert(~(~isempty(p.Results.forward) & (...
        ~isempty(p.Results.translation) | ~isempty(p.Results.rotation) ...
          | ~isempty(p.Results.crossterms))), ...
        ['Only forward or any of (forward|rotation|crossterms) ',
          'must be set, not both']);

      assert(~(~isempty(p.Results.inverse) & ( ...
        ~isempty(p.Results.itranslation) | ~isempty(p.Results.irotation) ...
          | ~isempty(p.Results.icrossterms))), ...
        ['Only inverse or any of (iforward|irotation|icrossterms) ',
          'must be set, not both']);

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
      obj.allowchanges = true;
      if p.Results.finalize
        if isempty(obj.forward) && ~isempty(obj.inverse)
          obj.forward = inv(obj.inverse);
          obj.allowchanges = false;
        elseif isempty(obj.inverse) && ~isempty(obj.forward)
          obj.inverse = inv(obj.forward);
          obj.allowchanges = false;
        end
      end
    end

    function val = get.translation(obj)
      val = obj.forward(1:3, 1:3);
    end
    function obj = set.translation(obj, val)
      assert(obj.nochanges, 'Change will invalidate inverse');
      assert(ismatrix(val) && size(val) == [3, 3], 'Must be 3x3 matrix');
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
      assert(ismatrix(val) && size(val) == [3, 3], 'Must be 3x3 matrix');
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
      assert(ismatrix(val) && size(val) == [3, 3], 'Must be 3x3 matrix');
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
      assert(ismatrix(val) && size(val) == [3, 3], 'Must be 3x3 matrix');
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
      assert(ismatrix(val) && size(val) == [3, 3], 'Must be 3x3 matrix');
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
      assert(ismatrix(val) && size(val) == [3, 3], 'Must be 3x3 matrix');
      if isempty(obj.inverse)
        obj.inverse = zeros(6, 6);
      end
      obj.inverse(1:3, 4:6) = val;
      obj.inverse(4:6, 1:3) = val.';
    end
  end

end
