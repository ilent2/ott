classdef NumericalVolume
% Calculate the volume numerically
%
% Uses ``integral3`` and ``insideXyz` to numerically compute the volume.
%
% Abstract properties
%   boundingBox
%
% Properties
%   volume
%
% Methods
%   computeVolume   -- Method called by volume property

  properties (Abstract)
    boundingBox
  end

  properties (Dependent)
    volume
  end

  methods
    function v = computeVolume(shape, varargin)
      % Numerically compute the shapes volume
      %
      % Usage
      %   v = shape.computeVolume(...)
      %   Additional parameters are passed to `integral3`.
      %
      % Optional named parameters
      %   - RelTol (numeric) -- Relative error tolerance.
      %     Default: ``1.0e-2`` (this differs from the default).

      p = inputParser;
      p.addParameter('RelTol', 1.0e-2);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      f = @(x, y, z) double(reshape(...
          shape.insideXyz([x(:), y(:), z(:)].', ...
          'origin', 'local'), size(x)));
      bb = shape.boundingBox();

      if shape.xySymmetry
        f = @(x, y, z) 2.*f(x, y, z);
        bb(3, 1) = 0.0;
      end

      if shape.zRotSymmetry == 0

        f = @(x, z) 2.*pi.*x.*f(x, zeros(size(x)), z);
        v = integral2(f, 0, bb(1, 2), bb(3, 1), bb(3, 2), ...
            'RelTol', p.Results.RelTol, unmatched{:});

      elseif shape.zRotSymmetry == 2

        v = integral3(f, bb(1, 1), bb(1, 2), ...
            0, bb(2, 2), bb(3, 1), bb(3, 2), ...
            'RelTol', p.Results.RelTol, unmatched{:});
        v = v .* 2;

      elseif shape.zRotSymmetry == 3

        ymin = @(x) max(0, -x.*tan(pi/3));
        v = integral3(f, bb(1, 1), bb(1, 2), ...
            ymin, bb(3, 2), bb(3, 1), bb(3, 2), ...
            'RelTol', p.Results.RelTol, unmatched{:});
        v = v .* 3;

      elseif shape.zRotSymmetry >= 4

        ymax = @(x) min(bb(2, 2), x.*tan(2*pi./shape.zRotSymmetry));
        v = integral3(f, 0, bb(1, 2), ...
            0, ymax, bb(3, 1), bb(3, 2), ...
            'RelTol', p.Results.RelTol, unmatched{:});
        v = v .* shape.zRotSymmetry;

      else

        v = integral3(f, bb(1, 1), bb(1, 2), ...
            bb(2, 1), bb(2, 2), bb(3, 1), bb(3, 2), ...
            'RelTol', p.Results.RelTol, unmatched{:});
      end
    end
  end

  methods % Getters/setters
    function v = get.volume(shape)
      v = shape.computeVolume();
    end
  end
end
