classdef FarfieldMapping
% Declares a far-field mapping property for beams
%
% Static methods
%   - paraxial2farfield     -- Convert paraxial coordinates to far-field
%   - farfield2paraxial     -- Convert farfield coordinates to paraxial
%
% Properties
%   - mapping       -- Far-field mapping (sintheta or tantheta)
%   - hemisphere    -- Far-field hemisphere (pos or neg)

  properties
    mapping
    hemisphere
  end

  methods (Static)
    function args = likeProperties(other, args)
      % Construct array of like-proeprties
      if isa(other, 'ott.beam.utils.FarfieldMapping')
        args = ott.utils.addDefaultParameter('mapping', other.mapping, args);
      end
    end

    function [nxyz, dA] = farfield2paraxial(rtp, varargin)
      % Convert form far-field coordinates to paraxial coordinates
      %
      % Usage
      %   [nxyz, dA] = farfield2paraxial(rtp, ...)
      %
      % Optional named arguments
      %   - mapping (enum) -- Mapping from theta-phi to far-field.
      %     Must be one of 'sin', 'tan' or 'theta'.
      %     Default: ``'sin'``.
      %
      %   - direction (enum) -- Beam hemisphere to calculate.
      %     Must be 'neg' or 'pos'.  Default: ``pos``.
      %
      %   - keepz (logical) -- If true, the output is a 3xN matrix.
      %     z-position is zeros.  Default: ``size(rtp, 1) == 3``.

      p = inputParser;
      p.addParameter('mapping', 'sin');
      p.addParameter('direction', 'pos');
      p.addParameter('keepz', size(rtp, 1) == 3);
      p.parse(varargin{:});

      assert(isnumeric(rtp) && any(size(rtp, 1) == [2, 3]), ...
        'rtp must be 2xN or 3xN numeric matrix');

      if size(rtp, 1) == 3
        theta = rtp(2, :);
        phi = rtp(3, :);
      else
        theta = rtp(1, :);
        phi = rtp(2, :);
      end

      % Apply mapping
      switch p.Results.mapping
        case 'sin'
          rr = sin(theta);
          dA = cos(theta);
        case 'tan'
          rr = tan(theta);
          dA = 1./cos(theta).^2;
        case 'theta'
          rr = theta;
          dA = 1.0;
        otherwise
          error('Unknown mapping argument value, must be sin, tan or theta');
      end

      y = rr .* sin(phi);
      x = rr .* cos(phi);

      % TODO: We don't use direction

      % Package output
      if p.Results.keepz
        nxyz = [x; y; zeros(size(x))];
      else
        nxyz = [x; y];
      end
    end

    function [rtp, dA] = paraxial2farfield(nxyz, varargin)
      % Convert from paraxial coordinates to far-field coordinates
      %
      % Usage
      %   rtp = paraxial2farfield(nxyz, ...)
      %
      % Optional named arguments
      %   - mapping (enum) -- Mapping from theta-phi to far-field.
      %     Must be one of 'sin', 'tan' or 'theta'.
      %     Default: ``'sin'``.
      %
      %   - direction (enum) -- Beam hemisphere to calculate.
      %     Must be 'neg' or 'pos'.  Default: ``pos``.
      %
      %   - keepr (logical) -- If true, the output is a 3xN matrix.
      %     Radial dimension is ones.  Default: ``size(nxyz, 1) == 3``.

      p = inputParser;
      p.addParameter('mapping', 'sin');
      p.addParameter('direction', 'pos');
      p.addParameter('keepr', size(nxyz, 1) == 3);
      p.parse(varargin{:});

      assert(isnumeric(nxyz) && any(size(nxyz, 1) == [2, 3]), ...
        'nxyz must be 2xN or 3xN numeric matrix');

      % Get phi/rr coordinates
      phi = atan2(nxyz(2, :), nxyz(1, :));
      rr = sqrt(nxyz(2, :).^2 + nxyz(1, :).^2);

      % Apply mapping
      switch p.Results.mapping
        case 'sin'
          theta = asin(rr);
          dA = 1./sqrt(1 - rr.^2);
        case 'tan'
          theta = atan(rr);
          dA = 1./(1 + rr.^2);
        case 'theta'
          theta = rr;
          dA = ones(size(rr));
        otherwise
          error('Unknown mapping argument value, must be sin, tan or theta');
      end

      % Filter out non-real theta
      mask = imag(theta) ~= 0;
      theta(mask) = nan;
      phi(mask) = nan;

      % Flip direction if needed
      switch p.Results.direction
        case 'neg'
          theta = pi - theta;
        case 'pos'
          % Nothing to do
        otherwise
          error('Unknown direction argument value, must be pos or neg');
      end

      % Package output
      if p.Results.keepr
        rtp = [ones(size(rr)); theta; phi];
      else
        rtp = [theta; phi];
      end
    end
  end

  methods
    function beam = FarfieldMapping(varargin)
      % Construct far-field mapping property
      %
      % Usage
      %   beam = beam@ott.beam.utils.FarfieldMapping(mapping, hemisphere, ...)
      %
      % Named arguments
      %   - mapping (enum) -- Initial mapping value.
      %   - hemisphere (enum) -- Hemisphere for paraxial field.

      p = inputParser;
      p.addOptional('mapping', [], ...
          @(x) any(strcmpi(x, {'sin', 'tan', 'theta'})));
      p.addOptional('hemisphere', [], ...
          @(x) any(strcmpi(x, {'pos', 'neg'})));
      p.parse(varargin{:});

      beam.mapping = p.Results.mapping;
      beam.hemisphere = p.Results.hemisphere;
    end
  end

  methods % Getters/setters
    function beam = set.mapping(beam, val)
      assert(any(strcmpi(val, {'sin', 'tan', 'theta'})), ...
          'mapping must be ''sin'', ''tan'' or ''theta''');
      beam.mapping = val;
    end

    function beam = set.hemisphere(beam, val)
      assert(any(strcmpi(val, {'pos', 'neg'})), ...
          'hemisphere must be pos or neg');
      beam.hemisphere = val;
    end
  end
end
