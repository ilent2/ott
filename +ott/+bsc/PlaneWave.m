classdef PlaneWave < ott.bsc.Bsc
% Bsc specialisation for plane waves
%
% Plane waves can be translated in any direction by simply applying
% a phase shift to the beam.  This class provides overloads for the
% beam translation functions implementing this optimisation.
%
% Properties
%   - direction     -- Propagation direction of plane wave
%
% Static methods
%   - FromDirection -- Construct a beam for the specified direction

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    direction       % Propagation direction of plane wave
  end

  methods (Static)
    function [beam, data] = FromDirection(Nmax, direction, ...
        polarisation, varargin)
      % Construct a beam for the specified direction and Nmax
      %
      % Usage
      %   [beam, data] = ott.bsc.PlaneWave.FromDirection(...
      %       Nmax, direction, polarisation)
      %
      % Parameters
      %   - Nmax (numeric) -- Size of beam shape coefficient data.
      %   - direction (3xN numeric) -- Directions for plane waves.
      %   - polarisation (3xN numeric) -- Directions for polarisation.
      %
      % Optional named arguments
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.

      p = inputParser;
      p.addParameter('data', ott.utils.VswfData(), ...
          @(x) isa(x, 'ott.utils.VswfData'));
      p.parse(varargin{:});

      assert(isnumeric(Nmax) && isscalar(Nmax) ...
          && Nmax >= 0 && round(Nmax) == Nmax, ...
          'Nmax should be an single positive integer');

      % Get spherical coordinates for plane wave vectors
      rtp = ott.utils.xyzv2rtpv(direction, zeros(size(direction)));
      Ertp = ott.utils.xyzv2rtpv(polarisation, zeros(size(polarisation)));

      ablength = ott.utils.combined_index(Nmax, Nmax);

      a = zeros(ablength, size(rtp, 2));
      b = zeros(ablength, size(rtp, 2));

      data = p.Results.data;

      for n = 1:Nmax
        iter= (n-1)*(n+1)+1:n*(n+2);
        leniter=2*n+1;

        %expand theta and phi components of field to match spherical harmonics
        ET=repmat(Ertp(2, :), [1,leniter]);
        EP=repmat(Ertp(3, :), [1,leniter]);

        %power normalisation.
        Nn = 1/sqrt(n*(n+1));

        %Generate the farfield components of the VSWFs
        ci = ott.utils.combined_index(n, -n:n);
        [~, dtY, dpY, data] = data.evaluateYtp(ci, rtp(2, :), rtp(3, :));

        %equivalent to dot((1i)^(n+1)*C,E);
        a(iter,:) = 4*pi*Nn*(-1i)^(n+1)*(dpY'.*ET - dtY'.*EP).';
        %equivalent to dot((1i)^(n)*B,E);
        b(iter,:) = 4*pi*Nn*(-1i)^(n)  *(dtY'.*ET + dpY'.*EP).';
      end

      beam = ott.bsc.PlaneWave(a, b, direction);
      beam = beam.makeSparse();
    end
  end

  methods
    function beam = PlaneWave(varargin)
      % Construct a new beam object
      %
      % Usage
      %   beam = PlaneWave() Constructs an empty Bsc beam.
      %
      %   beam = PlaneWave(a, b, direction, ...) constructs beam from
      %   a/b coefficients and direction vector.
      %   Does not verify that a/b coefficients match direction.
      %   See :meth:`FromDirection` for a constructor with only direction.
      %
      % Parameters
      %   - a,b (numeric) -- Vectors of VSWF coefficients
      %   - direction (3xN numeric) -- Direction vectors for each beam

      p = inputParser;
      p.addOptional('a', [], @isnumeric);
      p.addOptional('b', [], @isnumeric);
      p.addOptional('direction', [], @isnumeric);
      p.parse(varargin{:});

      beam = beam@ott.bsc.Bsc(p.Results.a, p.Results.b);
      beam.direction = p.Results.direction;
    end

    function beam = translateZ(beam, z, varargin)
      % Apply translation along z axis using phase shift.
      %
      % Usage
      %   beam = beam.translateZ(z)

      ott.utils.nargoutCheck(beam, nargout);

      xyz = [0*z(:), 0*z(:), z(:)].';
      beam = beam.translateXyzInternal(xyz);
    end
  end

  methods (Hidden)
    function beam = rotateInternal(beam, R, varargin)
      % Apply rotation to beam and beam data
      %
      % Calls the base class rotation function and also updates the
      % internal beam direction vector.

      rbeam = rotateInternal@ott.bsc.Bsc(beam, R, varargin{:});
      beam = ott.bsc.PlaneWave(rbeam.a, rbeam.b, R * beam.direction);
    end

    function beam = translateXyzInternal(beam, xyz, varargin)
      % Apply translation to beam using phase shift
      %
      % Usage
      %   beam = beam.translateXyzInternal(xyz)

      ott.utils.nargoutCheck(beam, nargout);

      assert(isnumeric(xyz) && ismatrix(xyz) && size(xyz, 1) == 3, ...
          'xyz must be 3xN numeric matrix');

      ibeam = beam;

      for ii = 1:size(xyz, 2)
        dz = dot(repmat(xyz(:, ii), 1, numel(ibeam)), [ibeam.direction]);
        dz = exp(1i.*dz.*2*pi);

        % Apply translation
        beam(ii) = ibeam .* dz;
      end
    end
  end

  methods % Getters/setters
    function beam = set.direction(beam, val)
      assert(isvector(val) && numel(val) == 3 && isnumeric(val), ...
          'Direction must be 3 element numeric vector');
      beam.direction = val(:);
    end
  end
end
