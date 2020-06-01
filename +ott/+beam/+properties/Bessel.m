classdef Bessel < ott.beam.properties.Material
% Properties of a Bessel beam.
%
% Properties
%   - angle       -- Far-field angle of Bessel beam (radians)
%   - field       -- Field in theta and phi directions
%   - lmode       -- Azimuthal angular momentum number

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Abstract, SetAccess=protected)
    angle       % Far-field angle of Bessel beam (radians)
    field       % Field in theta and phi directions
    lmode       % Azimuthal angular momentum number
  end

  methods (Static)
    function args = likeProperties(other, args)
      if isa(other, 'ott.beam.properties.Bessel')
        args = ott.utils.addDefaultParameter('angle', other.angle, args);
        args = ott.utils.addDefaultParameter('field', other.field, args);
        args = ott.utils.addDefaultParameter('lmode', other.lmode, args);
      end
      args = ott.beam.properties.Material.likeProperties(other, args);
    end

    function Etp = CartesianFieldWeights(angle, Exy)
      % Calculate weights for spherical to Cartesian polarisation conversion.
      %
      % The weights can be used to convert from spherical to Cartesian
      % Bessel beams.  For example, to construct a beam with [1;0]
      % polarisation, use::
      %
      %   beam = weights(1, :).*beam_theta + weights(2, :).*beam_phi
      %
      % where `beam_theta` and `beam_phi` are Bessel beams with
      % unitary polarisation in the theta/phi spherical directions.
      %
      % Usage
      %   weights = CartesianFieldWeights(angle, Exy)
      %
      % Parameters
      %   - angle (N numeric) -- Incoming Bessel angle (radians).
      %   - Exy (2xN numeric) -- Field in Cartesian coordinates [x;y].
      %
      % Returns
      %   - weights (2xN numeric) -- Conversion weights

      assert(isnumeric(Exy) && ismatrix(Exy) && size(Exy, 1) == 2, ...
          'Exy must be 2xN numeric');
      assert(numel(angle) == 1 || size(Exy, 2) == 1 ...
          || numel(angle) == size(Exy, 2), ...
          'angle must be scalar, Exy must be scalar, or size must match');

      Etp = [1,1i; 1,-1i] * Exy./2;
      Etp(1, :) = Etp(1, :) .* sign(cos(angle(:).'));
    end
  end

  methods
    function beam = Bessel(varargin)
      % Construct Bessel properties
      %
      % Usage
      %   beam = beam@Bessel(angle, field, lmode, ...)
      %   Parameters can also be passed as named inputs.
      %
      % Parameters
      %   - angle (N numeric) -- Incoming angle (radians).
      %
      %   - field (2xN numeric) -- Field in spherical coordinates.
      %
      %   - lmode (N numeric) -- Azimuthal mode number.
      %     Default: ``0``.

      p = inputParser;
      p.addOptional('angle', [], @isnumeric);
      p.addOptional('field', [], @isnumeric);
      p.addOptional('lmode', 0, @isnumeric);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.properties.Material(unmatched{:});
      beam = beam.setData(p.Results.angle, p.Results.field, p.Results.lmode);
    end

    function beam = setData(beam, angle, field, lmode)
      % Set beam data
      %
      % Usage
      %   beam = beam.setData(angle, field, lmode)
      %
      % Parameters
      %   - angle (N numeric) -- Incoming angle (radians).
      %
      %   - field (2xN numeric) -- Field in spherical coordinates.
      %
      %   - lmode (N numeric) -- Azimuthal mode number.

      ott.utils.nargoutCheck(beam, nargout);

      assert(isnumeric(angle) && isvector(angle), ...
          'angle must be numeric vector');
      assert(isnumeric(lmode) && isvector(lmode), ...
          'lmode must be numeric vector');
      assert(isnumeric(field) && ismatrix(field) && size(field, 1) == 2, ...
          'field must be 2xN numeric matrix');

      angle = angle(:).';
      lmode = lmode(:).';

      Ntheta = numel(angle);
      Npol = size(field, 2);
      Nlmode = numel(lmode);
      Nwork = max([Ntheta, Npol, Nlmode]);

      assert(Ntheta == 1 || Ntheta == Nwork, ...
          'angle must be scalar or inputs must have matching size');
      assert(Npol == 1 || Npol == Nwork, ...
          'field must be scalar or inputs must have matching size');
      assert(Nlmode == 1 || Nlmode == Nwork, ...
          'lmode must be scalar or inputs must have matching size');

      % Duplicate as needed
      if Ntheta == 1, angle = repmat(angle, Nwork, 1); end
      if Npol == 1, field = repmat(field, Nwork, 1); end
      if Nlmode == 1, lmode = repmat(lmode, Nwork, 1); end

      beam.angle = angle;
      beam.field = field;
      beam.lmode = lmode;
    end
  end
end
