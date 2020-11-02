classdef Mathieu < ott.beam.BscInfinite & ott.beam.properties.Mathieu
% Construct a VSWF representation of a Mathieu beam.
% Inherits from :class:`+ott.+beam.BscInfinite` and
% :class:`+ott.+beam.+properties.Mathieu`.
%
% Properties
%   - theta       -- Annular angle [radians]
%   - morder      -- Mathieu beam mode number
%   - ellipticity -- Ellipticity of Mathieu beam
%   - parity      -- Parity of beam ('even' or 'odd')

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function beam = Mathieu(varargin)
      % Construct a new VSWF Mathieu beam representation.
      %
      % Usage
      %   beam = Mathieu(theta, morder, ellipticity, parity, ...)
      %
      % Optional named arguments
      %   - theta (numeric) -- Annular beam angle from -z direction (radians).
      %     Default: ``pi/4``.
      %
      %   - morder (numeric) -- Mathieu beam order.  Default: ``0``.
      %
      %   - ellipticity (numeric) -- Ellipticity of Mathieu beam.
      %     Default: ``1.0``.
      %
      %   - parity (enum) -- Parity of beam.  Either 'even' or 'odd'
      %     Default: ``'even'``.
      %
      %   - initial_Nmax (numeric) -- Initial beam Nmax.  Default: ``0``.
      %     This parameter automatically grows when the beam is used.
      %     See also :meth:`recalculate` and :meth:`getData`.

      p = inputParser;
      p.addOptional('theta', pi/4, @isnumeric);
      p.addOptional('morder', 0, @isnumeric);
      p.addOptional('ellipticity', 1, @isnumeric);
      p.addOptional('parity', 'even', ...
          @(x) sum(strcmpi(x, {'even', 'odd'})) == 1);
      p.addParameter('initial_Nmax', 0, @isnumeric);
      p.parse(varargin{:});

      beam.theta = p.Results.theta;
      beam.morder = p.Results.morder;
      beam.ellipticity = p.Results.ellipticity;
      beam.parity = p.Results.parity;
      beam = beam.recalculate(p.Results.initial_Nmax);
    end
  end

  methods (Hidden)
    function [data, vswfData] = recalculateInternal(beam, Nmax, vswfData)
      % Re-calculate BSC data for specified Nmax.

      [data, vswfData] = ott.bsc.Annular.FromMathieu(Nmax, ...
          beam.theta, beam.morder, beam.ellipticity, ...
          beam.parity, 'data', vswfData);
    end
  end
end


