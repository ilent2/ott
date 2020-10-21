classdef VswfData
% A structure for storing VSWF data used for repeated field calculations.
%
% Properties
%   - sbesselj    -- Data for regular basis
%   - sbesselh1   -- Data for outgoing basis
%   - sbesselh2   -- Data for incoming basis
%   - spharm      -- Data for spherical harmonics
%
% Methods
%   - evaluateBessel -- Evaluate hn and dhn for near-field calculations
%   - evaluateYtp    -- Evaluate spherical harmonics

% Copyright 2020 Isaac Lenton (aka ilent2)
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    sbesselj      % Data for regular basis
    sbesselh1     % Data for outgoing basis
    sbesselh2     % Data for incoming basis
    spharm        % Data for spherical harmonics
  end

  methods
    function data = VswfData()
      % Construct a new VSWF data storage structure
      %
      % Usage
      %   data = VswfData()

      data.sbesselj = ott.utils.BesselData(@ott.utils.sbesselj);
      data.sbesselh1 = ott.utils.BesselData(@ott.utils.sbesselh1);
      data.sbesselh2 = ott.utils.BesselData(@ott.utils.sbesselh2);
      data.spharm = ott.utils.SpharmData();
    end

    function [hn, dhn, data] = evaluateBessel(data, n, kr, basis)
      % Evaluate Bessel function for specified basis
      %
      % Usage
      %   [hn, dhn, data] = data.evaluateBessel(kr, basis)

      switch basis
        case 'incoming'
          [hn, dhn, data.sbesselh2] = data.sbesselh2.evaluate(n, kr);
          hn = hn ./ 2;
          dhn = dhn ./ 2;
        case 'outgoing'
          [hn, dhn, data.sbesselh1] = data.sbesselh1.evaluate(n, kr);
          hn = hn ./ 2;
          dhn = dhn ./ 2;
        case 'regular'
          [hn, dhn, data.sbesselj] = data.sbesselj.evaluate(n, kr);
        otherwise
          error('ott:utils:VswfData:evaluateBessel:unknown_basis', ...
              'Unknown value for basis parameter');
      end
    end

    function [Y, Ytheta, Yphi, data] = evaluateYtp(data, ci, theta, phi)
      % Evaluate scalar spherical harmonics and update stored data
      %
      % Usage
      %   [Y, Ytheta, Yphi, data] = data.evaluateYtp(ci, theta, phi)

      [Y, Ytheta, Yphi, data.spharm] = data.spharm.evaluate(ci, theta, phi);
    end
  end
end

