classdef SpharmData
% Stores pre-computed VSWF field data, useful for repeated field calculations.
%
% Properties
%   - uphi            -- Unique list of phi for this data set
%   - utheta          -- Unique list of theta for this data set
%   - uci             -- Combined index for stored values
%
%   - expimphi        -- Pre-computed exp(1i*m.*phi)
%   - Y               -- Pre-computed scalar spherical harmonic
%   - Ytheta          -- Pre-computed theta angular derivative
%   - Yphi            -- Pre-computed phi angular derivative
%
% Methods
%   - VswfData        -- Constructor
%   - evaluate        -- Evaluate Y, Ytheta, and Yphi
%   - updateCi        -- Update mode indices in data set
%   - updatePhi       -- Update phi values in data set
%   - updateTheta     -- Update theta values in data set

% Copyright 2020 Isaac Lenton (aka ilent2)
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    uphi             % Unique list of phi for this data set
    utheta           % Unique list of theta for this data set
    um               % m index for stored values (expimphi values)
    uci              % Combined index for stored values (Y values)

    expimphi         % Pre-computed exp(1i*m.*phi)
    Y                % Pre-computed scalar spherical harmonic
    Ytheta           % Pre-computed theta angular derivative
    Yphi             % Pre-computed phi angular derivative
  end

  methods
    function data = SpharmData()
      % Initialise and empty data set

      data.uphi = [];
      data.utheta = [];
      data.uci = [];

      data.expimphi = [];
      data.Y = [];
      data.Ytheta = [];
      data.Yphi = [];
    end

    function [Y, Ytheta, Yphi, data] = evaluate(data, ci, theta, phi)
      % Evaluate scalar harmonics and update stored data
      %
      % Usage
      %   [Y, Ytheta, Yphi, data] = data.evaluateYtp(ci, theta, phi)
      %
      % Parameters
      %   - ci (N numeric) -- Combined indices for modes to evaluate.
      %   - theta, phi (M numeric) -- Angular coordinates for point.
      %
      % Returns
      %   - Y, Ytheta, Yphi (NxM numeric) -- Spherical harmonics/derivatives.

      assert(isvector(ci) && isnumeric(ci), ...
          'ci must be numeric vector');
      assert(isvector(theta) && isnumeric(theta), ...
          'theta must be numeric vector');
      assert(isvector(phi) && isnumeric(phi), ...
          'phi must be numeric vector');

      % Update internal data for theta/phi/ci
      data = data.updateCi(ci);
      data = data.updateTheta(theta);
      data = data.updatePhi(phi);
      
      [phi, theta] = ott.utils.matchsize(phi(:), theta(:));
      keep = ~isnan(phi) & ~isnan(theta);

      % Get indices for elements
      [~, dci] = ismember(ci, data.uci);
      [~, dtheta] = ismember(theta, data.utheta);
      [~, dphi] = ismember(phi, data.uphi);
      [~, m] = ott.utils.combined_index(ci);
      [~, dm] = ismember(m, data.um);
      
      % Allocate memory for output
      Y = zeros(numel(ci), numel(phi));
      Ytheta = Y;
      Yphi = Y;
      
      % Assign nan to invalid values
      Y(:, ~keep) = nan;
      Ytheta(:, ~keep) = nan;
      Yphi(:, ~keep) = nan;

      % Retrieve stored data
      oexpimphi = data.expimphi(dm, dphi(keep));
      Y(:, keep) = data.Y(dci, dtheta(keep)) .* oexpimphi;
      Ytheta(:, keep) = data.Ytheta(dci, dtheta(keep)) .* oexpimphi;
      Yphi(:, keep) = data.Yphi(dci, dtheta(keep)) .* oexpimphi;
    end

    function data = updateCi(data, ci)
      % Ensure specified ci are in data set, otherwise compute them
      %
      % Usage
      %   data = data.updateCi(ci)

      ott.utils.nargoutCheck(data, nargout);

      % Find which um and uci need to be updated
      ci = unique(ci);
      [~, dci] = ismember(ci, data.uci);
      ci_new = ci(dci == 0);
      [~, m] = ott.utils.combined_index(ci_new);
      m = unique(m);
      [~, dm] = ismember(m, data.um);
      m_new = m(dm == 0);

      % Update expimphi
      if ~isempty(m_new)
        data.um = [data.um, m_new(:).'];
        if ~isempty(data.uphi)
          data.expimphi = [data.expimphi; exp(1i*m_new(:) .* data.uphi)];
        end
      end

      % Update Y/dtY/dpY
      if ~isempty(ci_new)
        [nn, mm] = ott.utils.combined_index(ci_new);
        N = size(data.Y, 1);
        [unn, ~, vv] = unique(nn);

        data.uci = [data.uci, ci_new(:).'];
        data.Y = [data.Y; zeros(numel(ci_new), size(data.Y, 2))];
        data.Ytheta = [data.Ytheta; zeros(numel(ci_new), size(data.Y, 2))];
        data.Yphi = [data.Yphi; zeros(numel(ci_new), size(data.Y, 2))];
        for ii = unn(:).'
          ov = find(ii == nn);

          [oY, oYtheta, oYphi] = ott.utils.spharm(ii, mm(ov), ...
              data.utheta, zeros(size(data.utheta)));

          data.Y(N+ov, :) = oY.';
          data.Ytheta(N+ov, :) = oYtheta.';
          data.Yphi(N+ov, :) = oYphi.';
        end
      end
    end

    function data = updatePhi(data, phi)
      % Ensure specified phi are in data set, otherwise calculate them
      %
      % Usage
      %   data = data.updatePhi(phi)

      ott.utils.nargoutCheck(data, nargout);

      % Find phis to update
      phi(isnan(phi)) = [];
      phi = unique(phi);
      [~, dphi] = ismember(phi, data.uphi);
      new_phi = phi(dphi == 0);

      % Update expimphi
      if ~isempty(new_phi)
        data.uphi = [data.uphi, new_phi(:).'];
        data.expimphi = [data.expimphi, exp(1i*data.um.' .* new_phi(:).')];
      end

    end

    function data = updateTheta(data, theta)
      % Ensure specified theta are in data set, otherwise calculate them
      %
      % Usage
      %   data = data.updateTheta(theta)

      ott.utils.nargoutCheck(data, nargout);

      % Find phis to update
      theta(isnan(theta)) = [];
      theta = unique(theta);
      [~, dtheta] = ismember(theta, data.utheta);
      new_theta = theta(dtheta == 0);

      % Update spharm
      if ~isempty(new_theta)
        data.utheta = [data.utheta, new_theta(:).'];
        N = size(data.Y, 2);
        data.Y = [data.Y, zeros(size(data.Y, 1), numel(new_theta))];
        data.Ytheta = [data.Ytheta, zeros(size(data.Y, 1), numel(new_theta))];
        data.Yphi = [data.Yphi, zeros(size(data.Y, 1), numel(new_theta))];

        [n, m] = ott.utils.combined_index(data.uci);

        for nn = unique(n(:).')
          vv = find(nn == n);

          [oY, oYtheta, oYphi] = ott.utils.spharm(nn, m(vv), ...
              new_theta, zeros(size(new_theta)));

          data.Y(vv, N+1:end) = oY.';
          data.Ytheta(vv, N+1:end) = oYtheta.';
          data.Yphi(vv, N+1:end) = oYphi.';
        end
      end
    end
  end
end

